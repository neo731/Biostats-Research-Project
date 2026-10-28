library(tidyverse)
library(rstan)
library(ggplot2)
library(magick)
library(truncnorm)
library(profvis)
library(foreach)
library(doParallel)

######################
# construction of design matrix with cohort size of 5 and 10
######################

# all possible combinations of cohort structure
n_grid <- expand.grid(
  replicate(7, c(5,10), simplify = FALSE)
)

names(n_grid) <- paste0("n", 1:7)

n_grid <- n_grid %>%
  rowwise() %>%
  filter(
    !is.unsorted(c_across(n1:n7)) # check if the row is sorted in ascending order
  ) %>% 
  ungroup()

for (i in 1:nrow(n_grid)) {
  
  count = 0 
  max_sample = 35
  
  for (j in 1:ncol(n_grid)) {
    count = count + n_grid[i, j]
    
    if (count > max_sample) {
      n_grid[i, j] = NA
    }
  }
}

# design matrix for all the possible combinations of design structure with the given constraints
design_matrix <- n_grid

# CRM setup same as before
rstan_options(auto_write = TRUE)

crm_stan2 <- stan_model("code/CRM stan.stan")


# functions
p_inf <- function(D, a, b){
  # beta-poisson
  1 - (1 + D/b)^(-a)
}

id70_extraction <- function(fit) {
  post <- rstan::extract(fit)
  id70_samples <- post$id70
  
  # the overall mean for ID70 from posterior
  id70_mean <- mean(id70_samples)
  id70_lower <- quantile(id70_samples, 0.025)
  id70_upper <- quantile(id70_samples, 0.975)
  
  data.frame(
    mean=id70_mean,
    id70_lower = id70_lower,
    id70_upper = id70_upper
  )
}

# use this function to start cluster/worker at the beginning of each design in the loop
start_stan_cluster <- function(
    n_workers = min(parallel::detectCores() - 1, 8)
) {
  
  cl <- parallel::makeCluster(n_workers)
  
  doParallel::registerDoParallel(cl)
  
  parallel::clusterEvalQ(
    cl,
    {
      library(rstan)
      
      rstan_options(auto_write = TRUE)
      
      setwd("")
      
      crm_stan2 <- stan_model("code/CRM stan.stan")
      
      NULL
    }
  )
  
  cl
}


# use foreach for all possible combinations of dose and hypotheticcal outcome to speed things up
# i.e. if cohort size is 10 then there are 11 hypothetical outcome x 4 doses = 44 pairings
decision_rule <- function(doses, cohort_size, y, n, true_p, sim_id, cohort) {
  
  # the grid for dose and hypothetical outcome combinations
  job_grid <- expand.grid(
    dose_idx = seq_along(doses),
    hyp_y = 0:cohort_size
  )
  
  
  results <- foreach(
    j = 1:nrow(job_grid),
    .packages = "rstan",
    .combine = rbind
  ) %dopar% {
    
    dose_idx <- job_grid$dose_idx[j]
    hyp_y <- job_grid$hyp_y[j]
    
    # set seed for each outcome so it is reproducible
    seed = 1000*sim_id + 100*cohort + 10*dose_idx + hyp_y
    
    
    # ---------------------------------------------
    # Hypothetical dataset
    # ---------------------------------------------
    
    n_new <- n
    y_new <- y
    
    n_new[dose_idx] <- n_new[dose_idx] + cohort_size
    
    y_new[dose_idx] <- y_new[dose_idx] + hyp_y
    
    
    # ---------------------------------------------
    # Hypothetical Stan posterior
    # ---------------------------------------------
    
    fit_new <- rstan::sampling(
      crm_stan2,
      data = list(
        K = length(doses),
        dose = doses,
        n = n_new,
        y = y_new
      ),
      chains = 4,
      iter = 2000,
      warmup = 1000,
      refresh = 0,
      cores = 1,
      seed = seed
    )
    
    new_post <- rstan::extract(fit_new)
    
    # ---------------------------------------------
    # Information
    # ---------------------------------------------
    
    theta <- cbind(
      log_id70 = new_post$log_id70,
      a = new_post$a
    )
    
    information <- 1 / det(cov(theta))
    
    # ---------------------------------------------
    # Target probability
    # ---------------------------------------------
    
    target_prob <- mean(new_post$p[, dose_idx] >= 0.6 & new_post$p[, dose_idx] <= 0.8)
    
    # ---------------------------------------------
    # Utility
    # ---------------------------------------------
    
    utility <- exp(log(information) +log(target_prob))
    
    # True probability of this hypothetical outcome
    outcome_prob <- dbinom(
      hyp_y,
      size = cohort_size,
      prob = true_p[dose_idx]
    )
    
    
    # since foreach loop returns a list
    c(
      dose_idx = dose_idx,
      hyp_y = hyp_y,
      outcome_prob = outcome_prob,
      information = information,
      target_prob = target_prob,
      utility = utility
    )
  }
  
  # combining all results from the foreach loop to a data frame
  results <- as.data.frame(results)
  
  
  # =========================================================
  # Expected utility for each candidate dose
  # =========================================================
  
  expected_utility <- sapply(
    seq_along(doses),
    function(d) {
      
      # now we map information for each dose together
      # extracting all utility and outcome prob that corresponds to the dose, then sum them together
      temp <- results[
        results$dose_idx == d,
      ]
      
      sum(
        temp$outcome_prob *
          temp$utility
      )
    }
  )
  
  
  # Select dose with maximum expected utility
  next_idx <- which.max(expected_utility)
  
  
  list(
    next_idx = next_idx,
    expected_utility = expected_utility,
    details = results
  )
}


# actual trial doses
doses <- 10^(3:7)


#n_cohorts <- length(doses)
cohort_size <- 5
max_sample <- 35
max_cohorts <- max_sample / cohort_size


# initially use 500 simulations per design
n_sim_design <- 500

set.seed(42)
id70_draws <- rtruncnorm(n_sim_design, a = 3, b = 7, mean = 4, sd = 0.8)
log_a_draws <- runif(n_sim_design, log(0.3), log(0.8))


# Generate scenarios once before the design loop so same set of true scenarios across every design
true_scenarios <- data.frame(
  sim_id = 1:n_sim_design,
  true_log10_id70 = id70_draws,
  true_log_a = log_a_draws
)

true_scenarios <- true_scenarios %>%
  mutate(
    true_a = exp(true_log_a),
    true_id70 = 10^true_log10_id70,
    true_b = true_id70 / (0.3^(-1 / true_a) - 1)
  )

# store the outcome sequence for each simulation so each design comes with the same
# infection sequence for a fair comparison (all design use the same random infection sequence)
potential_outcomes <- vector(
  "list",
  n_sim_design
)

set.seed(42)

# generating the binary sequence (length = 35) of infection for 100 simulation using the true_p
for (sim_id in 1:n_sim_design) {
  
  true_p_sim <- p_inf(
    doses,
    true_scenarios$true_a[sim_id],
    true_scenarios$true_b[sim_id]
  )
  
  potential_outcomes[[sim_id]] <- lapply(
    true_p_sim,
    function(prob) {
      
      rbinom(max_sample, size = 1, prob = prob)
    }
  )
}


# Main object for comparing designs
design_results <- list()
design_counter <- 1

final_id70_results <- list()
id70_counter <- 1


# ============================================================
# LOOP 1: EACH DESIGN
# ============================================================

for (design_id in 1:nrow((design_matrix))) {
  
  cl <- start_stan_cluster()
  
  cat(
    "Running design",
    design_id,
    "\n"
  )
  
  # Extract this design
  design_row <- design_matrix[design_id, ]
  
  # ==========================================================
  # LOOP 2: REPLICATE SIMULATIONS
  # ==========================================================
  
  for (sim_id in 1:n_sim_design) {
    
    # reset within each simulation
    n <- rep(0L, length(doses))
    y <- rep(0L, length(doses))
    
    # picking out random values from the proposed prior distribution to act as a true prior,
    # so that we can produce a 'true' probability of infection at each dose
    true_log10_id70 <- true_scenarios$true_log10_id70[sim_id]
    true_log_a <- true_scenarios$true_log_a[sim_id]
    
    true_a <- true_scenarios$true_a[sim_id]
    true_id70 <- true_scenarios$true_id70[sim_id]
    true_b <- true_scenarios$true_b[sim_id]
    
    true_p <- p_inf(doses, true_a, true_b)
    
    
    # starting dose
    current_idx <- 2
    
    # starting info
    current_expected_utility = 0
    
    # ========================================================
    # LOOP 3: COHORT
    # ========================================================
    
    for(cohort in 1:max_cohorts){
      
      # get design parameters for this cohort
      cohort_size <- as.numeric(
        design_row[[paste0("n", cohort)]]
      )
      
      
      # NA means this design cannot enrol another cohort
      # without exceeding the maximum N = 35
      if (is.na(cohort_size)) {
        break
      }
      # if existing accumulated patients plus next cohort size is greater than 35, break
      if (sum(n) + cohort_size > max_sample) {
        break
      }
      
      next_idx <- NA_integer_
      selected_expected_utility <- NA_real_
      
      # Treat current cohort
      current_dose <- doses[current_idx]
      
      # Select the next unused potential patients at this dose
      start_id <- n[current_idx] + 1L
      end_id   <- n[current_idx] + cohort_size
      
      # number of infected patients at the current cohort from the pre-generated sequence
      infected <- sum(
        potential_outcomes[[sim_id]][[current_idx]][start_id:end_id]
      )
      
      # updating the data for stan later
      n[current_idx] <- n[current_idx] + cohort_size
      y[current_idx] <- y[current_idx] + infected
      
      
      # stopping criteria
      # determine if the next cohort needs to be introduced
      total_enrolled <- sum(n)
      
      # if there is a next cohort
      if (cohort < max_cohorts) {
        
        next_cohort_size <- as.numeric(
          design_row[[paste0("n", cohort + 1)]]
        )
        
        # Check that next cohort exists under this design
        if (
          !is.na(next_cohort_size) &&
          total_enrolled + next_cohort_size <= max_sample) {
          
          # ----------------------------------------------
          # Decision for NEXT cohort
          # 1. Calculate expected utility for all doses
          # 2. Pick dose with highest utility
          # ----------------------------------------------
          
          decision <- decision_rule(
            doses = doses,
            cohort_size = next_cohort_size,
            y = y, # y and n are updated from above
            n = n,
            true_p = true_p,
            sim_id = sim_id,
            cohort = cohort
          )
          
          next_idx <- decision$next_idx
          
          selected_expected_utility <- decision$expected_utility[next_idx]
          
        } 
      }
      
      
      
      # store actual path / decision information
      design_results[[design_counter]] <- data.frame(
        design_id = design_id,
        sim_id = sim_id,
        cohort = cohort,
        cohort_size = cohort_size,
        total_enrolled = total_enrolled,
        current_dose = current_dose,
        infected = infected,
        next_dose = ifelse(
          is.na(next_idx),
          NA,
          doses[next_idx]
        ),
        current_expected_utility = current_expected_utility,
        selected_expected_utility = selected_expected_utility
      )
      
      design_counter <- design_counter + 1
      
      # End trial when this design has no further valid cohort
      if (is.na(next_idx)) {
        break
      }
      
      # otherwise update for the next cohort
      current_idx <- next_idx
      current_expected_utility = selected_expected_utility
    }
    
    final_fit <- rstan::sampling(
      crm_stan2,
      data = list(
        K = length(doses),
        dose = doses,
        n = n,
        y = y
      ),
      chains = 4,
      iter = 2000,
      warmup = 1000,
      refresh = 0,
      cores = 1)
    
    final_id70 <- id70_extraction(final_fit)
    final_id70$design_id <- design_id
    final_id70$sim_id <- sim_id
    final_id70$current_dose <- current_dose
    
    id70_results[[id70_counter]] <- final_id70
    id70_counter = id70_counter + 1
    
    
  }
  
  parallel::stopCluster(cl)
  
}


design_results_df <- bind_rows(design_results)



write.csv(design_results_df, file = "simulation result/expected_info_2/optimisation_results_design.csv", row.names = FALSE)
write.csv(true_scenarios, file = "simulation result/expected_info_2/true_priors.csv", row.names = FALSE)
