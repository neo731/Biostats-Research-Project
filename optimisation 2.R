library(tidyverse)
library(rstan)
library(ggplot2)
library(magick)
library(truncnorm)
library(profvis)

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
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


setwd('C:/Users/Neo/OneDrive/Documents/Master of Biostats/Research project')
crm_stan2 <- stan_model("code/CRM stan.stan")


p_inf <- function(D, a, b){
  # beta-poisson
  1 - (1 + D/b)^(-a)
}



expected_info_one_dose <- function(cohort_size, doses, dose_idx, y, n, true_p) {
  
  # All possible numbers infected in the next cohort
  y_values <- 0:cohort_size
  
  # Probability of each possible outcomes using True P
  prob_outcome <- dbinom(
    x = y_values,
    size = cohort_size,
    prob = true_p[dose_idx]
  )
  
  # storing the values
  target_prob_by_y <- numeric(length(y_values))
  info_by_y <- numeric(length(y_values))
  utility_by_y <- numeric(length(y_values))
  
  
  # Evaluate every possible hypothetical outcome
  for (i in seq_along(y_values)) {
    
    hyp_y <- y_values[i]
    
    #
    # 1. Add hypothetical outcome to current data and re-fit Stan
    #
    
    # reset within each simulation
    n_new <- n
    y_new <- y
    
    # update posterior
    y_new[dose_idx] <- y_new[dose_idx] + hyp_y
    n_new[dose_idx] <- n_new[dose_idx] + cohort_size
    
    # refitting stan for posterior for each unique outcome
    fit_new <- sampling(
      crm_stan2,
      data=list(
        K=length(doses),
        dose=doses,
        n=n_new,
        y=y_new
      ),
      chains=4,
      iter=2000,
      warmup=1000,
      refresh=0
    )
    
    # extract new posterior
    new_post <- rstan::extract(fit_new)
    
    #
    # 2. Posterior covariance after hypothetical outcome (Information component)
    #
    theta <- cbind(
      log_id70 = new_post$log_id70, 
      log_a = new_post$log_a
    )
    
    # information conditional on this outcome (want to maximise this)
    info_by_y[i] <- 1/det(cov(theta))
    
    #
    # 3. Target prob from hypothetical posterior (want to maximise this)
    #
    target_prob_by_y[i] <- mean(new_post$p[, dose_idx] >= 0.6 & new_post$p[, dose_idx] <= 0.8)
    
    # Utility conditional on hypothetical outcome
    utility_by_y[i] <- exp(log(info_by_y[i]) + log(target_prob_by_y[i]))
    
  }
  
  # summing up all information gained by each possible infection outcome with their respective weights
  expected_utility <- sum(prob_outcome * utility_by_y)
  
  return(expected_utility)
  
}



# Just look at which candidate dose has the highest expected information gain

decision_rule <- function(doses, cohort_size, y, n, true_p) {
  
  # store values
  expected_utility <- rep(NA_real_, length(doses))
  
  # 1. calculate expected utility for all doses
  
  for (i in seq_along(doses)) {
    
    expected_utility[i] <- expected_info_one_dose(
      cohort_size = cohort_size, 
      doses = doses, 
      dose_idx = i,
      y = y, 
      n = n,
      true_p = true_p
    )
    
  }

  # 2. choose the dose with highest expected utility
  next_idx <- which.max(expected_utility)

  
  list(
    next_idx = next_idx,
    expected_utility = expected_utility
  )
}

# simulation begins here
set.seed(42)

# actual trial doses
doses <- 10^(4:7)

# smooth dose grid for plotting only
dose_grid <- 10^seq(1, 9, length.out = 101)

#n_cohorts <- length(doses)
cohort_size <- 5
max_sample <- 35
max_cohorts <- max_sample / cohort_size


# initially use 1 simulations per design
n_sim_design <- 1

# Generate scenarios once before the design loop so same set of true scenarios across every design
true_scenarios <- data.frame(
  sim_id = 1:n_sim_design,
  true_log10_id70 = rtruncnorm(
    n_sim_design,
    a = 3,
    b = 7,
    mean = 4,
    sd = 0.8
  ),
  true_log_a = runif(
    n_sim_design,
    log(0.3),
    log(0.8)
  )
)

true_scenarios <- true_scenarios %>%
  mutate(
    true_a = exp(true_log_a),
    true_id70 = 10^true_log10_id70,
    true_b = true_id70 / (0.3^(-1 / true_a) - 1)
  )

# store the outcome sequence for each simulation so each design comes with the same
# infection sequence for a fairer comparison (all design use the same random infection numbers)
potential_outcomes <- vector(
  "list",
  n_sim_design
)


# Main object for comparing designs
design_results <- list()
design_counter <- 1



# ============================================================
# LOOP 1: EACH DESIGN
# ============================================================

for (design_id in 1:nrow(design_matrix)) {
  
  cat(
    "Running design",
    design_id,
    "\n"
  )
  
  # Extract this design
  design_row <- design_matrix[design_id, ]
  
  set.seed(42)
  
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
    
    potential_outcomes[[sim_id]] <- lapply(
      true_p,
      function(prob) {
        # generate a 35 binary sequence for all 4 doses, so later can extract how many are infected
        # in total at that cohort at that dose. This makes every design to have same number of infected 
        # people in each cohort so it's fair for performance comparison
        rbinom(max_sample, size = 1, prob = prob)
        
      }
    )
    
    
    # starting dose
    current_idx <- 1
    
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
      
      # Treat current cohort
      current_dose <- doses[current_idx]
      
      # Select the next unused potential patients at this dose
      start_id <- n[current_idx] + 1
      end_id   <- n[current_idx] + cohort_size
      
      # number of infected patients at the current cohort
      infected <- sum(
        potential_outcomes[[sim_id]][[current_idx]][start_id:end_id]
      )
      
      # infected <- rbinom(1, cohort_size, true_p[current_idx])
    
      
      n[current_idx] <- n[current_idx] + cohort_size
      y[current_idx] <- y[current_idx] + infected
      
      # update posterior
      stan_data <- list(
        K=length(doses),
        dose=doses,
        n=n,
        y=y
      )
      
      fit <- sampling(
        crm_stan2,
        data=stan_data,
        chains=4,
        iter=2000,
        warmup=1000,
        refresh=0
      )
      
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
            y = y,
            n = n,
            true_p = true_p
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
  }
}

design_results_df <- bind_rows(design_results)


write.csv(design_results_df, file = "simulation result/expected info/optimisation_results_2.csv", row.names = FALSE)
write.csv(true_scenarios, file = "simulation result/expected info/true_priors.csv", row.names = FALSE)