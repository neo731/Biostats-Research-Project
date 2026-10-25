library(tidyverse)
library(rstan)
library(ggplot2)
library(magick)
library(truncnorm)

# all possible combinations of cohort structure
n_grid <- expand.grid(
  replicate(7, c(5,10), simplify = FALSE)
)

#q_grid <- expand.grid(
#  replicate(7, c(0.6, 0.7, 0.8), simplify = FALSE)
#)

names(n_grid) <- paste0("n", 1:7)
#names(q_grid) <- paste0("q", 1:7)

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

#q_grid <- q_grid %>%
#  rowwise() %>%
#  filter(
#    !is.unsorted(c_across(q1:q7)) # check if the row is sorted in ascending order
#  ) %>% 
#  ungroup()

# design matrix for all the possible combinations of design structure with the given constraints
# design_matrix <- merge(n_grid, q_grid)
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

calc_ID <- function(a, b) {
  b * ((0.3^-1/a) - 1)
}


posterior_extraction <- function(fit, dose_grid){
  
  post <- rstan::extract(fit)
  
  a_samples <- post$a
  b_samples <- post$b
  
  p_samples <- sapply(dose_grid, function(d) p_inf(d, a_samples, b_samples))
  
  # getting the posterior mean for probability of infection at each dose level
  p_mean <- colMeans(p_samples)
  
  # getting the quantiles to construct 95% credible intervals later for the graph
  p_lower <- apply(p_samples, 2, quantile, 0.025)
  p_upper <- apply(p_samples, 2, quantile, 0.975)
  q1 <- apply(p_samples, 2, quantile, 0.25)
  q3 <- apply(p_samples, 2, quantile, 0.75)
  
  
  data.frame(
    dose=dose_grid,
    mean=p_mean,
    lower=p_lower,
    upper=p_upper,
    q1=q1,
    q3=q3
  )
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



expected_info_one_dose <- function(post, cohort_size, doses, dose_idx, n_sim=500, y, n, seed=42) {
  
  set.seed(42)
  
  # extracting posterior prob for the dose
  p_dose <- post$p[, dose_idx]
  
  # simulate 500 times to observe infection outcome
  simulated_y <- rbinom(
    n_sim,
    size = cohort_size,
    prob = sample(
      p_dose, n_sim, replace = TRUE
    )
  )
  
  # averaging frequency of y with total number of simulation to be used as weights
  y_freq <- table(simulated_y) / n_sim
  
  info_by_y <- numeric(length(y_freq))
  
  y_values <- as.numeric(names(y_freq))
  
  for (i in seq_along(y_values)) {
    
    # reset within each simulation
    n_new <- n
    y_new <- y
    
    # update posterior
    y_new[dose_idx] <- y_new[dose_idx] + y_values[i]
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
    
    theta <- cbind(
      log_id70 = new_post$log_id70, 
      a = new_post$a
    )
    
    info_by_y[i] <- 1/det(cov(theta))
    
  }
  
  # summing up all information gained by each possible infection outcome with their respective weights, divided by cohort size
  # to scale the expected info
  return(sum(as.numeric(y_freq) * info_by_y))
  
}



# Just look at which candidate dose has the highest expected information gain

decision_rule <- function(fit, doses, cohort_size, y, n, current_idx) {
  
  post <- rstan::extract(fit)
  
  
  # 1. P(0.6 < p_inf < 0.8) for every dose and only retain doses that satisfy q
  target_prob <- sapply(seq_along(doses), function(d) {
    mean(post$p[, d] >= 0.6 & post$p[, d] <= 0.8)
  })
  
  expected_info <- rep(NA_real_, length(doses))

  # 2. calculate expected information only for eligible doses
  
  for (i in seq_along(doses)) {
    
    expected_info[i] <- expected_info_one_dose(
      post = post, 
      cohort_size = cohort_size, 
      doses = doses, 
      dose_idx = i, 
      n_sim=500, 
      y = y, 
      n = n
    )
  }
  
  # 3. Combine utility by multiplying target probability and expected info
  # expected info gain is potentially a large value which will dominate target prob. So log is taken
  utility <- log(expected_info) + log(target_prob)
  
  # choose the dose with highest utility
  next_idx <- which.max(utility)
  
  
  
  list(
    next_idx = next_idx,
    target_prob = target_prob,
    expected_info = expected_info,
    utility = utility
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


# initially use 10 simulations per design
n_sim_design <- 10

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


# store actual cohort path data from all simulations (from actual dose level, 10^(4:7))
all_design_paths <- list()
path_counter <- 1


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
    current_idx <- 1
    
    # starting info
    current_expected_info = 0
    
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
      
      if (sum(n) + cohort_size > max_sample) {
        break
      }
      
      
      # Treat current cohort
      current_dose <- doses[current_idx]
      
      # number of observed infected individuals at the current dose/cohort
      infected <- rbinom(1, cohort_size, true_p[current_idx])
      
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
      
      stop_reason <- NA_character_
      
      if (n[current_idx] >= 20) {
        stop_reason <- "20 participants at target dose"
      } else if (total_enrolled == max_sample) {
        stop_reason <- "35 participants enrolled"
      }
      
      # if continues, determine the next n and q for next cohort
      if (is.na(stop_reason) && cohort < max_cohorts) {
        
        next_cohort_size <- as.numeric(
          design_row[[paste0("n", cohort + 1)]]
        )
        
        # Check that next cohort exists under this design
        if (
          !is.na(next_cohort_size) &&
          total_enrolled + next_cohort_size <= max_sample) {
          
          # ----------------------------------------------
          # Decision for NEXT cohort
          #
          # 1. Calculate P(0.6 < p < 0.8) for all doses
          # 2. Calculate expected information for all doses
          # 3. Pick dose with highest expected information
          # ----------------------------------------------
          
          decision <- decision_rule(
            fit = fit,
            doses = doses,
            cohort_size = next_cohort_size,
            y = y,
            n = n,
            current_idx = current_idx
          )
          
          next_idx <- decision$next_idx
          
          selected_expected_info <- decision$expected_info[next_idx]
          
          selected_target_prob <- decision$target_prob[next_idx]
          
          selected_utility <- decision$utility
          
          
        } else {
          
          next_idx <- NA_integer_
          selected_expected_info <- NA_real_
          
          stop_reason <- "Next cohort exceeds maximum sample"
        }
        
      } else {
        
        next_idx <- NA_integer_
        selected_expected_info <- NA_real_
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
        current_expected_info = current_expected_info,
        selected_expected_info = selected_expected_info,
        target_prob = selected_target_prob,
        utility = selected_utility
      )
      design_counter <- design_counter + 1
      
      # stop trial early if stopping rule reached
      if (!is.na(stop_reason)) {
        break
      }
      
      # otherwise update for the next cohort
      current_idx <- next_idx
      current_expected_info = selected_expected_info
    }
  }
}

design_results_df <- bind_rows(design_results)


write.csv(design_results_df, file = "simulation result/expected info/optimisation_results.csv", row.names = FALSE)
write.csv(true_scenarios, file = "simulation result/expected info/true_priors.csv", row.names = FALSE)