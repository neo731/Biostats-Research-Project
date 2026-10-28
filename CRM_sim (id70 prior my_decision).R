library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(magick)
library(truncnorm)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd('')
crm_stan2 <- stan_model("code/CRM stan 2.stan")


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

# The function below is to implement the decision rules. This is to determine the dose level that the next cohort should be inoculating in the simulation, based on the current observation for infection from the current dose/cohort. It uses posterior p at ACTUAL trial doses only. Here are the concepts of how decisions are being updated:
  
# 1. If ID70 is more than 2 log10 units above the current dose, we escalate the next dose by 2 log10 units. Move on to the next check if this is not met.

# 2. We look at the proportion of infection probability at the current dose that exceeds 0.6, if the percentage is lower than 75%, escalate by one log10 unit, otherwise, move on to the next check.

# 3. Now look at the proportion of infection probability at the current dose that exceeds 0.8, if the percentage is higher than 95%, de-escalate by one log10 unit.

# 4. If any of the conditions are not met above, maintain the current dose for the next cohrot.

decision_rule <- function(current_idx, fit, doses, id70) {
  
  current_dose <- doses[current_idx]
  
  post <- rstan::extract(fit)
  
  # extracting all posterior probabilities
  p_samples <- post$p
  
  if (log10(id70) >= (2 + log10(current_dose))) {
    next_idx <- min(current_idx + 2, length(doses))
  }
  else {
    # If there is less than a 75% posterior probability of infection at the
    # current dose exceeds 0.6, escalate by one log10 unit
    if (mean(p_samples[, current_idx] > 0.6) < 0.75) {
      next_idx <- min(current_idx + 1, length(doses))
    }
    # If there is greater than a 95% posterior probability of infection at the
    # current dose exceeds 0.8, de-escalate by one log10 unit.
    else if (mean(p_samples[, current_idx] > 0.8) > 0.95) {
      next_idx <- max(current_idx - 1, 1)
    }
    else {
      next_idx <- current_idx
    }
  }
  
  return(next_idx)
}

# Number of CPU cores available
available_cores <- parallel::detectCores(logical = TRUE)

# Leave two cores available for other tasks
n_workers <- max(1, available_cores - 2)

# Configure parallel workers
plan(multisession, workers = n_workers)



# simulation begins here
set.seed(42)
n_sim <- 1000

# actual trial doses
doses <- 10^(3:7)

# smooth dose grid for plotting only
dose_grid <- 10^seq(1, 9, length.out = 101)

#n_cohorts <- length(doses)
cohort_size <- 5
max_cohorts <- 35 / cohort_size

# store all smooth posterior summaries from all simulations/cohorts
all_results <- list()
row_counter <- 1

# store actual cohort path data from all simulations (from actual dose level, 10^(4:8))
all_paths <- list()
path_counter <- 1

id70_results <- list()

# Store true parameters for each simulation
true_parameters_results <- vector("list", n_sim)


for (i in 1:n_sim) {
  
  # reset within each simulation
  n <- rep(0L, length(doses))
  y <- rep(0L, length(doses))
  
  # picking out random values from the proposed prior distribution to act as a true prior,
  # so that we can produce a 'true' probability of infection at each dose
  true_log10_id70 <- rtruncnorm(1, a = 3, b = 7, mean = 4, sd = 0.8)
  true_log_a <- runif(1, log(0.3), log(0.8))
  
  true_a <- exp(true_log_a)
  true_id70 <- 10^true_log10_id70
  true_b <- true_id70 / (0.3^(-1/true_a) - 1)
  
  true_p <- p_inf(doses, true_a, true_b)
  
  true_parameters_results[[i]] <- data.frame(
    sim_id = i,
    true_log10_id70 = true_log10_id70,
    true_id70 = true_id70,
    true_log_a = true_log_a,
    true_a = true_a,
    true_b = true_b
  )
  
  # starting dose
  current_idx <- 2
  
  for(cohort in 1:max_cohorts){
    
    current_dose <- doses[current_idx]
    
    # number of observed infected individuals at the current dose/cohort
    infected <- rbinom(1, cohort_size, true_p[current_idx])
    
    n[current_idx] <- n[current_idx] + cohort_size
    y[current_idx] <- y[current_idx] + infected
    
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
    
    # extracting id70 results from stan
    id70 <- id70_extraction(fit)
    id70$sim_id <- i
    id70$cohort <- cohort
    id70$current_dose <- current_dose
    
    id70_results[[row_counter]] <- id70
    
    
    # extracting results to construct smooth plot
    df <- posterior_extraction(fit, dose_grid)
    
    # add identifiers and current trial state to the df
    df$sim_id <- i
    df$cohort <- cohort
    df$current_dose <- current_dose
    
    all_results[[row_counter]] <- df
    row_counter <- row_counter + 1
    
    # this part is to add into the df later, to check if the implemented decision rules were correctly executed
    post <- rstan::extract(fit)
    p_samples <- post$p
    
    p_gt_06 <- mean(p_samples[, current_idx] > 0.6)
    p_gt_08 <- mean(p_samples[, current_idx] > 0.8)
    
    # work out next dose
    next_idx <- decision_rule(current_idx, fit, doses, id70$mean)
    
    # stopping criteria
    #target_idx <- which.min(abs(doses - id70$mean)) # target dose = the dose that is closest to ID70?
    total_enrolled <- sum(n)
    
    stop_reason <- NA_character_
    
    if (n[current_idx] == 20) {
      stop_reason <- "20 participants at target dose"
    } else if (total_enrolled == 35) {
      stop_reason <- "35 participants enrolled"
    }
    
    # store actual path / decision information
    all_paths[[path_counter]] <- data.frame(
      sim_id = i,
      cohort = cohort,
      current_dose = current_dose,
      infected = infected,
      cum_n_at_dose = n[current_idx],
      cum_y_at_dose = y[current_idx],
      p_gt_06 = p_gt_06,
      p_gt_08 = p_gt_08,
      next_dose = doses[next_idx]
    )
    path_counter <- path_counter + 1
    
    # stop trial early if stopping rule reached
    if (!is.na(stop_reason)) {
      break
    }
    
    # otherwise update for the next cohort
    current_idx <- next_idx
    
  }
  
  cat(
    "\nSim",
    i,
    "completed successfully.\n"
  )
  
}

# combining all results together into a single dataframe
all_posteriors <- bind_rows(all_results)
all_paths_df <- bind_rows(all_paths)
all_paths_df$prob <- all_paths_df$infected / cohort_size
id70_df <- bind_rows(id70_results)
true_parameters_df <- bind_rows(true_parameters_results)

# saving all the outputs into csv files for markdown files later
write.csv(all_posteriors, file = "simulation result/old decision rule/posterior_sim_result_3.csv", row.names = FALSE)
write.csv(all_paths_df, file = "simulation result/old decision rule/all_paths_3.csv", row.names = FALSE)
write.csv(id70_df, file = "simulation result/old decision rule/id70_3.csv", row.names = FALSE)
write.csv(true_parameters_df, "simulation result/old decision rule/1000_true_params.csv", row.names = FALSE)