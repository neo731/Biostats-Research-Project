library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(magick)
library(truncnorm)

setwd('C:/Users/Neo/OneDrive/Documents/Master of Biostats/Research project')
crm_stan2 <- stan_model("code/CRM stan 2.stan")


set.seed(42)

# beta-poisson dose-response
p_inf <- function(D, a, b) {
  1 - (1 + D / b)^(-a)
}

# draw from prior
draw_prior_parameters <- function(n_draws = 20000) {
  log10_id70 <- rtruncnorm(
    n_draws, a = 3, b = 7, mean = 4, sd = 0.8
  )
  log_a <- runif(n_draws, log(0.3), log(0.8))
  
  a <- exp(log_a)
  id70 <- 10^log10_id70
  b <- id70 / (0.3^(-1 / a) - 1)
  
  data.frame(
    log10_id70 = log10_id70,
    id70 = id70,
    a = a,
    b = b
  )
}

prior_draws <- draw_prior_parameters(20000)

# target percentile locations
target_probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
target_id70 <- quantile(prior_draws$id70, probs = target_probs)

target_id70 # if we convert them to log10 scale it matches the tutorial paper's true ID70 values


# obtaining the set of parameters for each of the 5 scenarios
scenario_params <- lapply(seq_along(target_probs), function(j) {
  idx <- which.min(abs(prior_draws$id70 - target_id70[j]))
  prior_draws[idx, ]
})


scenario_params <- bind_rows(scenario_params) %>%
  mutate(
    scenario = paste0("Scenario_", c("05", "25", "50", "75", "95")),
    percentile = target_probs
  )


# functions to extract posterior summary from CRM simulation
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


# function that runs simulation in one iteration/trial, and saves the results
simulate_one_trial <- function(sim_id, scenario_name, true_a, true_b,
                               doses = 10^(4:7),
                               dose_grid = 10^seq(1, 9, length.out = 101),
                               cohort_size = 5,
                               max_cohorts = 35 / 5) {
  
  n <- rep(0L, length(doses))
  y <- rep(0L, length(doses))
  true_p <- p_inf(doses, true_a, true_b)
  
  current_idx <- 1
  
  all_results <- list()
  all_paths <- list()
  all_id70 <- list()
  
  for (cohort in 1:max_cohorts) {
    
    current_dose <- doses[current_idx]
    infected <- rbinom(1, cohort_size, true_p[current_idx])
    
    n[current_idx] <- n[current_idx] + cohort_size
    y[current_idx] <- y[current_idx] + infected
    
    stan_data <- list(
      K = length(doses),
      dose = doses,
      n = n,
      y = y
    )
    
    fit <- sampling(
      crm_stan2,
      data = stan_data,
      chains = 4,
      iter = 2000,
      warmup = 1000,
      refresh = 0
    )
    
    # ID70 summary
    id70 <- id70_extraction(fit)
    id70$sim_id <- sim_id
    id70$scenario <- scenario_name
    id70$cohort <- cohort
    id70$current_dose <- current_dose
    all_id70[[cohort]] <- id70
    
    # posterior curve summary on smooth grid
    df <- posterior_extraction(fit, dose_grid)
    df$sim_id <- sim_id
    df$scenario <- scenario_name
    df$cohort <- cohort
    df$current_dose <- current_dose
    all_results[[cohort]] <- df
    
    # posterior probabilities at actual doses
    post <- rstan::extract(fit)
    p_samples <- post$p
    
    p_gt_06 <- mean(p_samples[, current_idx] > 0.6)
    p_gt_08 <- mean(p_samples[, current_idx] > 0.8)
    
    next_idx <- decision_rule(current_idx, fit, doses, id70$mean)
    
    total_enrolled <- sum(n)
    
    stop_reason <- NA_character_
    if (n[current_idx] == 20) {
      stop_reason <- "20 participants at target dose"
    } else if (total_enrolled == 35) {
      stop_reason <- "35 participants enrolled"
    }
    
    all_paths[[cohort]] <- data.frame(
      sim_id = sim_id,
      scenario = scenario_name,
      cohort = cohort,
      current_dose = current_dose,
      infected = infected,
      cum_n_at_target_dose = n[current_idx],
      cum_y_at_target_dose = y[current_idx],
      total_enrolled = total_enrolled,
      next_dose = doses[next_idx],
      stop_reason = stop_reason
    )
    
    if (!is.na(stop_reason)) break
    
    current_idx <- next_idx
  }
  
  list(
    posterior = bind_rows(all_results),
    path = bind_rows(all_paths),
    id70 = bind_rows(all_id70)
  )
}


# the simulation begins here, where loop over the sequence of 5 scenarios, with 100 simulations per scenario
n_per_scenario <- 100

scenario_posteriors <- list()
scenario_paths <- list()
scenario_id70 <- list()

counter <- 1

for (s in seq_len(nrow(scenario_params))) {
  
  scen_name <- scenario_params$scenario[s]
  a <- scenario_params$a[s]
  b <- scenario_params$b[s]
  
  message("Running ", scen_name)
  
  for (rep in 1:n_per_scenario) {
    out <- simulate_one_trial(
      sim_id = rep,
      scenario_name = scen_name,
      true_a = a,
      true_b = b,
      doses = 10^(4:7),
      dose_grid = 10^seq(1, 9, length.out = 101),
      cohort_size = 5,
      max_cohorts = 35 / 5
    )
    
    scenario_posteriors[[counter]] <- out$posterior
    scenario_paths[[counter]] <- out$path
    scenario_id70[[counter]] <- out$id70
    counter <- counter + 1
  }
}

all_posteriors <- bind_rows(scenario_posteriors)
all_paths_df <- bind_rows(scenario_paths)
id70_df <- bind_rows(scenario_id70)

write.csv(all_posteriors, file = "simulation result/old decision rule/posterior_scenarios_sim.csv", row.names = FALSE)
write.csv(all_paths_df, file = "simulation result/old decision rule/all_paths_scenarios.csv", row.names = FALSE)
write.csv(id70_df, file = "simulation result/old decision rule/id70_scenarios.csv", row.names = FALSE)