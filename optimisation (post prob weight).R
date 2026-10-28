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
setwd('')
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
  id70_lower <- unname(quantile(id70_samples, 0.025))
  id70_upper <- unname(quantile(id70_samples, 0.975))
  
  data.frame(
    mean=id70_mean,
    id70_lower = id70_lower,
    id70_upper = id70_upper
  )
}

# ============================================================
# CHECKPOINT SETUP
# ============================================================

checkpoint_dir <- file.path(
  "simulation result",
  "adaptive_design_checkpoints2"
)

if (!dir.exists(checkpoint_dir)) {
  
  dir.create(
    checkpoint_dir,
    recursive = TRUE
  )
}


checkpoint_file <- function(design_id, sim_id) {
  
  file.path(
    checkpoint_dir,
    sprintf(
      "design_%02d_sim_%03d.rds",
      design_id,
      sim_id
    )
  )
}

# use this function to start cluster/worker at the beginning of each design in the loop
start_stan_cluster <- function(n_workers = min(parallel::detectCores() - 1, 8)) {
  
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

stop_stan_cluster <- function(cl) {
  
  if (!is.null(cl)) {
    
    try(
      parallel::stopCluster(cl),
      silent = TRUE
    )
  }
  
  # Switch foreach back to sequential execution
  foreach::registerDoSEQ()
  
  # Clean up any unused connections left by dead workers
  closeAllConnections()
  
  invisible(gc())
}


# use foreach for all possible combinations of dose and hypotheticcal outcome to speed things up
# i.e. if cohort size is 10 then there are 11 hypothetical outcome x 4 doses = 44 pairings
decision_rule <- function(doses, cohort_size, y, n, current_post, sim_id, cohort) {
  
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
    
    # -------------------------------------------------------
    # POSTERIOR PREDICTIVE PROBABILITY OF THIS OUTCOME
    # -------------------------------------------------------
    
    # Posterior draws of infection probability at candidate dose
    p_draws <- current_post$p[, dose_idx]
    
    # For every posterior draw:
    # P(Y = hyp_y | p_draw)
    outcome_likelihood <- dbinom(
      hyp_y,
      size = cohort_size,
      prob = p_draws
    )
    
    # Integrate over current posterior:
    # P(Y = hyp_y | current data)
    outcome_prob <- mean(outcome_likelihood)
    
    
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
    
    if (information > 0 && target_prob > 0) {
      
      utility <- exp(log(information) + log(target_prob))
      
    } else {
      
      utility <- 0
    }
    
    
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

#write.csv(true_scenarios, file = "simulation result/expected_info_2/true_priors.csv", row.names = FALSE)

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



run_one_simulation <- function(design_id, sim_id, design_row) {
  
  # ==========================================================
  # INITIALISE TRIAL
  # ==========================================================
  
  n <- rep(0L, length(doses))
  y <- rep(0L, length(doses))
  
  
  # True scenario
  true_log10_id70 <- true_scenarios$true_log10_id70[sim_id]
  
  true_log_a <- true_scenarios$true_log_a[sim_id]
  
  true_a <- true_scenarios$true_a[sim_id]
  
  true_id70 <- true_scenarios$true_id70[sim_id]
  
  true_b <- true_scenarios$true_b[sim_id]
  
  
  # Starting dose: index 2 = 10^4 with doses <- 10^(3:7)
  current_idx <- 2
  
  current_expected_utility <- 0
  
  
  # Storage for THIS simulation only
  sim_results <- list()
  selected_dose_info <- list()
  sim_counter <- 1
  post_id70 <- list()
  row_counter <- 1
  
  # ==========================================================
  # COHORT LOOP
  # ==========================================================
  
  for (cohort in 1:max_cohorts) {
    
    cohort_size <- as.numeric(
      design_row[[paste0("n", cohort)]]
    )
    
    
    # Design has no further cohort
    if (is.na(cohort_size)) {
      break
    }
    
    
    # Maximum sample size
    if (sum(n) + cohort_size > max_sample) {
      break
    }
    
    
    next_idx <- NA_integer_
    selected_expected_utility <- NA_real_
    
    
    # --------------------------------------------------------
    # Treat current cohort
    # --------------------------------------------------------
    
    current_dose <- doses[current_idx]
    
    
    start_id <- n[current_idx] + 1L
    end_id   <- n[current_idx] + cohort_size
    
    
    infected <- sum(
      potential_outcomes[[sim_id]][[current_idx]][
        start_id:end_id
      ]
    )
    
    
    # Update observed trial data
    n[current_idx] <- n[current_idx] + cohort_size
    
    y[current_idx] <- y[current_idx] + infected
    
    
    total_enrolled <- sum(n)
    
    
    # ========================================================
    # CURRENT POSTERIOR
    # ========================================================
    
    current_seed <- 500000 + 1000 * sim_id + 10 * cohort
    
    current_fit <- rstan::sampling(
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
      cores = 1,
      seed = current_seed
    )
    
    
    current_post <- rstan::extract(current_fit)
    
    current_id70 <- id70_extraction(current_fit)
    current_id70$design_id <- design_id
    current_id70$sim_id <- sim_id
    current_id70$cohort <- cohort
    current_id70$current_dose <- current_dose
    current_id70$true_log10_id70 = true_log10_id70
    current_id70$true_id70 = true_id70
    
    post_id70[[row_counter]] <- current_id70
    row_counter = row_counter + 1
    
    
    # ========================================================
    # DECISION FOR NEXT COHORT
    # ========================================================
    
    if (cohort < max_cohorts) {
      
      next_cohort_size <- as.numeric(
        design_row[[paste0("n", cohort + 1)]]
      )
      
      
      if (
        !is.na(next_cohort_size) &&
        total_enrolled + next_cohort_size <= max_sample
      ) {
        
        decision <- decision_rule(
          doses = doses,
          cohort_size = next_cohort_size,
          y = y,
          n = n,
          current_post = current_post,
          sim_id = sim_id,
          cohort = cohort
        )
        
        
        next_idx <- decision$next_idx
        
        selected_expected_utility <- decision$expected_utility[next_idx]
        
        hyp_details <- decision$details[decision$details$dose_idx == next_idx, , drop = FALSE]
      }
    }
    
    
    # ========================================================
    # SAVE CURRENT COHORT PATH
    # ========================================================
    
    sim_results[[sim_counter]] <- data.frame(
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
    
    selected_dose_info[[sim_counter]] <- data.frame(
      design_id = design_id,
      sim_id = sim_id,
      cohort = cohort,
      cohort_size = cohort_size,
      total_enrolled = total_enrolled,
      current_dose = current_dose,
      selected_outcome_prob = hyp_details$outcome_prob,
      
      selected_information = hyp_details$information,
      
      selected_target_prob = hyp_details$target_prob
      
    )
    
    
    sim_counter <- sim_counter + 1
    
    
    # ========================================================
    # END OR CONTINUE
    # ========================================================
    
    if (is.na(next_idx)) {
      break
    }
    
    
    current_idx <- next_idx
    
    current_expected_utility <- selected_expected_utility
    
    
    # Free the Stan fit itself.
    # current_post remains available for this iteration.
    rm(current_fit)
    gc()
  }
  
  
  # ==========================================================
  # COMPLETE PATH
  # ==========================================================
  
  sim_results_df <- bind_rows(sim_results)
  post_id70_df <- bind_rows(post_id70)
  selected_dose_info_df <- bind_rows(selected_dose_info)
  
  
  # Return everything from this simulated trial
  list(
    path = sim_results_df,
    id70 = post_id70_df,
    all_selected_dose_info = selected_dose_info_df
  )
}

# ============================================================
# AUTOMATIC CHECKPOINT / RECOVERY SIMULATION
# ============================================================

max_retries <- 20


for (design_id in 1:nrow(design_matrix)) {
  
  cat(
    "\n========================================\n",
    "DESIGN ", design_id, "\n",
    "========================================\n",
    sep = ""
  )
  
  
  design_row <- design_matrix[design_id, ]
  
  design_complete <- FALSE
  retry_count <- 0
  
  
  # ==========================================================
  # KEEP RESTARTING THIS DESIGN UNTIL ALL SIMULATIONS COMPLETE
  # ==========================================================
  
  while (!design_complete) {
    
    cl <- NULL
    
    
    design_complete <- tryCatch(
      {
        
        # ----------------------------------------------------
        # Fresh cluster
        # ----------------------------------------------------
        
        cl <- start_stan_cluster()
        
        # ----------------------------------------------------
        # Simulation loop
        # ----------------------------------------------------
        
        for (sim_id in 1:n_sim_design) {
          
          checkpoint <-checkpoint_file(design_id, sim_id)
          
          
          # ==================================================
          # ALREADY COMPLETED -> SKIP
          # ==================================================
          
          if (file.exists(checkpoint)) {
            
            cat(
              "Design", design_id,
              "| Simulation", sim_id,
              "| already complete - skipping\n"
            )
            
            next
          }
          
          
          cat(
            "Design", design_id,
            "| Simulation", sim_id,
            "| running...\n"
          )
          
          
          # ==================================================
          # RUN ONE COMPLETE SIMULATION
          # ==================================================
          
          result <- run_one_simulation(
            design_id = design_id,
            sim_id = sim_id,
            design_row = design_row
          )
          
          
          # ==================================================
          # SAVE ONLY AFTER SUCCESSFUL COMPLETION
          # ==================================================
          
          saveRDS(
            result,
            checkpoint
          )
          
          
          cat(
            "Design", design_id,
            "| Simulation", sim_id,
            "| saved successfully\n"
          )
          
          
          rm(result)
          gc()
        }
        
        
        # If loop reaches here, every simulation has either
        # completed now or already had a checkpoint.
        
        TRUE
      },
      
      
      # ======================================================
      # ERROR HANDLER
      # ======================================================
      
      error = function(e) {
        
        retry_count <<- retry_count + 1
        
        
        message(
          "\n----------------------------------------"
        )
        
        message(
          "Design ", design_id,
          " crashed."
        )
        
        message(
          "Error: ",
          conditionMessage(e)
        )
        
        message(
          "Retry number: ",
          retry_count
        )
        
        message(
          "The cluster will be recreated and the design ",
          "will resume from the first unfinished simulation."
        )
        
        message(
          "----------------------------------------\n"
        )
        
        
        FALSE
      },
      
      
      # ======================================================
      # ALWAYS DESTROY CURRENT CLUSTER
      # ======================================================
      
      finally = {
        
        stop_stan_cluster(cl)
      }
    )
    
    
    # ========================================================
    # PROTECT AGAINST INFINITE RETRY LOOP
    # ========================================================
    
    if (!design_complete) {
      
      if (retry_count >= max_retries) {
        
        stop(
          paste0(
            "Design ",
            design_id,
            " failed ",
            max_retries,
            " consecutive times. ",
            "Check stan_workers.log."
          )
        )
      }
      
      
      message(
        "Waiting 15 seconds before restarting..."
      )
      
      Sys.sleep(15)
    }
  }
  
  
  cat(
    "\nDesign",
    design_id,
    "completed successfully.\n"
  )
  
  
  gc()
}


# ============================================================
# COMBINE ALL COMPLETED CHECKPOINTS
# ============================================================

checkpoint_files <- list.files(
  checkpoint_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)


saved_results <- lapply(
  checkpoint_files,
  readRDS
)


# Complete cohort-level trajectories
design_results_df <- dplyr::bind_rows(
  lapply(
    saved_results,
    function(x) x$path
  )
)


# Final posterior ID70 summaries
id70_df <- dplyr::bind_rows(
  lapply(
    saved_results,
    function(x) x$id70
  )
)

selected_dose_info_df <- dplyr::bind_rows(
  lapply(
    saved_results,
    function(x) x$all_selected_dose_info
  )
)


# Sort
design_results_df <- design_results_df %>%
  arrange(
    design_id,
    sim_id,
    cohort
  )


id70_df <- id70_df %>%
  arrange(
    design_id,
    sim_id
  )

selected_dose_info_df <- selected_dose_info_df %>%
  arrange(
    design_id,
    sim_id
  )

write.csv(
  design_results_df,
  "simulation result/expected_info_2/adaptive_design_results.csv",
  row.names = FALSE
)


write.csv(
  id70_df,
  "simulation result/expected_info_2/adaptive_post_id70.csv",
  row.names = FALSE
)

write.csv(
  selected_dose_info_df,
  "simulation result/expected_info_2/adaptive_selected_dose_info.csv",
  row.names = FALSE
)

