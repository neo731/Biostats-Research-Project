library(dplyr)
library(ggplot2)
library(truncnorm)

# dose-response model
p_inf <- function(D, a, b) {
  1 - (1 + D / b)^(-a)
}

# smooth dose grid
dose_grid <- 10^seq(1, 9, length.out = 200)

# prior predictive draws
n_draws <- 20000
set.seed(42)

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

# prior predictive probability matrix
p_prior <- sapply(dose_grid, function(d) {
  p_inf(d, prior_draws$a, prior_draws$b)
})

prior_df <- data.frame(
  dose = dose_grid,
  mean = colMeans(p_prior),
  lower95 = apply(p_prior, 2, quantile, 0.025),
  upper95 = apply(p_prior, 2, quantile, 0.975),
  lower50 = apply(p_prior, 2, quantile, 0.25),
  upper50 = apply(p_prior, 2, quantile, 0.75)
)

# scenario curves
scenario_curves <- bind_rows(
  lapply(seq_len(nrow(scenario_params)), function(i) {
    data.frame(
      dose = dose_grid,
      p = p_inf(dose_grid, scenario_params$a[i], scenario_params$b[i]),
      scenario_label = as.character(i)
    )
  })
)

label_positions <- c(3.15, 3.7, 4.2, 4.9, 5.65)  # on log10 scale

label_df <- scenario_curves %>%
  mutate(log_dose = log10(dose)) %>%
  group_by(scenario_label) %>%
  slice(which.min(abs(log_dose - label_positions[as.numeric(first(scenario_label))]))) %>%
  ungroup()

# plot
p <- ggplot() +
  geom_ribbon(
    data = prior_df,
    aes(x = dose, ymin = lower95, ymax = upper95),
    fill = "steelblue",
    alpha = 0.20
  ) +
  geom_ribbon(
    data = prior_df,
    aes(x = dose, ymin = lower50, ymax = upper50),
    fill = "blue",
    alpha = 0.25
  ) +
  geom_line(
    data = scenario_curves,
    aes(x = dose, y = p, group = scenario_label),
    linetype = "dotted",
    linewidth = 1
  ) +
  geom_text(
    data = label_df,
    aes(x = dose, y = 0.75, label = scenario_label),
    #hjust = -0.2,
    size = 4
  ) +
  geom_hline(
    yintercept = 0.70, color = "red", linewidth = 0.4
  ) +
  geom_hline(
    yintercept = c(0.60, 0.80),
    color = "red", 
    linetype = "dashed", 
    linewidth = 0.2
  ) +
  scale_x_log10(
    breaks = 10^(1:9),               # Sets ticks at each power of 10
    labels = function(x) log10(x)    # Converts 10, 100, 1000 -> 1, 2, 3
  ) +
  theme_bw() +
  labs(
    title = paste("Prior Predictive Dose-response Relationship"),
    x = expression("Dose (" * log[10] ~ "CFU)"),
    y = "Probability of infection"
  ) 

p
ggsave("plots/prior_predictive_scenarios.png", 
       plot = p,
       width = 7, height = 5)

# now put 20 random samples of priors on the plot to showcase this as an example of how the 1000 simulations work
n_prior_samples <- 20

prior_sample_draws <- prior_draws %>%
  slice_sample(n = n_prior_samples) %>%
  mutate(draw_id = row_number())

prior_sample_curves <- bind_rows(
  lapply(seq_len(nrow(prior_sample_draws)), function(i) {
    data.frame(
      draw_id = prior_sample_draws$draw_id[i],
      dose = dose_grid,
      p = p_inf(
        dose_grid,
        prior_sample_draws$a[i],
        prior_sample_draws$b[i]
      )
    )
  })
)

p2 <- ggplot() +
  
  geom_ribbon(
    data = prior_df,
    aes(x = dose, ymin = lower95, ymax = upper95),
    fill = "steelblue",
    alpha = 0.20
  ) +
  geom_ribbon(
    data = prior_df,
    aes(x = dose, ymin = lower50, ymax = upper50),
    fill = "blue",
    alpha = 0.25
  ) +
  # 20 random prior predictive curves
  geom_line(
    data = prior_sample_curves,
    aes(x = dose, y = p, group = draw_id),
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 0.70, color = "red", linewidth = 0.4
  ) +
  geom_hline(
    yintercept = c(0.60, 0.80),
    color = "red", 
    linetype = "dashed", 
    linewidth = 0.2
  ) +
  scale_x_log10(
    breaks = 10^(1:9),
    labels = function(x) log10(x)
  ) +
  theme_bw() +
  labs(
    title = "Prior Predictive Dose-response Relationship",
    x = expression("Dose (" * log[10] ~ "CFU)"),
    y = "Probability of infection"
  )

p2
ggsave("plots/prior_predictive_sample_priors.png", 
       plot = p2,
       width = 7, height = 5)