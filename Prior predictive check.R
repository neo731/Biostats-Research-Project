library(ggplot2)
library(dplyr)
library(tidyr)
library(truncnorm)

setwd('')

# prior predictive check
# beta-poisson function
beta_poisson <- function(dose, a, b) {
  1 - (1 + dose/b)^-a
}

calc_ID <- function(a, b) {
  b * ((0.3^-1/a) - 1)
}

# we will do 1000 simulations to check the shape of dose-response model with our choice of prior
set.seed(42)
n_sims <- 1000
doses <- seq(10^3, 10^9, by=10000)
dose_grid <- 10^seq(1, 9, length.out = 101)

# sampling prior distribution
# when b is small relative to the dose, the curve becomes flat at 1
# I discovered that the value of b cannot be too far off the dose level
# Thus, I attempted to use log uniform distribution for b, such that log10(b) ~ Unif(4.5,5.5)
# The value of b is therefore 10^runif(4.5, 5.5)
prior_data <- data.frame(
  sim = 1:n_sims,
  a = runif(n_sims, 0.2, 0.7),
  b = 10^runif(n_sims, 4.5, 5.5)
)

# creating a data frame for the plot later, which includes the dose level, and their respective 
# probabilities from the dose-response function
plot_data <- prior_data %>%
  group_by(sim) %>%
  do(data.frame(dose=doses, prob = beta_poisson(doses, .$a, .$b)))

plot1 <- ggplot() +
  geom_line(data = plot_data, aes(x=dose, y=prob, group=sim), alpha=0.15, colour = 'steelblue') +
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "darkgreen") +
  scale_x_log10(limits = c(10^3, 10^9), breaks = 10^(3:9), labels = scales::label_log())

plot1
hist(prior_data$b)


# Attempting with different prior distribution for b, as this parameter seemed to be more difficult
# to adjust. Distribution of a remains unchanged, but now log(b) ~ normal(5, 0.05)

prior_data2 <- data.frame(
  sim = 1:n_sims,
  a = runif(n_sims, 0.2, 0.7),
  b = 10^rnorm(n_sims, 5, 0.05)
)

# generate probabilities for each dose level
plot_data2 <- prior_data2 %>%
  group_by(sim) %>%
  do(data.frame(dose=doses, prob = beta_poisson(doses, .$a, .$b)))

# Calculate the 95% Credible Interval (Ribbon Data)
# group by dose to find the range of 'prob' across all 1000 sims
cri_data <- plot_data2 %>%
  group_by(dose) %>%
  summarise(
    lower = quantile(prob, 0.025),
    median = quantile(prob, 0.50),
    upper = quantile(prob, 0.975),
    Q1 = quantile(prob, 0.25),
    Q3 = quantile(prob, 0.75)
  )

plot2 <- ggplot() +
  # The 95% Credible Interval Ribbon
  geom_ribbon(data = cri_data, aes(x = dose, ymin = lower, ymax = upper), 
              fill = "steelblue", alpha = 0.2) +
  # The 50% Credible Interval Ribbon
  geom_ribbon(data = cri_data, aes(x = dose, ymin = Q1, ymax = Q3), 
              fill = "blue", alpha = 0.2) +
  # geom_line(data = plot_data2, aes(x=dose, y=prob, group=sim), alpha=0.15, colour = 'steelblue') +
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "darkgreen") +
  scale_x_log10(limits = c(10^3, 10^9), breaks = 10^(3:9), labels = scales::label_log())

plot2
ggsave("prior_check1.png", plot = plot2)
hist(prior_data2$b)


## using the priors from the tutorial paper
log_a <- runif(n_sims, log(0.3), log(0.8))
log10_id70 <- rtruncnorm(
  n_sims, a = 3, b = 7, mean = 4, sd = 0.8
)
a = exp(log_a)
id70 <- 10^log10_id70

prior_data_3 <- data.frame(
  sim = 1:n_sims,
  a = a,
  b = id70 / (0.3^(-1 / a) - 1)
)

# generate probabilities for each dose level
plot_data_3 <- prior_data_3 %>%
  group_by(sim) %>%
  do(data.frame(dose=dose_grid, prob = beta_poisson(dose_grid, .$a, .$b)))

# Calculate the 95% Credible Interval (Ribbon Data)
# group by dose to find the range of 'prob' across all 1000 sims
cri_data <- plot_data_3 %>%
  group_by(dose) %>%
  summarise(
    lower = quantile(prob, 0.025),
    median = quantile(prob, 0.50),
    upper = quantile(prob, 0.975),
    Q1 = quantile(prob, 0.25),
    Q3 = quantile(prob, 0.75)
  )

plot3 <- ggplot() +
  # The 95% Credible Interval Ribbon
  geom_ribbon(data = cri_data, aes(x = dose, ymin = lower, ymax = upper), 
              fill = "steelblue", alpha = 0.2) +
  # The 50% Credible Interval Ribbon
  geom_ribbon(data = cri_data, aes(x = dose, ymin = Q1, ymax = Q3), 
              fill = "blue", alpha = 0.2) +
  geom_hline(yintercept = 0.70, color = "red", linewidth = 0.4) +
  geom_hline(yintercept = c(0.60, 0.80),
             color = "red", linetype = "dashed", linewidth = 0.2) +
  scale_x_log10(
    breaks = 10^(1:9),               # Sets ticks at each power of 10
    labels = function(x) log10(x)    # Converts 10, 100, 1000 -> 1, 2, 3
  ) +
  theme_bw() +
  labs(
    title = paste("Prior Predictive Distribution"),
    x = expression("Dose (" * log[10] ~ "CFU)"),
    y = "Probability of infection"
  )

plot3
ggsave("prior_check.png", plot = plot3)


# getting the curve to demonstrate the dose-response shape
# getting example a and b values for the demonstration curve

p_inf <- function(x, a, b) {
  1 - (1 + x/b)^-a
}

set.seed(42)
a <- exp(runif(1, log(0.3), log(0.8)))
id70 <- 10^(rtruncnorm(1, a = 3, b = 7, mean = 4, sd = 0.8))
b = id70 / (0.3^(-1 / a) - 1)

example_dose_response_curve <- ggplot() +
  geom_function(fun = p_inf, args = list(a = a, b = b), linewidth = 1, colour = "black") +
  # Vertical segment: (id70, 0) to (id70, 0.7)
  annotate("segment", 
           x = id70, xend = id70, 
           y = 0, yend = 0.7, 
           color = "red", linetype = "dotted", linewidth = 0.8) +
  # Horizontal segment: (10^1, 0.7) to (id70, 0.7)
  annotate("segment", 
    x = 10^1, xend = id70, 
    y = 0.7, yend = 0.7, 
    color = "red", linetype = "dotted", linewidth = 0.8
  ) +
  # Set log10 scale for x-axis with ticks at 10^1 to 10^9 labeled as 1 to 9
  scale_x_log10(
    limits = c(10^1, 10^9),
    breaks = 10^(1:9),
    labels = 1:9,
    expand = c(0, 0)
  ) +
  # Constrain y-axis limits
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "Example Beta-Poisson Dose-Response Curve",
    x = expression("Dose (" * log[10] ~ "CFU)"),
    y = "Probability of Infection"
  ) +
  theme_classic() 

# Save high-resolution PNG output matching original dimensions/resolution
ggsave(
  "example_dose_response_curve.png", 
  plot = example_dose_response_curve, 
  width = 6, 
  height = 4, 
  dpi = 300
)


# getting 5 curves to demonstrate the dose-response shape under the same ID70
# getting example a and b values for the demonstration curve
set.seed(42)

# One common ID70 for all curves
id70 <- 10^(rtruncnorm(
  1, a = 3, b = 7, mean = 4, sd = 0.8
))

# Five different shape parameters
a_values <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

# Calculate b to ensure all curves share the same ID70
curve_params <- data.frame(
  a = a_values
) %>%
  mutate(
    b = id70 / (0.3^(-1 / a) - 1),
    shape = factor(
      paste0("a = ", a),
      levels = paste0("a = ", a_values)
    )
  )

# Generate curves over the log-dose range
curve_data <- expand_grid(
  dose = 10^seq(1, 9, length.out = 500),
  curve_params
) %>%
  mutate(
    p = p_inf(dose, a, b)
  )

example_dose_response_curve <- ggplot(curve_data, aes(x = dose, y = p, colour = shape)) +
  geom_line(linewidth = 1.1) +
  
  # Shared ID70 vertical reference
  annotate(
    "segment",
    x = id70, xend = id70,
    y = 0, yend = 0.7,
    colour = "red",
    linetype = "dotted",
    linewidth = 0.9
  ) +
  
  # Horizontal reference at 70% infection
  annotate(
    "segment",
    x = 10^1, xend = id70,
    y = 0.7, yend = 0.7,
    colour = "red",
    linetype = "dotted",
    linewidth = 0.9
  ) +
  
  # Highlight common ID70
  annotate(
    "point",
    x = id70, y = 0.7,
    colour = "red",
    size = 2.8
  ) +
  
  scale_x_log10(
    limits = c(10^1, 10^9),
    breaks = 10^(1:9),
    labels = 1:9,
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  
  scale_colour_viridis_d(
    option = "D",
    end = 0.9,
    name = "Shape parameter"
  ) +
  
  labs(
    title = "Beta-Poisson Dose-Response Curves",
    subtitle = "Different curve shapes with a common ID70",
    x = expression("Dose (" * log[10] ~ "CFU)"),
    y = "Probability of Infection"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right"
  )

# Save high-resolution PNG output matching original dimensions/resolution
ggsave(
  "example_dose_response_curve2.png", 
  plot = example_dose_response_curve, 
  width = 6, 
  height = 4, 
  dpi = 300
)


# Beta-Poisson dose-response function
p_inf <- function(dose, a, b) {
  1 - (1 + dose / b)^(-a)
}

# Parameter values
a_values <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
b_values <- 10^seq(2, 5, length.out = 6)

# Fixed reference parameters
a_fixed <- 0.5
b_fixed <- 1000

# Dose range
dose_values <- 10^seq(1, 9, length.out = 500)

# Panel A: Vary a, keep b fixed
panel_a <- expand_grid(
  dose = dose_values,
  a = a_values
) %>%
  mutate(
    b = b_fixed,
    p = p_inf(dose, a, b),
    parameter = factor(
      paste0("a = ", a),
      levels = paste0("a = ", a_values)
    ),
    panel = "A. Varying a (fixed b = 1000)"
  )

# Panel B: Vary b, keep a fixed
panel_b <- expand_grid(
  dose = dose_values,
  b = b_values
) %>%
  mutate(
    a = a_fixed,
    p = p_inf(dose, a, b),
    parameter = factor(
      paste0("b = ", signif(b, 3)),
      levels = paste0("b = ", signif(b_values, 3))
    ),
    panel = "B. Varying b (fixed a = 0.5)"
  )

# Combine data
curve_data <- bind_rows(panel_a, panel_b) %>%
  mutate(
    panel = factor(
      panel,
      levels = c(
        "A. Varying a (fixed b = 1000)",
        "B. Varying b (fixed a = 0.5)"
      )
    )
  )

# Two-panel figure
example_dose_response_curve <- ggplot(
  curve_data,
  aes(x = dose, y = p, colour = parameter)
) +
  geom_line(linewidth = 1.1) +
  
  geom_hline(
    yintercept = 0.7,
    linetype = "dotted",
    colour = "red",
    linewidth = 0.7
  ) +
  
  facet_wrap(
    ~ panel,
    ncol = 2
  ) +
  
  scale_x_log10(
    limits = c(10^1, 10^9),
    breaks = 10^(1:9),
    labels = 1:9,
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  
  scale_colour_viridis_d(
    option = "D",
    end = 0.9,
    name = "Parameter value"
  ) +
  
  labs(
    title = "Effects of Beta-Poisson Parameters on Dose-Response Curves",
    x = expression("Dose (" * log[10] ~ "CFU)"),
    y = "Probability of Infection"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5)
  )

example_dose_response_curve
  
