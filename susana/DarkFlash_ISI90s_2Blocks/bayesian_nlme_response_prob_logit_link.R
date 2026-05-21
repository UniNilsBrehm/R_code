# ==============================================================================
# Bayesian nonlinear mixed model for zebrafish dark-flash habituation
# ==============================================================================
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(brms)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(loo)
library(DHARMa)
library(bayesplot)

# ==============================================================================
# HELPER
# Helper to save consistently
save_plot <- function(p, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(diag_dir, filename),
    plot     = p,
    width    = width,
    height   = height,
    dpi      = 150,
    bg       = "white"
  )
}
# ==============================================================================

# source("C:/UniFreiburg/Code/R_code/susana/nlme_utils.R")
# source("C:/UniFreiburg/Code/R_code/susana/nlme_plot_utils.R")

source("D:/Behavior_Data/R_code/susana/nlme_utils.R")
source("D:/Behavior_Data/R_code/susana/plot_utils.R")


# base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <- "D:/WorkingData/Susana/DarkFlash_ISI90s_2Blocks"
base_dir <- "D:/Behavior_Data/DarkFlash_ISI90s_2Blocks"

file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

var_name = 'response_prob'
col_name = 'move'

save_fig_dir = file.path(base_dir, "figs", "nlme", var_name)
save_results_dir = file.path(base_dir, "results", "nlme", var_name)

# Create directories if they do not exist
dir.create(save_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Load and prepare data
# ==============================================================================
res <- load_data_darkflash_60s(file_dir, move_th = 0, take_peak = 0)

df_final     <- res$df_final

df_resp <- df_final %>%
  mutate(
    stimulus = as.numeric(stimulus),
    stimulus0 = stimulus - 1,
    Block = factor(Block),
    Genotype = factor(Genotype),
    Video = factor(Video),
    Well = factor(Well),
    animal = interaction(Video, Well, drop = TRUE)
    
  )


# ==============================================================================
# The Model
# ==============================================================================
model <- bf(
  move ~ A + (R0 - A) * exp(-exp(logk) * stimulus0),
  
  A     ~ 0 + Genotype * Block + (1 | animal), # initial response level on logit scale
  R0    ~ 0 + Genotype * Block + (1 | animal), # final/asymptotic response level on logit scale
  logk ~ 0 + Genotype * Block + (1 | animal),  # log decay rate: k = exp(logk), half_life = log(2) / exp(logk)
  
  nl = TRUE
)

# ==============================================================================
# The Priors
# ==============================================================================
priors <- c(
  prior(normal(0, 1.0), class = "b", nlpar = "A"),
  prior(normal(1.5, 1), class = "b", nlpar = "R0"),
  prior(normal(-3, 1), class = "b", nlpar = "logk"),
  
  prior(exponential(2), class = "sd", nlpar = "A"),
  prior(exponential(2), class = "sd", nlpar = "R0"),
  prior(exponential(2), class = "sd", nlpar = "logk")
)

# ==============================================================================
# Validate Priors
# ==============================================================================
fit_prior <- brm(
  formula = model,
  data = df_resp,
  family = bernoulli(link = "logit"),
  prior = priors,
  sample_prior = "only",          # <-- key change
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  iter = 2000,                    # fewer iters are fine for prior checks
  warmup = 500,
  seed = 42
)


diag_dir <- file.path(base_dir, "models",  "diagnostics", "nlme", "priors", var_name)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

p1 <- pp_check(fit_prior, ndraws = 100) +
  labs(title = "Prior predictive: response density")
save_plot(p1, "01_prior_ppc_density.png")

p2 <- pp_check(fit_prior, type = "stat_grouped", group = "stimulus0",
               stat = "mean", ndraws = 200) +
  labs(title = "Prior predictive: mean response per stimulus")
save_plot(p2, "02_prior_ppc_mean_per_stimulus.png", width = 10, height = 7)

new_data <- expand_grid(
  stimulus0 = 0:max(df_resp$stimulus0),
  Genotype  = levels(df_resp$Genotype),
  Block     = levels(df_resp$Block)
) %>%
  mutate(animal = NA)

draws <- posterior_epred(
  fit_prior,
  newdata    = new_data,
  re_formula = NA,
  ndraws     = 200
)

plot_df <- as.data.frame(t(draws)) %>%
  bind_cols(new_data) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to  = "draw",
    values_to = "p_response"
  )

p3 <- ggplot(plot_df, aes(x = stimulus0, y = p_response, group = draw)) +
  geom_line(alpha = 0.1, color = "steelblue") +
  facet_grid(Block ~ Genotype) +
  labs(
    title = "Prior predictive habituation curves",
    x = "Stimulus number (0-indexed)",
    y = "P(response)"
  ) +
  ylim(0, 1) +
  theme_minimal()
save_plot(p3, "03_prior_curves_by_group.png", width = 10, height = 7)

prior_samples <- as_draws_df(fit_prior)

p4 <- prior_samples %>%
  mutate(k = exp(b_lk_Intercept)) %>%
  ggplot(aes(x = k)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  scale_x_log10() +
  labs(title = "Prior on decay rate k (log scale)", x = "k")
save_plot(p4, "04_prior_k_density.png")

p5 <- prior_samples %>%
  transmute(
    p0_prob   = plogis(b_p0_Intercept),
    pinf_prob = plogis(b_pinf_Intercept)
  ) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value, fill = name)) +
  geom_density(alpha = 0.5) +
  labs(title = "Prior on initial and asymptotic response probability",
       x = "Probability")
save_plot(p5, "05_prior_p0_pinf_density.png")

param_summary <- prior_samples %>%
  transmute(
    p0_prob    = plogis(b_p0_Intercept),
    pinf_prob  = plogis(b_pinf_Intercept),
    k          = exp(b_lk_Intercept),
    half_life  = log(2) / exp(b_lk_Intercept)  # stimuli to halve gap
  ) %>%
  pivot_longer(everything(), names_to = "param") %>%
  group_by(param) %>%
  summarise(
    mean    = mean(value),
    median  = median(value),
    q025    = quantile(value, 0.025),
    q25     = quantile(value, 0.25),
    q75     = quantile(value, 0.75),
    q975    = quantile(value, 0.975),
    .groups = "drop"
  )

print(param_summary)
write.csv(param_summary,
          file.path(diag_dir, "prior_param_summary.csv"),
          row.names = FALSE)

curve_summary <- plot_df %>%
  group_by(Genotype, Block, stimulus0) %>%
  summarise(
    mean   = mean(p_response),
    median = median(p_response),
    q025   = quantile(p_response, 0.025),
    q25    = quantile(p_response, 0.25),
    q75    = quantile(p_response, 0.75),
    q975   = quantile(p_response, 0.975),
    .groups = "drop"
  )

# How wide is the 95% prior predictive interval at each stimulus?
interval_widths <- curve_summary %>%
  mutate(width_95 = q975 - q025) %>%
  group_by(stimulus0) %>%
  summarise(mean_width_95 = mean(width_95))

print(interval_widths)
write.csv(curve_summary,
          file.path(diag_dir, "prior_curve_summary.csv"),
          row.names = FALSE)

# For each draw, compute key features per group
per_draw_features <- plot_df %>%
  group_by(draw, Genotype, Block) %>%
  summarise(
    p_first      = p_response[stimulus0 == 0],
    p_last       = p_response[stimulus0 == max(stimulus0)],
    decrease     = p_first - p_last,
    decreasing   = p_first > p_last,
    .groups = "drop"
  )

feature_summary <- per_draw_features %>%
  summarise(
    pct_decreasing       = mean(decreasing),
    pct_strong_habituate = mean(decrease > 0.2),
    pct_no_change        = mean(abs(decrease) < 0.05),
    pct_increasing       = mean(decrease < -0.05),
    median_decrease      = median(decrease)
  )

print(feature_summary)
write.csv(feature_summary,
          file.path(diag_dir, "prior_feature_summary.csv"),
          row.names = FALSE)

# ==============================================================================
# Fit Fast Test Model
# ==============================================================================
# Create a tiny xx% subset of your data for rapid prototyping
df_test_sub <- df_resp %>% 
  group_by(Genotype, Block) %>% 
  slice_sample(prop = 0.90) %>% 
  ungroup()

# Fit using Meanfield Variational Inference (VI)
fit_vi_test <- brm(
  formula = model,
  data = df_resp,               
  family = bernoulli(link = "logit"),
  prior = priors,
  backend = "cmdstanr",
  algorithm = "meanfield",          # <--- Reliable, fast VI method
  # algorithm = "fullrank",
  iter = 10000                      # VI likes higher iterations (it's still lightning fast)
)

fit_model <- fit_vi_test

# ==============================================================================
# Fit the model
# ==============================================================================
fit_model <- brm(
  formula = model,
  data = df_resp,
  family = bernoulli(link = "logit"),
  prior = priors,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  threads = threading(6),
  iter = 3000,
  warmup = 1000,
  seed = 42,
  control = list(adapt_delta = 0.90, max_treedepth = 10)
)

# identity because our formula yields the logit directly, but brms needs to 
# know it's a probability mapping if we don't wrap it, OR use link="logit" if 
# we wrap the right side. 
# Alternative cleaner approach for brms logit link:
# family = bernoulli(link = "logit"), but then the formula predicts the 
# logit directly:
# move ~ A + B * exp(-exp(c) * stimulus0)

# ==============================================================================
# Save fitted model to HDD
# ==============================================================================
saveRDS(
  fit_model,
  file = file.path(base_dir, "models", paste0("bayesian_nlme_", var_name,"_results.rds"))
)

# ==============================================================================
# Load fitted model if available
# ==============================================================================
fit_model <- readRDS(
  file.path(base_dir, "models", paste0("bayesian_nlme_", var_name,"_results.rds"))
)

# ==============================================================================
# Get summary
# ==============================================================================
# --- Summary table: R-hat, ESS ---
sink(file.path(save_results_dir, paste0(var_name, '_model_summary.txt')))
fit_summary <- summary(fit_model)
print(fit_summary)
sink()

posterior_summary(fit_model, pars = "^b_A")
posterior_summary(fit_model, pars = "^b_R0")
posterior_summary(fit_model, pars = "^b_logk")
# ==============================================================================
# Plot Habituation curves
# ==============================================================================
save_fig_dir <- NULL
p_hab <- plot_habituation_probability(
  df_resp = df_resp,
  fit_model = fit_model,
  save_fig_dir = save_fig_dir,
  var_name = var_name
)

p_hab

p_latent <- plot_habituation_latent(
  df_resp = df_resp,
  fit_model = fit_model,
  save_fig_dir = save_fig_dir,
  var_name = var_name
)

p_latent

# ==============================================================================
# Compare nonlinear model parameters
# ==============================================================================
# ==============================================================================
# 1. Compare habituation rate k = exp(lk)
# ==============================================================================

k_draws_df <- make_nlpar_draws(
  fit_model = fit_model,
  df_resp = df_resp,
  nlpar = "lk",
  transform_fun = exp
) %>%
  rename(k = value, lk = value_raw)

k_summary <- k_draws_df %>%
  group_by(Genotype, Block) %>%
  summarise(
    k_median = median(k),
    k_low = quantile(k, 0.025),
    k_high = quantile(k, 0.975),
    half_life_median = median(log(2) / k),
    half_life_low = quantile(log(2) / k, 0.025),
    half_life_high = quantile(log(2) / k, 0.975),
    .groups = "drop"
  )

print(k_summary)

write.csv(
  k_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_k_summary.csv")),
  row.names = FALSE
)

k_comparison <- compare_nlpar(
  k_draws_df %>% rename(value = k),
  value_name = "k",
  ratio = TRUE
)

k_comparison_summary <- k_comparison$summary
print(k_comparison_summary)

write.csv(
  k_comparison_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_k_all_pairwise_comparisons.csv")),
  row.names = FALSE
)


# ==============================================================================
# 2. Compare asymptote pinf on probability scale
# ==============================================================================

pinf_draws_df <- make_nlpar_draws(
  fit_model = fit_model,
  df_resp = df_resp,
  nlpar = "pinf",
  transform_fun = plogis
) %>%
  rename(pinf = value, pinf_logit = value_raw)

pinf_summary <- summarise_nlpar(
  pinf_draws_df %>% rename(value = pinf),
  "pinf"
)

print(pinf_summary)

write.csv(
  pinf_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_pinf_summary.csv")),
  row.names = FALSE
)

pinf_comparison <- compare_nlpar(
  pinf_draws_df %>% rename(value = pinf),
  value_name = "pinf",
  rope = 0.02,
  ratio = FALSE
)

pinf_comparison_summary <- pinf_comparison$summary
print(pinf_comparison_summary)

write.csv(
  pinf_comparison_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_pinf_all_pairwise_comparisons.csv")),
  row.names = FALSE
)


# ==============================================================================
# 3. Compare initial response probability p0
# ==============================================================================

p0_draws_df <- make_nlpar_draws(
  fit_model = fit_model,
  df_resp = df_resp,
  nlpar = "p0",
  transform_fun = plogis
) %>%
  rename(p0 = value, p0_logit = value_raw)

p0_summary <- summarise_nlpar(
  p0_draws_df %>% rename(value = p0),
  "p0"
)

print(p0_summary)

write.csv(
  p0_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_p0_summary.csv")),
  row.names = FALSE
)

p0_comparison <- compare_nlpar(
  p0_draws_df %>% rename(value = p0),
  value_name = "p0",
  rope = 0.02,
  ratio = FALSE
)

p0_comparison_summary <- p0_comparison$summary
print(p0_comparison_summary)

write.csv(
  p0_comparison_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_p0_all_pairwise_comparisons.csv")),
  row.names = FALSE
)


# ==============================================================================
# 4. Plots
# ==============================================================================

# k posterior density

p_k <- plot_posterior_density(
  k_draws_df,
  k,
  "Posterior distributions of habituation rate k",
  "Habituation rate k",
  paste0("nlme_", var_name, "_k_posteriors.png"),
  save_fig_dir = save_fig_dir,
  block_limits = list(
    Block1 = c(0, 2),
    Block2 = c(0, 0.10)
  )
)


# k pairwise ratios

p_k_compare <- k_comparison_summary %>%
  ggplot(
    aes(
      x = comparison,
      y = median_ratio,
      ymin = ratio_low,
      ymax = ratio_high,
      color = Genotype_1
    )
  ) +
  facet_wrap(~Block, scales = "free_x") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(linewidth = 0.8) +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x = "Comparison",
    y = "Habituation rate ratio",
    title = "All pairwise habituation rate comparisons",
    color = "Numerator genotype"
  )

print(p_k_compare)

ggsave(
  file.path(save_fig_dir, paste0("nlme_", var_name, "_k_all_pairwise_ratios.png")),
  p_k_compare,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)


# pinf posterior density

p_pinf <- plot_posterior_density(
  pinf_draws_df,
  pinf,
  "Posterior distributions of asymptote pinf",
  "Asymptote response probability pinf",
  paste0("nlme_", var_name, "_pinf_posteriors.png"),
  save_fig_dir = save_fig_dir
)


# pinf pairwise differences

p_pinf_compare <- plot_pairwise_differences(
  pinf_comparison_summary,
  "pinf difference",
  "All pairwise pinf comparisons",
  paste0("nlme_", var_name, "_pinf_all_pairwise_differences.png"),
  save_fig_dir = save_fig_dir
)


# p0 posterior density

p_p0 <- plot_posterior_density(
  p0_draws_df,
  p0,
  "Posterior distributions of initial response probability p0",
  "Initial response probability p0",
  paste0("nlme_", var_name, "_p0_posteriors.png"),
  save_fig_dir = save_fig_dir
)


# p0 pairwise differences

p_p0_compare <- plot_pairwise_differences(
  p0_comparison_summary,
  "p0 difference",
  "All pairwise p0 comparisons",
  paste0("nlme_", var_name, "_p0_all_pairwise_differences.png"),
  save_fig_dir = save_fig_dir
)