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

source("C:/UniFreiburg/Code/R_code/susana/nlme_utils.R")
source("C:/UniFreiburg/Code/R_code/susana/nlme_plot_utils.R")

# source("D:/Behavior_Data/R_code/susana/nlme_utils.R")
# source("D:/Behavior_Data/R_code/susana/nlme_plot_utils.R")


# base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_ISI90s_2Blocks"
base_dir <- "D:/WorkingData/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <- "D:/Behavior_Data/DarkFlash_ISI90s_2Blocks"

file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

var_name = 'response_prob'
col_name = 'move'

save_fig_dir = file.path(base_dir, "figs", "nlme", var_name)
save_results_dir = file.path(base_dir, "results", "nlme", var_name)

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
  # logit(p_t) = pinf + (p0 - pinf) * exp(-k * stimulus0)
  move ~ pinf + (p0 - pinf) * exp(-exp(lk) * stimulus0),
  
  pinf ~ Genotype * Block + (1 | animal),
  p0   ~ Genotype * Block + (1 | animal),
  lk   ~ Genotype * Block + (1 | animal),
  
  nl = TRUE
)


# ==============================================================================
# The Priors
# ==============================================================================
# Comments on Priors:
# normal(0, 5) on logit-scale coefficients: 95% prior mass roughly between 
# logit values of -10 and +10, which covers probabilities 
# from ~0.00005 to ~0.99995. Effectively flat over the plausible range but rules 
# out absurd values that would destabilize sampling.

# normal(0, 2) on lk: allows k to range from 
# exp(-4) ≈ 0.018 to exp(4) ≈ 54.6

# student_t(3, 0, 2.5) on SDs: the brms/Stan default for variance components 
# — heavy tails allow large group-level variation if the data support it, 
# but pulls toward zero in the absence of evidence (partial pooling).

priors <- c(
  # Group coefficients (Genotype, Block, and their interactions)
  prior(normal(0, 5), nlpar = "p0",   class = "b"),
  prior(normal(0, 5), nlpar = "pinf", class = "b"),
  prior(normal(0, 2), nlpar = "lk",   class = "b"),
  
  # Random effect SDs (half-t due to positivity constraint)
  prior(student_t(3, 0, 2.5), nlpar = "p0",   class = "sd"),
  prior(student_t(3, 0, 2.5), nlpar = "pinf", class = "sd"),
  prior(student_t(3, 0, 2.5), nlpar = "lk",   class = "sd")
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
# directly on prob scale
fit_model <- brm(
  formula = model,
  data = df_resp,
  family = bernoulli(link = "logit"),
  prior = priors,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  threads = threading(6),
  iter = 5000,
  warmup = 1500,
  seed = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

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

# ==============================================================================
# Get diagnostics
# ==============================================================================
# ==== Convergence Diagnostics ====

# --- Extract diagnostics as a tidy table ---
draws_array <- as_draws_array(fit_model)
diag_df <- summarise_draws(draws_array, default_convergence_measures())
print(diag_df, n = Inf)
write.csv(diag_df,
          file.path(save_results_dir, "convergence_diagnostics.csv"),
          row.names = FALSE)

# --- Flag any problems ---
bad_rhat <- diag_df %>% filter(rhat > 1.01)
bad_ess  <- diag_df %>% filter(ess_bulk < 400 | ess_tail < 400)

cat("\nParameters with R-hat > 1.01:", nrow(bad_rhat), "\n")
cat("Parameters with ESS < 400:    ", nrow(bad_ess),  "\n")

# --- HMC diagnostics: divergences, treedepth, energy ---
np <- nuts_params(fit_model)
cat("\nDivergent transitions:", sum(subset(np, Parameter == "divergent__")$Value), "\n")
cat("Max treedepth hits:    ", sum(subset(np, Parameter == "treedepth__")$Value >= 12), "\n")

# ==== Trace Plots ====
# --- Trace plots for key parameters ---
trace_pars <- variables(fit_model)[grepl("^b_|^sd_", variables(fit_model))]

p_trace <- mcmc_trace(fit_model, pars = trace_pars[1:min(12, length(trace_pars))]) +
  ggtitle("Trace plots (first 12 fixed/SD parameters)")
save_plot(p_trace, "10_trace_plots.png", width = 12, height = 8)

# --- Rank plots: alternative to trace, often more sensitive ---
p_rank <- mcmc_rank_overlay(fit_model, pars = trace_pars[1:min(6, length(trace_pars))]) +
  ggtitle("Rank plots")
save_plot(p_rank, "11_rank_plots.png", width = 12, height = 8)

# --- Pairs plot for the intercepts (catches funnel-y geometry) ---
p_pairs <- mcmc_pairs(
  fit_model,
  pars = c("b_p0_Intercept", "b_pinf_Intercept", "b_lk_Intercept"),
  off_diag_args = list(size = 0.5, alpha = 0.3)
)
save_plot(p_pairs, "12_pairs_intercepts.png", width = 9, height = 9)

# ==== Posterior Predictive Checks ====
# --- Overall fit: predicted vs observed response density ---
p_ppc_dens <- pp_check(fit_model, ndraws = 100) +
  ggtitle("Posterior predictive: response density")
save_plot(p_ppc_dens, "20_ppc_density.png")

# --- Mean response per stimulus: KEY check for habituation models ---
p_ppc_stim <- pp_check(fit_model, type = "stat_grouped",
                       group = "stimulus0", stat = "mean", ndraws = 200) +
  ggtitle("Posterior predictive: mean response per stimulus")
save_plot(p_ppc_stim, "21_ppc_mean_per_stimulus.png", width = 12, height = 8)

# --- Mean response per genotype × block ---
df_resp$geno_block <- interaction(df_resp$Genotype, df_resp$Block, drop = TRUE)
p_ppc_grp <- pp_check(fit_model, type = "stat_grouped",
                      group = "geno_block", stat = "mean", ndraws = 200) +
  ggtitle("Posterior predictive: mean response per group")
save_plot(p_ppc_grp, "22_ppc_mean_per_group.png", width = 12, height = 8)

# --- Per-animal check (subset for readability) ---
animals_subset <- sample(levels(df_resp$animal), min(20, nlevels(df_resp$animal)))
df_subset_idx  <- df_resp$animal %in% animals_subset
p_ppc_anim <- pp_check(fit_model, type = "stat_grouped",
                       group = "animal", stat = "mean", ndraws = 100,
                       newdata = df_resp[df_subset_idx, ]) +
  ggtitle("Posterior predictive: per-animal mean (subset)")
save_plot(p_ppc_anim, "23_ppc_per_animal.png", width = 14, height = 10)

# ==== Bayesian R² and LOO ====
# --- Bayesian R² ---
r2 <- bayes_R2(fit_model)
print(r2)
write.csv(r2, file.path(save_results_dir, "bayes_R2.csv"))

# --- LOO cross-validation ---
loo_result <- loo(fit_model)
print(loo_result)

# Save LOO summary
sink(file.path(save_results_dir, "loo_summary.txt"))
print(loo_result)
sink()

# Check for problematic observations
pareto_k <- loo_result$diagnostics$pareto_k
cat("\nPareto k > 0.7 (problematic):", sum(pareto_k > 0.7), "\n")
cat("Pareto k > 1.0 (very bad):    ", sum(pareto_k > 1.0), "\n")

# ==============================================================================
# --- Build prediction grid ---
new_data <- expand_grid(
  stimulus0 = 0:max(df_resp$stimulus0),
  Genotype  = levels(df_resp$Genotype),
  Block     = levels(df_resp$Block)
) %>%
  mutate(
    stimulus = stimulus0 + 1,
    animal   = NA  # population-level
  )

# --- Get posterior expected values (probabilities), marginalizing over animals ---
pred <- new_data %>%
  add_epred_draws(fit_model, re_formula = NA, ndraws = 1000) %>%
  group_by(Genotype, Block, stimulus0) %>%
  summarise(
    mean   = mean(.epred),
    median = median(.epred),
    lower  = quantile(.epred, 0.025),
    upper  = quantile(.epred, 0.975),
    lower50 = quantile(.epred, 0.25),
    upper50 = quantile(.epred, 0.75),
    .groups = "drop"
  )

# --- Compute empirical means per stimulus × group for overlay ---
emp <- df_resp %>%
  group_by(Genotype, Block, stimulus0) %>%
  summarise(
    p_obs = mean(move),
    n     = n(),
    .groups = "drop"
  )

# --- Main habituation curve plot ---
p_curves <- ggplot(pred, aes(x = stimulus0)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Genotype), alpha = 0.2) +
  geom_ribbon(aes(ymin = lower50, ymax = upper50, fill = Genotype), alpha = 0.35) +
  geom_line(aes(y = mean, color = Genotype), linewidth = 1) +
  geom_point(data = emp, aes(y = p_obs, color = Genotype, size = n),
             alpha = 0.5, show.legend = TRUE) +
  scale_size_continuous(range = c(0.5, 2.5), name = "n trials") +
  facet_grid(Block ~ Genotype) +
  labs(
    title    = "Posterior habituation curves",
    subtitle = "Ribbons: 50% and 95% credible intervals. Points: observed proportions.",
    x = "Stimulus number (0-indexed)",
    y = "P(response)"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position = "bottom")

save_plot(p_curves, "30_habituation_curves.png", width = 14, height = 7)


# ==============================================================================
# Plot Habituation curves
# ==============================================================================
save_fig_dir = NULL

p_prob <- plot_habituation_probability(
  df_resp = df_resp,
  fit_model = fit_model,
  save_dir = save_fig_dir,
  var_name = var_name
)

p_prob_raw <- plot_habituation_probability_raw(
  df_resp = df_resp,
  fit_model = fit_model,
  save_dir = save_fig_dir,
  var_name = var_name
)

p_logit_raw <- plot_habituation_logit_raw(
  df_resp = df_resp,
  fit_model = fit_model,
  save_dir = save_fig_dir,
  var_name = var_name
)

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