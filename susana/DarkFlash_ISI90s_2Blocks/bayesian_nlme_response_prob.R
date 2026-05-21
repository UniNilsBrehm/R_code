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

# source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/utils.R")
# source("C:/UniFreiburg/Code/R_code/susana/utils.R")
source("D:/Behavior_Data/R_code/susana/utils.R")
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
# Directly on prob scale
model <- bf(
  # R0 (dp0) must be larger than Asymptote (pinf)
  move ~ inv_logit(pinf) + (1 - inv_logit(pinf)) * inv_logit(dp0) * exp(-exp(lk) * stimulus0),
  
  # More general version: But can produce rising curves
  # move ~ inv_logit(p_inf_raw) + (inv_logit(p0_raw) - inv_logit(p_inf_raw)) * exp(-exp(lk) * stimulus0),
  
  # inv_logit(pinf) keeps the asymptote between 0 and 1.
  # inv_logit(dp0) keeps the starting offset fraction between 0 and 1.
  # exp(lk) keeps the decay rate positive.
  # Gurantees: 0≤p(move=1)≤1
  
  pinf ~ Genotype * Block + (1 | animal),  # on logit scale, transform with inv_logit(pinf)
  dp0  ~ Genotype * Block + (1 | animal),  # on logit scale, fraction of the remaining distance from the asymptote to 1.
  lk   ~ Genotype * Block + (1 | animal),  # log rate (k=exp(lk))
  
  nl = TRUE
)


# ==============================================================================
# The Priors
# ==============================================================================

# Directly on prob scale
priors <- c(
  # pinf: asymptotic/lower probability
  # inv_logit(0) = 0.50
  set_prior("normal(0, 1)", nlpar = "pinf", class = "b"),
  
  # dp0: fraction of remaining distance from pinf to 1
  # inv_logit(1.4) ≈ 0.80
  set_prior("normal(1.4, 1)", nlpar = "dp0", class = "b"),
  
  # lk: log rate
  # exp(-1.5) ≈ 0.22
  set_prior("normal(-1.5, 1)", nlpar = "lk", class = "b"),
  
  # animal-level variation
  set_prior("exponential(2)", nlpar = "pinf", class = "sd"),
  set_prior("exponential(2)", nlpar = "dp0", class = "sd"),
  set_prior("exponential(2)", nlpar = "lk", class = "sd")
)

# ==============================================================================
# Validate Priors
# ==============================================================================
fit_prior_only <- brm(
  formula = model,
  data = df_resp,
  family = bernoulli(link = "logit"),
  prior = priors,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  iter = 2000,          
  warmup = 1000,
  seed = 42,
  sample_prior = "only" 
)

validate_response_prob_priors(fit_prior_only)

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
  data = df_test_sub,               # Using the xx% subset
  # family = bernoulli(link = "logit"),
  family = bernoulli(link = "identity"),
  prior = priors,
  backend = "cmdstanr",
  algorithm = "meanfield",          # <--- Reliable, fast VI method
  iter = 10000                      # VI likes higher iterations (it's still lightning fast)
)

# fit_model <- fit_vi_test

# ==============================================================================
# Fit the model
# ==============================================================================
# directly on prob scale
fit_model <- brm(
  formula = model,
  data = df_resp,
  family = bernoulli(link = "identity"),
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
sink(file.path(save_results_dir, paste0(var_name, '_model_summary.txt')))
summary(fit_model)
sink()

# ==============================================================================
# Get diagnostics
# ==============================================================================
diag_dir <- file.path(base_dir, "models", "diagnostics", "nlme", var_name)
compute_diagnostics(fit_model, diag_dir, var_name)

# Check the correlations between the main non-linear parameters
# identifiability of nonlinear parameters
save_plot_as_png(
  paste0("nlme_", var_name, "_pairs_plot.png"),
  quote(
    pairs(
      fit_model,
      variable = c(
        "b_pinf_Intercept",
        "b_dp0_Intercept",
        "b_lk_Intercept"
      )
    )
  ),
  width = 2200,
  height = 2200
)


# Random effects
re_df <- ranef(fit_model)$animal

re_long <- bind_rows(
  as.data.frame(re_df[, , "pinf_Intercept"]) %>%
    mutate(animal = rownames(re_df), nlpar = "Asym"),
  as.data.frame(re_df[, , "dp0_Intercept"]) %>%
    mutate(animal = rownames(re_df), nlpar = "R0"),
  as.data.frame(re_df[, , "lk_Intercept"]) %>%
    mutate(animal = rownames(re_df), nlpar = "lrc")
)

p_re <- ggplot(re_long,
               aes(x = reorder(animal, Estimate),
                   y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
                width = 0) +
  coord_flip() +
  facet_wrap(~nlpar, scales = "free_x") +
  theme_bw() +
  labs(x = "Animal", y = "Random effect estimate")

ggsave(
  file.path(diag_dir, paste0("nlme_", var_name, "_random_effects.png")),
  p_re,
  width = 12,
  height = 10,
  bg = "white",
  dpi = 300
)


# ==============================================================================
# Plot Habituation curves
# ==============================================================================
p_prob <- plot_habituation_probability(
  df_resp,
  fit_model,
  save_fig_dir,
  var_name
)

print(p_prob)

p_prob_raw <- plot_habituation_probability_raw(
  df_resp,
  fit_model,
  save_fig_dir,
  var_name
)

print(p_prob_raw)

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
    .groups = "drop"
  )

print(k_summary)

k_comparison <- compare_nlpar(
  k_draws_df %>% rename(value = k),
  value_name = "k",
  ratio = TRUE
)

comparison_summary <- k_comparison$summary
print(comparison_summary)

write.csv(
  comparison_summary,
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
# 3. Compare dp0 on probability scale
# ==============================================================================

dp0_draws_df <- make_nlpar_draws(
  fit_model = fit_model,
  df_resp = df_resp,
  nlpar = "dp0",
  transform_fun = plogis
) %>%
  rename(dp0 = value, dp0_logit = value_raw)

dp0_summary <- summarise_nlpar(
  dp0_draws_df %>% rename(value = dp0),
  "dp0"
)

print(dp0_summary)

dp0_comparison <- compare_nlpar(
  dp0_draws_df %>% rename(value = dp0),
  value_name = "dp0",
  rope = 0.02,
  ratio = FALSE
)

dp0_comparison_summary <- dp0_comparison$summary
print(dp0_comparison_summary)

write.csv(
  dp0_comparison_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_dp0_all_pairwise_comparisons.csv")),
  row.names = FALSE
)

# ==============================================================================
# 4. Optional: compare actual starting probability p0
# p0 = pinf + (1 - pinf) * dp0
# ==============================================================================

p0_draws_df <- pinf_draws_df %>%
  select(draw, Genotype, Block, pinf) %>%
  inner_join(
    dp0_draws_df %>%
      select(draw, Genotype, Block, dp0),
    by = c("draw", "Genotype", "Block")
  ) %>%
  mutate(
    p0 = pinf + (1 - pinf) * dp0
  )

p0_summary <- p0_draws_df %>%
  group_by(Genotype, Block) %>%
  summarise(
    p0_median = median(p0),
    p0_low = quantile(p0, 0.025),
    p0_high = quantile(p0, 0.975),
    .groups = "drop"
  )

print(p0_summary)

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
# 5. Plots
# ==============================================================================
# Plot Functions
plot_posterior_density <- function(
    draws_df,
    xvar,
    title,
    xlab,
    filename,
    block_limits = NULL
) {
  
  p <- ggplot(
    draws_df,
    aes(x = {{ xvar }}, fill = Genotype)
  ) +
    facet_wrap(~Block, scales = "free") +
    geom_density(alpha = 0.35) +
    theme_pubr(base_size = 14) +
    labs(
      x = xlab,
      y = "Posterior density",
      title = title,
      fill = "Genotype"
    )
  
  # Apply separate x-axis limits per Block
  if (!is.null(block_limits)) {
    
    scale_list <- lapply(
      names(block_limits),
      function(b) {
        as.formula(
          paste0(
            'Block == "', b,
            '" ~ scale_x_continuous(limits = c(',
            block_limits[[b]][1], ", ",
            block_limits[[b]][2], "))"
          )
        )
      }
    )
    
    p <- p +
      ggh4x::facetted_pos_scales(
        x = scale_list
      )
  }
  
  print(p)
  
  ggsave(
    file.path(save_fig_dir, filename),
    p,
    width = 10,
    height = 5,
    dpi = 300
  )
  
  p
}

plot_pairwise_differences <- function(summary_df, ylab, title, filename) {
  
  p <- summary_df %>%
    ggplot(
      aes(
        x = comparison,
        y = median_difference,
        ymin = diff_low,
        ymax = diff_high,
        color = Genotype_1
      )
    ) +
    facet_wrap(~Block, scales = "free_x") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(linewidth = 0.8) +
    coord_flip() +
    theme_pubr(base_size = 14) +
    labs(
      x = "Comparison",
      y = ylab,
      title = title,
      color = "Genotype 1"
    )
  
  print(p)
  
  ggsave(
    file.path(save_fig_dir, filename),
    p,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  p
}

# k posterior and ratio plot

p_k <- plot_posterior_density(
  k_draws_df,
  k,
  "Posterior distributions of habituation rate k",
  "Habituation rate k",
  paste0("nlme_", var_name, "_k_posteriors.png"),
  
  block_limits = list(
    Block1 = c(0, 2),
    Block2 = c(0, 0.10)
  )
)

p_k_compare <- comparison_summary %>%
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
  dpi = 300
)

# pinf plots

p_pinf <- plot_posterior_density(
  pinf_draws_df,
  pinf,
  "Posterior distributions of asymptote pinf",
  "Asymptote probability pinf",
  paste0("nlme_", var_name, "_pinf_posteriors.png")
)

p_pinf_compare <- plot_pairwise_differences(
  pinf_comparison_summary,
  "pinf difference",
  "All pairwise pinf comparisons",
  paste0("nlme_", var_name, "_pinf_all_pairwise_differences.png")
)

# dp0 plots

p_dp0 <- plot_posterior_density(
  dp0_draws_df,
  dp0,
  "Posterior distributions of starting offset fraction dp0",
  "Starting offset fraction dp0",
  paste0("nlme_", var_name, "_dp0_posteriors.png")
)

p_dp0_compare <- plot_pairwise_differences(
  dp0_comparison_summary,
  "dp0 difference",
  "All pairwise dp0 comparisons",
  paste0("nlme_", var_name, "_dp0_all_pairwise_differences.png")
)

# optional p0 plots

p_p0 <- plot_posterior_density(
  p0_draws_df,
  p0,
  "Posterior distributions of starting probability p0",
  "Starting probability p0",
  paste0("nlme_", var_name, "_p0_posteriors.png")
)

p_p0_compare <- plot_pairwise_differences(
  p0_comparison_summary,
  "p0 difference",
  "All pairwise p0 comparisons",
  paste0("nlme_", var_name, "_p0_all_pairwise_differences.png")
)
