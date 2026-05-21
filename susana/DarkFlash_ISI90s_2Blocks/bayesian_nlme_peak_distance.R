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

source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/utils.R")
source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/plot_utils.R")

# source("C:/UniFreiburg/Code/R_code/susana/nlme_utils.R")
# source("C:/UniFreiburg/Code/R_code/susana/plot_utils.R")

# source("D:/Behavior_Data/R_code/susana/nlme_utils.R")
# source("D:/Behavior_Data/R_code/susana/plot_utils.R")

base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <- "D:/WorkingData/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <-

file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

var_name = 'peak_distance'
col_name = 'max_peak'

# ==============================================================================
# Load and prepare data
# ==============================================================================
res <- load_data_darkflash_60s(file_dir, move_th = 0, take_peak = 0)

df_final     <- res$df_final
df_final_sub <- res$df_final_sub

df_resp <- df_final_sub %>%
  mutate(
    stimulus = as.numeric(stimulus),
    stimulus0 = stimulus - 1,
    Block = factor(Block),
    Genotype = factor(Genotype),
    Video = factor(Video),
    Well = factor(Well),
    animal = interaction(Video, Well, drop = TRUE)
    
  )

# Check if all values are larger than zero (for Gamma)
summary(df_resp[[col_name]])
any(df_resp[[col_name]] <= 0, na.rm = TRUE)

# Histogram
ggplot(df_resp, aes(x = .data[[col_name]])) +
  geom_histogram(
    bins = 100,
    fill = "lightblue",
    color = "white"
  ) +
  labs(
    title = paste("Histogram of", col_name),
    x = col_name,
    y = "Count"
  ) +
  theme_minimal()

# ==============================================================================
# The Model
# ==============================================================================
model <- bf(
  max_peak ~ exp(A) + (exp(R0) - exp(A)) * exp(-exp(logk) * stimulus0),
  A    ~ 1 + Genotype * Block + (1 | animal),   # <-- 1 +
  R0   ~ 1 + Genotype * Block + (1 | animal),
  logk ~ 1 + Genotype * Block + (1 | animal),
  nl = TRUE
)

# ==============================================================================
# The Priors
# ==============================================================================
priors <- c(
  # Reference cell baseline
  prior(normal(2.0, 0.3),  nlpar = "R0",   class = "b", coef = "Intercept"),
  prior(normal(0.5, 0.3),  nlpar = "A",    class = "b", coef = "Intercept"),
  prior(normal(-2.3, 0.4), nlpar = "logk", class = "b", coef = "Intercept"),
  # Deviations from reference (all genotype contrasts, block, interactions)
  prior(normal(0, 0.2), nlpar = "R0",   class = "b"),
  prior(normal(0, 0.2), nlpar = "A",    class = "b"),
  prior(normal(0, 0.2), nlpar = "logk", class = "b"),
  # Animal REs
  prior(student_t(3, 0, 0.2), nlpar = "R0",   class = "sd"),
  prior(student_t(3, 0, 0.2), nlpar = "A",    class = "sd"),
  prior(student_t(3, 0, 0.2), nlpar = "logk", class = "sd"),
  prior(gamma(2, 0.1), class = "shape")
)

# ==============================================================================
# Validate Priors
# ==============================================================================
fit_prior_only <- brm(
  formula = model,
  data = df_resp,
  family = Gamma(link = "identity"),
  prior = priors,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  seed = 42,
  sample_prior = "only"
)

ep <- posterior_epred(fit_prior_only, ndraws = 500)

quantile(as.vector(ep), c(.001,.01,.05,.5,.95,.99,.999), na.rm = TRUE)
quantile(df_resp[[col_name]], c(.001,.01,.05,.5,.95,.99,.999), na.rm = TRUE)

prior_draws <- as_draws_df(fit_prior_only)
summary(prior_draws)


# ==============================================================================
# Fit Fast Test Model
# ==============================================================================
# Fit using Meanfield Variational Inference (VI)
fit_vi_test <- brm(
  formula = model,
  data = df_resp,               
  family  = Gamma(link = "identity"),
  prior = priors,
  backend = "cmdstanr",
  algorithm = "meanfield",          # <--- Reliable, fast VI method
  # algorithm = "fullrank",
  iter = 10000                      # VI likes higher iterations (it's still lightning fast)
)

fit_nuts_test <- brm(
  formula = model,
  data = df_resp,
  family  = Gamma(link = "identity"),
  prior = priors,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr",
  chains = 2,
  cores = 2,
  threads = threading(6),
  iter = 1000,
  warmup = 500,
  seed = 42,
  control = list(adapt_delta = 0.90, max_treedepth = 10),
  init    = 0,                        # important for nonlinear models
)

fit_model <- fit_vi_test
fit_model <- fit_nuts_test

# ==============================================================================
# Fit the model
# ==============================================================================
fit_model <- brm(
  model,
  data    = df_resp,
  family  = Gamma(link = "identity"),
  prior   = priors,
  save_pars = save_pars(all = TRUE),
  recompile = TRUE,
  chains  = 4, 
  cores = 4,
  threads = threading(6),
  iter    = 4000, 
  warmup = 1500,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  init    = 0,
  backend = "cmdstanr"
)

# ==============================================================================
# Get summary and diagnostics
# ==============================================================================

diag_dir <- file.path(base_dir, "models", "diagnostics", var_name)

# Sampler diagnostics
# Top-level summary
summary(fit_model)             # all Rhat < 1.01, Bulk_ESS > 400, Tail_ESS > 400

# Divergences
sum(subset(nuts_params(fit_model), Parameter == "divergent__")$Value)
# 0 ideal, single digits tolerable, more is concerning

# Energy diagnostics
bayesplot::mcmc_nuts_energy(nuts_params(fit_model))

# Posterior predictive check
pp_check(fit_model, type = "dens_overlay", ndraws = 100)  # marginal fit
pp_check(fit_model, type = "stat_grouped", group = "stimulus0",
         stat = "median", ndraws = 100)             # does decay shape fit?
pp_check(fit_model, type = "stat_grouped", group = "Genotype",
         stat = "mean", ndraws = 100)               # do genotype means fit?

df_resp |>
  add_predicted_draws(fit_model, ndraws = 200) |>
  group_by(Genotype, Block, stimulus0) |>
  median_qi(.prediction, .width = c(0.5, 0.9)) |>
  ggplot(aes(x = stimulus0, y = .prediction)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  geom_point(
    data = df_resp |>
      group_by(Genotype, Block, stimulus0) |>
      summarise(observed = median(max_peak), .groups = "drop"),
    aes(x = stimulus0, y = observed),    # <-- added x = stimulus0
    color = "red", size = 0.8, inherit.aes = FALSE
  ) +
  scale_fill_brewer() +
  facet_grid(Block ~ Genotype) +
  labs(y = "max_peak", x = "stimulus")

# Conditional curves per genotype × block
conditional_effects(fit_model, effects = "stimulus0:Genotype",
                    conditions = data.frame(Block = c("Block1", "Block2")))

bayesplot::mcmc_pairs(
  fit_model, np = nuts_params(fit_model),
  pars = c("b_logk_Intercept", "sd_animal__logk_Intercept"),
  off_diag_args = list(size = 0.5)
)

# Residuals
yrep <- posterior_predict(fit_model, ndraws = 200)

sim_res <- createDHARMa(
  simulatedResponse = t(yrep),
  observedResponse = df_resp$max_peak,
  fittedPredictedResponse = fitted(fit_model)[, "Estimate"]
)

save_plot_as_png(
  paste0("nlme_", var_name, "_DHARMa_residuals.png"),
  quote(plot(sim_res))
)

# Compute leave-one-out cross-validation:
loo_var <- loo(fit_model)
print(loo_var)
save_plot_as_png(
  paste0("nlme_", var_name, "_loo_plot.png"),
  quote(plot(loo_var))
)
# loo_compare(loo1, loo2)  # compare models

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
# Plot Habituation curves
# ==============================================================================
new_peak <- df_resp %>%
  group_by(Genotype, Block) %>%
  summarise(
    stim_min = min(stimulus, na.rm = TRUE),
    stim_max = max(stimulus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Genotype, Block) %>%
  reframe(
    stimulus = seq(stim_min, stim_max, length.out = 100)
  ) %>%
  mutate(
    stimulus0 = stimulus - 1,
    Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
    Block = factor(Block, levels = levels(df_resp$Block))
  )

pred_peak <- fitted(
  fit_model,
  newdata = new_peak,
  re_formula = NA,
  summary = TRUE
)

pred_peak_data <- bind_cols(new_peak, as.data.frame(pred_peak)) %>%
  rename(
    fit = Estimate,
    CI_low = Q2.5,
    CI_high = Q97.5
  )

raw_peak <- df_resp %>%
  group_by(Genotype, Block, stimulus) %>%
  summarise(
    mean_peak = mean(max_peak, na.rm = TRUE),
    .groups = "drop"
  )

p_peak_exp <- ggplot(pred_peak_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  geom_point(
    data = raw_peak,
    aes(x = stimulus, y = mean_peak),
    alpha = 0.35,
    size = 1,
    inherit.aes = FALSE
  ) +
  
  geom_ribbon(
    aes(ymin = CI_low, ymax = CI_high),
    alpha = 0.15,
    color = NA
  ) +
  
  geom_line(
    aes(y = fit),
    linewidth = 1.3
  ) +
  
  theme_pubr(base_size = 14) +
  labs(
    x = "Stimulus number within block",
    y = "Peak distance moved",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_peak_exp)

ggsave(
  filename = file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_habituation_curves.png")
  ),
  plot = p_peak_exp,
  width = 14,
  height = 7,
  dpi = 300
)
# ==============================================================================
# Plot latent-scale (log-scale) exponential curves
# ==============================================================================
# ------------------------------------------------------------------------------
# 1. Create prediction grid
# ------------------------------------------------------------------------------
new_peak_log <- df_resp %>%
  group_by(Genotype, Block) %>%
  summarise(
    stim_min = min(stimulus, na.rm = TRUE),
    stim_max = max(stimulus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Genotype, Block) %>%
  reframe(
    stimulus = seq(stim_min, stim_max, length.out = 200)
  ) %>%
  mutate(
    stimulus0 = stimulus - 1,
    Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
    Block = factor(Block, levels = levels(df_resp$Block))
  )

# ------------------------------------------------------------------------------
# 2. Get latent-scale predictions: log(expected max_peak)
# ------------------------------------------------------------------------------

pred_log_peak <- fitted(
  fit_model,
  newdata = new_peak_log,
  re_formula = NA,
  scale = "linear",
  summary = TRUE
)

pred_log_peak_data <- bind_cols(
  new_peak_log,
  as.data.frame(pred_log_peak)
) %>%
  rename(
    fit = Estimate,
    CI_low = Q2.5,
    CI_high = Q97.5
  )

# ------------------------------------------------------------------------------
# 3. Convert raw peak distances to log scale
# ------------------------------------------------------------------------------

raw_peak_log <- df_resp %>%
  filter(max_peak > 0) %>%
  group_by(Genotype, Block, stimulus) %>%
  summarise(
    mean_peak = mean(max_peak, na.rm = TRUE),
    median_peak = median(max_peak, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    log_mean_peak = log(mean_peak),
    log_median_peak = log(median_peak)
  )

# ------------------------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------------------------

p_log_peak <- ggplot(
  pred_log_peak_data,
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  # Raw grouped means on log scale
  geom_point(
    data = raw_peak_log,
    aes(
      x = stimulus,
      y = log_mean_peak,
      color = Genotype
    ),
    inherit.aes = FALSE,
    alpha = 0.5,
    size = 1
  ) +
  
  # Credible intervals
  geom_ribbon(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    alpha = 0.15,
    color = NA
  ) +
  
  # Model curve
  geom_line(
    aes(y = fit),
    linewidth = 1.3
  ) +
  
  theme_pubr(base_size = 14) +
  
  labs(
    x = "Stimulus number within block",
    y = "Latent peak distance moved (log scale)",
    title = "Bayesian nonlinear habituation curves on latent log scale",
    color = "Genotype",
    fill = "Genotype"
  ) +
  
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_log_peak)

ggsave(
  filename = file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_habituation_curves_logit_scale_with_raw.png")
  ),
  plot = p_log_peak,
  width = 14,
  height = 7,
  dpi = 300
)

# ==============================================================================
# Plot Habituation curves WITH true raw peak-distance data
# ==============================================================================

new_peak <- df_resp %>%
  group_by(Genotype, Block) %>%
  summarise(
    stim_min = min(stimulus, na.rm = TRUE),
    stim_max = max(stimulus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Genotype, Block) %>%
  reframe(
    stimulus = seq(stim_min, stim_max, length.out = 100)
  ) %>%
  mutate(
    stimulus0 = stimulus - 1,
    Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
    Block = factor(Block, levels = levels(df_resp$Block))
  )

pred_peak <- fitted(
  fit_model,
  newdata = new_peak,
  re_formula = NA,
  summary = TRUE
)

pred_peak_data <- bind_cols(
  new_peak,
  as.data.frame(pred_peak)
) %>%
  rename(
    fit = Estimate,
    CI_low = Q2.5,
    CI_high = Q97.5
  )

p_peak_exp <- ggplot(
  pred_peak_data,
  aes(
    x = stimulus,
    color = Genotype,
    fill = Genotype
  )
) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  geom_jitter(
    data = df_resp,
    aes(
      x = stimulus,
      y = max_peak,
      color = Genotype
    ),
    inherit.aes = FALSE,
    width = 0.15,
    height = 0,
    alpha = 0.12,
    size = 0.7
  ) +
  
  geom_ribbon(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    alpha = 0.15,
    color = NA
  ) +
  
  geom_line(
    aes(y = fit),
    linewidth = 1.3
  ) +
  
  theme_pubr(base_size = 14) +
  labs(
    x = "Stimulus number within block",
    y = "Peak distance moved",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_peak_exp)

ggsave(
  filename = file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_habituation_curves_true_raw_data.png")
  ),
  plot = p_peak_exp,
  width = 14,
  height = 7,
  dpi = 300
)
# ==============================================================================
# Plot latent-scale (log-scale) exponential curves WITH true raw peak-distance data
# ==============================================================================

new_peak_log <- df_resp %>%
  group_by(Genotype, Block) %>%
  summarise(
    stim_min = min(stimulus, na.rm = TRUE),
    stim_max = max(stimulus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Genotype, Block) %>%
  reframe(
    stimulus = seq(stim_min, stim_max, length.out = 200)
  ) %>%
  mutate(
    stimulus0 = stimulus - 1,
    Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
    Block = factor(Block, levels = levels(df_resp$Block))
  )

pred_log_peak <- fitted(
  fit_model,
  newdata = new_peak_log,
  re_formula = NA,
  scale = "linear",
  summary = TRUE
)

pred_log_peak_data <- bind_cols(
  new_peak_log,
  as.data.frame(pred_log_peak)
) %>%
  rename(
    fit = Estimate,
    CI_low = Q2.5,
    CI_high = Q97.5
  )

raw_peak_log <- df_resp %>%
  filter(max_peak > 0) %>%
  mutate(
    log_peak = log(max_peak)
  )

p_log_peak <- ggplot(
  pred_log_peak_data,
  aes(
    x = stimulus,
    color = Genotype,
    fill = Genotype
  )
) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  geom_jitter(
    data = raw_peak_log,
    aes(
      x = stimulus,
      y = log_peak,
      color = Genotype
    ),
    inherit.aes = FALSE,
    width = 0.15,
    height = 0,
    alpha = 0.12,
    size = 0.7
  ) +
  
  geom_ribbon(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    alpha = 0.15,
    color = NA
  ) +
  
  geom_line(
    aes(y = fit),
    linewidth = 1.3
  ) +
  
  theme_pubr(base_size = 14) +
  labs(
    x = "Stimulus number within block",
    y = "Latent peak distance moved (log scale)",
    title = "Bayesian nonlinear habituation curves on latent log scale",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_log_peak)

ggsave(
  filename = file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_habituation_curves_logit_scale_true_raw_binary.png")
  ),
  plot = p_log_peak,
  width = 14,
  height = 7,
  dpi = 300
)
# ==============================================================================
# Compare habituation rate k between all Genotypes
# ==============================================================================
# ------------------------------------------------------------------------------
# 1. Create grid of Genotype × Block combinations
# ------------------------------------------------------------------------------
k_grid <- expand.grid(
  Genotype = levels(df_resp$Genotype),
  Block = levels(df_resp$Block)
) %>%
  mutate(
    stimulus0 = 0,
    Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
    Block = factor(Block, levels = levels(df_resp$Block))
  )

# ------------------------------------------------------------------------------
# 2. Extract posterior draws for lrc
# ------------------------------------------------------------------------------
lrc_draws <- posterior_linpred(
  fit_model,
  newdata = k_grid,
  nlpar = "lrc",
  re_formula = NA,
  transform = FALSE
)


# ------------------------------------------------------------------------------
# 3. Convert to k = exp(lrc)
# ------------------------------------------------------------------------------

k_draws_df <- bind_rows(
  lapply(seq_len(nrow(k_grid)), function(i) {
    
    tibble(
      draw = seq_len(nrow(lrc_draws)),
      Genotype = k_grid$Genotype[i],
      Block = k_grid$Block[i],
      lrc = lrc_draws[, i],
      k = exp(lrc_draws[, i])
    )
    
  })
)

# ------------------------------------------------------------------------------
# 4. Summary table
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# 5. All pairwise genotype comparisons within each Block
# ------------------------------------------------------------------------------

comparison_df <- k_draws_df %>%
  rename(
    Genotype_1 = Genotype,
    k_1 = k
  ) %>%
  inner_join(
    k_draws_df %>%
      rename(
        Genotype_2 = Genotype,
        k_2 = k
      ),
    by = c("draw", "Block"),
    relationship = "many-to-many"
  ) %>%
  filter(as.character(Genotype_1) < as.character(Genotype_2)) %>%
  mutate(
    comparison = paste(Genotype_1, "vs", Genotype_2),
    k_difference = k_1 - k_2,
    k_ratio = k_1 / k_2
  )

# ------------------------------------------------------------------------------
# 6. Summarise all comparisons
# ------------------------------------------------------------------------------

comparison_summary <- comparison_df %>%
  group_by(Block, Genotype_1, Genotype_2, comparison) %>%
  summarise(
    median_difference = median(k_difference),
    diff_low = quantile(k_difference, 0.025),
    diff_high = quantile(k_difference, 0.975),
    
    median_ratio = median(k_ratio),
    ratio_low = quantile(k_ratio, 0.025),
    ratio_high = quantile(k_ratio, 0.975),
    
    prob_Genotype_1_faster = mean(k_1 > k_2),
    prob_Genotype_2_faster = mean(k_1 < k_2),
    
    .groups = "drop"
  )

print(comparison_summary)

write.csv(
  comparison_summary,
  file.path(
    base_dir,
    "results",
    paste0("nlme_", var_name, "_habituation_rate_all_pairwise_comparisons.csv")
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7. Plot posterior distributions of k
# ------------------------------------------------------------------------------
p_k <- ggplot(
  k_draws_df,
  aes(x = k, fill = Genotype)
) +
  
  facet_wrap(~Block, scales = "free") +
  
  geom_density(
    alpha = 0.35
  ) +
  
  theme_pubr(base_size = 14) +
  
  labs(
    x = "Habituation rate k",
    y = "Posterior density",
    title = "Posterior distributions of habituation rate",
    fill = "Genotype"
  )

print(p_k)

ggsave(
  file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_habituation_rate_posteriors.png")
  ),
  p_k,
  width = 10,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 8. Plot all pairwise comparisons
# ------------------------------------------------------------------------------

p_compare <- comparison_summary %>%
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
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_pointrange(
    linewidth = 0.8
  ) +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x = "Comparison",
    y = "Habituation rate ratio",
    title = "All pairwise habituation rate comparisons",
    color = "Numerator genotype"
  )

print(p_compare)

ggsave(
  file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_habituation_rate_all_pairwise_ratios.png")
  ),
  p_compare,
  width = 10,
  height = 8,
  dpi = 300
)

# ==============================================================================
# Compare Asym between all Genotypes
# ==============================================================================

asym_grid <- expand.grid(
  Genotype = levels(df_resp$Genotype),
  Block = levels(df_resp$Block)
) %>%
  mutate(
    stimulus0 = 0,
    Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
    Block = factor(Block, levels = levels(df_resp$Block))
  )

asym_draws <- posterior_linpred(
  fit_model,
  newdata = asym_grid,
  nlpar = "Asym",
  re_formula = NA,
  transform = FALSE
)

asym_draws_df <- bind_rows(
  lapply(seq_len(nrow(asym_grid)), function(i) {
    tibble(
      draw = seq_len(nrow(asym_draws)),
      Genotype = asym_grid$Genotype[i],
      Block = asym_grid$Block[i],
      Asym = asym_draws[, i]
    )
  })
)

asym_summary <- asym_draws_df %>%
  group_by(Genotype, Block) %>%
  summarise(
    Asym_median = median(Asym),
    Asym_low = quantile(Asym, 0.025),
    Asym_high = quantile(Asym, 0.975),
    .groups = "drop"
  )

print(asym_summary)

asym_comparison_df <- asym_draws_df %>%
  rename(
    Genotype_1 = Genotype,
    Asym_1 = Asym
  ) %>%
  inner_join(
    asym_draws_df %>%
      rename(
        Genotype_2 = Genotype,
        Asym_2 = Asym
      ),
    by = c("draw", "Block"),
    relationship = "many-to-many"
  ) %>%
  filter(as.character(Genotype_1) < as.character(Genotype_2)) %>%
  mutate(
    comparison = paste(Genotype_1, "vs", Genotype_2),
    Asym_difference = Asym_1 - Asym_2
  )

asym_comparison_summary <- asym_comparison_df %>%
  group_by(Block, Genotype_1, Genotype_2, comparison) %>%
  summarise(
    median_difference = median(Asym_difference),
    diff_low = quantile(Asym_difference, 0.025),
    diff_high = quantile(Asym_difference, 0.975),
    prob_Genotype_1_higher = mean(Asym_1 > Asym_2),
    prob_Genotype_2_higher = mean(Asym_1 < Asym_2),
    .groups = "drop"
  )

print(asym_comparison_summary)

write.csv(
  asym_comparison_summary,
  file.path(
    base_dir,
    "results",
    paste0("nlme_", var_name, "_Asym_all_pairwise_comparisons.csv")
  ),
  row.names = FALSE
)
# ------------------------------------------------------------------------------
# Plot Asym posterior distributions
# ------------------------------------------------------------------------------

p_asym <- ggplot(
  asym_draws_df,
  aes(x = Asym, fill = Genotype)
) +
  facet_wrap(~Block, scales = "free") +
  geom_density(alpha = 0.35) +
  theme_pubr(base_size = 14) +
  labs(
    x = "Asym",
    y = "Posterior density",
    title = "Posterior distributions of Asym",
    fill = "Genotype"
  )

print(p_asym)

ggsave(
  file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_Asym_posteriors.png")
  ),
  p_asym,
  width = 10,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------------------------
# Plot all pairwise Asym comparisons
# ------------------------------------------------------------------------------

p_asym_compare <- asym_comparison_summary %>%
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
    y = "Asym difference",
    title = "All pairwise Asym comparisons",
    color = "Genotype 1"
  )

print(p_asym_compare)


ggsave(
  file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_Asym_all_pairwise_differences.png")
  ),
  p_asym_compare,
  width = 10,
  height = 8,
  dpi = 300
)

# ==============================================================================
# Compare R0 between all Genotypes
# ==============================================================================

r0_grid <- expand.grid(
  Genotype = levels(df_resp$Genotype),
  Block = levels(df_resp$Block)
) %>%
  mutate(
    stimulus0 = 0,
    Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
    Block = factor(Block, levels = levels(df_resp$Block))
  )

r0_draws <- posterior_linpred(
  fit_model,
  newdata = r0_grid,
  nlpar = "R0",
  re_formula = NA,
  transform = FALSE
)

r0_draws_df <- bind_rows(
  lapply(seq_len(nrow(r0_grid)), function(i) {
    tibble(
      draw = seq_len(nrow(r0_draws)),
      Genotype = r0_grid$Genotype[i],
      Block = r0_grid$Block[i],
      R0 = r0_draws[, i]
    )
  })
)

r0_summary <- r0_draws_df %>%
  group_by(Genotype, Block) %>%
  summarise(
    R0_median = median(R0),
    R0_low = quantile(R0, 0.025),
    R0_high = quantile(R0, 0.975),
    .groups = "drop"
  )

print(r0_summary)

r0_comparison_df <- r0_draws_df %>%
  rename(
    Genotype_1 = Genotype,
    R0_1 = R0
  ) %>%
  inner_join(
    r0_draws_df %>%
      rename(
        Genotype_2 = Genotype,
        R0_2 = R0
      ),
    by = c("draw", "Block"),
    relationship = "many-to-many"
  ) %>%
  filter(as.character(Genotype_1) < as.character(Genotype_2)) %>%
  mutate(
    comparison = paste(Genotype_1, "vs", Genotype_2),
    R0_difference = R0_1 - R0_2
  )

r0_comparison_summary <- r0_comparison_df %>%
  group_by(Block, Genotype_1, Genotype_2, comparison) %>%
  summarise(
    median_difference = median(R0_difference),
    diff_low = quantile(R0_difference, 0.025),
    diff_high = quantile(R0_difference, 0.975),
    prob_Genotype_1_higher = mean(R0_1 > R0_2),
    prob_Genotype_2_higher = mean(R0_1 < R0_2),
    .groups = "drop"
  )

print(r0_comparison_summary)

write.csv(
  r0_comparison_summary,
  file.path(
    base_dir,
    "results",
    paste0("nlme_", var_name, "_R0_all_pairwise_comparisons.csv")
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Plot R0 posterior distributions
# ------------------------------------------------------------------------------

p_r0 <- ggplot(
  r0_draws_df,
  aes(x = R0, fill = Genotype)
) +
  facet_wrap(~Block, scales = "free") +
  geom_density(alpha = 0.35) +
  theme_pubr(base_size = 14) +
  labs(
    x = "R0",
    y = "Posterior density",
    title = "Posterior distributions of R0",
    fill = "Genotype"
  )

print(p_r0)

ggsave(
  file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_R0_posteriors.png")
  ),
  p_r0,
  width = 10,
  height = 5,
  dpi = 300
)
# ------------------------------------------------------------------------------
# Plot all pairwise R0 comparisons
# ------------------------------------------------------------------------------

p_r0_compare <- r0_comparison_summary %>%
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
    y = "R0 difference",
    title = "All pairwise R0 comparisons",
    color = "Genotype 1"
  )

print(p_r0_compare)

ggsave(
  file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_R0_all_pairwise_differences.png")
  ),
  p_r0_compare,
  width = 10,
  height = 8,
  dpi = 300
)
