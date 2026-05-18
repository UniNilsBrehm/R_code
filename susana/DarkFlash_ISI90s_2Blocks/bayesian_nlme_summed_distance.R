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


source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/utils.R")
# source("C:/UniFreiburg/Code/R_code/susana/utils.R")

base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <- "D:/WorkingData/Susana/DarkFlash_ISI90s_2Blocks"
file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

var_name = 'summed_distance'

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
summary(df_resp$max_peak)
any(df_resp$max_peak <= 0, na.rm = TRUE)

# ==============================================================================
# The Model
# ==============================================================================
model <- bf(
  max_cumsum ~ Asym + (R0 - Asym) * exp(-exp(lrc) * stimulus0),
  
  Asym ~ Genotype * Block + (1 | animal),
  R0   ~ Genotype * Block + (1 | animal),
  lrc  ~ Genotype * Block + (1 | animal),
  
  nl = TRUE
)

# ==============================================================================
# The Priors
# ==============================================================================
priors<- c(
  # Population-level nonlinear parameters
  # log-scale because family = Gamma(link = "log")
  set_prior("normal(1.6, 0.5)", nlpar = "Asym", class = "b"),
  set_prior("normal(2.2, 0.5)", nlpar = "R0",   class = "b"),
  
  # Habituation rate: exp(lrc)
  # lrc around -2.5 means rate ≈ 0.08 per stimulus
  set_prior("normal(-2.5, 0.5)", nlpar = "lrc", class = "b"),
  
  # Animal-level variation
  set_prior("exponential(5)", nlpar = "Asym", class = "sd"),
  set_prior("exponential(5)", nlpar = "R0",   class = "sd"),
  set_prior("exponential(5)", nlpar = "lrc",  class = "sd"),
  
  # Gamma shape
  set_prior("exponential(1)", class = "shape")
)

# ==============================================================================
# Fit the model
# ==============================================================================
fit_model <- brm(
  formula = model,
  data = df_resp,
  family = Gamma(link = "log"),
  prior = priors,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  threads = threading(4),
  iter = 5000,
  warmup = 1500,
  seed = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
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

