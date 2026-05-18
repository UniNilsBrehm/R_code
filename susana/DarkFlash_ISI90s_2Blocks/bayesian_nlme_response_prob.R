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

source("C:/UniFreiburg/Code/R_code/susana/utils.R")

base_dir <- "D:/WorkingData/Susana/DarkFlash_ISI90s_2Blocks"
file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

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
bf_resp_exp <- bf(
  move ~ Asym + (R0 - Asym) * exp(-exp(lrc) * stimulus0),
  
  Asym ~ Genotype * Block + (1 | animal),
  R0   ~ Genotype * Block + (1 | animal),
  lrc  ~ Genotype * Block + (1 | animal),
  
  nl = TRUE
)

# ==============================================================================
# The Priors
# ==============================================================================
priors_resp_exp <- c(
  set_prior("normal(0, 1.5)", nlpar = "Asym", class = "b"),
  set_prior("normal(0, 2)",   nlpar = "R0", class = "b"),
  set_prior("normal(-1.5, 1)", nlpar = "lrc", class = "b"),
  
  set_prior("exponential(2)", nlpar = "Asym", class = "sd"),
  set_prior("exponential(2)", nlpar = "R0", class = "sd"),
  set_prior("exponential(2)", nlpar = "lrc", class = "sd")
)

# ==============================================================================
# Fit the model
# ==============================================================================
fit_resp_exp <- brm(
  formula = bf_resp_exp,
  data = df_resp,
  family = bernoulli(link = "logit"),
  prior = priors_resp_exp,
  
  backend = "cmdstanr",
  
  chains = 4,
  cores = 4,
  threads = threading(4),
  
  iter = 5000,
  warmup = 1500,
  seed = 42,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

# ==============================================================================
# Get summary and diagnostics
# ==============================================================================
summary(fit_resp_exp)
plot(fit_resp_exp)
pp_check(fit_resp_exp)
# pp_check(fit_resp_exp, type = "bars", ndraws = 100)
# Check the correlations between the main non-linear parameters
pairs(fit_resp_exp, variable = c("b_Asym_Intercept", "b_R0_Intercept", "b_lrc_Intercept"))

# ==============================================================================
# Save fitted model to HDD
# ==============================================================================
saveRDS(
  fit_resp_exp,
  file = file.path(base_dir, "data_files", "bayesian_nlme_response_prob_results.rds")
)

# ==============================================================================
# Load fitted model if available
# ==============================================================================
fit_resp_exp <- readRDS(
  file.path(base_dir, "data_files", "bayesian_nlme_response_prob_results.rds")
)

# ==============================================================================
# Plot Habituation curves
# ==============================================================================
new_resp <- df_resp %>%
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

pred_resp <- fitted(
  fit_resp_exp,
  newdata = new_resp,
  re_formula = NA,
  summary = TRUE
)

pred_resp_data <- bind_cols(new_resp, as.data.frame(pred_resp)) %>%
  rename(
    fit = Estimate,
    CI_low = Q2.5,
    CI_high = Q97.5
  )

raw_resp <- df_resp %>%
  group_by(Genotype, Block, stimulus) %>%
  summarise(
    response_prob = mean(move, na.rm = TRUE),
    .groups = "drop"
  )

p_resp_exp <- ggplot(pred_resp_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  geom_point(
    data = raw_resp,
    aes(x = stimulus, y = response_prob),
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
    y = "Response probability",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_resp_exp)

ggsave(
  filename = file.path(
    base_dir,
    "figs",
    "response_prob_nlme_habituation_curves.png"
  ),
  plot = p_resp_exp,
  width = 14,
  height = 7,
  dpi = 300
)


# ==============================================================================
# Plot latent-scale (logit-scale) exponential curves WITH raw data
# ==============================================================================
# ------------------------------------------------------------------------------
# 1. Create prediction grid
# ------------------------------------------------------------------------------
new_resp_logit <- df_resp %>%
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
# 2. Get latent-scale predictions
# ------------------------------------------------------------------------------

pred_logit <- fitted(
  fit_resp_exp,
  newdata = new_resp_logit,
  re_formula = NA,
  scale = "linear",
  summary = TRUE
)

pred_logit_data <- bind_cols(
  new_resp_logit,
  as.data.frame(pred_logit)
) %>%
  rename(
    fit = Estimate,
    CI_low = Q2.5,
    CI_high = Q97.5
  )

# ------------------------------------------------------------------------------
# 3. Convert raw probabilities to logit scale
# ------------------------------------------------------------------------------

raw_prob <- df_resp %>%
  group_by(Genotype, Block, stimulus) %>%
  summarise(
    response_prob = mean(move, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  
  # avoid Inf logits
  mutate(
    response_prob_clipped = pmin(
      pmax(response_prob, 0.001),
      0.999
    ),
    
    logit_prob = qlogis(response_prob_clipped)
  )

# ------------------------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------------------------

p_logit <- ggplot(
  pred_logit_data,
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  # Raw data points on logit scale
  geom_point(
    data = raw_prob,
    aes(
      x = stimulus,
      y = logit_prob,
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
    y = "Latent response tendency (logit scale)",
    title = "Bayesian nonlinear habituation curves on latent logit scale",
    color = "Genotype",
    fill = "Genotype"
  ) +
  
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_logit)

ggsave(
  filename = file.path(
    base_dir,
    "figs",
    "response_prob_nlme_habituation_curves_logit_scale_with_raw.png"
  ),
  plot = p_logit,
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
  fit_resp_exp,
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
    "response_prob_nlme_habituation_rate_all_pairwise_comparisons.csv"
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
    "response_prob_nlme_habituation_rate_posteriors.png"
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
    "response_prob_nlme_habituation_rate_all_pairwise_ratios.png"
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
  fit_resp_exp,
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
    "response_prob_nlme_Asym_all_pairwise_comparisons.csv"
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
    "response_prob_nlme_Asym_posteriors.png"
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
    "response_prob_nlme_Asym_all_pairwise_differences.png"
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
  fit_resp_exp,
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
    "response_prob_nlme_R0_all_pairwise_comparisons.csv"
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
    "response_prob_nlme_R0_posteriors.png"
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
    "response_prob_nlme_R0_all_pairwise_differences.png"
  ),
  p_r0_compare,
  width = 10,
  height = 8,
  dpi = 300
)

