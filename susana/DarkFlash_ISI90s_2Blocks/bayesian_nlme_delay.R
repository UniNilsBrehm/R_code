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
# source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/plot_utils.R")

# source("C:/UniFreiburg/Code/R_code/susana/nlme_utils.R")
# source("C:/UniFreiburg/Code/R_code/susana/plot_utils.R")

source("D:/Behavior_Data/R_code/susana/nlme_utils.R")
source("D:/Behavior_Data/R_code/susana/plot_utils.R")

base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <- "D:/WorkingData/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <-

file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

var_name = 'delay'
col_name = 'delay'

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
    delay = ordered(delay, levels = c(0, 1, 2, 3, 4)),
    animal = interaction(Video, Well, drop = TRUE)
  )

# Check ordinal delay distribution
summary(df_resp[[col_name]])
table(df_resp[[col_name]], useNA = "ifany")
is.ordered(df_resp[[col_name]])

# ==============================================================================
# The Model
# ==============================================================================
  model <- bf(
    as.formula(paste0(col_name, " ~ R0 + (4 - R0) * stimulus0 / (exp(logK) + stimulus0)")),
    
    R0   ~ Genotype * Block + (1 | animal),
    logK ~ Genotype * Block + (1 | animal),
    
    nl = TRUE
  )

# ==============================================================================
# The Priors
# ==============================================================================
priors <- c(
  set_prior("normal(0, 1)", class = "b", nlpar = "R0"),
  set_prior("normal(log(3), 0.7)", class = "b", nlpar = "logK"),
  
  set_prior("exponential(2)", class = "sd", nlpar = "R0"),
  set_prior("exponential(2)", class = "sd", nlpar = "logK"),
  
  set_prior("normal(0, 2)", class = "Intercept")
)

# ==============================================================================
# Validate Priors
# ==============================================================================
fit_prior_only <- brm(
  formula = model,
  data = df_resp,
  family = cumulative(link = "logit"),
  prior = priors,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  iter = 2000,          
  warmup = 1000,
  seed = 42,
  sample_prior = "only" 
)

# Prior predictive check for ordinal delay response
# ==============================================================================
yrep_prior <- posterior_predict(fit_prior_only, ndraws = 50)

# Make observed response numeric 0–4 for bayesplot
y_delay <- as.numeric(as.character(df_resp[[col_name]]))

ppc_bars(
  y = y_delay,
  yrep = yrep_prior
) +
  ggtitle("Prior Predictive Check: Ordinal Delay")

# Analytical prior validation for ordinal nonlinear delay model
# ==============================================================================

rm(list = intersect(c("sim_asym", "sim_r0", "sim_lrc", "sim_k"), ls()))

set.seed(42)

# Match priors for latent-scale nonlinear parameters
sim_asym <- rnorm(10000, mean = 0, sd = 1)
sim_r0   <- rnorm(10000, mean = 0, sd = 1)

# lrc is log-rate; k = exp(lrc)
sim_lrc <- rnorm(10000, mean = -1.5, sd = 1)
sim_k   <- exp(sim_lrc)

print("Prior latent-scale quantiles for Asym:")
quantile(sim_asym, probs = c(0.025, 0.5, 0.975))

print("Prior latent-scale quantiles for R0:")
quantile(sim_r0, probs = c(0.025, 0.5, 0.975))

print("Prior habituation-rate quantiles for k = exp(lrc):")
quantile(sim_k, probs = c(0.025, 0.5, 0.975))

print("Prior half-life quantiles, log(2) / k:")
quantile(log(2) / sim_k, probs = c(0.025, 0.5, 0.975))

# ==============================================================================
# Fit Fast Test Model
# ==============================================================================
# Step 1: Create a tiny xx% subset of your data for rapid prototyping
df_test_sub <- df_resp %>% 
  group_by(Genotype, Block) %>% 
  slice_sample(prop = 0.99) %>% 
  ungroup()

# Step 2: Fit using Meanfield Variational Inference (VI)
fit_vi_test <- brm(
  formula = model,
  data = df_resp,
  family = cumulative(link = "logit"),
  prior = priors,
  backend = "cmdstanr",
  algorithm = "meanfield",          # <--- Reliable, fast VI method
  iter = 10000                      # VI likes higher iterations (it's still lightning fast)
)

fit_model <- fit_vi_test

# # COMPARE TWO MODELS
# # 1. Compute LOO for both models
# loo_1 <- loo(fit_1, cores = 4)
# loo_2    <- loo(fit_2, cores = 4)
# 
# # 2. Compare them
# loo_compare(loo_1, loo_2)
# 
# # WAIC
# waic_1 <- waic(fit_1)
# waic_2    <- waic(fit_2)
# 
# loo_compare(waic_1, waic_2)

# ==============================================================================
# Fit the model
# ==============================================================================
fit_delay <- brm(
  formula = model,
  data = df_resp,
  
  family = cumulative(link = "logit"),
  
  prior = priors,
  
  backend = "cmdstanr",
  
  chains = 4,
  cores = 4,
  threads = threading(6),
  
  iter = 4000,
  warmup = 2000,
  
  seed = 42,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  )
)

# ==============================================================================
# Get summary and diagnostics
# ==============================================================================
summary(fit_model)

diag_dir <- file.path(base_dir, "models", "diagnostics", var_name)

# Trace plots
trace_plots <- plot(fit_model, ask = FALSE)
for (i in seq_along(trace_plots)) {
  ggsave(
    filename = file.path(
      diag_dir,
      paste0("nlme_", var_name, "_traceplot_", i, ".png")
    ),
    plot = trace_plots[[i]],
    width = 12,
    height = 8,
    dpi = 300
  )
}

# Posterior predictive checks
p1 <- pp_check(fit_model, ndraws = 100)

ggsave(
  file.path(diag_dir, paste0("nlme_", var_name,"_ppcheck_default.png")),
  p1,
  width = 10,
  height = 8,
  dpi = 300
)


p2 <- pp_check(
  fit_model,
  type = "dens_overlay",
  ndraws = 100
)

ggsave(
  file.path(diag_dir, paste0("nlme_", var_name,"_densed_overlay.png")),
  p2,
  width = 10,
  height = 8,
  dpi = 300
)


p3 <- pp_check(
  fit_model,
  type = "hist",
  ndraws = 100
)

ggsave(
  file.path(diag_dir,paste0("nlme_", var_name,"hist.png")),
  p3,
  width = 10,
  height = 8,
  dpi = 300
)

p4 <- pp_check(
  fit_model,
  type = "ecdf_overlay",
  ndraws = 100
)

ggsave(
  file.path(diag_dir, paste0("nlme_", var_name, "_ppcheck_ecdf.png")),
  p4,
  width = 10,
  height = 8,
  dpi = 300
)


p5 <- pp_check(
  fit_model,
  type = "stat_grouped",
  group = "stimulus0",
  ndraws = 100
)

ggsave(
  file.path(diag_dir, paste0("nlme_", var_name,"_ppcheck_stimulus.png")),
  p5,
  width = 12,
  height = 8,
  dpi = 300
)


p6 <- pp_check(
  fit_model,
  type = "stat_grouped",
  group = "Genotype",
  ndraws = 100
)

ggsave(
  file.path(diag_dir, paste0("nlme_", var_name, "_ppcheck_genotype.png")),
  p6,
  width = 10,
  height = 8,
  dpi = 300
)

p7 <- pp_check(
  fit_model,
  type = "stat_grouped",
  group = "Block",
  ndraws = 100
)

ggsave(
  file.path(diag_dir, paste0("nlme_", var_name, "_ppcheck_block.png")),
  p7,
  width = 10,
  height = 8,
  dpi = 300
)

# Check the correlations between the main non-linear parameters
# identifiability of nonlinear parameters
save_plot_as_png(
  paste0("nlme_", var_name, "_pairs_plot.png"),
  quote(
    pairs(
      fit_model,
      variable = c(
        "b_Asym_Intercept",
        "b_R0_Intercept",
        "b_lrc_Intercept"
      )
    )
  ),
  width = 2200,
  height = 2200
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

# Random effects
re_df <- ranef(fit_model)$animal

re_long <- bind_rows(
  as.data.frame(re_df[, , "Asym_Intercept"]) %>%
    mutate(animal = rownames(re_df), nlpar = "Asym"),
  as.data.frame(re_df[, , "R0_Intercept"]) %>%
    mutate(animal = rownames(re_df), nlpar = "R0"),
  as.data.frame(re_df[, , "lrc_Intercept"]) %>%
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
  dpi = 300
)

# Conditional effects plots
ce <- conditional_effects(
  fit_model,
  effects = "stimulus0:Genotype",
  re_formula = NA
)
save_plot_as_png(
  paste0("nlme_", var_name, "_conditional_effects.png"),
  quote(plot(ce)),
  width = 2200,
  height = 1800
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
# Plot delay habituation curves: expected delay scale
# ==============================================================================
delay_levels <- c(0, 1, 2, 3, 4)

new_delay <- df_resp %>%
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

# posterior_epred gives category probabilities:
# draws x observations x response categories
epred_delay <- posterior_epred(
  fit_model,
  newdata = new_delay,
  re_formula = NA
)

# Convert category probabilities to expected delay values
delay_draws <- apply(epred_delay, c(1, 2), function(p) {
  sum(p * delay_levels)
})

pred_delay_data <- bind_cols(
  new_delay,
  tibble(
    fit = apply(delay_draws, 2, median),
    CI_low = apply(delay_draws, 2, quantile, probs = 0.025),
    CI_high = apply(delay_draws, 2, quantile, probs = 0.975)
  )
)

raw_delay <- df_resp %>%
  mutate(
    delay_num = as.numeric(as.character(delay))
  ) %>%
  group_by(Genotype, Block, stimulus) %>%
  summarise(
    mean_delay = mean(delay_num, na.rm = TRUE),
    median_delay = median(delay_num, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_delay_exp <- ggplot(
  pred_delay_data,
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  geom_point(
    data = raw_delay,
    aes(x = stimulus, y = mean_delay),
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
  
  scale_y_continuous(
    limits = c(0, 4),
    breaks = delay_levels
  ) +
  
  theme_pubr(base_size = 14) +
  labs(
    x = "Stimulus number within block",
    y = "Expected delay",
    title = "Bayesian nonlinear ordinal delay habituation curves",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_delay_exp)

ggsave(
  filename = file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_delay_habituation_curves_expected.png")
  ),
  plot = p_delay_exp,
  width = 14,
  height = 7,
  dpi = 300
)

# ==============================================================================
# Plot delay habituation curves WITH true raw delay data
# ==============================================================================

p_delay_raw <- ggplot(
  pred_delay_data,
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  geom_jitter(
    data = df_resp %>%
      mutate(delay_num = as.numeric(as.character(delay))),
    aes(
      x = stimulus,
      y = delay_num,
      color = Genotype
    ),
    inherit.aes = FALSE,
    width = 0.15,
    height = 0.08,
    alpha = 0.12,
    size = 0.7
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
  
  scale_y_continuous(
    limits = c(0, 4),
    breaks = delay_levels
  ) +
  
  theme_pubr(base_size = 14) +
  labs(
    x = "Stimulus number within block",
    y = "Delay",
    title = "Bayesian nonlinear ordinal delay curves with raw data",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_delay_raw)

ggsave(
  filename = file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_delay_habituation_curves_true_raw_data.png")
  ),
  plot = p_delay_raw,
  width = 14,
  height = 7,
  dpi = 300
)

# ==============================================================================
# Plot latent-scale nonlinear delay curves
# ==============================================================================
new_delay_latent <- df_resp %>%
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

latent_draws <- posterior_linpred(
  fit_model,
  newdata = new_delay_latent,
  re_formula = NA,
  transform = FALSE
)

pred_latent_delay_data <- bind_cols(
  new_delay_latent,
  tibble(
    fit = apply(latent_draws, 2, median),
    CI_low = apply(latent_draws, 2, quantile, probs = 0.025),
    CI_high = apply(latent_draws, 2, quantile, probs = 0.975)
  )
)

p_delay_latent <- ggplot(
  pred_latent_delay_data,
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
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
    y = "Latent delay tendency",
    title = "Bayesian nonlinear delay curves on latent ordinal scale",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_delay_latent)

ggsave(
  filename = file.path(
    base_dir,
    "figs",
    paste0("nlme_", var_name, "_delay_habituation_curves_latent_scale.png")
  ),
  plot = p_delay_latent,
  width = 14,
  height = 7,
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

