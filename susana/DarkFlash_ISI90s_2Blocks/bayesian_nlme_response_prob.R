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
source("C:/UniFreiburg/Code/R_code/susana/utils.R")
# source("D:/Behavior_Data/R_code/susana/utils.R")

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
  as.formula(paste0(col_name, " ~ Asym + (R0 - Asym) * exp(-exp(lrc) * stimulus0)")),
  
  Asym ~ Genotype * Block + (1 | animal),
  R0   ~ Genotype * Block + (1 | animal),
  lrc  ~ Genotype * Block + (1 | animal),

  nl = TRUE
)

# Directly on prob scale
model <- bf(
  # R0 (dp0) must be larger than Asymptote (pinf)
  move ~ inv_logit(pinf) +
    (1 - inv_logit(pinf)) * inv_logit(dp0) * exp(-exp(lk) * stimulus0),
  
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
priors <- c(
  # Shifted to center around a ~75% baseline response probability
  set_prior("normal(1.1, 1.0)", nlpar = "Asym", class = "b"),
  set_prior("normal(1.1, 1.0)",  nlpar = "R0", class = "b"),
  set_prior("normal(-1.5, 1)",  nlpar = "lrc", class = "b"), # Rate parameter stays the same
  
  # Group-level variations
  set_prior("exponential(2)", nlpar = "Asym", class = "sd"),
  set_prior("exponential(2)", nlpar = "R0", class = "sd"),
  set_prior("exponential(2)", nlpar = "lrc", class = "sd")
)

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

# 2. Visual Prior Predictive Checks for Binary Data
# BAR PLOT instead of Density: Compares the counts of 0s and 1s
ppc_bars(
  y = df_resp[[col_name]], 
  yrep = posterior_predict(fit_prior_only, ndraws = 50)
) + 
  ggtitle("Prior Predictive Check: Counts of 0s and 1s")

# PROPORTION STAT instead of Mean/Max: Checks if the predicted probability 
# of a response (y = 1) matches reality.
# We define a custom function to calculate the mean (which equals the proportion of 1s)
prop_one <- function(x) mean(x == 1)

ppc_stat(
  y = df_resp[[col_name]], 
  yrep = posterior_predict(fit_prior_only, ndraws = 200), 
  stat = "prop_one"
) +
  ggtitle("Prior Predictive Check: Proportion of Responses (1s)")

# 3. Analytically Prior Validation (Logit Scale Simulation)
# Clear the old simulations just to be safe
rm(sim_asym, sim_r0)

# Simulate 10,000 draws using the updated priors: normal(1.1, 1.0)
sim_asym <- plogis(rnorm(10000, mean = 1.1, sd = 1.0))
sim_r0   <- plogis(rnorm(10000, mean = 1.1, sd = 1.0))

# Check the new probability quantiles
print("UPDATED Prior Probability Quantiles for Asym (0.025, 0.5, 0.975):")
quantile(sim_asym, probs = c(0.025, 0.5, 0.975))

print("UPDATED Prior Probability Quantiles for R0 (0.025, 0.5, 0.975):")
quantile(sim_r0, probs = c(0.025, 0.5, 0.975))


# 4. Diagnostic & Effect Plots
# MCMC density plot of the prior distributions
mcmc_plot(fit_prior_only, type = "dens")

# Plot the conditional effects implied ONLY by your priors
plot(conditional_effects(fit_prior_only))

# ==============================================================================
# Fit Fast Test Model
# ==============================================================================
# Step 1: Create a tiny xx% subset of your data for rapid prototyping
df_test_sub <- df_resp %>% 
  group_by(Genotype, Block) %>% 
  slice_sample(prop = 0.90) %>% 
  ungroup()

# Step 2: Fit using Meanfield Variational Inference (VI)
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
  iter = 5000,
  warmup = 1500,
  seed = 42,
  control = list(adapt_delta = 0.90, max_treedepth = 10)
)

# directly on prob scale
fit <- brm(
  formula = model,
  data = df_resp,
  family = bernoulli(link = "identity"),
  prior = priors_prob,
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
# Get summary and diagnostics
# ==============================================================================
sink(file.path(save_results_dir, paste0(var_name, '_model_summary.txt')))
summary(fit_model)
sink()

diag_dir <- file.path(base_dir, "models", "diagnostics", "nlme", var_name)

bayes_R2(fit_model)

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
  observedResponse = df_resp$move,
  fittedPredictedResponse = fitted(fit_model)[, "Estimate"]
)

save_plot_as_png(
  paste0("nlme_", var_name, "_DHARMa_residuals.png"),
  quote(plot(sim_res))
)

# Compute leave-one-out cross-validation:
loo_sum <- loo(fit_model, cores=4)
print(loo_sum)
save_plot_as_png(
  paste0("nlme_", var_name, "_loo_plot.png"),
  quote(plot(loo_sum))
)

# fit_model 1: with corr
# fit_model 2: no corr
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
  fit_model,
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
  filename = file.path(save_fig_dir, paste0("nlme_", var_name, "_habituation_curves.png")),
  plot = p_resp_exp,
  width = 14,
  height = 7,
  dpi = 300
)


# ==============================================================================
# Plot latent-scale (logit-scale) exponential curves
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
  fit_model,
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
    save_fig_dir,
    paste0("nlme_", var_name, "_habituation_curves_logit_scale_with_raw.png")
  ),
  plot = p_logit,
  width = 14,
  height = 7,
  dpi = 300
)

# ==============================================================================
# Plot Habituation curves WITH true raw binary data
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
    Genotype = factor(
      Genotype,
      levels = levels(df_resp$Genotype)
    ),
    Block = factor(
      Block,
      levels = levels(df_resp$Block)
    )
  )

# ------------------------------------------------------------------------------
# Get fitted response-scale predictions
# ------------------------------------------------------------------------------

pred_resp <- fitted(
  fit_model,
  newdata = new_resp,
  re_formula = NA,
  summary = TRUE
)

pred_resp_data <- bind_cols(
  new_resp,
  as.data.frame(pred_resp)
) %>%
  rename(
    fit = Estimate,
    CI_low = Q2.5,
    CI_high = Q97.5
  )

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

p_resp_exp <- ggplot(
  pred_resp_data,
  aes(
    x = stimulus,
    color = Genotype,
    fill = Genotype
  )
) +
  
  facet_grid(
    Block ~ Genotype,
    scales = "fixed"
  ) +
  
  # TRUE raw Bernoulli observations
  geom_jitter(
    data = df_resp,
    aes(
      x = stimulus,
      y = move,
      color = Genotype
    ),
    inherit.aes = FALSE,
    width = 0.15,
    height = 0.03,
    alpha = 0.12,
    size = 0.7
  ) +
  
  # Credible interval
  geom_ribbon(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    alpha = 0.15,
    color = NA
  ) +
  
  # Model fit
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
    save_fig_dir,
    paste0("nlme_", var_name, "_habituation_curves_true_raw_data.png")
  ),
  plot = p_resp_exp,
  width = 14,
  height = 7,
  dpi = 300
)

# ==============================================================================
# Plot latent-scale (logit-scale) exponential curves WITH true raw binary data
# ==============================================================================

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

pred_logit <- fitted(
  fit_model,
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

# place true binary observations at finite logit reference values
raw_binary_logit <- df_resp %>%
  mutate(
    logit_binary_display = ifelse(move == 1, 4, -4)
  )

p_logit <- ggplot(
  pred_logit_data,
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  geom_jitter(
    data = raw_binary_logit,
    aes(
      x = stimulus,
      y = logit_binary_display,
      color = Genotype
    ),
    inherit.aes = FALSE,
    width = 0.15,
    height = 0.15,
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
    save_fig_dir,
    paste0("nlme_", var_name, "_habituation_curves_logit_scale_true_raw_binary.png")
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
    
    ROPE_prob_10percent = mean(abs(log(k_ratio)) < log(1.10)),
    
    evidence_strength = case_when(
      prob_Genotype_1_faster > 0.995 ~ "extreme",
      prob_Genotype_1_faster > 0.97  ~ "very strong",
      prob_Genotype_1_faster > 0.90  ~ "strong",
      prob_Genotype_1_faster > 0.75  ~ "moderate",
      TRUE ~ "weak"
    ),
    
    .groups = "drop"
  )

print(comparison_summary)

write.csv(
  comparison_summary,
  file.path(
    save_results_dir,
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
    save_fig_dir,
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
    save_fig_dir,
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
    
    # ROPE on logit scale
    ROPE_prob_small_effect = mean(abs(Asym_difference) < 0.2),
    
    evidence_strength = case_when(
      prob_Genotype_1_higher > 0.995 ~ "extreme",
      prob_Genotype_1_higher > 0.97  ~ "very strong",
      prob_Genotype_1_higher > 0.90  ~ "strong",
      prob_Genotype_1_higher > 0.75  ~ "moderate",
      TRUE ~ "weak"
    ),
    
    .groups = "drop"
  )

print(asym_comparison_summary)

write.csv(
  asym_comparison_summary,
  file.path(
    save_results_dir,
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
    save_fig_dir,
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
    save_fig_dir,
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
    
    # ROPE on logit scale
    ROPE_prob_small_effect = mean(abs(R0_difference) < 0.2),
    
    evidence_strength = case_when(
      prob_Genotype_1_higher > 0.995 ~ "extreme",
      prob_Genotype_1_higher > 0.97  ~ "very strong",
      prob_Genotype_1_higher > 0.90  ~ "strong",
      prob_Genotype_1_higher > 0.75  ~ "moderate",
      TRUE ~ "weak"
    ),
    
    .groups = "drop"
  )
print(r0_comparison_summary)

write.csv(
  r0_comparison_summary,
  file.path(
    save_results_dir,
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
    save_fig_dir,
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
    save_fig_dir,
    paste0("nlme_", var_name, "_R0_all_pairwise_differences.png")
  ),
  p_r0_compare,
  width = 10,
  height = 8,
  dpi = 300
)

