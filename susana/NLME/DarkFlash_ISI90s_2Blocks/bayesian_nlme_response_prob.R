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

# source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/nlme_utils.R")
# source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/plot_utils.R")

source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/nlme_utils.R")
source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/plot_utils.R")

# source("D:/Behavior_Data/R_code/susana/nlme_utils.R")
# source("D:/Behavior_Data/R_code/susana/plot_utils.R")

# base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/NLME//DarkFlash_ISI90s_2Blocks"
base_dir <- "D:/WorkingData/Susana/NLME/DarkFlash_ISI90s_2Blocks"
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

df_resp %>%
  distinct(animal, Genotype) %>%
  count(Genotype)

# ==============================================================================
# The Model
# ==============================================================================
model <- bf(
  move ~ inv_logit(A) + (inv_logit(R0) - inv_logit(A)) * exp(-exp(logk) * stimulus0),
  
  A    ~ 1 + Genotype * Block + (1 | animal),
  R0   ~ 1 + Genotype * Block + (1 | animal),
  logk ~ 1 + Genotype * Block + (1 | animal),
  
  nl = TRUE
)

# ==============================================================================
# The Priors
# ==============================================================================
priors <- c(
  prior(normal(0, 1),     class = "b", coef = "Intercept", nlpar = "A"),
  prior(normal(1.5, 1),   class = "b", coef = "Intercept", nlpar = "R0"),
  prior(normal(-3, 1),    class = "b", coef = "Intercept", nlpar = "logk"),
  
  prior(normal(0, 0.75), class = "b", nlpar = "A"),
  prior(normal(0, 0.75), class = "b", nlpar = "R0"),
  prior(normal(0, 0.75), class = "b", nlpar = "logk"),
  
  prior(exponential(4), class = "sd", nlpar = "A"),
  prior(exponential(4), class = "sd", nlpar = "R0"),
  prior(exponential(4), class = "sd", nlpar = "logk")
)
# validate_prior(priors, formula = model, data = df_resp, family = bernoulli(link = "identity"))

# ==============================================================================
# Validate Priors
# ==============================================================================
fit_prior <- brm(
  formula = model,
  data = df_resp,
  family = bernoulli(link = "identity"),
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

ep <- posterior_epred(fit_prior, ndraws = 500)

quantile(as.vector(ep), c(.001,.01,.05,.5,.95,.99,.999), na.rm = TRUE)
quantile(df_resp[[col_name]], c(.001,.01,.05,.5,.95,.99,.999), na.rm = TRUE)

prior_draws <- as_draws_df(fit_prior)
summary(prior_draws)


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
  family = bernoulli(link = "identity"),
  prior = priors,
  backend = "cmdstanr",
  algorithm = "meanfield",          # <--- Reliable, fast VI method
  # algorithm = "fullrank",
  iter = 10000                      # VI likes higher iterations (it's still lightning fast)
)

fit_vi_fullrank <- brm(
  formula = model,
  data = df_resp,
  family = bernoulli(link = "identity"),
  prior = priors,
  backend = "cmdstanr",
  algorithm = "fullrank",
  iter = 10000,
  seed = 42
)

fit_nuts_test <- brm(
  formula = model,
  data = df_resp,
  family = bernoulli(link = "identity"),
  prior = priors,
  save_pars = save_pars(all = TRUE),
  init    = 0,                        # important for nonlinear models
  backend = "cmdstanr",
  chains = 2,
  cores = 2,
  threads = threading(12),
  iter = 1000,
  warmup = 500,
  seed = 42,
  control = list(
    adapt_delta = 0.90,
    max_treedepth = 10
  )
)

fit_model <- fit_vi_test

fit_model <- fit_vi_fullrank

fit_model <- fit_nuts_test

# ==============================================================================
# Fit the model
# ==============================================================================
fit_model <- brm(
  formula = model,
  data = df_resp,
  family = bernoulli(link = "identity"),
  prior = priors,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  threads = threading(6),
  iter = 6000,
  warmup = 3000,
  seed = 42,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  init    = 0,                        # important for nonlinear models
  file = file.path(base_dir, "models", paste0("bayesian_nlme_", var_name,"_results.rds")),
  file_refit = "on_change"
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

# Divergences
sum(subset(nuts_params(fit_model), Parameter == "divergent__")$Value)
# 0 ideal, single digits tolerable, more is concerning

# Energy diagnostics
bayesplot::mcmc_nuts_energy(nuts_params(fit_model))

# Posterior predictive check
pp_check(fit_model, type = "dens_overlay", ndraws = 100)  # marginal fit
pp_check(fit_model, type = "stat_grouped", group = "stimulus0",
         stat = "median", ndraws = 100)             # does decay shape fit?


bayesplot::mcmc_pairs(
  fit_model, np = nuts_params(fit_model),
  off_diag_args = list(size = 0.5)
)

# Residuals
yrep <- posterior_predict(fit_model, ndraws = 1000)

fit_mu <- fitted(fit_model)[, "Estimate"]

sim_res <- createDHARMa(
  simulatedResponse = t(yrep),
  observedResponse = df_resp$move,
  fittedPredictedResponse = fit_mu,
  integerResponse = TRUE
)

plot(sim_res)
plotResiduals(sim_res, df_resp$stimulus)
plotResiduals(sim_res, df_resp$stimulus0)
plotResiduals(sim_res, df_resp$Genotype)
plotResiduals(sim_res, df_resp$Block)
plotResiduals(sim_res, df_resp$animal)

# Population-level residuals
yrep <- posterior_predict(
  fit_model,
  ndraws = 500,
  re_formula = NA
)

fit_mu <- fitted(
  fit_model,
  re_formula = NA
)[, "Estimate"]

sim_res_pop <- createDHARMa(
  simulatedResponse = t(yrep),
  observedResponse = df_resp$move,
  fittedPredictedResponse = fit_mu,
  integerResponse = TRUE
)

plot(sim_res_pop)

# Compute leave-one-out cross-validation:
loo_var <- loo(fit_model, moment_match = TRUE)
print(loo_var)
save_plot_as_png(
  paste0("nlme_", var_name, "_loo_plot.png"),
  quote(plot(loo_var))
)
# loo_compare(loo1, loo2)  # compare models


diag_plots <- plot(fit_model)
out_dir <- file.path(base_dir, "models", "diagnostics", "nlme", var_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(diag_plots)) {
  ggsave(
    filename = file.path(out_dir, paste0("fit_model_diagnostic_", i, ".png")),
    plot = diag_plots[[i]],
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

# ==============================================================================
# Plot Habituation curves
# ==============================================================================
p_hab_ind <- plot_habituation_probability(
  df_resp = df_resp,
  fit_model = fit_model,
  save_fig_dir = save_fig_dir,
  var_name = var_name,
  raw_data = "individual",
  ndraws = 500
)

p_hab_ind

p_hab_raw <- plot_habituation_probability(
  df_resp = df_resp,
  fit_model = fit_model,
  save_fig_dir = save_fig_dir,
  var_name = var_name,
  raw_data = "binary"
)

p_hab_raw

p_hab_agg <- plot_habituation_probability(
  df_resp = df_resp,
  fit_model = fit_model,
  save_fig_dir = save_fig_dir,
  var_name = var_name,
  # raw_data = "binary"
  raw_data = "aggregate"
)

p_hab_agg

# Only for diagnostics
p_animal_avg <- plot_habituation_probability_animal_averaged(
  df_resp = df_resp,
  fit_model = fit_model,
  save_fig_dir = save_fig_dir,
  var_name = "response_prob",
  ndraws = 500
)

print(p_animal_avg)

# ==============================================================================
# Compare nonlinear model parameters
# ==============================================================================

# ==============================================================================
# 1. Compare habituation rate k = exp(logk)
# ==============================================================================
k_draws_df <- make_nlpar_draws(
  fit_model    = fit_model,
  df_resp      = df_resp,
  nlpar        = "logk",          # brms nlpar name in the model
  transform_fun = exp
) %>%
  rename(k = value, logk = value_raw)

k_summary <- k_draws_df %>%
  group_by(Genotype, Block) %>%
  summarise(
    k_median         = median(k),
    k_low            = quantile(k, 0.025),
    k_high           = quantile(k, 0.975),
    half_life_median = median(log(2) / k),
    half_life_low    = quantile(log(2) / k, 0.025),
    half_life_high   = quantile(log(2) / k, 0.975),
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
# 2. Compare asymptote A on probability scale (pinf = inv_logit(A))
# ==============================================================================
pinf_draws_df <- make_nlpar_draws(
  fit_model     = fit_model,
  df_resp       = df_resp,
  nlpar         = "A",            # brms nlpar name in the model
  transform_fun = plogis          # inv_logit
) %>%
  rename(pinf = value, A_logit = value_raw)

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
# 3. Compare initial response probability R0 on probability scale (p0 = inv_logit(R0))
# ==============================================================================
p0_draws_df <- make_nlpar_draws(
  fit_model     = fit_model,
  df_resp       = df_resp,
  nlpar         = "R0",           # brms nlpar name in the model
  transform_fun = plogis          # inv_logit
) %>%
  rename(p0 = value, R0_logit = value_raw)

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
    Block1 = c(0, 1),
    Block2 = c(0, 0.05)
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
  save_fig_dir = save_fig_dir,
  block_limits = list(
    Block1 = c(0, 1),
    Block2 = c(0, 1)
  )
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
  save_fig_dir = save_fig_dir,
  block_limits = list(
    Block1 = c(0.7, 1),
    Block2 = c(0.7, 1)
  )
)


# p0 pairwise differences

p_p0_compare <- plot_pairwise_differences(
  p0_comparison_summary,
  "p0 difference",
  "All pairwise p0 comparisons",
  paste0("nlme_", var_name, "_p0_all_pairwise_differences.png"),
  save_fig_dir = save_fig_dir
)

# ==============================================================================




# ==============================================================================
# Endpoint contrasts: average of last 3 stimuli (stimuli 58, 59, 60)
# ==============================================================================

last_n <- 5
last_stims <- (max(df_resp$stimulus) - last_n + 1):max(df_resp$stimulus)
# e.g. stimuli 58, 59, 60  →  stimulus0 = 57, 58, 59

grid_end <- expand.grid(
  Genotype  = levels(df_resp$Genotype),
  Block     = levels(df_resp$Block),
  stimulus  = last_stims
) %>%
  mutate(
    stimulus0 = stimulus - 1,
    Genotype  = factor(Genotype, levels = levels(df_resp$Genotype)),
    Block     = factor(Block,    levels = levels(df_resp$Block))
  )

# Full posterior draws of the population-level prediction
pred_end_draws <- fitted(
  fit_model,
  newdata    = grid_end,
  re_formula = NA,
  summary    = FALSE   # draws × rows  (ndraws × nrows)
)

# Average over the last 3 stimuli per genotype × block
# grid_end rows: for each (Genotype, Block) combo, 3 consecutive stimulus rows
# We average the 3 stimulus draws for each genotype × block combination
n_combos  <- nrow(grid_end) / last_n   # = n_genotypes × n_blocks
n_draws   <- nrow(pred_end_draws)

# Build a draws data frame averaged over the last 3 stimuli
draws_end <- bind_rows(
  lapply(levels(df_resp$Genotype), function(g) {
    bind_rows(
      lapply(levels(df_resp$Block), function(b) {
        row_idx <- which(grid_end$Genotype == g & grid_end$Block == b)
        tibble(
          draw     = seq_len(n_draws),
          Genotype = factor(g, levels = levels(df_resp$Genotype)),
          Block    = factor(b, levels = levels(df_resp$Block)),
          value    = rowMeans(pred_end_draws[, row_idx])
        )
      })
    )
  })
)

# Verify before proceeding
glimpse(draws_end)
print(distinct(draws_end, Genotype, Block))
# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
end_summary <- draws_end %>%
  group_by(Genotype, Block) %>%
  summarise(
    median = median(value),
    low    = quantile(value, 0.025),
    high   = quantile(value, 0.975),
    .groups = "drop"
  )

print(end_summary)

write.csv(
  end_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_endpoint_summary.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Pairwise comparisons
# ------------------------------------------------------------------------------
end_comparison <- compare_nlpar(
  draws_end,
  value_name = "p_end",
  rope       = 0.02,
  ratio      = FALSE
)

end_comparison_summary <- end_comparison$summary
print(end_comparison_summary)

write.csv(
  end_comparison_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_endpoint_all_pairwise_comparisons.csv")),
  row.names = FALSE
)

# ==============================================================================
# Plots
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Posterior density of endpoint response probability
# ------------------------------------------------------------------------------
p_end_density <- ggplot(
  draws_end,
  aes(x = value, fill = Genotype)
) +
  facet_wrap(~Block, scales = "free_y") +
  geom_density(alpha = 0.35) +
  geom_vline(
    data = end_summary,
    aes(xintercept = median, color = Genotype),
    linewidth = 0.8, linetype = "dashed"
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_pubr(base_size = 14) +
  labs(
    x     = paste0("Response probability (mean of last ", last_n, " stimuli)"),
    y     = "Posterior density",
    title = paste0("Posterior distributions of endpoint response probability\n(mean of last ", last_n, " stimuli)"),
    fill  = "Genotype",
    color = "Genotype"
  )

print(p_end_density)

ggsave(
  file.path(save_fig_dir, paste0("nlme_", var_name, "_endpoint_posteriors.png")),
  p_end_density,
  width = 10, height = 5, dpi = 300, bg = "white"
)

# ------------------------------------------------------------------------------
# 2. Pointrange summary plot
# ------------------------------------------------------------------------------
p_end_summary <- end_summary %>%
  ggplot(aes(
    x     = Genotype,
    y     = median,
    ymin  = low,
    ymax  = high,
    color = Genotype
  )) +
  facet_wrap(~Block) +
  geom_pointrange(linewidth = 0.8, size = 0.8) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x     = "Genotype",
    y     = paste0("Response probability (last ", last_n, " stimuli)"),
    title = paste0("Endpoint response probability (last ", last_n, " stimuli)")
  )

print(p_end_summary)

ggsave(
  file.path(save_fig_dir, paste0("nlme_", var_name, "_endpoint_summary_pointrange.png")),
  p_end_summary,
  width = 8, height = 5, dpi = 300, bg = "white"
)

# ------------------------------------------------------------------------------
# 3. Pairwise differences plot
# ------------------------------------------------------------------------------
p_end_compare <- plot_pairwise_differences(
  end_comparison_summary,
  "Endpoint probability difference",
  paste0("All pairwise endpoint comparisons (last ", last_n, " stimuli)"),
  paste0("nlme_", var_name, "_endpoint_all_pairwise_differences.png"),
  save_fig_dir = save_fig_dir
)

print(p_end_compare)

# ------------------------------------------------------------------------------
# 4. Habituation curves with endpoint marked
# ------------------------------------------------------------------------------
# Add a vertical line and the endpoint summaries overlaid on the habituation plot
p_hab_with_endpoint <- plot_habituation_probability(
  df_resp    = df_resp,
  fit_model  = fit_model,
  raw_data   = "aggregate"
) +
  # shaded region for last 3 stimuli
  geom_vline(
    xintercept = min(last_stims) - 0.5,
    linetype   = "dashed",
    color      = "grey40",
    linewidth  = 0.6
  ) +
  # endpoint posterior median + 95% CI as a cross
  geom_pointrange(
    data = end_summary,
    aes(
      x     = max(last_stims) + 1.5,  # just to the right of the last stimulus
      y     = median,
      ymin  = low,
      ymax  = high,
      color = Genotype
    ),
    inherit.aes = FALSE,
    size      = 0.8,
    linewidth = 0.8
  ) +
  labs(caption = paste0(
    "Dashed line marks start of last ", last_n,
    " stimuli window. Pointrange to the right shows posterior median + 95% CI of endpoint."
  ))

print(p_hab_with_endpoint)

ggsave(
  file.path(save_fig_dir, paste0("nlme_", var_name, "_habituation_curves_with_endpoint.png")),
  p_hab_with_endpoint,
  width = 14, height = 7, dpi = 300, bg = "white"
)




# ==============================================================================
# Extrapolation plot (set stimulus range here)
# ==============================================================================
extrap_to <- 600

new_extrap <- df_resp %>%
  group_by(Genotype, Block) %>%
  summarise(.groups = "drop") %>%
  reframe(
    stimulus = seq(1, extrap_to, length.out = extrap_to * 2),
    .by = c(Genotype, Block)
  ) %>%
  mutate(
    stimulus0 = stimulus - 1,
    Genotype  = factor(Genotype, levels = levels(df_resp$Genotype)),
    Block     = factor(Block,    levels = levels(df_resp$Block))
  )

pred_extrap <- fitted(
  fit_model,
  newdata    = new_extrap,
  re_formula = NA,
  summary    = TRUE
)

pred_extrap_data <- bind_cols(new_extrap, as.data.frame(pred_extrap)) %>%
  rename(fit = Estimate, CI_low = Q2.5, CI_high = Q97.5)

pinf_lines <- pinf_draws_df %>%
  group_by(Genotype, Block) %>%
  summarise(
    pinf_median = median(pinf),
    .groups = "drop"
  )

p_extrap <- ggplot(
  pred_extrap_data,
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  geom_ribbon(
    aes(ymin = CI_low, ymax = CI_high),
    alpha = 0.15, color = NA
  ) +
  
  geom_line(aes(y = fit), linewidth = 1.0) +
  
  geom_vline(
    xintercept = max(df_resp$stimulus),
    linetype   = "dashed",
    color      = "grey30",
    linewidth  = 0.6
  ) +
  
  geom_hline(
    data      = pinf_lines,
    aes(yintercept = pinf_median, color = Genotype),
    linetype  = "dotted",
    linewidth = 0.8
  ) +
  
  geom_text(
    data  = pinf_lines,
    aes(
      x     = extrap_to * 0.95,
      y     = pinf_median + 0.03,
      label = paste0("p\u221e=", round(pinf_median, 2)),
      color = Genotype
    ),
    inherit.aes = FALSE,
    size        = 3,
    hjust       = 1
  ) +
  
  annotate(
    "text",
    x     = max(df_resp$stimulus) + extrap_to * 0.01,
    y     = 0.5,
    label = "data ends",
    size  = 3,
    color = "grey30",
    hjust = 0,
    angle = 90
  ) +
  
  coord_cartesian(ylim = c(0, 1), xlim = c(1, extrap_to + 10)) +
  
  theme_pubr(base_size = 14) +
  labs(
    x       = "Stimulus number",
    y       = "Response probability",
    title   = paste0(
      "Extrapolation to stimulus ", extrap_to,
      " (dashed = data end, dotted = p\u221e)"
    ),
    color   = "Genotype",
    fill    = "Genotype",
    caption = "Dotted horizontal line = posterior median asymptote (p\u221e). Dashed vertical line = end of observed data."
  ) +
  theme(
    legend.position = "top",
    panel.spacing   = unit(1.2, "lines")
  )

print(p_extrap)

ggsave(
  file.path(
    save_fig_dir,
    paste0("nlme_", var_name, "_habituation_curves_extrapolated_", extrap_to, ".png")
  ),
  p_extrap,
  width  = 14,
  height = 7,
  dpi    = 300,
  bg     = "white"
)