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

source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/nlme_utils.R")
source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/plot_utils.R")

# source("C:/UniFreiburg/Code/R_code/susana/nlme_utils.R")
# source("C:/UniFreiburg/Code/R_code/susana/plot_utils.R")

# source("D:/Behavior_Data/R_code/susana/nlme_utils.R")
# source("D:/Behavior_Data/R_code/susana/plot_utils.R")

base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <- "D:/WorkingData/Susana/DarkFlash_ISI90s_2Blocks"
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
save_fig_dir <- NULL

p_hab <- plot_habituation_probability(
  df_resp = df_resp,
  fit_model = fit_model,
  save_fig_dir = save_fig_dir,
  var_name = var_name,
  raw_data = "binary"
)

p_hab

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
  ndraws = 300
)

print(p_animal_avg)

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


# ==============================================================================
# Simple aggregate response-probability exponential fits
# ==============================================================================

library(purrr)
library(broom)

# ------------------------------------------------------------------------------
# 1. Compute response probability per Genotype × Block × stimulus
# ------------------------------------------------------------------------------

df_prob_agg <- df_final %>%
  mutate(
    stimulus = as.numeric(stimulus),
    stimulus0 = stimulus - 1,
    Genotype = factor(Genotype),
    Block = factor(Block),
    move = as.integer(move)
  ) %>%
  group_by(Genotype, Block, stimulus, stimulus0) %>%
  summarise(
    n_animals = n(),
    n_response = sum(move, na.rm = TRUE),
    response_prob = mean(move, na.rm = TRUE),
    .groups = "drop"
  )

print(df_prob_agg)

# ------------------------------------------------------------------------------
# 2. Fit simple exponential per Genotype × Block
# response_prob = A + (R0 - A) * exp(-exp(logk) * stimulus0)
# ------------------------------------------------------------------------------

fit_simple_exp <- function(dat) {
  tryCatch(
    {
      nls(
        response_prob ~ A + (R0 - A) * exp(-exp(logk) * stimulus0),
        data = dat,
        weights = n_animals,
        start = list(
          A    = max(0.01, min(dat$response_prob, na.rm = TRUE)),
          R0   = min(0.99, max(dat$response_prob, na.rm = TRUE)),
          logk = log(0.1)
        ),
        algorithm = "port",
        lower = c(
          A    = 0.001,
          R0   = 0.001,
          logk = log(0.0001)
        ),
        upper = c(
          A    = 0.999,
          R0   = 0.999,
          logk = log(10)
        ),
        control = nls.control(
          maxiter = 500,
          warnOnly = TRUE
        )
      )
    },
    error = function(e) {
      message(
        "NLS failed for ",
        unique(dat$Genotype), " ",
        unique(dat$Block), ": ",
        e$message
      )
      return(NULL)
    }
  )
}

simple_fits <- df_prob_agg %>%
  group_by(Genotype, Block) %>%
  group_split() %>%
  map(fit_simple_exp)

fit_keys <- df_prob_agg %>%
  distinct(Genotype, Block) %>%
  arrange(Genotype, Block)

# ------------------------------------------------------------------------------
# 3. Extract fitted parameters
# ------------------------------------------------------------------------------

simple_exp_params <- map2_dfr(
  simple_fits,
  seq_along(simple_fits),
  function(mod, i) {
    if (is.null(mod)) {
      tibble(
        Genotype = fit_keys$Genotype[i],
        Block = fit_keys$Block[i],
        A = NA_real_,
        R0 = NA_real_,
        logk = NA_real_,
        k = NA_real_
      )
    } else {
      cc <- coef(mod)
      
      tibble(
        Genotype = fit_keys$Genotype[i],
        Block = fit_keys$Block[i],
        A = unname(cc["A"]),
        R0 = unname(cc["R0"]),
        logk = unname(cc["logk"]),
        k = exp(unname(cc["logk"]))
      )
    }
  }
) %>%
  mutate(
    half_life_stimuli = log(2) / k
  )

print(simple_exp_params)

write.csv(
  simple_exp_params,
  file = file.path(base_dir, "results", "simple_aggregate_exp_response_prob_parameters.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Generate smooth fitted curves
# ------------------------------------------------------------------------------

simple_exp_pred <- df_prob_agg %>%
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
    stimulus0 = stimulus - 1
  ) %>%
  left_join(simple_exp_params, by = c("Genotype", "Block")) %>%
  mutate(
    fit = A + (R0 - A) * exp(-exp(logk) * stimulus0)
  )

# ------------------------------------------------------------------------------
# 5. Plot simple aggregate exponential curves
# ------------------------------------------------------------------------------

p_simple_exp <- ggplot(simple_exp_pred, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  geom_point(
    data = df_prob_agg,
    aes(x = stimulus, y = response_prob, color = Genotype),
    alpha = 0.45,
    size = 1,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    aes(y = fit),
    linewidth = 1.3
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)
  ) +
  
  theme_pubr(base_size = 14) +
  labs(
    x = "Stimulus number within block",
    y = "Response probability",
    color = "Genotype",
    fill = "Genotype",
    title = "Simple aggregate exponential fits"
  ) +
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_simple_exp)

ggsave(
  filename = file.path(base_dir, "figs", "Simple_aggregate_exp_response_prob_curves.png"),
  plot = p_simple_exp,
  width = 14,
  height = 7,
  dpi = 300
)

# ==============================================================================
# Create brms priors from aggregate simple-fit parameters
# ==============================================================================

library(dplyr)
library(brms)

# Priors on logit scale
simple_exp_priors <- simple_exp_params %>%
  mutate(
    A_logit  = qlogis(A),
    R0_logit = qlogis(R0),
    logk_prior = logk
  )

# priors on prob. scale
# simple_exp_priors <- simple_exp_params %>%
#   mutate(
#     A_logit  = A,
#     R0_logit = R0,
#     logk_prior = logk
#   )

prior_summary <- simple_exp_priors %>%
  summarise(
    A_mean = mean(A_logit, na.rm = TRUE),
    A_sd   = sd(A_logit, na.rm = TRUE),
    
    R0_mean = mean(R0_logit, na.rm = TRUE),
    R0_sd   = sd(R0_logit, na.rm = TRUE),
    
    logk_mean = mean(logk_prior, na.rm = TRUE),
    logk_sd   = sd(logk_prior, na.rm = TRUE)
  ) %>%
  mutate(
    A_sd_prior    = max(A_sd, 1.0),
    R0_sd_prior   = max(R0_sd, 1.0),
    logk_sd_prior = max(logk_sd, 1.0)
  )

print(prior_summary)

priors <- c(
  prior_string(
    paste0(
      "normal(",
      round(prior_summary$A_mean, 3), ", ",
      round(prior_summary$A_sd_prior, 3),
      ")"
    ),
    class = "b",
    nlpar = "A"
  ),
  
  prior_string(
    paste0(
      "normal(",
      round(prior_summary$R0_mean, 3), ", ",
      round(prior_summary$R0_sd_prior, 3),
      ")"
    ),
    class = "b",
    nlpar = "R0"
  ),
  
  prior_string(
    paste0(
      "normal(",
      round(prior_summary$logk_mean, 3), ", ",
      round(prior_summary$logk_sd_prior, 3),
      ")"
    ),
    class = "b",
    nlpar = "logk"
  ),
  
  prior(exponential(2), class = "sd", nlpar = "A"),
  prior(exponential(2), class = "sd", nlpar = "R0"),
  prior(exponential(2), class = "sd", nlpar = "logk")
)

print(priors)

