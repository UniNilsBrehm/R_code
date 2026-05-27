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

source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/utils.R")
source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/plot_utils.R")

# source("D:/Behavior_Data/R_code/susana/nlme_utils.R")
# source("D:/Behavior_Data/R_code/susana/plot_utils.R")

# base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_ISI90s_2Blocks"
base_dir <- "D:/WorkingData/Susana/NLME/DarkFlash_ISI90s_2Blocks"
# base_dir <-

file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

var_name = 'peak_distance'
col_name = 'max_peak'

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
  A    ~ 1 + Genotype * Block + (1 | animal),
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
  chains  = 4, 
  cores = 4,
  threads = threading(6),
  iter    = 4000, 
  warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  init    = 0,
  backend = "cmdstanr"
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

diag_dir <- file.path(base_dir, "models", "diagnostics", var_name)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

# Sampler diagnostics
# Top-level summary
summary(fit_model)             # all Rhat < 1.01, Bulk_ESS > 400, Tail_ESS > 400

# Divergences
sum(subset(nuts_params(fit_model), Parameter == "divergent__")$Value)
# 0 ideal, single digits tolerable, more is concerning

# Energy diagnostics
energy_plot <- bayesplot::mcmc_nuts_energy(nuts_params(fit_model))
ggsave(
  filename = file.path(diag_dir, "nuts_energy.png"),
  plot     = energy_plot,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# Posterior predictive check
ppc_marginal <- pp_check(fit_model, type = "dens_overlay", ndraws = 200)
ppc_decay <- pp_check(fit_model, type = "stat_grouped", group = "stimulus0",
                      stat = "median", ndraws = 200)
ppc_genotype <- pp_check(fit_model, type = "stat_grouped", group = "Genotype",
                         stat = "mean", ndraws = 200)
ppc_list <- list(
  ppc_marginal  = ppc_marginal,
  ppc_decay     = ppc_decay,
  ppc_genotype  = ppc_genotype
)
purrr::iwalk(ppc_list, ~ ggsave(
  filename = file.path(diag_dir, paste0(.y, ".png")),
  plot     = .x,
  width    = 8,
  height   = 6,
  dpi      = 300
))


mcmc_pairs_plot <- bayesplot::mcmc_pairs(
  fit_model, np = nuts_params(fit_model),
  pars = c("b_logk_Intercept", "sd_animal__logk_Intercept"),
  off_diag_args = list(size = 0.5)
)
ggsave(
  filename = file.path(diag_dir, "mcmc_pairs.png"),
  plot     = mcmc_pairs_plot,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# Residuals
yrep <- posterior_predict(fit_model, ndraws = 1000)

sim_res <- createDHARMa(
  simulatedResponse = t(yrep),
  observedResponse = df_resp$max_peak,
  fittedPredictedResponse = fitted(fit_model)[, "Estimate"]
)

print(plot(sim_res))
plotResiduals(sim_res, df_resp$stimulus)
plotResiduals(sim_res, df_resp$stimulus0)
plotResiduals(sim_res, df_resp$Genotype)
plotResiduals(sim_res, df_resp$Block)
plotResiduals(sim_res, df_resp$animal)

save_plot_as_png(
  paste0("nlme_", var_name, "_DHARMa_residuals.png"),
  quote(plot(sim_res))
)


# Compute leave-one-out cross-validation:
loo_var <- loo(fit_model)

print(loo_var)
plot(loo_var)
save_plot_as_png(
  paste0("nlme_", var_name, "_loo_plot.png"),
  quote(plot(loo_var))
)
# Identify which observations they are
plot(loo_var, diagnostic = "k")

# Or extract directly
problematic <- which(loo_var$diagnostics$pareto_k > 0.7)
df_resp[problematic, ]

# loo_compare(loo1, loo2)  # compare models


# ==============================================================================
# Plot Habituation curves
# ==============================================================================
p_hab_ind <- plot_habituation_response(
  df_resp, fit_model,
  response_var = "max_peak",
  y_label = "Peak distance moved",
  y_limits = c(0, 10),
  raw_data = "individual",
  ndraws        = 200, 
  save_fig_dir = save_fig_dir
)

p_hab_ind

p_hab <- plot_habituation_response(
  df_resp, fit_model,
  response_var = "max_peak",
  y_label = "Peak distance moved",
  y_limits = c(0, 10),
  raw_data = "trials",
  save_fig_dir = save_fig_dir
)

p_hab

p_hab_agg <- plot_habituation_response(
  df_resp, fit_model,
  response_var = "max_peak",
  y_label = "Peak distance moved",
  y_limits = c(0, 10),
  raw_data = "aggregate",
  save_fig_dir = save_fig_dir
)

p_hab_agg


p_hab_averaged <- plot_habituation_response_animal_averaged(
  df_resp, fit_model,
  response_var = "max_peak",
  ndraws = 500,
  y_label = "Peak distance moved",
  y_limits = c(0, 10),
  save_fig_dir = save_fig_dir
)

p_hab_averaged