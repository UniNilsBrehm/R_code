# ==============================================================================
# Bayesian nonlinear mixed model for zebrafish MASSED TRAINING habituation
# ------------------------------------------------------------------------------
# Mirrors the structure of the DarkFlash ISI60 NLME script, with priors
# adjusted for the much longer Block 1 (479 stimuli vs 9 in dark-flash).
#
# Model (per animal):
#   p(move) = inv_logit(A) + (inv_logit(R0) - inv_logit(A)) * exp(-exp(logk)*stim0)
#
#   A    : asymptotic response probability (logit scale)
#   R0   : initial response probability    (logit scale)
#   logk : log habituation rate            (per stimulus)
#
# All three nonlinear parameters get a full Genotype * Block fixed effect
# and a per-animal random intercept (1 | animal).
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
library(purrr)
library(broom)


# ==============================================================================
# HELPER
# ==============================================================================
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
# Paths
# ==============================================================================
source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/DarkFlash_Massed/utils.R")
source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/DarkFlash_Massed/plot_utils.R")
base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_Massed/"

var_name <- "response_prob"
col_name <- "move"

save_fig_dir     <- file.path(base_dir, "figs",    "nlme", var_name)
save_results_dir <- file.path(base_dir, "results", "nlme", var_name)
save_model_dir   <- file.path(base_dir, "models")

dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_model_dir,   recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# Load and prepare data
# ==============================================================================
file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_Massed_Training_7Nov2025.csv"
)


res <- load_data(
  # "D:/WorkingData/Susana/SPZ_Massed_Training_7Nov2025.csv",
  file_dir,
  move_th = 1,
  # keep = c("ABTL", "th, th2, tyr", "th, tyr"),
  # keep = c("ABTL", "th, th2, tyr")
)

df_final     <- res$df_final
df_final_sub <- res$df_final_sub

df_resp <- df_final %>%
  mutate(
    stimulus  = as.numeric(stimulus),
    stimulus0 = stimulus + 1,
    Block     = factor(Block),
    Genotype  = factor(Genotype),
    Video     = factor(Video),
    Well      = factor(Well),
    animal    = interaction(Video, Well, drop = TRUE)
  )

# Sanity check: stimulus range per block (Block 1 should hit ~478, Block 2 ~9)
print(
  df_resp %>%
    group_by(Block) %>%
    summarise(
      n_stim_unique = n_distinct(stimulus0),
      min_stim      = min(stimulus0),
      max_stim      = max(stimulus0),
      n_animals     = n_distinct(animal),
      .groups = "drop"
    )
)

print(
  df_resp %>%
    distinct(animal, Genotype) %>%
    count(Genotype)
)

cat("Observed marginal P(move):", round(mean(df_resp$move), 3), "\n")
# Should print ~0.69

ggplot(df_resp, aes(x = move)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Distribution of Binary Movement Responses (Massed)",
    x = "Response (0 = No, 1 = Yes)",
    y = "Count"
  )

ggplot(df_resp, aes(x = max_peak)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Distribution of Binary Movement Responses (Massed)",
    x = "Response (0 = No, 1 = Yes)",
    y = "Count"
  )


# ==============================================================================
# The Model
# ==============================================================================
# Identical structure to the dark-flash NLME so all downstream helpers
# (plot_habituation_probability, make_nlpar_draws, compare_nlpar, ...) work.
model <- bf(
  move ~ inv_logit(A) + (inv_logit(R0) - inv_logit(A)) * exp(-exp(logk) * stimulus0),
  
  A    ~ 1 + Genotype * Block + (1 | animal),
  R0   ~ 1 + Genotype * Block + (1 | animal),
  logk ~ 1 + Genotype * Block + (1 | animal),
  
  nl = TRUE
)


# ==============================================================================
# The Priors  -- ADJUSTED for the long Block 1 (479 stimuli)
# ==============================================================================
# Dark-flash used normal(-3, 1) on logk because blocks were 9 stimuli long.
# Here Block 1 is ~479 stimuli, so habituation rates are MUCH smaller.
# Reasonable half-lives are ~30-200 stimuli -> k in [0.0035, 0.023]
#                                          -> logk in [-5.6, -3.8]
# Centering logk at -4.5 with sd=1 covers half-lives ~10 to ~500 stimuli,
# which is wide enough to also accommodate the shorter Block 2.
#
# A and R0 priors are unchanged from the dark-flash script (they live on
# the logit probability scale; the stimulus axis does not affect them).
priors <- c(
  prior(normal(0,    1),   class = "b", coef = "Intercept", nlpar = "A"),
  prior(normal(1.5,  1),   class = "b", coef = "Intercept", nlpar = "R0"),
  prior(normal(-4.5, 1),   class = "b", coef = "Intercept", nlpar = "logk"),
  
  prior(normal(0, 0.75), class = "b", nlpar = "A"),
  prior(normal(0, 0.75), class = "b", nlpar = "R0"),
  prior(normal(0, 0.75), class = "b", nlpar = "logk"),
  
  prior(exponential(4), class = "sd", nlpar = "A"),
  prior(exponential(4), class = "sd", nlpar = "R0"),
  prior(exponential(4), class = "sd", nlpar = "logk")
)


# ==============================================================================
# Validate Priors (prior predictive check)
# ==============================================================================
fit_prior <- brm(
  formula      = model,
  data         = df_resp,
  family       = bernoulli(link = "identity"),
  prior        = priors,
  sample_prior = "only",
  backend      = "cmdstanr",
  chains       = 4,
  cores        = 4,
  iter         = 2000,
  warmup       = 500,
  seed         = 42
)


diag_dir <- file.path(save_model_dir, "diagnostics", "nlme", "priors", var_name)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

ep <- posterior_epred(fit_prior, ndraws = 500)

cat("\nPrior predictive quantiles (expected p(move)):\n")
print(quantile(as.vector(ep), c(.001, .01, .05, .5, .95, .99, .999), na.rm = TRUE))

cat("\nObserved data quantiles:\n")
print(quantile(df_resp[[col_name]], c(.001, .01, .05, .5, .95, .99, .999), na.rm = TRUE))

# look at predicted trajectories
prior_curves <- df_resp %>%
  distinct(Genotype, Block) %>%
  group_by(Genotype, Block) %>%
  reframe(stimulus0 = seq(0, ifelse(Block == "1", 477, 8), length.out = 100)) %>%
  ungroup() %>%
  mutate(
    stimulus = stimulus0,
    animal   = factor(NA, levels = levels(df_resp$animal))
  ) %>%
  add_epred_draws(fit_prior, ndraws = 50, re_formula = NA)

ggplot(prior_curves, aes(x = stimulus0, y = .epred, group = interaction(.draw, Genotype, Block))) +
  geom_line(alpha = 0.2) +
  facet_grid(Block ~ Genotype, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr() +
  labs(title = "Prior predictive habituation curves",
       y = "P(move)", x = "Stimulus within block")

# ==============================================================================
# Fast test fits (variational + short NUTS) for iteration
# ==============================================================================
fit_vi_test <- brm(
  formula   = model,
  data      = df_resp,
  family    = bernoulli(link = "identity"),
  prior     = priors,
  backend   = "cmdstanr",
  algorithm = "meanfield",
  iter      = 10000
)

fit_vi_fullrank <- brm(
  formula   = model,
  data      = df_resp,
  family    = bernoulli(link = "identity"),
  prior     = priors,
  backend   = "cmdstanr",
  algorithm = "fullrank",
  iter      = 10000,
  seed      = 42
)

fit_nuts_test <- brm(
  formula   = model,
  data      = df_resp,
  family    = bernoulli(link = "identity"),
  prior     = priors,
  save_pars = save_pars(all = TRUE),
  init      = 0,
  backend   = "cmdstanr",
  chains    = 1,
  cores     = 1,
  threads   = threading(12),
  iter      = 600,
  warmup    = 300,
  seed      = 42,
  control   = list(
    adapt_delta   = 0.90,
    max_treedepth = 10
  )
)

fit_model <- fit_vi_test
fit_model <- fit_nuts_test

# ==============================================================================
# Fit the full model
# ==============================================================================
# NOTE: with 52 animals and ~489 stimuli/animal the dataset is bigger than
# the dark-flash one. Expect longer compile + sampling times. If divergences
# appear, bump adapt_delta to 0.995.
fit_model <- brm(
  formula   = model,
  data      = df_resp,
  family    = bernoulli(link = "identity"),
  prior     = priors,
  save_pars = save_pars(all = TRUE),
  backend   = "cmdstanr",
  chains    = 4,
  cores     = 4,
  threads   = threading(6),
  iter      = 6000,
  warmup    = 3000,
  seed      = 42,
  control   = list(adapt_delta = 0.99, max_treedepth = 15),
  init      = 0,
  file      = file.path(save_model_dir, paste0("bayesian_nlme_massed_", var_name, "_results.rds")),
  file_refit = "on_change"
)


# ==============================================================================
# Save / load fitted model
# ==============================================================================
saveRDS(
  fit_model,
  file = file.path(save_model_dir, paste0("bayesian_nlme_massed_", var_name, "_results.rds"))
)

# fit_model <- readRDS(
#   file.path(save_model_dir, paste0("bayesian_nlme_massed_", var_name, "_results.rds"))
# )


# ==============================================================================
# Summary and diagnostics
# ==============================================================================
sink(file.path(save_results_dir, paste0(var_name, "_model_summary.txt")))
fit_summary <- summary(fit_model)
print(fit_summary)
sink()

# Divergences (0 ideal, single digits tolerable)
cat("\nN divergent transitions:\n")
print(sum(subset(nuts_params(fit_model), Parameter == "divergent__")$Value))

# Energy diagnostics
bayesplot::mcmc_nuts_energy(nuts_params(fit_model))

# Posterior predictive checks
pp_check(fit_model, type = "dens_overlay", ndraws = 100)
pp_check(fit_model, type = "stat_grouped", group = "stimulus0",
         stat = "median", ndraws = 100)

bayesplot::mcmc_pairs(
  fit_model, np = nuts_params(fit_model),
  off_diag_args = list(size = 0.5)
)


# ------------------------------------------------------------------------------
# Conditional residuals (with random effects)
# ------------------------------------------------------------------------------
yrep    <- posterior_predict(fit_model, ndraws = 1000)
fit_mu  <- fitted(fit_model)[, "Estimate"]

sim_res <- createDHARMa(
  simulatedResponse       = t(yrep),
  observedResponse        = df_resp$move,
  fittedPredictedResponse = fit_mu,
  integerResponse         = TRUE
)

plot(sim_res)
plotResiduals(sim_res, df_resp$stimulus)
plotResiduals(sim_res, df_resp$stimulus0)
plotResiduals(sim_res, df_resp$Genotype)
plotResiduals(sim_res, df_resp$Block)
plotResiduals(sim_res, df_resp$animal)


# ------------------------------------------------------------------------------
# Population-level residuals
# ------------------------------------------------------------------------------
yrep_pop <- posterior_predict(fit_model, ndraws = 500, re_formula = NA)
fit_mu_pop <- fitted(fit_model, re_formula = NA)[, "Estimate"]

sim_res_pop <- createDHARMa(
  simulatedResponse       = t(yrep_pop),
  observedResponse        = df_resp$move,
  fittedPredictedResponse = fit_mu_pop,
  integerResponse         = TRUE
)

plot(sim_res_pop)


# ------------------------------------------------------------------------------
# LOO cross-validation
# ------------------------------------------------------------------------------
loo_var <- loo(fit_model, moment_match = TRUE)
print(loo_var)

save_plot_as_png(
  paste0("nlme_massed_", var_name, "_loo_plot.png"),
  quote(plot(loo_var))
)


# ------------------------------------------------------------------------------
# Parameter trace / density diagnostic plots
# ------------------------------------------------------------------------------
diag_plots <- plot(fit_model)
out_dir <- file.path(save_model_dir, "diagnostics", "nlme", var_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(diag_plots)) {
  ggsave(
    filename = file.path(out_dir, paste0("fit_model_diagnostic_", i, ".png")),
    plot     = diag_plots[[i]],
    width    = 10,
    height   = 7,
    dpi      = 300,
    bg       = "white"
  )
}


# ==============================================================================
# Plot habituation curves
# ==============================================================================
p_hab <- plot_habituation_probability(
  df_resp, fit_model,
  raw_data     = "aggregate",
  stim0_offset = 0,
  facet_layout = "genotype_rows"
)

print(p_hab)


# ==============================================================================
# Simple aggregate response-probability exponential fits
# (sanity check + prior derivation, unchanged in spirit from dark-flash script)
# ==============================================================================
df_prob_agg <- df_resp %>%
  mutate(
    stimulus  = as.numeric(stimulus),
    stimulus0 = stimulus,
    Genotype  = factor(Genotype),
    Block     = factor(Block),
    move      = as.integer(move)
  ) %>%
  group_by(Genotype, Block, stimulus, stimulus0) %>%
  summarise(
    n_animals     = n(),
    n_response    = sum(move, na.rm = TRUE),
    response_prob = mean(move, na.rm = TRUE),
    .groups = "drop"
  )

print(head(df_prob_agg, 20))


fit_simple_exp <- function(dat) {
  tryCatch(
    {
      nls(
        response_prob ~ A + (R0 - A) * exp(-exp(logk) * stimulus0),
        data    = dat,
        weights = n_animals,
        start = list(
          A    = max(0.01, min(dat$response_prob, na.rm = TRUE)),
          R0   = min(0.99, max(dat$response_prob, na.rm = TRUE)),
          # start logk much smaller now (Block 1 is long)
          logk = log(0.01)
        ),
        algorithm = "port",
        lower = c(A = 0.001, R0 = 0.001, logk = log(1e-5)),
        upper = c(A = 0.999, R0 = 0.999, logk = log(10)),
        control = nls.control(maxiter = 500, warnOnly = TRUE)
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

simple_exp_params <- map2_dfr(
  simple_fits,
  seq_along(simple_fits),
  function(mod, i) {
    if (is.null(mod)) {
      tibble(
        Genotype = fit_keys$Genotype[i],
        Block    = fit_keys$Block[i],
        A = NA_real_, R0 = NA_real_, logk = NA_real_, k = NA_real_
      )
    } else {
      cc <- coef(mod)
      tibble(
        Genotype = fit_keys$Genotype[i],
        Block    = fit_keys$Block[i],
        A    = unname(cc["A"]),
        R0   = unname(cc["R0"]),
        logk = unname(cc["logk"]),
        k    = exp(unname(cc["logk"]))
      )
    }
  }
) %>%
  mutate(half_life_stimuli = log(2) / k)

print(simple_exp_params)

write.csv(
  simple_exp_params,
  file = file.path(save_results_dir, "simple_aggregate_exp_response_prob_parameters.csv"),
  row.names = FALSE
)

simple_exp_pred <- df_prob_agg %>%
  group_by(Genotype, Block) %>%
  summarise(
    stim_min = min(stimulus, na.rm = TRUE),
    stim_max = max(stimulus, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  group_by(Genotype, Block) %>%
  reframe(stimulus = seq(stim_min, stim_max, length.out = 200)) %>%
  mutate(stimulus0 = stimulus) %>%
  left_join(simple_exp_params, by = c("Genotype", "Block")) %>%
  mutate(fit = A + (R0 - A) * exp(-exp(logk) * stimulus0))


p_simple_exp <- ggplot(simple_exp_pred, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_grid(Block ~ Genotype, scales = "free_x") +
  geom_point(
    data = df_prob_agg,
    aes(x = stimulus, y = response_prob, color = Genotype),
    alpha = 0.45, size = 1, inherit.aes = FALSE
  ) +
  geom_line(aes(y = fit), linewidth = 1.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  theme_pubr(base_size = 14) +
  labs(
    x     = "Stimulus number within block",
    y     = "Response probability",
    color = "Genotype",
    fill  = "Genotype",
    title = "Simple aggregate exponential fits (massed)"
  ) +
  theme(legend.position = "top", panel.spacing = unit(1.2, "lines"))

print(p_simple_exp)

ggsave(
  filename = file.path(save_fig_dir, "Simple_aggregate_exp_response_prob_curves.png"),
  plot     = p_simple_exp,
  width    = 14, height = 7, dpi = 300, bg = "white"
)
