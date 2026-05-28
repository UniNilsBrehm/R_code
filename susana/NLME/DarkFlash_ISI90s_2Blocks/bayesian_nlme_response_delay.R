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

source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/nlme_utils.R")
source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/plot_utils.R")

# source("D:/Behavior_Data/R_code/susana/nlme_utils.R")
# source("D:/Behavior_Data/R_code/susana/plot_utils.R")

# base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_ISI90s_2Blocks"
base_dir <- "D:/WorkingData/Susana/NLME/DarkFlash_ISI90s_2Blocks"

file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

var_name = 'response_delay'
col_name = 'delay'

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

df_resp <- df_resp %>%
  mutate(delay = factor(delay, levels = 0:4, ordered = TRUE))

# Check if all values are larger than zero (for Gamma)
summary(df_resp[[col_name]])

# ==============================================================================
# The Model
# ==============================================================================
# In ordinal regression, you don't model the response directly.
# You model a latent continuous variable that gets thresholded into categories. 
# The link function controls how that latent variable maps to category probabilities.

model_mm <- bf(
  delay ~ R0 + (A - R0) * stimulus0 / (exp(logk) + stimulus0),  # Michaelis-Menten

  A    ~ 1 + Genotype * Block + (1 | animal),
  R0   ~ 1 + Genotype * Block + (1 | animal),
  logk ~ 1 + Genotype * Block + (1 | animal),
  nl = TRUE
  
)

model_exp <- bf(
  # delay ~ R0 + (A - R0) * stimulus0 / (exp(logk) + stimulus0),  # Michaelis-Menten
  delay ~ A + (R0 - A) * exp(-exp(logk) * stimulus0),           # Exponential
  
  A    ~ 1 + Genotype * Block + (1 | animal),
  R0   ~ 1 + Genotype * Block + (1 | animal),
  logk ~ 1 + Genotype * Block + (1 | animal),
  nl = TRUE
  
)

# ==============================================================================
# The Priors
# ==============================================================================
priors_mm <- c(
  # Reference-cell intercepts on latent (logit) scale
  prior(normal(-3.0, 0.4), nlpar = "R0", class = "b", coef = "Intercept"),
  prior(normal(-1.5, 0.4), nlpar = "A",  class = "b", coef = "Intercept"),
  prior(normal( 2.0, 0.5), nlpar = "logk", class = "b", coef = "Intercept"), 
  
  # Deviations from reference
  prior(normal(0, 0.3), nlpar = "R0",   class = "b"),
  prior(normal(0, 0.3), nlpar = "A",    class = "b"),
  prior(normal(0, 0.3), nlpar = "logk", class = "b"),
  
  # Animal random effects
  prior(student_t(3, 0, 0.5), nlpar = "R0",   class = "sd"),
  prior(student_t(3, 0, 0.5), nlpar = "A",    class = "sd"),
  prior(student_t(3, 0, 0.3), nlpar = "logk", class = "sd")
)

# Exponential Model Priors
priors_exp <- c(
  # Reference-cell intercepts on latent (logit) scale
  prior(normal(-3.0, 0.4), nlpar = "R0",   class = "b", coef = "Intercept"),
  prior(normal(-1.5, 0.4), nlpar = "A",    class = "b", coef = "Intercept"),
  prior(normal(-2.5, 0.5), nlpar = "logk", class = "b", coef = "Intercept"),
  
  # Deviations from reference
  prior(normal(0, 0.3), nlpar = "R0",   class = "b"),
  prior(normal(0, 0.3), nlpar = "A",    class = "b"),
  prior(normal(0, 0.3), nlpar = "logk", class = "b"),
  
  # Animal random effects
  prior(student_t(3, 0, 0.5), nlpar = "R0",   class = "sd"),
  prior(student_t(3, 0, 0.5), nlpar = "A",    class = "sd"),
  prior(student_t(3, 0, 0.3), nlpar = "logk", class = "sd")
)

# ==============================================================================
# Validate Priors
# ==============================================================================
fit_delay_prior <- brm(
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

# Predicted category distribution
pp <- posterior_predict(fit_delay_prior, ndraws = 500)
round(prop.table(table(pp)), 3)
round(prop.table(table(df_resp$delay)), 3)

# Predicted latent eta
# ep <- posterior_epred(fit_delay_prior, ndraws = 500)  
# Note: epred returns category probabilities, shape draws × N × 5
# Or use fitted on linear predictor:
fit_lin <- fitted(fit_delay_prior, scale = "linear", summary = TRUE)
summary(fit_lin[, "Estimate"])

# ==============================================================================
# Fit Fast Test Model
# ==============================================================================
# Fit using Meanfield Variational Inference (VI)
fit_vi_test <- brm(
  formula = model,
  data = df_resp,               
  family = cumulative(link = "logit"),
  prior = priors,
  backend = "cmdstanr",
  algorithm = "meanfield",          # <--- Reliable, fast VI method
  # algorithm = "fullrank",
  iter = 10000                      # VI likes higher iterations (it's still lightning fast)
)

fit_nuts_test <- brm(
  formula = model,
  data = df_resp,
  family = cumulative(link = "logit"),
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
# ---- Setup ------------------------------------------------------------------
models_dir <- file.path(base_dir, "models")
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

file_mm  <- file.path(models_dir, paste0("bayesian_nlme_", var_name, "_mm_results.rds"))
file_exp <- file.path(models_dir, paste0("bayesian_nlme_", var_name, "_exp_results.rds"))

# ---- MM fit -----------------------------------------------------------------
message("Starting MM fit at ", Sys.time())
t0 <- Sys.time()

fit_model_mm <- brm(
  model_mm,
  data      = df_resp,
  family    = cumulative(link = "logit"),
  prior     = priors_mm,
  save_pars = save_pars(all = TRUE),
  chains    = 4,
  cores     = 4,
  threads   = threading(6),
  iter      = 4000,
  warmup    = 2000,
  control   = list(adapt_delta = 0.99, max_treedepth = 15),
  init      = 0,
  backend   = "cmdstanr"
)

# Save IMMEDIATELY after first fit
saveRDS(fit_model_mm, file = file_mm)
message("MM fit done and saved at ", Sys.time(), " (took ", 
        round(difftime(Sys.time(), t0, units = "mins"), 1), " min)")

# ---- Exponential fit --------------------------------------------------------
message("Starting Exp fit at ", Sys.time())
t1 <- Sys.time()

fit_model_exp <- brm(
  model_exp,
  data      = df_resp,
  family    = cumulative(link = "logit"),
  prior     = priors_exp,
  save_pars = save_pars(all = TRUE),
  chains    = 4,
  cores     = 4,
  threads   = threading(6),
  iter      = 4000,
  warmup    = 2000,
  control   = list(adapt_delta = 0.99, max_treedepth = 15),
  init      = 0,
  backend   = "cmdstanr"
)

saveRDS(fit_model_exp, file = file_exp)
message("Exp fit done and saved at ", Sys.time(), " (took ", 
        round(difftime(Sys.time(), t1, units = "mins"), 1), " min)")

message("All done. Total time: ", 
        round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

# ==============================================================================
# Load fitted model if available
# ==============================================================================
fit_model_exp <- readRDS(
  file.path(base_dir, "models", paste0("bayesian_nlme_", var_name,"_exp_results.rds"))
)
fit_model_mm <- readRDS(
  file.path(base_dir, "models", paste0("bayesian_nlme_", var_name,"_mm_results.rds"))
)

# ==============================================================================
# Plot Habituation curves
# ==============================================================================
fit_model <- fit_model_exp

# Population-level
p_hab_ind <- plot_habituation_ordinal(
  df_resp, fit_model,
  response_var = "delay",
  raw_data = "individual",
  y_label = "Delay (Exp model)",
  y_limits = c(0, 4),
  save_fig_dir = save_fig_dir,
  ndraws = 200
)

print(p_hab_ind)

# Population-level
p_hab <- plot_habituation_ordinal(
  df_resp, fit_model,
  response_var = "delay",
  raw_data = "trials",
  y_label = "Delay (Exp model)",
  y_limits = c(0, 4),
  save_fig_dir = save_fig_dir
)

print(p_hab)

p_hab_agg <- plot_habituation_ordinal(
  df_resp, fit_model,
  response_var = "delay",
  raw_data = "aggregate",
  y_label = "Delay (Exp model)",
  y_limits = c(0, 4),
  save_fig_dir = save_fig_dir
)

print(p_hab_agg)

# Animal-averaged
p_hab_averaged <-plot_habituation_ordinal_animal_averaged(
  df_resp, fit_model,
  response_var = "delay",
  ndraws = 500,
  y_limits = c(0, 4),
  y_breaks = 0:4,
  save_fig_dir = save_fig_dir
)

print(p_hab_averaged)

# ==============================================================================
# COMPARE MODELS
loo_mm  <- loo(fit_model_mm)
loo_exp <- loo(fit_model_exp)

print(summary(fit_model_exp))
print(summary(fit_model_mm))
print(loo_mm)
print(loo_exp)
loo_compare(loo_mm, loo_exp)

# ==============================================================================
# Habituation Curves Plot
# Population-level
p_hab_exp_ind <- plot_habituation_ordinal(
  df_resp, fit_model_exp,
  response_var = "delay",
  raw_data = "individual",
  y_label = "Delay (Exp model)",
  y_limits = c(0, 4)
)


p_hab_exp <- plot_habituation_ordinal(
  df_resp, fit_model_exp,
  response_var = "delay",
  raw_data = "trials",
  y_label = "Delay (Exp model)",
  y_limits = c(0, 4)
)
p_hab_mm <- plot_habituation_ordinal(
  df_resp, fit_model_mm,
  response_var = "delay",
  raw_data = "trials",
  y_label = "Delay (MM model)",
  y_limits = c(0, 4)
)


print(p_hab_exp)
print(p_hab_mm)


p_hab_exp <- plot_habituation_ordinal(
  df_resp, fit_model_exp,
  response_var = "delay",
  raw_data = "trials",
  y_label = "Delay (Exp model)",
  y_limits = c(0, 4)
)

p_hab_exp_agg <- plot_habituation_ordinal(
  df_resp, fit_model_exp,
  response_var = "delay",
  raw_data = "aggregate",
  y_label = "Delay (Exp model)",
  y_limits = c(0, 4)
)
print(p_hab_exp)
print(p_hab_exp_agg)

# ==============================================================================
# Compute Comparisons of parameters (for exp model)
# ==============================================================================
# ==============================================================================
# Compute hybrid-scale contrasts for ordinal delay model:
#   - R0: expected delay in seconds at stimulus = stimulus_for_R0
#   - A:  expected delay in seconds at stimulus = stimulus_for_A
#   - logk: rate ratio vs reference (multiplicative)
# ==============================================================================
compute_contrasts_hybrid <- function(
    fit_model,
    df_resp,
    reference_level = "ABTL",
    category_values = 0:4,
    stimulus_for_R0 = 1,
    stimulus_for_A  = 60,
    ndraws          = NULL
) {
  
  genotypes <- levels(df_resp$Genotype)
  blocks    <- levels(df_resp$Block)
  
  # --------------------------------------------------------------------------
  # Part 1: response-scale contrasts for R0 and A
  # --------------------------------------------------------------------------
  
  # Build a prediction grid for each (Genotype × Block) at two stimulus points:
  # one near the start of the block (R0-dominated) and one at the end (A-dominated).
  grid <- expand.grid(
    Genotype = genotypes,
    Block    = blocks,
    stimulus = c(stimulus_for_R0, stimulus_for_A),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      stimulus0 = stimulus - 1,
      Genotype  = factor(Genotype, levels = genotypes),
      Block     = factor(Block,    levels = blocks),
      which_par = ifelse(stimulus == stimulus_for_R0, "R0", "A")
    )
  
  # Posterior epred at every grid point (returns array [draws × N × K])
  epred_args <- list(
    object     = fit_model,
    newdata    = grid,
    re_formula = NA,
    summary    = FALSE
  )
  if (!is.null(ndraws)) epred_args$ndraws <- ndraws
  
  ep_array <- do.call(posterior_epred, epred_args)
  # ep_array dimensions: [draws, N_rows, K_categories]
  
  # Collapse categories to expected delay: weighted sum
  exp_delay <- array(0, dim = dim(ep_array)[1:2])
  for (k in seq_along(category_values)) {
    exp_delay <- exp_delay + ep_array[, , k] * category_values[k]
  }
  # exp_delay: [draws × N_rows], each value is "expected delay in seconds"
  
  # Build a lookup: column index -> (Genotype, Block, which_par)
  grid_idx <- grid %>%
    mutate(.col = row_number())
  
  # For each (Block, which_par), compute contrasts: each Genotype vs reference
  results_seconds <- list()
  
  for (par in c("R0", "A")) {
    for (blk in blocks) {
      ref_idx <- grid_idx %>%
        filter(which_par == par, Block == blk, Genotype == reference_level) %>%
        pull(.col)
      
      ref_draws <- exp_delay[, ref_idx]
      
      for (geno in setdiff(genotypes, reference_level)) {
        comp_idx <- grid_idx %>%
          filter(which_par == par, Block == blk, Genotype == geno) %>%
          pull(.col)
        
        diff_draws <- exp_delay[, comp_idx] - ref_draws
        
        results_seconds[[length(results_seconds) + 1]] <- tibble(
          parameter = par,
          Block     = blk,
          contrast  = paste0(geno, " - ", reference_level),
          estimate  = mean(diff_draws),
          lower     = quantile(diff_draws, 0.025, names = FALSE),
          upper     = quantile(diff_draws, 0.975, names = FALSE),
          pd        = max(mean(diff_draws > 0), mean(diff_draws < 0)),
          scale     = "seconds"
        )
      }
    }
  }
  
  seconds_tbl <- bind_rows(results_seconds) %>%
    mutate(Block = factor(Block, levels = blocks))
  
  # --------------------------------------------------------------------------
  # Part 2: rate-ratio contrasts for logk
  # --------------------------------------------------------------------------
  
  emm_logk <- emmeans(fit_model, ~ Genotype * Block, nlpar = "logk", re_formula = NA)
  
  logk_ctr <- contrast(emm_logk, method = "trt.vs.ctrl", ref = reference_level, by = "Block") %>%
    gather_emmeans_draws() %>%
    group_by(contrast, Block) %>%
    summarise(
      # Exponentiate the draws to get rate ratios, then summarize
      estimate = mean(exp(.value)),
      lower    = quantile(exp(.value), 0.025, names = FALSE),
      upper    = quantile(exp(.value), 0.975, names = FALSE),
      pd       = max(mean(.value > 0), mean(.value < 0)),
      .groups  = "drop"
    ) %>%
    mutate(
      parameter = "logk",
      scale     = "rate_ratio",
      Block     = factor(Block, levels = blocks)
    ) %>%
    select(parameter, Block, contrast, estimate, lower, upper, pd, scale)
  
  # --------------------------------------------------------------------------
  # Combine
  # --------------------------------------------------------------------------
  bind_rows(seconds_tbl, logk_ctr)
}


# ==============================================================================
# Forest plot for hybrid-scale contrasts
# ==============================================================================
plot_contrasts_forest_hybrid <- function(
    contrasts_tbl,
    pd_threshold     = 0.95,
    genotype_order   = NULL,
    show_pd_labels   = TRUE,
    point_size       = 3,
    line_size        = 0.7,
    seconds_unit     = "s"   # label for the seconds unit (e.g., "s", "sec", "delay units")
) {
  
  # ---- Parameter labels with biology -----------------------------------------
  contrasts_tbl <- contrasts_tbl %>%
    mutate(
      parameter_label = case_when(
        parameter == "R0"   ~ paste0("R0 — initial delay\n(", seconds_unit, " vs ABTL)"),
        parameter == "A"    ~ paste0("A — asymptotic delay\n(", seconds_unit, " vs ABTL)"),
        parameter == "logk" ~ "logk — rate ratio\n(× vs ABTL; >1 = faster)",
        TRUE ~ parameter
      ),
      parameter_label = factor(
        parameter_label,
        levels = c(
          paste0("R0 — initial delay\n(", seconds_unit, " vs ABTL)"),
          paste0("A — asymptotic delay\n(", seconds_unit, " vs ABTL)"),
          "logk — rate ratio\n(× vs ABTL; >1 = faster)"
        )
      ),
      contrast_clean = stringr::str_remove_all(contrast, "\\s*-\\s*ABTL.*$") %>%
        stringr::str_remove_all("[()]") %>%
        stringr::str_trim(),
      is_significant = pd >= pd_threshold,
      # pd label formatting
      pd_label = ifelse(pd >= 0.999, "pd > 0.999", sprintf("pd = %.2f", pd)),
      # Display value labels (formatted differently for seconds vs ratio)
      value_label = case_when(
        parameter == "logk" ~ sprintf("%.2fx", estimate),
        TRUE                ~ sprintf("%+.2f", estimate)
      ),
      # Reference line (0 for seconds, 1 for rate ratio)
      ref_line = ifelse(parameter == "logk", 1, 0)
    )
  
  # ---- Genotype ordering -----------------------------------------------------
  if (is.null(genotype_order)) {
    genotype_order <- c("tyr", "th, tyr", "th2, tyr", "th, th2, tyr")
  }
  
  contrasts_tbl <- contrasts_tbl %>%
    filter(contrast_clean %in% genotype_order) %>%
    mutate(contrast_clean = factor(contrast_clean, levels = rev(genotype_order)))
  
  # ---- Reference lines: data frame for facet-specific vertical lines --------
  ref_lines <- contrasts_tbl %>%
    distinct(parameter_label, ref_line)
  
  # ---- Build plot ------------------------------------------------------------
  p <- ggplot(
    contrasts_tbl,
    aes(x = estimate, y = contrast_clean, color = Block, group = Block)
  ) +
    geom_vline(
      data = ref_lines,
      aes(xintercept = ref_line),
      linetype = "dashed",
      color = "grey50",
      inherit.aes = FALSE
    ) +
    geom_linerange(
      aes(xmin = lower, xmax = upper),
      position = position_dodge(width = 0.6, reverse = TRUE),
      linewidth = line_size,
      alpha = 0.85
    ) +
    geom_point(
      aes(shape = is_significant),
      position = position_dodge(width = 0.6, reverse = TRUE),
      size = point_size,
      stroke = 1.2
    ) +
    scale_shape_manual(
      name   = paste0("pd ≥ ", pd_threshold),
      values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("TRUE" = "Yes", "FALSE" = "No")
    )
  
  if (show_pd_labels) {
    p <- p +
      geom_text(
        aes(label = pd_label, x = upper),
        position = position_dodge(width = 0.6, reverse = TRUE),
        hjust = -0.15,
        size = 3,
        color = "grey30",
        show.legend = FALSE
      )
  }
  
  p <- p +
    facet_wrap(~ parameter_label, scales = "free_x", ncol = 3) +
    scale_color_manual(values = c("Block1" = "#E69F00", "Block2" = "#56B4E9")) +
    labs(
      x        = NULL,
      y        = NULL,
      title    = "Genotype contrasts vs ABTL across habituation parameters",
      subtitle = sprintf("Solid points: pd ≥ %.2f (high posterior certainty); open points: lower certainty",
                         pd_threshold),
      caption  = "th = tyrosine hydroxylase   |   th2 = tyrosine hydroxylase 2   |   tyr = tyrosinase"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      legend.box      = "horizontal",
      panel.spacing.x = unit(1.5, "lines"),
      plot.caption    = element_text(hjust = 0, color = "grey30",
                                     size = 9, face = "italic"),
      plot.subtitle   = element_text(color = "grey30", size = 10),
      strip.text      = element_text(face = "bold", size = 11),
      axis.text.y     = element_text(face = "italic"),
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor   = element_blank()
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.25)))
  
  p
}

library(tidybayes)

# ==============================================================================
# Posterior densities of nonlinear parameters on response scale
# ==============================================================================
# - For delay model: R0/A as expected seconds, logk as rate per stimulus
# - For distance models: R0/A as exp(R0)/exp(A), logk as exp(logk)

extract_param_posteriors_delay <- function(
    fit_model,
    df_resp,
    category_values = 0:4,
    stimulus_for_R0 = 1,
    stimulus_for_A  = 60,
    ndraws          = NULL
) {
  
  genotypes <- levels(df_resp$Genotype)
  blocks    <- levels(df_resp$Block)
  
  # ---- R0 and A in seconds ---------------------------------------------------
  grid <- expand.grid(
    Genotype = genotypes,
    Block    = blocks,
    stimulus = c(stimulus_for_R0, stimulus_for_A),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      stimulus0 = stimulus - 1,
      Genotype  = factor(Genotype, levels = genotypes),
      Block     = factor(Block,    levels = blocks),
      which_par = ifelse(stimulus == stimulus_for_R0, "R0", "A")
    )
  
  epred_args <- list(object = fit_model, newdata = grid,
                     re_formula = NA, summary = FALSE)
  if (!is.null(ndraws)) epred_args$ndraws <- ndraws
  
  ep <- do.call(posterior_epred, epred_args)
  exp_delay <- array(0, dim = dim(ep)[1:2])
  for (k in seq_along(category_values)) {
    exp_delay <- exp_delay + ep[, , k] * category_values[k]
  }
  
  draws_R0_A <- as.data.frame(exp_delay) %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(-.draw, names_to = ".col", values_to = "value") %>%
    mutate(.col = as.integer(gsub("V", "", .col))) %>%
    left_join(grid %>% mutate(.col = row_number()), by = ".col") %>%
    rename(parameter = which_par) %>%
    select(.draw, parameter, Genotype, Block, value)
  
  # ---- logk as rate per stimulus --------------------------------------------
  emm_logk <- emmeans(fit_model, ~ Genotype * Block,
                      nlpar = "logk", re_formula = NA)
  
  draws_logk <- emm_logk %>%
    gather_emmeans_draws() %>%
    mutate(
      parameter = "logk",
      value     = exp(.value)
    ) %>%
    select(.draw, parameter, Genotype, Block, value)
  
  bind_rows(draws_R0_A, draws_logk) %>%
    mutate(
      parameter = factor(parameter, levels = c("R0", "A", "logk")),
      Block     = factor(Block,     levels = levels(df_resp$Block)),
      Genotype  = factor(Genotype,  levels = levels(df_resp$Genotype))
    )
}


# ==============================================================================
# Plot: overlaid posterior densities per Block × Parameter
# ==============================================================================
plot_posterior_densities <- function(
    posterior_df,
    parameter_labels  = NULL,
    genotype_order    = NULL,
    reference_level   = "ABTL",
    palette           = NULL
) {
  
  if (is.null(parameter_labels)) {
    parameter_labels <- c(
      R0   = "R0 — initial delay (s)",
      A    = "A — asymptotic delay (s)",
      logk = "logk — rate per stimulus"
    )
  }
  
  if (is.null(genotype_order)) {
    genotype_order <- c("ABTL", "tyr", "th, tyr", "th2, tyr", "th, th2, tyr")
  }
  
  posterior_df <- posterior_df %>%
    mutate(
      Genotype = factor(Genotype, levels = genotype_order),
      parameter_label = factor(
        parameter_labels[as.character(parameter)],
        levels = parameter_labels
      )
    )
  
  # Default colorblind-safe palette
  if (is.null(palette)) {
    palette <- c(
      "ABTL"          = "#000000",
      "tyr"           = "#E69F00",
      "th, tyr"       = "#56B4E9",
      "th2, tyr"      = "#009E73",
      "th, th2, tyr"  = "#CC79A7"
    )
  }
  
  ggplot(posterior_df, aes(x = value, color = Genotype, fill = Genotype)) +
    geom_density(alpha = 0.2, linewidth = 0.8) +
    facet_grid(Block ~ parameter_label, scales = "free", switch = "y") +
    scale_color_manual(values = palette) +
    scale_fill_manual(values  = palette) +
    labs(
      x = "Parameter value",
      y = "Posterior density",
      title = "Posterior distributions of habituation parameters by genotype",
      caption = "th = tyrosine hydroxylase   |   th2 = tyrosine hydroxylase 2   |   tyr = tyrosinase"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.spacing   = unit(1.5, "lines"),
      strip.text      = element_text(face = "bold", size = 11),
      strip.placement = "outside",
      plot.caption    = element_text(hjust = 0, color = "grey30",
                                     size = 9, face = "italic"),
      panel.grid.minor = element_blank(),
      legend.text     = element_text(face = "italic")
    )
}

##
# Compute hybrid contrasts
contrasts_hybrid <- compute_contrasts_hybrid(
  fit_model_exp,
  df_resp,
  reference_level = "ABTL",
  category_values = 0:4,
  stimulus_for_R0 = 1,
  stimulus_for_A  = 60,
  ndraws          = 2000      # subsample for speed; use NULL for all draws
)

print(contrasts_hybrid, n = Inf)

# Plot
p_forest_hybrid <- plot_contrasts_forest_hybrid(
  contrasts_hybrid,
  pd_threshold = 0.95
)
print(p_forest_hybrid)

post_df <- extract_param_posteriors_delay(
  fit_model_exp,
  df_resp,
  category_values = 0:4,
  stimulus_for_R0 = 1,
  stimulus_for_A  = 60,
  ndraws          = 2000
)

p_post <- plot_posterior_densities(post_df)
print(p_post)



# -------------------------------------------------------------------------
# pd      Frequentist analog   Verbal phrasing
# -------------------------------------------------------------------------
# 0.50    p ≈ 1.0              no evidence about direction
# 0.75    p ≈ 0.50             weakly suggestive
# 0.85    p ≈ 0.30             leans toward an effect
# 0.90    p ≈ 0.20             suggestive but not conclusive
# 0.95    p ≈ 0.10             moderate evidence
# 0.975   p ≈ 0.05             credibly different (standard threshold)
# 0.99    p ≈ 0.02             strong evidence
# 0.999   p ≈ 0.002            very strong evidence
# -------------------------------------------------------------------------