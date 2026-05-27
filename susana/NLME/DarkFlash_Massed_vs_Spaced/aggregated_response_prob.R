# ==============================================================================
# Aggregated response probability — power law / exponential fit
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(readr)
library(minpack.lm)

# ==============================================================================
# Paths
# ==============================================================================
source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_Massed_vs_Spaced//utils.R")
# source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/NLME/DarkFlash_Massed_vs_Spaced//utils.R")

base_dir <- "D:/WorkingData/Susana/NLME/DarkFlash_Joint_SpacedVsMassed/"
# base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/NLME/DarkFlash_Joint_SpacedVsMassed/"

select_first_responders <- FALSE
th_for_move <- 1

file_massed <- file.path(base_dir, "data_files", "SPZ_Massed_Training_7Nov2025.csv")
file_spaced <- file.path(base_dir, "data_files", "SPZ_Spaced_Training_Nov2025.csv")

save_fig_dir     <- file.path(base_dir, "figs",    "aggregated", "response_prob")
save_results_dir <- file.path(base_dir, "results", "aggregated", "response_prob")
save_model_dir   <- file.path(base_dir, "models")

dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_model_dir,   recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Load and prepare both datasets
# ==============================================================================
res_massed <- load_data(file_massed, move_th = th_for_move, drop = c('th2, tyr', 'tyr'))
res_spaced <- load_data(file_spaced, move_th = th_for_move, drop = c('th2, tyr', 'tyr'))

df_massed <- res_massed$df_final
df_spaced <- res_spaced$df_final

massed_block_test <- "2"
spaced_block_test <- "5"

df_massed_tagged <- df_massed %>%
  mutate(
    Training  = "massed",
    BlockRole = ifelse(as.character(Block) == massed_block_test, "test", "training")
  )

df_spaced_tagged <- df_spaced %>%
  mutate(
    Training  = "spaced",
    BlockRole = ifelse(as.character(Block) == spaced_block_test, "test", "training")
  )

df_all <- bind_rows(df_massed_tagged, df_spaced_tagged) %>%
  mutate(
    stimulus     = as.numeric(stimulus),
    stimulus0    = stimulus - 1,
    stimulus_log = log(stimulus),
    Training     = factor(Training,  levels = c("massed", "spaced")),
    Block        = factor(Block),
    BlockRole    = factor(BlockRole, levels = c("training", "test")),
    Genotype     = factor(Genotype),
    Video        = factor(Video),
    Well         = factor(Well),
    animal       = factor(paste0(Training, "_", Video, ".", Well))
  )

if (select_first_responders) {
  message("Removing all non-responders to first stimulus")
  df_filtered <- df_all %>%
    group_by(animal) %>%
    filter(!any(Block == 1 & stimulus == 1 & move == 0)) %>%
    ungroup()
  df_all <- df_filtered
}

# ==============================================================================
# Switches — change these only
# ==============================================================================
use_log_stimulus  <- TRUE   # TRUE = log(stimulus) = power law, FALSE = raw stimulus
use_asymptote     <- FALSE  # TRUE = add asymptote, FALSE = no asymptote
use_double_exp    <- FALSE  # TRUE = double exponential, FALSE = single
fit_max_stimulus  <- Inf    # Inf = use all; e.g. 100 = fit only stimulus 1:100
extrapolate_curve <- TRUE   # TRUE = draw curve over full range, FALSE = clip to fit range

# ==============================================================================
# Fitting function — response probability version
# Key differences from continuous version:
#   * amplitude upper bound = 1 (probability cannot exceed 1)
#   * asymptote upper bound = 1
#   * predictions clamped to [0, 1]
# ==============================================================================
fit_resp_prob <- function(dat, use_log, use_asym, use_double) {
  
  dat <- dat %>%
    mutate(x = if (use_log) stimulus_log else stimulus) %>%
    filter(is.finite(x), is.finite(p_move), !is.na(p_move)) %>%
    arrange(x)
  
  n_stim <- length(unique(dat$x))
  if (nrow(dat) < 3 || n_stim < 3) {
    return(tibble(
      amplitude_fast = NA_real_, rate_fast = NA_real_,
      amplitude_slow = NA_real_, rate_slow = NA_real_,
      asymptote      = NA_real_,
      fit_type       = "not_fitted",
      converged      = FALSE,
      AIC            = NA_real_, BIC = NA_real_,
      RSS            = NA_real_, n_obs = nrow(dat)
    ))
  }
  
  # --- starting values ---
  lm_init <- tryCatch(
    lm(log(p_move + 1e-6) ~ x, data = dat),
    error = function(e) NULL
  )
  if (!is.null(lm_init)) {
    start_amp  <- min(max(exp(coef(lm_init)[1]), 1e-6), 1)
    start_rate <- max(-coef(lm_init)[2], 1e-6)
  } else {
    start_amp  <- min(max(dat$p_move[which.min(dat$x)], 1e-6), 1)
    start_rate <- 0.5
  }
  start_asym <- max(min(dat$p_move, na.rm = TRUE), 0)
  
  rate_candidates <- unique(c(start_rate, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0))
  
  # --- inner fitter ---
  try_fit <- function(formula_nls, starts, lower_bounds, upper_bounds, is_double) {
    fit <- NULL
    for (rate_start in rate_candidates) {
      starts$rate_fast <- rate_start * ifelse(is_double, 5, 1)
      if (is_double && "rate_slow" %in% names(starts)) {
        starts$rate_slow <- rate_start * 0.2
      }
      fit <- tryCatch(
        withCallingHandlers(
          nlsLM(
            formula_nls,
            data    = dat,
            start   = starts,
            lower   = lower_bounds,
            upper   = upper_bounds,
            control = nls.lm.control(maxiter = 500, ftol = 1e-6, ptol = 1e-6)
          ),
          warning = function(w) invokeRestart("muffleWarning")
        ),
        error = function(e) NULL
      )
      if (!is.null(fit) && any(!is.finite(coef(fit)))) fit <- NULL
      if (!is.null(fit)) break
    }
    fit
  }
  
  fit      <- NULL
  fit_used <- NULL
  
  # --- attempt double exponential ---
  if (use_double) {
    if (use_asym) {
      formula_nls  <- p_move ~ asymptote +
        amplitude_fast * exp(-rate_fast * x) +
        amplitude_slow * exp(-rate_slow * x)
      starts       <- list(amplitude_fast = start_amp * 0.5, rate_fast = start_rate * 5,
                           amplitude_slow = start_amp * 0.5, rate_slow = start_rate * 0.2,
                           asymptote      = start_asym)
      lower_bounds <- c(amplitude_fast = 0, rate_fast = 1e-6,
                        amplitude_slow = 0, rate_slow = 1e-6, asymptote = 0)
      upper_bounds <- c(amplitude_fast = 1, rate_fast = 20,         # capped at 1
                        amplitude_slow = 1, rate_slow = 20, asymptote = 1)
    } else {
      formula_nls  <- p_move ~ amplitude_fast * exp(-rate_fast * x) +
        amplitude_slow * exp(-rate_slow * x)
      starts       <- list(amplitude_fast = start_amp * 0.5, rate_fast = start_rate * 5,
                           amplitude_slow = start_amp * 0.5, rate_slow = start_rate * 0.2)
      lower_bounds <- c(amplitude_fast = 0, rate_fast = 1e-6,
                        amplitude_slow = 0, rate_slow = 1e-6)
      upper_bounds <- c(amplitude_fast = 1, rate_fast = 20,         # capped at 1
                        amplitude_slow = 1, rate_slow = 20)
    }
    fit      <- try_fit(formula_nls, starts, lower_bounds, upper_bounds, is_double = TRUE)
    fit_used <- if (!is.null(fit)) "double" else NULL
  }
  
  # --- fallback to single exponential ---
  if (is.null(fit)) {
    if (use_asym) {
      formula_nls  <- p_move ~ asymptote + amplitude_fast * exp(-rate_fast * x)
      starts       <- list(amplitude_fast = start_amp, rate_fast = start_rate,
                           asymptote = start_asym)
      lower_bounds <- c(amplitude_fast = 0,   rate_fast = 1e-6, asymptote = 0)
      upper_bounds <- c(amplitude_fast = 1,   rate_fast = 20,   asymptote = 1)
    } else {
      formula_nls  <- p_move ~ amplitude_fast * exp(-rate_fast * x)
      starts       <- list(amplitude_fast = start_amp, rate_fast = start_rate)
      lower_bounds <- c(amplitude_fast = 0,   rate_fast = 1e-6)
      upper_bounds <- c(amplitude_fast = 1,   rate_fast = 20)
    }
    fit      <- try_fit(formula_nls, starts, lower_bounds, upper_bounds, is_double = FALSE)
    fit_used <- if (!is.null(fit)) "single_fallback" else NULL
  }
  
  if (is.null(fit)) {
    return(tibble(
      amplitude_fast = NA_real_, rate_fast = NA_real_,
      amplitude_slow = NA_real_, rate_slow = NA_real_,
      asymptote      = NA_real_,
      fit_type       = "failed",
      converged      = FALSE,
      AIC            = NA_real_, BIC = NA_real_,
      RSS            = NA_real_, n_obs = nrow(dat)
    ))
  }
  
  coefs    <- coef(fit)
  get_coef <- function(name) if (name %in% names(coefs)) unname(coefs[name]) else NA_real_
  
  tibble(
    amplitude_fast = get_coef("amplitude_fast"),
    rate_fast      = get_coef("rate_fast"),
    amplitude_slow = get_coef("amplitude_slow"),
    rate_slow      = get_coef("rate_slow"),
    asymptote      = get_coef("asymptote"),
    fit_type       = paste0(
      fit_used, "_",
      if (use_asym)   "asym"    else "no_asym", "_",
      if (use_log)    "logstim" else "rawstim"
    ),
    converged = TRUE,
    AIC       = AIC(fit),
    BIC       = BIC(fit),
    RSS       = sum(residuals(fit)^2),
    n_obs     = nrow(dat)
  )
}

# ==============================================================================
# Aggregate response probability
# ==============================================================================
raw_summary_joint <- df_all %>%
  group_by(Training, Block, BlockRole, Genotype, stimulus, stimulus0, stimulus_log) %>%
  summarise(
    n_trials = sum(!is.na(move)),
    n_move   = sum(move, na.rm = TRUE),
    p_move   = mean(move, na.rm = TRUE),
    .groups  = "drop"
  )

# ==============================================================================
# Fit
# ==============================================================================
exp_params <- raw_summary_joint %>%
  filter(stimulus <= fit_max_stimulus) %>%
  group_by(Training, Block, BlockRole, Genotype) %>%
  group_modify(~ fit_resp_prob(
    dat        = .x,
    use_log    = use_log_stimulus,
    use_asym   = use_asymptote,
    use_double = use_double_exp
  )) %>%
  ungroup()

print(exp_params)
readr::write_csv(
  exp_params,
  file.path(save_results_dir, "exp_params_response_prob.csv")
)

# ==============================================================================
# Prediction grid
# ==============================================================================
new_data_exp <- raw_summary_joint %>%
  group_by(Training, Block, BlockRole, Genotype) %>%
  summarise(
    stim_min = min(stimulus),
    stim_max = if (extrapolate_curve) max(stimulus)
    else min(max(stimulus), fit_max_stimulus),
    .groups  = "drop"
  ) %>%
  left_join(exp_params, by = c("Training", "Block", "BlockRole", "Genotype")) %>%
  rowwise() %>%
  mutate(grid = list(tibble(stimulus = seq(stim_min, stim_max, length.out = 200)))) %>%
  unnest(grid) %>%
  ungroup() %>%
  mutate(
    stimulus0    = stimulus - 1,
    stimulus_log = log(stimulus),
    x_pred       = if (use_log_stimulus) stimulus_log else stimulus,
    fit = if (use_double_exp & use_asymptote) {
      asymptote + amplitude_fast * exp(-rate_fast * x_pred) +
        amplitude_slow * exp(-rate_slow * x_pred)
    } else if (use_double_exp & !use_asymptote) {
      amplitude_fast * exp(-rate_fast * x_pred) +
        amplitude_slow * exp(-rate_slow * x_pred)
    } else if (!use_double_exp & use_asymptote) {
      asymptote + amplitude_fast * exp(-rate_fast * x_pred)
    } else {
      amplitude_fast * exp(-rate_fast * x_pred)
    },
    fit      = pmin(pmax(fit, 0), 1),           # clamp to [0, 1]
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  ) %>%
  select(-x_pred)

# ==============================================================================
# Plot helper
# ==============================================================================
make_plot <- function(training_label, legend_pos = "none") {
  ggplot(
    new_data_exp %>% filter(Training == training_label),
    aes(x = stimulus, color = Genotype)
  ) +
    facet_grid(Genotype ~ Block, scales = "free_x") +
    geom_point(
      data = raw_summary_joint %>% filter(Training == training_label),
      aes(x = stimulus, y = p_move, color = Genotype),
      inherit.aes = FALSE, alpha = 0.45, size = 1.0
    ) +
    geom_line(aes(y = fit), linewidth = 1.1, na.rm = TRUE) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_pubr(base_size = 12) +
    labs(
      x        = "Stimulus number within block",
      y        = "Response probability",
      title    = paste0("Aggregated response probability: ",
                        tools::toTitleCase(training_label), " training"),
      subtitle = paste0(
        "Points = group means; lines = ",
        if (use_asymptote) "asymptote + " else "",
        if (use_double_exp) "double " else "single ",
        "exp decay on ",
        if (use_log_stimulus) "log(stimulus)" else "stimulus"
      )
    ) +
    theme(legend.position = legend_pos)
}

p_massed <- make_plot("massed", legend_pos = "none")
p_spaced <- make_plot("spaced", legend_pos = "top")

print(p_massed)
print(p_spaced)

ggsave(
  file.path(save_fig_dir, paste0('response_prob', "aggregated_massed.png")),
  p_massed , width = 10, height = 8, dpi = 300, bg = "white"
)
ggsave(
  file.path(save_fig_dir, paste0('response_prob', "aggregated_spaced.png")),
  p_spaced , width = 10, height = 8, dpi = 300, bg = "white"
)


# ==============================================================================
# Optional: model comparison (uncomment to run)
# ==============================================================================
# exp_params_log <- raw_summary_joint %>%
#   filter(stimulus <= fit_max_stimulus) %>%
#   group_by(Training, Block, BlockRole, Genotype) %>%
#   group_modify(~ fit_resp_prob(.x, use_log = TRUE,  use_asym = use_asymptote, use_double = use_double_exp)) %>%
#   ungroup()
#
# exp_params_raw <- raw_summary_joint %>%
#   filter(stimulus <= fit_max_stimulus) %>%
#   group_by(Training, Block, BlockRole, Genotype) %>%
#   group_modify(~ fit_resp_prob(.x, use_log = FALSE, use_asym = use_asymptote, use_double = use_double_exp)) %>%
#   ungroup()
#
# model_comparison <- bind_rows(exp_params_log, exp_params_raw) %>%
#   filter(converged) %>%
#   select(Training, Block, BlockRole, Genotype, fit_type, AIC, BIC, RSS, n_obs) %>%
#   arrange(Training, Block, BlockRole, Genotype, AIC)
#
# print(model_comparison)
# readr::write_csv(model_comparison, file.path(save_results_dir, "model_comparison_response_prob.csv"))