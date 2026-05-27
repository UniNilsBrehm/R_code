# ==============================================================================
# Aggregated response delay — saturating / logistic growth fit
# delay ∈ {0,1,2,3,4} seconds; rises quickly then plateaus
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(readr)
library(minpack.lm)

# ==============================================================================
# Paths  (unchanged — adjust as needed)
# ==============================================================================
source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_Massed_vs_Spaced//utils.R")

base_dir <- "D:/WorkingData/Susana/NLME/DarkFlash_Joint_SpacedVsMassed/"

select_first_responders <- FALSE
conditional_on_response <- TRUE
th_for_move <- 1

file_massed <- file.path(base_dir, "data_files", "SPZ_Massed_Training_7Nov2025.csv")
file_spaced <- file.path(base_dir, "data_files", "SPZ_Spaced_Training_Nov2025.csv")

save_fig_dir     <- file.path(base_dir, "figs",    "aggregated", "response_delay")
save_results_dir <- file.path(base_dir, "results", "aggregated", "response_delay")
save_model_dir   <- file.path(base_dir, "models")

dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_model_dir,   recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Load and prepare both datasets  (identical to response_prob script)
# ==============================================================================
# res_massed <- load_data(file_massed, move_th = th_for_move, drop = c('th2, tyr', 'tyr'))
# res_spaced <- load_data(file_spaced, move_th = th_for_move, drop = c('th2, tyr', 'tyr'))

res_massed <- load_data(file_massed, move_th = th_for_move)
res_spaced <- load_data(file_spaced, move_th = th_for_move)

if (conditional_on_response) {
  df_massed <- res_massed$df_final_sub
  df_spaced <- res_spaced$df_final_sub
  
}else {
  df_massed <- res_massed$df_final
  df_spaced <- res_spaced$df_final
}


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
  df_all <- df_all %>%
    group_by(animal) %>%
    filter(!any(Block == 1 & stimulus == 1 & move == 0)) %>%
    ungroup()
}

# ==============================================================================
# Switches
# ==============================================================================
use_log_stimulus  <- TRUE   # TRUE = power-law x-axis (log stimulus)
use_logistic      <- FALSE    # TRUE = logistic (3-param); FALSE = saturating exponential + baseline
fit_max_stimulus  <- Inf
extrapolate_curve <- TRUE

DELAY_MAX         <- 4      # hard ceiling for predictions
rho_threshold     <- 0.15   # |Spearman rho| below this → flat
range_threshold   <- 1.0    # delay range below this → flat

# ==============================================================================
# Fitting function
# ==============================================================================
fit_resp_delay <- function(dat, use_log, use_logistic_model) {
  
  dat <- dat %>%
    mutate(x = if (use_log) stimulus_log else stimulus) %>%
    filter(is.finite(x), is.finite(mean_delay), !is.na(mean_delay)) %>%
    arrange(x)
  
  n_stim <- length(unique(dat$x))
  if (nrow(dat) < 3 || n_stim < 3) {
    return(tibble(
      baseline  = NA_real_, plateau = NA_real_,
      rate      = NA_real_, x0      = NA_real_,
      fit_type  = "not_fitted", converged = FALSE,
      AIC = NA_real_, BIC = NA_real_, RSS = NA_real_, n_obs = nrow(dat)
    ))
  }
  
  # --- starting values, robust to flat or declining trends ---
  baseline_start <- max(min(dat$mean_delay, na.rm = TRUE), 0)
  max_delay      <- max(dat$mean_delay, na.rm = TRUE)
  plateau_start  <- max(max_delay - baseline_start, 0.1)
  
  if (mean(diff(dat$mean_delay)) <= 0) {
    baseline_start <- max(mean(dat$mean_delay, na.rm = TRUE) - 0.1, 0)
    plateau_start  <- 0.1
  }
  
  rate_candidates <- c(0.5, 1, 2, 5, 0.1, 0.05, 3, 10, 0.01, 0.005, 0.001, 20)
  
  try_fit <- function(formula_nls, starts, lower_bounds, upper_bounds) {
    fit <- NULL
    for (r in rate_candidates) {
      starts$rate <- r
      if ("x0" %in% names(starts)) starts$x0 <- median(dat$x, na.rm = TRUE)
      fit <- tryCatch(
        withCallingHandlers(
          nlsLM(
            formula_nls,
            data    = dat,
            start   = starts,
            lower   = lower_bounds,
            upper   = upper_bounds,
            control = nls.lm.control(maxiter = 500, ftol = 1e-7, ptol = 1e-7)
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
  
  if (use_logistic_model) {
    formula_nls  <- mean_delay ~ plateau / (1 + exp(-rate * (x - x0)))
    starts       <- list(plateau = baseline_start + plateau_start, rate = 2,
                         x0 = median(dat$x, na.rm = TRUE))
    lower_bounds <- c(plateau = 0,         rate = 1e-4, x0 = min(dat$x))
    upper_bounds <- c(plateau = DELAY_MAX, rate = 20,   x0 = max(dat$x))
    fit_label    <- paste0("logistic_", if (use_log) "logstim" else "rawstim")
  } else {
    formula_nls  <- mean_delay ~ baseline + plateau * (1 - exp(-rate * x))
    starts       <- list(baseline = baseline_start, plateau = plateau_start, rate = 1)
    lower_bounds <- c(baseline = 0,         plateau = 0,         rate = 1e-4)
    upper_bounds <- c(baseline = DELAY_MAX, plateau = DELAY_MAX, rate = 50)
    fit_label    <- paste0("sat_exp_baseline_", if (use_log) "logstim" else "rawstim")
  }
  
  fit <- try_fit(formula_nls, starts, lower_bounds, upper_bounds)
  
  # --- fallback 1: fix rate, fit only baseline + plateau ---
  if (is.null(fit) && !use_logistic_model) {
    fit <- tryCatch(
      withCallingHandlers(
        nlsLM(
          mean_delay ~ baseline + plateau * (1 - exp(-0.01 * x)),
          data    = dat,
          start   = list(baseline = baseline_start, plateau = plateau_start),
          lower   = c(baseline = 0,         plateau = 0),
          upper   = c(baseline = DELAY_MAX, plateau = DELAY_MAX),
          control = nls.lm.control(maxiter = 500)
        ),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) NULL
    )
    if (!is.null(fit)) fit_label <- paste0(fit_label, "_rate_fixed")
  }
  
  # --- fallback 2: horizontal grand mean line ---
  if (is.null(fit)) {
    return(tibble(
      baseline  = mean(dat$mean_delay, na.rm = TRUE),
      plateau   = 0,
      rate      = 0,
      x0        = NA_real_,
      fit_type  = paste0(fit_label, "_mean_fallback"),
      converged = FALSE,
      AIC       = NA_real_,
      BIC       = NA_real_,
      RSS       = sum((dat$mean_delay - mean(dat$mean_delay, na.rm = TRUE))^2),
      n_obs     = nrow(dat)
    ))
  }
  
  coefs    <- coef(fit)
  get_coef <- function(name) if (name %in% names(coefs)) unname(coefs[name]) else NA_real_
  
  tibble(
    baseline  = get_coef("baseline"),
    plateau   = get_coef("plateau"),
    rate      = get_coef("rate"),
    x0        = get_coef("x0"),
    fit_type  = fit_label,
    converged = TRUE,
    AIC       = AIC(fit),
    BIC       = BIC(fit),
    RSS       = sum(residuals(fit)^2),
    n_obs     = nrow(dat)
  )
}

# ==============================================================================
# Aggregate mean delay (df_all, responders already filtered upstream)
# ==============================================================================
raw_summary_delay <- df_all %>%
  group_by(Training, Block, BlockRole, Genotype, stimulus, stimulus0, stimulus_log) %>%
  summarise(
    n_trials   = sum(!is.na(delay)),
    mean_delay = mean(delay, na.rm = TRUE),
    se_delay   = sd(delay,   na.rm = TRUE) / sqrt(n_trials),
    .groups    = "drop"
  )

# ==============================================================================
# Flatness detection
# ==============================================================================
flat_flags <- raw_summary_delay %>%
  group_by(Training, Block, BlockRole, Genotype) %>%
  summarise(
    spearman_rho = cor(stimulus, mean_delay, method = "spearman", use = "complete.obs"),
    delay_range  = max(mean_delay, na.rm = TRUE) - min(mean_delay, na.rm = TRUE),
    grand_mean   = mean(mean_delay, na.rm = TRUE),
    grand_se     = sd(mean_delay,   na.rm = TRUE) / sqrt(sum(!is.na(mean_delay))),
    is_flat      = delay_range < range_threshold,
    .groups      = "drop"
  )

message("Flat groups (no curve fitted, dashed mean line):")
print(flat_flags %>% filter(is_flat) %>%
        select(Training, Block, Genotype, spearman_rho, delay_range, grand_mean))

# ==============================================================================
# Fit — non-flat groups only
# ==============================================================================
delay_params <- raw_summary_delay %>%
  left_join(
    flat_flags %>% select(Training, Block, BlockRole, Genotype, is_flat),
    by = c("Training", "Block", "BlockRole", "Genotype")
  ) %>%
  filter(!is_flat, stimulus <= fit_max_stimulus) %>%
  group_by(Training, Block, BlockRole, Genotype) %>%
  group_modify(~ fit_resp_delay(
    dat                = .x,
    use_log            = use_log_stimulus,
    use_logistic_model = use_logistic
  )) %>%
  ungroup()

print(delay_params)
readr::write_csv(delay_params,
                 file.path(save_results_dir, "delay_params_response_delay.csv"))
readr::write_csv(flat_flags,
                 file.path(save_results_dir, "flat_flags_response_delay.csv"))

# ==============================================================================
# Prediction grid — non-flat groups only
# ==============================================================================
new_data_delay <- raw_summary_delay %>%
  left_join(
    flat_flags %>% select(Training, Block, BlockRole, Genotype, is_flat),
    by = c("Training", "Block", "BlockRole", "Genotype")
  ) %>%
  filter(!is_flat) %>%
  group_by(Training, Block, BlockRole, Genotype) %>%
  summarise(
    stim_min = min(stimulus),
    stim_max = if (extrapolate_curve) max(stimulus)
    else min(max(stimulus), fit_max_stimulus),
    .groups  = "drop"
  ) %>%
  left_join(delay_params, by = c("Training", "Block", "BlockRole", "Genotype")) %>%
  rowwise() %>%
  mutate(grid = list(tibble(stimulus = seq(stim_min, stim_max, length.out = 200)))) %>%
  unnest(grid) %>%
  ungroup() %>%
  mutate(
    stimulus0    = stimulus - 1,
    stimulus_log = log(stimulus),
    x_pred       = if (use_log_stimulus) stimulus_log else stimulus,
    fit = if (use_logistic) {
      plateau / (1 + exp(-rate * (x_pred - x0)))
    } else {
      baseline + plateau * (1 - exp(-rate * x_pred))
    },
    fit      = pmin(pmax(fit, 0), DELAY_MAX),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  ) %>%
  select(-x_pred)

# ==============================================================================
# Plot helper
# ==============================================================================
make_delay_plot <- function(training_label, legend_pos = "none") {
  
  flat_bands <- flat_flags %>%
    filter(Training == training_label, is_flat)
  
  # linetype: solid for converged fits, dashed for mean_fallback
  line_data <- new_data_delay %>%
    filter(Training == training_label) %>%
    mutate(ltype = if_else(converged, "solid", "dashed"))
  
  ggplot(
    line_data,
    aes(x = stimulus, color = Genotype)
  ) +
    facet_grid(Genotype ~ Block, scales = "free_x") +
    # shaded ±SE band for genuinely flat groups
    geom_rect(
      data        = flat_bands,
      aes(xmin = -Inf, xmax = Inf,
          ymin = grand_mean - grand_se,
          ymax = grand_mean + grand_se,
          fill = Genotype),
      inherit.aes = FALSE, alpha = 0.15, color = NA
    ) +
    # dashed mean line for genuinely flat groups
    geom_hline(
      data        = flat_bands,
      aes(yintercept = grand_mean, color = Genotype),
      linetype = "dashed", linewidth = 0.9
    ) +
    # points — all groups
    geom_point(
      data = raw_summary_delay %>% filter(Training == training_label),
      aes(x = stimulus, y = mean_delay, color = Genotype),
      inherit.aes = FALSE, alpha = 0.55, size = 1.2
    ) +
    # fitted curve — solid if converged, dashed if mean_fallback
    geom_line(aes(y = fit, linetype = ltype), linewidth = 1.1, na.rm = TRUE) +
    scale_linetype_identity() +
    coord_cartesian(ylim = c(0, DELAY_MAX)) +
    scale_y_continuous(breaks = 0:DELAY_MAX) +
    theme_pubr(base_size = 12) +
    labs(
      x        = "Stimulus number within block",
      y        = "Mean response delay (s)",
      title    = paste0("Aggregated response delay: ",
                        tools::toTitleCase(training_label), " training"),
      subtitle = paste0(
        "Solid: sat-exp + baseline on log(stimulus); ",
        "dashed: mean fallback or flat block"
      )
    ) +
    theme(legend.position = legend_pos)
}

p_massed_delay <- make_delay_plot("massed", legend_pos = "none")
p_spaced_delay <- make_delay_plot("spaced", legend_pos = "top")

print(p_massed_delay)
print(p_spaced_delay)

ggsave(
  file.path(save_fig_dir, paste0('delay', "_aggregated_massed.png")),
  p_massed_delay , width = 10, height = 8, dpi = 300, bg = "white"
)
ggsave(
  file.path(save_fig_dir, paste0('delay', "_aggregated_spaced.png")),
  p_spaced_delay , width = 10, height = 8, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# Per-animal raw traces — response delay
# --------------------------------------------------------------------------
raw_per_animal_delay <- df_all %>%
  group_by(Training, Block, BlockRole, Genotype, animal, stimulus) %>%
  summarise(
    mean_delay = mean(delay, na.rm = TRUE),
    .groups    = "drop"
  )

make_individual_delay_plot <- function(training_label) {
  ggplot(
    raw_per_animal_delay %>% filter(Training == training_label),
    aes(x = stimulus, y = mean_delay, group = animal, color = Genotype)
  ) +
    facet_grid(Genotype ~ Block, scales = "free_x") +
    geom_line(alpha = 0.5, linewidth = 0.5) +
    geom_point(alpha = 0.4, size = 0.8) +
    # group-level fit overlaid in black
    geom_line(
      data        = new_data_delay %>% filter(Training == training_label),
      aes(x = stimulus, y = fit, group = Genotype),
      inherit.aes = FALSE,
      color       = "black",
      linewidth   = 1.2,
      na.rm       = TRUE
    ) +
    coord_cartesian(ylim = c(0, DELAY_MAX)) +
    scale_y_continuous(breaks = 0:DELAY_MAX) +
    theme_pubr(base_size = 11) +
    labs(
      x        = "Stimulus number within block",
      y        = "Response delay (s)",
      title    = paste0("Individual traces — response delay: ",
                        tools::toTitleCase(training_label), " training"),
      subtitle = "One line per fish; black = group fit"
    ) +
    theme(legend.position = "none")
}

p_ind_massed_delay <- make_individual_delay_plot("massed")
p_ind_spaced_delay <- make_individual_delay_plot("spaced")

print(p_ind_massed_delay)
print(p_ind_spaced_delay)

ggsave(
  file.path(save_fig_dir, "delay_individual_massed.png"),
  p_ind_massed_delay, width = 10, height = 12, dpi = 300, bg = "white"
)
ggsave(
  file.path(save_fig_dir, "delay_individual_spaced.png"),
  p_ind_spaced_delay, width = 10, height = 12, dpi = 300, bg = "white"
)