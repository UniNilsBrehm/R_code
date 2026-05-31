###############################################################################
# GLMM Analysis: Dark Flash Blocks
# for response probability, response strength and response delay
# Author: Nils Brehm
# Date: 05/2026
#
# Description:
#   This script analyzes habituation behavior of larval zebrafish to dark flash 
#   experiments. It fits a GLMM model predicting the response strength across 
#   stimulus blocks, validates the model, visualizes habituation curves, and 
#   computes estimated marginal means (EMMs) and contrasts between genotypes 
#   and blocks.
#
# Experimental Design:
#   In each block dark flash (DF: brief period of darkness) is presented every
#   60 seconds. There are two blocks with a inter-block pause of 1 hour. Each
#   block gas 60 DF stimuli. The analysis is based on the "distance moved" in
#   response to each DF.
#
###############################################################################

# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
library(readr)        # Reading CSV files
library(DHARMa)       # Residual diagnostics for (G)LMMs
library(emmeans)      # Estimated marginal means (EMMs) and contrasts
library(glmmTMB)      # Generalized linear mixed models
library(ggplot2)      # Visualization
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(performance)  # Model diagnostics (AIC, R², etc.)
library(ggpubr)       # Publication-ready plots
library(ordinal)

# ==============================================================================
machine <- Sys.info()[["nodename"]]

paths <- switch(
  machine,
  "DESKTOP-5N5AJ0U"     = list(
    code = "C:/Users/NilsPC/Desktop/Susana/R_code/susana",
    data = "C:/Users/NilsPC/Desktop/Susana/Susana"),

  "UNIFREIBURG" = list(
    code = "C:/UniFreiburg/Code/R_code/susana",
    data = "D:/WorkingData/Susana"),
  
  stop("Unknown machine: ", machine, " — add it to the switch block")
)

source(file.path(paths$code, "GLMM/DarkFlash_ISI90s_2Blocks/utils.R"))
base_dir <- file.path(paths$data, "GLMM/DarkFlash_ISI90s_2Blocks")


file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
  # "SPZ_ISI60_ALL_GENOTYPES_Block1_Block2_baseline_subtracted_metrics.csv"
)
# ==============================================================================
# Load Data
# ==============================================================================
message("Loading data...")
res <- load_data_darkflash_60s(file_dir, move_th = 0, , take_peak = 0)

df_final <- res$df_final
df_final_sub <- res$df_final_sub
df_final_sub$delay_ord <- factor(
  df_final_sub$delay,
  ordered = TRUE
)

# Number of animals per genotype
n_per_genotype <- df_final %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

print(n_per_genotype)

df_final <- df_final %>%
  mutate(fish_id = paste(Video, Well, sep = "_"))
df_final_sub <- df_final_sub %>%
  mutate(fish_id = paste(Video, Well, sep = "_"))

df_final_sub <- df_final_sub %>%
  mutate(peak_log = log(max_peak))

# base_dir <- "D:/WorkingData/Susana/GLMM/DarkFlash_ISI90s_2Blocks/TEST"
save_fig_dir = file.path(base_dir, "figs")
save_results_dir = file.path(base_dir, "results")
models_dir = file.path(base_dir, "models")

# Create directories if they do not exist
dir.create(save_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)


geno_colors <- c(
  "ABTL" = "#440154",
  "tyr" = "#3B528B",
  "th2, tyr" = "#21908C",
  "th, tyr" = "#5DC863",
  "th, th2, tyr" = "#FDE725"
)

# ==============================================================================
# Explore Distributions (Histograms)
# ==============================================================================

# # Peak Distance Distribution
# h1 <- ggplot(df_final_sub, aes(x = max_peak)) +
#   geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
#   facet_grid(cols = vars(Genotype)) +
#   theme_minimal() +
#   labs(
#     title = "Peak Distance Distribution",
#     x = "Peak Distance Moved (mm)",
#     y = "Count"
#   )
# h1
# ggsave(
#   filename = file.path(base_dir, "figs", "peak_distance_dist.png"),
#   plot = h1,
#   width = 6,
#   height = 4,
#   bg = "white",
#   dpi = 300
# )
# 
# # Summed Distance Distribution
# h2 <- ggplot(df_final_sub, aes(x = max_cumsum)) +
#   geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
#   facet_grid(cols = vars(Genotype)) +
#   theme_minimal() +
#   labs(
#     title = "Summed Distance Distribution",
#     x = "Summed Distance Moved (mm)",
#     y = "Count"
#   )
# h2
# ggsave(
#   filename = file.path(base_dir, "figs", "summed_distance_dist.png"),
#   plot = h2,
#   width = 6,
#   height = 4,
#   bg = "white",
#   dpi = 300
# )
# 
# # Summed Distance Distribution
# h3 <- ggplot(df_final_sub, aes(x = delay)) +
#   geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
#   facet_grid(cols = vars(Genotype)) +
#   theme_minimal() +
#   labs(
#     title = "Delay Distribution",
#     x = "Delay (s)",
#     y = "Count"
#   )
# h3
# ggsave(
#   filename = file.path(base_dir, "figs", "delay_dist.png"),
#   plot = h3,
#   width = 6,
#   height = 4,
#   bg = "white",
#   dpi = 300
# )


# ==============================================================================
# Fit Models
# ==============================================================================
message("Fitting GLMM models...")

df_final_sub$stimulus_inv <- 1 / df_final_sub$stimulus

# --- Model 1: Peak Movement (Gamma GLMM) -------------------------------------
m_peak <- glmmTMB(
  max_peak ~ Genotype * stimulus_log * Block + (stimulus_log | Video/Well),
  family = Gamma(link = "log"),
  # family = t_family(link = "log"),
  data = df_final_sub
)

m_peak_v2 <- glmmTMB(
  max_peak ~ Genotype * stimulus_log * Block + stimulus_inv + (stimulus_log | Video/Well),
  family = Gamma(link = "log"),
  data = df_final_sub
)

# --- Model 2: Summed Distance (Gamma GLMM) -----------------------------------
m_sum <- glmmTMB(
  max_cumsum ~ Genotype * stimulus_log * Block + (stimulus_log | Video/Well),
  family = Gamma(link = "log"),
  data = df_final_sub
)

m_sum_v2 <- glmmTMB(
  max_cumsum ~ Genotype * stimulus_log * Block + stimulus_inv + (stimulus_log | Video/Well),
  family = Gamma(link = "log"),
  data = df_final_sub
)

# --- Model 3: Response Delay (Gaussian GLMM) ----------------------------------
df_final_sub$delay_non_zero <- df_final_sub$delay + 0.001
m_delay <- glmmTMB(
  delay_non_zero ~ Genotype * stimulus_log * Block + (stimulus_log || Video/Well),
  family = gaussian(link = "identity"),
  data = df_final_sub
)

m_delay_ordinal <- clmm(
  delay_ord ~ Genotype * stimulus_log * Block +
    (1 | Video) +
    (1 | Video:Well),
  data = df_final_sub,
  link = "logit"
)

# --- Model 4: Response Prob (Binomial GLMM) -----------------------------------
m_prob <- glmmTMB(
  move ~ Genotype * Block * stimulus_log + (stimulus_log | Video/Well),
  family = binomial(link = "logit"),
  data = df_final
)

# ==============================================================================
# Save Model Fits to HDD
# ==============================================================================
saveRDS(m_peak, file = file.path(base_dir, "models", "m_peak.rds"))
saveRDS(m_sum, file = file.path(base_dir, "models", "m_sum.rds"))
saveRDS(m_delay, file = file.path(base_dir, "models", "m_delay.rds"))
saveRDS(m_delay_ordinal, file = file.path(base_dir, "models", "m_delay_ordinal.rds"))
saveRDS(m_prob, file = file.path(base_dir, "models", "m_prob.rds"))

# ==============================================================================
# Load Model Fits to HDD
# ==============================================================================
m_peak <- readRDS(file.path(base_dir, "models", "m_peak.rds"))
m_sum <- readRDS(file.path(base_dir, "models", "m_sum.rds"))
m_delay <- readRDS(file.path(base_dir, "models", "m_delay.rds"))
m_delay_ordinal <- readRDS(file.path(base_dir, "models", "m_delay_ordinal.rds"))
m_prob <- readRDS(file.path(base_dir, "models", "m_prob.rds"))

# ==============================================================================
# Model Validation
# ==============================================================================
message("Validating models...")

# Residual check1
model_residuals_check(m_peak, df_final_sub)
model_residuals_check(m_sum, df_final_sub)
model_residuals_check(m_delay, df_final_sub)
model_residuals_check(m_prob, df_final)

# Full validation
validate_model(m_peak, df_final_sub)
validate_model(m_sum, df_final_sub)
validate_model(m_delay, df_final_sub)
validate_model(m_prob, df_final)

# ==============================================================================
# COMPARISON TESTS
# ==============================================================================
# ------------------------------------------------------------------------------
# Response Prob.
# ------------------------------------------------------------------------------
get_all_comparisons(
  m_prob, df_final, n_stim=60, label_name='response prob.', 
  save_dir=file.path(base_dir, "results", "glmm_response_prob_comparisons.txt")
  )

# ------------------------------------------------------------------------------
# Peak Distance
# ------------------------------------------------------------------------------
get_all_comparisons(
  m_peak, df_final_sub, n_stim=60, label_name='peak distance', 
  save_dir=file.path(base_dir, "results", "glmm_peak_distance_comparisons.txt")
)

# ------------------------------------------------------------------------------
# Summed Distance
# ------------------------------------------------------------------------------
get_all_comparisons(
  m_sum, df_final_sub, n_stim=60, label_name='summed distance', 
  save_dir=file.path(base_dir, "results", "glmm_summed_distance_comparisons.txt")
)

# ------------------------------------------------------------------------------
# Response Delay (Ordinal Model)
# ------------------------------------------------------------------------------
get_all_comparisons(
  m_delay_ordinal, df_final_sub, n_stim=60, label_name='delay (ordinal)', 
  save_dir=file.path(base_dir, "results", "glmm_delay_ordinal_comparisons.txt")
)


# ==============================================================================
# Plot Habituation Curves Separated - loop over raw_points styles
# ==============================================================================
plot_habituation_glmm_by_genotype_block(
  df_final  = df_final_sub,
  model     = m_sum,
  label     = "Summed distance moved (mm)",
  raw_var   = "max_cumsum",
  transform = "exp",
  raw_points = "mean",
  y_limits  = c(0, 20),
  y_break   = 2,
  fish_id        = "fish_id",
  colors         = geno_colors,
  genotype_order = names(geno_colors),
)

plot_habituation_glmm_by_genotype_block_random_effects(
  df_final  = df_final_sub,
  model     = m_sum_v2,
  label     = "Peak distance moved (mm)",
  raw_var   = "max_cumsum",
  transform = "exp",
  raw_points = "mean",
  y_limits  = c(0, 20),
  y_break   = 2,
  fish_id        = "fish_id",
  colors         = geno_colors,
  genotype_order = names(geno_colors),
)

plot_configs <- list(
  list(
    fn        = plot_habituation_glmm_by_genotype_block,
    model     = m_peak,
    label     = "Peak distance moved (mm)",
    raw_var   = "max_peak",
    transform = "exp",
    filename  = "GLMM_peak_distance",
    y_limits  = c(0, 10),
    y_break   = 2
  ),
  list(
    fn        = plot_habituation_glmm_by_genotype_block,
    model     = m_sum,
    label     = "Summed distance moved (mm)",
    raw_var   = "max_cumsum",
    transform = "exp",
    filename  = "GLMM_summed_distance",
    y_limits  = c(0, 20),
    y_break   = 5
  ),
  list(
    fn        = plot_habituation_glmm_by_genotype_block,
    model     = m_delay,
    label     = "Response delay (s)",
    raw_var   = "delay",
    transform = "none",
    filename  = "GLMM_delay_gaussian",
    y_limits  = c(0, 3),
    y_break   = 1
  ),
  list(
    fn        = plot_habituation_clmm_by_genotype_block,
    model     = m_delay_ordinal,
    label     = "Response delay (s)",
    raw_var   = "delay",
    transform = NULL,   # clmm function has no transform arg
    filename  = "GLMM_delay_ordinal",
    y_limits  = c(0, 3),
    y_break   = 1
  ),
  list(
    fn        = plot_habituation_glmm_by_genotype_block,
    model     = m_prob,
    label     = "Response probability",
    raw_var   = "move",
    transform = "plogis",
    filename  = "GLMM_response_prob",
    df        = df_final,   # uses full df_final, not df_final_sub
    y_limits  = c(0, 1),
    y_break   = 0.2
  )
)

for (rp in c("raw", "mean")) {
  for (cfg in plot_configs) {
    
    df <- if (!is.null(cfg$df)) cfg$df else df_final_sub
    
    # Build argument list shared by both functions
    args <- list(
      df_final       = df,
      model          = cfg$model,
      label          = cfg$label,
      raw_var        = cfg$raw_var,
      raw_points     = rp,
      fish_id        = "fish_id",
      colors         = geno_colors,
      genotype_order = names(geno_colors),
      y_limits       = cfg$y_limits,
      y_break        = cfg$y_break
    )
    
    # glmm function also takes transform; clmm does not
    if (!is.null(cfg$transform)) {
      args$transform <- cfg$transform
    }
    
    p <- do.call(cfg$fn, args)
    
    out_file <- file.path(base_dir, "figs",
                          paste0(cfg$filename, "_", rp, ".png"))
    print(p)
    ggsave(filename = out_file, plot = p, width = 14, height = 7, dpi = 300)
    message("Saved: ", out_file)
  }
}

# ==============================================================================
# GLMM Contrast Plots
# ==============================================================================
# ==============================================================================
# 0. Helpers
# ==============================================================================

p_to_stars <- function(p) {
  case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE        ~ "ns"
  )
}

normalise_emm <- function(df) {
  # Response column → y_val
  if ("response" %in% names(df)) {
    df$y_val <- df$response
  } else if ("prob" %in% names(df)) {
    df$y_val <- df$prob
  } else if ("emmean" %in% names(df)) {
    df$y_val <- df$emmean
  } else {
    stop("Cannot find response column. Available: ",
         paste(names(df), collapse = ", "))
  }
  
  # CI columns → y_lwr, y_upr
  if ("asymp.LCL" %in% names(df)) {
    df$y_lwr <- df$asymp.LCL
    df$y_upr <- df$asymp.UCL
  } else if ("lower.CL" %in% names(df)) {
    df$y_lwr <- df$lower.CL
    df$y_upr <- df$upper.CL
  } else {
    df$y_lwr <- df$y_val - 1.96 * df$SE
    df$y_upr <- df$y_val + 1.96 * df$SE
  }
  df
}

normalise_contrast <- function(df) {
  # Value column
  if ("odds.ratio" %in% names(df)) {
    df$value    <- df$odds.ratio
    df$is_ratio <- TRUE
  } else if ("ratio" %in% names(df)) {
    df$value    <- df$ratio
    df$is_ratio <- TRUE
  } else if ("estimate" %in% names(df)) {
    df$value    <- df$estimate
    df$is_ratio <- FALSE
  } else {
    stop("Cannot find value column. Available: ",
         paste(names(df), collapse = ", "))
  }
  
  df %>% mutate(
    lwr      = if (unique(is_ratio)) value / exp(1.96 * SE)
    else                  value - 1.96 * SE,
    upr      = if (unique(is_ratio)) value * exp(1.96 * SE)
    else                  value + 1.96 * SE,
    stars    = p_to_stars(p.value),
    sig      = p.value < 0.05,
    null_val = if (unique(is_ratio)) 1    else 0,
    y_label  = if (unique(is_ratio)) "Ratio" else "Difference"
  )
}

# ==============================================================================
# 1. Extract comparisons from model
# ==============================================================================

extract_comparisons <- function(model, df, n_stim = 60) {
  
  stim_vals <- log(c(
    mean(1:round(n_stim * 0.1)),
    mean(round(n_stim * 0.4):round(n_stim * 0.6)),
    mean(round(n_stim * 0.9):n_stim)
  ))
  stim_full <- mean(log(1:n_stim))
  
  results <- list()
  
  for (blk in c("Block1", "Block2")) {
    for (win in c("full", "first", "middle", "last")) {
      
      sv <- switch(win,
                   full   = stim_full,
                   first  = stim_vals[1],
                   middle = stim_vals[2],
                   last   = stim_vals[3]
      )
      
      emm <- emmeans(
        model,
        ~ Genotype,
        at   = list(stimulus_log = sv, Block = blk),
        type = "response"
      )
      
      con <- contrast(emm, method = "pairwise") %>%
        as.data.frame() %>%
        mutate(Block = blk, Window = win)
      
      results[[paste0(blk, "_", win, "_emm")]] <- as.data.frame(emm) %>%
        mutate(Block = blk, Window = win)
      results[[paste0(blk, "_", win, "_con")]] <- con
    }
  }
  
  # Slopes
  slopes <- emtrends(
    model, ~ Genotype | Block,
    var = "stimulus_log"
  ) %>% as.data.frame()
  
  slope_con <- emtrends(
    model, pairwise ~ Genotype | Block,
    var = "stimulus_log"
  )$contrasts %>% as.data.frame()
  
  # Between-block EMMs
  bb_emm <- emmeans(
    model,
    ~ Genotype * Block,
    at   = list(stimulus_log = stim_full),
    type = "response"
  ) %>% as.data.frame()
  
  # Between-block contrasts
  bb_con <- emmeans(
    model,
    ~ Block | Genotype,
    at   = list(stimulus_log = stim_full),
    type = "response"
  ) %>%
    contrast(method = "pairwise") %>%
    as.data.frame()
  
  list(
    emm       = normalise_emm(
      bind_rows(results[grep("_emm$", names(results))])),
    contrasts = normalise_contrast(
      bind_rows(results[grep("_con$", names(results))])),
    slopes    = slopes,
    slope_con = slope_con,
    bb_emm    = normalise_emm(bb_emm),
    bb_con    = normalise_contrast(bb_con)
  )
}

# ==============================================================================
# 2. Plot: pairwise contrasts at each stimulus window
# ==============================================================================

plot_contrasts <- function(contrast_list, outcome_label) {
  
  df_con <- bind_rows(
    lapply(names(contrast_list), function(nm) {
      contrast_list[[nm]] %>% mutate(window = nm)
    })
  ) %>%
    tidyr::separate(window, into = c("Block", "Window"), sep = "_") %>%
    mutate(
      Window = factor(Window, levels = c("first", "middle", "last", "full")),
      Block  = factor(Block)
    )
  
  if (!"value" %in% names(df_con)) df_con <- normalise_contrast(df_con)
  
  null_val <- unique(df_con$null_val)
  y_label  <- unique(df_con$y_label)
  
  ggplot(
    df_con,
    aes(x = contrast, y = value,
        ymin = lwr, ymax = upr,
        color = sig)
  ) +
    geom_hline(yintercept = null_val,
               linetype = "dashed", color = "grey50") +
    geom_pointrange(linewidth = 0.7, size = 0.5) +
    geom_text(
      aes(label = stars, y = upr),
      vjust = -0.5, size = 3.5, color = "black"
    ) +
    coord_flip() +
    facet_grid(Block ~ Window) +
    scale_color_manual(
      values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
      labels = c("TRUE" = "p < 0.05", "FALSE" = "ns"),
      name   = NULL
    ) +
    labs(
      title    = paste0("Pairwise contrasts — ", outcome_label),
      subtitle = "Red = significant (p < 0.05, Tukey-adjusted)  |  windows: first / middle / last / full block",
      x        = NULL,
      y        = y_label
    ) +
    theme_pubr(base_size = 12) +
    theme(
      legend.position = "top",
      strip.text      = element_text(size = 10)
    )
}

# ==============================================================================
# 3. Plot: habituation slopes
# ==============================================================================

plot_slopes <- function(slope_df, slope_contrast_df, outcome_label) {
  
  # p-value for each slope vs zero
  slope_df <- slope_df %>%
    mutate(
      stars = p_to_stars(
        2 * pnorm(abs(stimulus_log.trend / SE), lower.tail = FALSE)
      )
    )
  
  p_slopes <- ggplot(
    slope_df,
    aes(x = Genotype, y = stimulus_log.trend,
        ymin = asymp.LCL, ymax = asymp.UCL,
        color = Genotype)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(linewidth = 0.8, size = 0.6) +
    geom_text(
      aes(label = stars, y = asymp.UCL),
      vjust = -0.5, size = 3.5, color = "black"
    ) +
    facet_wrap(~Block) +
    scale_color_manual(values = geno_colors) +
    labs(
      title    = paste0("Habituation slopes — ", outcome_label),
      subtitle = "Slope on log(stimulus) | negative = habituation",
      x        = NULL,
      y        = "Slope (log-stimulus scale)"
    ) +
    theme_pubr(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 30, hjust = 1)
    )
  
  slope_contrast_df <- slope_contrast_df %>%
    mutate(
      stars = p_to_stars(p.value),
      sig   = p.value < 0.05,
      lwr   = estimate - 1.96 * SE,
      upr   = estimate + 1.96 * SE
    )
  
  p_slope_con <- ggplot(
    slope_contrast_df,
    aes(x = contrast, y = estimate,
        ymin = lwr, ymax = upr,
        color = sig)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(linewidth = 0.7, size = 0.5) +
    geom_text(
      aes(label = stars, y = upr),
      vjust = -0.5, size = 3.5, color = "black"
    ) +
    coord_flip() +
    facet_wrap(~Block) +
    scale_color_manual(
      values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
      name   = NULL
    ) +
    labs(
      title    = paste0("Slope pairwise contrasts — ", outcome_label),
      subtitle = "Differences in habituation rate between genotypes",
      x        = NULL,
      y        = "Slope difference"
    ) +
    theme_pubr(base_size = 12) +
    theme(legend.position = "top")
  
  p_slopes / p_slope_con +
    plot_annotation(
      title = paste0("Habituation rate analysis — ", outcome_label),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
}

# ==============================================================================
# 4. Plot: between-block comparison
# ==============================================================================

plot_between_blocks <- function(bb_emm_df, bb_contrast_df, outcome_label) {
  
  p_emm <- ggplot(
    bb_emm_df,
    aes(x = Genotype, y = y_val,
        ymin = y_lwr, ymax = y_upr,
        color = Genotype, shape = Block)
  ) +
    geom_pointrange(
      linewidth = 0.8, size = 0.6,
      position  = position_dodge(width = 0.5)
    ) +
    scale_color_manual(values = geno_colors) +
    scale_shape_manual(values = c("Block1" = 16, "Block2" = 17)) +
    labs(
      title    = paste0("Block1 vs Block2 EMMs — ", outcome_label),
      subtitle = "Circles = Block1  |  Triangles = Block2",
      x        = NULL,
      y        = outcome_label
    ) +
    theme_pubr(base_size = 13) +
    theme(
      legend.position = "top",
      axis.text.x     = element_text(angle = 30, hjust = 1)
    )
  
  null_val <- unique(bb_contrast_df$null_val)
  y_label  <- unique(bb_contrast_df$y_label)
  subtitle <- if (unique(bb_contrast_df$is_ratio))
    "Ratio < 1 = lower response in Block2 = retention"
  else
    "Negative = lower response in Block2 = retention"
  
  p_con <- ggplot(
    bb_contrast_df,
    aes(x = Genotype, y = value,
        ymin = lwr, ymax = upr,
        color = sig)
  ) +
    geom_hline(yintercept = null_val,
               linetype = "dashed", color = "grey50") +
    geom_pointrange(linewidth = 0.8, size = 0.6) +
    geom_text(
      aes(label = stars, y = upr),
      vjust = -0.5, size = 4, color = "black"
    ) +
    scale_color_manual(
      values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
      name   = NULL
    ) +
    labs(
      title    = paste0("Block2 / Block1 — ", outcome_label),
      subtitle = subtitle,
      x        = NULL,
      y        = y_label
    ) +
    theme_pubr(base_size = 13) +
    theme(
      legend.position = "top",
      axis.text.x     = element_text(angle = 30, hjust = 1)
    )
  
  p_emm | p_con
}

# ==============================================================================
# 5. Run extraction
# ==============================================================================

comp_peak <- extract_comparisons(m_peak, df_final_sub, n_stim = 60)
comp_prob <- extract_comparisons(m_prob, df_final,     n_stim = 60)
comp_sum  <- extract_comparisons(m_sum,  df_final_sub, n_stim = 60)

# ==============================================================================
# 6. Generate all plots
# ==============================================================================

# --- Pairwise contrasts ---
p_con_peak <- plot_contrasts(
  split(comp_peak$contrasts,
        paste0(comp_peak$contrasts$Block, "_",
               comp_peak$contrasts$Window)),
  "Peak distance"
)

p_con_prob <- plot_contrasts(
  split(comp_prob$contrasts,
        paste0(comp_prob$contrasts$Block, "_",
               comp_prob$contrasts$Window)),
  "Response probability"
)

p_con_sum <- plot_contrasts(
  split(comp_sum$contrasts,
        paste0(comp_sum$contrasts$Block, "_",
               comp_sum$contrasts$Window)),
  "Summed distance"
)

# --- Slopes ---
p_slope_peak <- plot_slopes(comp_peak$slopes, comp_peak$slope_con,
                            "Peak distance")
p_slope_prob <- plot_slopes(comp_prob$slopes, comp_prob$slope_con,
                            "Response probability")
p_slope_sum  <- plot_slopes(comp_sum$slopes,  comp_sum$slope_con,
                            "Summed distance")

# --- Between blocks ---
p_bb_peak <- plot_between_blocks(comp_peak$bb_emm, comp_peak$bb_con,
                                 "Peak distance")
p_bb_prob <- plot_between_blocks(comp_prob$bb_emm, comp_prob$bb_con,
                                 "Response probability")
p_bb_sum  <- plot_between_blocks(comp_sum$bb_emm,  comp_sum$bb_con,
                                 "Summed distance")

# ==============================================================================
# 7. Save
# ==============================================================================

ggsave(file.path(save_fig_dir, "GLMM_contrasts_peak_pairwise.png"),
       p_con_peak,   width = 16, height = 10, dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "GLMM_contrasts_prob_pairwise.png"),
       p_con_prob,   width = 16, height = 10, dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "GLMM_contrasts_sum_pairwise.png"),
       p_con_sum,    width = 16, height = 10, dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "GLMM_slopes_peak.png"),
       p_slope_peak, width = 12, height = 8,  dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "GLMM_slopes_prob.png"),
       p_slope_prob, width = 12, height = 8,  dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "GLMM_slopes_sum.png"),
       p_slope_sum,  width = 12, height = 8,  dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "GLMM_between_blocks_peak.png"),
       p_bb_peak,    width = 12, height = 6,  dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "GLMM_between_blocks_prob.png"),
       p_bb_prob,    width = 12, height = 6,  dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "GLMM_between_blocks_sum.png"),
       p_bb_sum,     width = 12, height = 6,  dpi = 300, bg = "white")
