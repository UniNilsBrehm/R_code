###############################################################################
# Joint GLMM Analysis: Spaced vs Massed Training - Peak Distance Moved
# Author: Nils Brehm
# Date: 2026
#
# Description:
#   Joint GLMM combining the spaced and massed training datasets into one
#   model to test the spaced-vs-massed memory contrast directly for:
#
#     max_peak = peak distance moved
#
#   Model:
#     max_peak ~ Genotype * Training * Block * stimulus_log + (1 | animal)
#
#   Model family:
#     Gamma(link = "log")
#
#   Memory contrasts per genotype:
#     (A) max_peak at first test stimulus, spaced vs massed
#     (B) mean max_peak over first 3 test stimuli, spaced vs massed
#     (C) mean max_peak over all shared test stimuli, spaced vs massed
#     (D) inter-block recovery in max_peak, spaced vs massed
#
#   Negative spaced-massed difference => spaced fish move less in the test block
#   => better memory retention, assuming smaller peak movement reflects stronger
#   retained habituation.
###############################################################################

# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
library(readr)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(tidyr)
library(performance)
library(ggpubr)
library(stringr)

# Load helper functions
source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/GLMM/DarkFlash_Massed_vs_Spaced/utils.R")
from_file <- TRUE  # load pre-fitted model from file

# ==============================================================================
# 1. Paths
# ==============================================================================
base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/GLMM/DarkFlash_Massed_vs_Spaced"

file_massed <- file.path(
  base_dir,
  "data_files",
  "SPZ_Massed_Training_7Nov2025.csv"
)

file_spaced <- file.path(
  base_dir,
  "data_files",
  "SPZ_Spaced_Training_Nov2025.csv"
)

save_results_dir <- file.path(base_dir, "results", "glmm_joint_max_peak")
save_fig_dir     <- file.path(base_dir, "figs",    "glmm_joint_max_peak")

dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. Load and prepare both datasets
# ==============================================================================
# For Gamma models, use the positive-response subset from load_data().
# This matches the old massed max_peak analysis, which used df_final_sub.
res_massed <- load_data(file_massed, move_th = 1, drop = c("th2, tyr", "tyr"))
res_spaced <- load_data(file_spaced, move_th = 1, drop = c("th2, tyr", "tyr"))

df_massed <- res_massed$df_final_sub
df_spaced <- res_spaced$df_final_sub

# Fallback in case your current utils.R does not return df_final_sub.
# Gamma models require strictly positive responses.
if (is.null(df_massed)) {
  df_massed <- res_massed$df_final %>%
    filter(!is.na(max_peak), max_peak > 0)
}
if (is.null(df_spaced)) {
  df_spaced <- res_spaced$df_final %>%
    filter(!is.na(max_peak), max_peak > 0)
}

# Safety filter for Gamma model
df_massed <- df_massed %>% filter(!is.na(max_peak), max_peak > 0)
df_spaced <- df_spaced %>% filter(!is.na(max_peak), max_peak > 0)


# Tag each row with Training, and define BlockRole.
# In each experiment the LAST block is the test block.
massed_test_block <- "2"
spaced_test_block <- "5"

df_massed_tagged <- df_massed %>%
  mutate(
    Training  = "massed",
    BlockRole = ifelse(as.character(Block) == massed_test_block, "test", "training")
  )

df_spaced_tagged <- df_spaced %>%
  mutate(
    Training  = "spaced",
    BlockRole = ifelse(as.character(Block) == spaced_test_block, "test", "training")
  )


# CRITICAL: make animal IDs unique across experiments.
# Same Video x Well combo from different experiments must not collide.
df_all <- bind_rows(df_massed_tagged, df_spaced_tagged) %>%
  mutate(
    Training  = factor(Training,  levels = c("massed", "spaced")),
    Block     = factor(Block),
    BlockRole = factor(BlockRole, levels = c("training", "test")),
    Genotype  = factor(Genotype),
    Video     = factor(Video),
    Well      = factor(Well),
    animal    = factor(paste0(Training, "_", Video, ".", Well))
  )


# Sanity checks ----------------------------------------------------------------
cat("\n--- Animals per Genotype x Training ---\n")
print(
  df_all %>%
    distinct(animal, Genotype, Training) %>%
    count(Genotype, Training)
)

cat("\n--- Rows per Training x Block ---\n")
print(df_all %>% count(Training, Block, BlockRole))

cat("\n--- Stimulus range per Training x Block ---\n")
print(
  df_all %>%
    group_by(Training, Block) %>%
    summarise(min_stim = min(stimulus), max_stim = max(stimulus), .groups = "drop")
)

cat("\n--- max_peak summary by Training x Block ---\n")
print(
  df_all %>%
    group_by(Training, Block, BlockRole) %>%
    summarise(
      n = n(),
      mean_max_peak = mean(max_peak, na.rm = TRUE),
      median_max_peak = median(max_peak, na.rm = TRUE),
      min_max_peak = min(max_peak, na.rm = TRUE),
      max_max_peak = max(max_peak, na.rm = TRUE),
      .groups = "drop"
    )
)


# ==============================================================================
# 3. Exploratory distribution plot
# ==============================================================================
message("Plotting max_peak distributions...")
p_dist <- ggplot(df_all, aes(x = max_peak)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_grid(Training ~ Genotype, scales = "free_y") +
  theme_pubr(base_size = 12) +
  labs(
    title = "Peak distance distribution: joint spaced + massed dataset",
    x = "Peak distance moved, max_peak",
    y = "Count"
  )

ggsave(
  file.path(save_fig_dir, "joint_max_peak_distribution.png"),
  p_dist, width = 12, height = 6, dpi = 300, bg = "white"
)


# ==============================================================================
# 4. Fit the joint GLMM
# ==============================================================================
message("Fitting joint Gamma GLMM for max_peak...")

# Same joint structure as the response-probability model,
# but with Gamma(log) because max_peak is positive continuous.
if (from_file){
  # Load fitted model if available
  m_joint_peak  <- readRDS(
    file.path(save_results_dir, "joint_glmm_spaced_vs_massed_max_peak.rds")
  )

} else{
  m_joint_peak <- glmmTMB(
    max_peak ~ Genotype * Training * Block * stimulus_log + (1 | animal),
    family = Gamma(link = "log"),
    data   = df_all
  )
  
  saveRDS(
    m_joint_peak,
    file.path(save_results_dir, "joint_glmm_spaced_vs_massed_max_peak.rds")
  )
}

capture.output(
  summary(m_joint_peak),
  file = file.path(save_results_dir, "summary_results_max_peak.txt")
)

print(summary(m_joint_peak))

# ==============================================================================
# 5. Validate model
# ==============================================================================
message("Validating joint max_peak model residuals...")
validate_model(m_joint_peak, df_all)


# ==============================================================================
# 6. Plot habituation curves: predicted max_peak across stimuli
# ==============================================================================
message("Plotting joint max_peak habituation curves...")

# Build prediction grid: per (Training, Block, Genotype) at each observed stimulus
new_data_joint <- df_all %>%
  group_by(Training, Block, Genotype) %>%
  summarise(
    stim_min = min(stimulus, na.rm = TRUE),
    stim_max = max(stimulus, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    grid = list(tibble(stimulus = seq(max(stim_min, 1), stim_max, length.out = 200)))
  ) %>%
  unnest(grid) %>%
  select(-stim_min, -stim_max) %>%
  mutate(
    stimulus_log = log(stimulus),
    Training     = factor(Training, levels = levels(df_all$Training)),
    Block        = factor(Block,    levels = levels(df_all$Block)),
    Genotype     = factor(Genotype, levels = levels(df_all$Genotype))
  )

pred_joint <- predict(
  m_joint_peak,
  newdata = new_data_joint,
  re.form = NA,
  se.fit  = TRUE,
  type    = "link"
)

new_data_joint <- new_data_joint %>%
  mutate(
    fit     = exp(pred_joint$fit),
    CI_low  = exp(pred_joint$fit - 1.96 * pred_joint$se.fit),
    CI_high = exp(pred_joint$fit + 1.96 * pred_joint$se.fit)
  )

raw_summary_joint <- df_all %>%
  group_by(Training, Block, Genotype, stimulus) %>%
  summarise(
    mean_max_peak = mean(max_peak, na.rm = TRUE),
    .groups = "drop"
  )


# Plot massed and spaced separately so panels remain readable
p_massed_curves <- ggplot(
  new_data_joint %>% filter(Training == "massed"),
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_joint %>% filter(Training == "massed"),
    aes(x = stimulus, y = mean_max_peak, color = Genotype),
    inherit.aes = FALSE, alpha = 0.25, size = 0.6
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Peak distance moved, max_peak",
    title = "Joint Gamma GLMM: Massed training"
  ) +
  theme(legend.position = "none")

p_spaced_curves <- ggplot(
  new_data_joint %>% filter(Training == "spaced"),
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_joint %>% filter(Training == "spaced"),
    aes(x = stimulus, y = mean_max_peak, color = Genotype),
    inherit.aes = FALSE, alpha = 0.25, size = 0.6
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Peak distance moved, max_peak",
    title = "Joint Gamma GLMM: Spaced training"
  ) +
  theme(legend.position = "top")

print(p_spaced_curves)
print(p_massed_curves)

ggsave(
  file.path(save_fig_dir, "joint_glmm_max_peak_curves_massed.png"),
  p_massed_curves, width = 10, height = 12, dpi = 300, bg = "white"
)
ggsave(
  file.path(save_fig_dir, "joint_glmm_max_peak_curves_spaced.png"),
  p_spaced_curves, width = 14, height = 12, dpi = 300, bg = "white"
)

# Test-block curves only: massed vs spaced side by side for max_peak
test_curve_data <- new_data_joint %>%
  filter(
    (Training == "massed" & as.character(Block) == massed_test_block) |
      (Training == "spaced" & as.character(Block) == spaced_test_block)
  ) %>%
  mutate(
    TestBlock = ifelse(Training == "massed", "Massed test block", "Spaced test block"),
    TestBlock = factor(TestBlock, levels = c("Massed test block", "Spaced test block"))
  )

test_raw_data <- raw_summary_joint %>%
  filter(
    (Training == "massed" & as.character(Block) == massed_test_block) |
      (Training == "spaced" & as.character(Block) == spaced_test_block)
  ) %>%
  mutate(
    TestBlock = ifelse(Training == "massed", "Massed test block", "Spaced test block"),
    TestBlock = factor(TestBlock, levels = c("Massed test block", "Spaced test block"))
  )

p_test_massed_spaced_side_by_side <- ggplot(
  test_curve_data,
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ TestBlock, scales = "free_x") +
  geom_point(
    data = test_raw_data,
    aes(x = stimulus, y = mean_max_peak, color = Genotype),
    inherit.aes = FALSE,
    alpha = 0.25,
    size = 0.7
  ) +
  geom_ribbon(
    aes(ymin = CI_low, ymax = CI_high),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  scale_color_manual(values = genotype_colors, drop = FALSE) +
  scale_fill_manual(values = genotype_colors, drop = FALSE) +
  theme_pub(base_size = 13) +
  labs(
    x = "Stimulus number within test block",
    y = "Peak distance moved, max_peak",
    title = "Joint Gamma GLMM: test-block max_peak curves",
    subtitle = "Massed and spaced test blocks shown side by side"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

print(p_test_massed_spaced_side_by_side)

ggsave(
  file.path(save_fig_dir, "joint_glmm_test_block_massed_vs_spaced_side_by_side_max_peak.png"),
  p_test_massed_spaced_side_by_side,
  width = 9,
  height = 12,
  dpi = 300,
  bg = "white"
)

# ==============================================================================
# 7. Helper functions for max_peak contrasts with proper joint-covariance SEs
# ==============================================================================
#
# For Gamma(log), emmeans(..., type = "response") gives response-scale summaries,
# but contrasts may still be handled on the link scale unless we regrid.
#
# Therefore:
#   emm_resp <- regrid(emm, transform = "response")
#
# This makes max_peak the working linear scale, so contrast() estimates
# differences in predicted max_peak, with SEs based on the full covariance matrix.
# ==============================================================================

genotypes <- levels(df_all$Genotype)

get_response_est_col <- function(df) {
  candidate_cols <- c("response", "rate", "prob", "emmean", "estimate")
  found <- candidate_cols[candidate_cols %in% names(df)]
  if (length(found) == 0) {
    stop("Could not find response estimate column. Columns are: ",
         paste(names(df), collapse = ", "))
  }
  found[1]
}

get_lcl_col <- function(df) {
  candidate_cols <- c("asymp.LCL", "lower.CL", "LCL")
  found <- candidate_cols[candidate_cols %in% names(df)]
  if (length(found) == 0) {
    stop("Could not find lower CI column. Columns are: ",
         paste(names(df), collapse = ", "))
  }
  found[1]
}

get_ucl_col <- function(df) {
  candidate_cols <- c("asymp.UCL", "upper.CL", "UCL")
  found <- candidate_cols[candidate_cols %in% names(df)]
  if (length(found) == 0) {
    stop("Could not find upper CI column. Columns are: ",
         paste(names(df), collapse = ", "))
  }
  found[1]
}


get_test_emm_cells_peak <- function(model, stim_log_values) {
  
  emm <- emmeans(
    model,
    specs = ~ Genotype * Training * Block * stimulus_log,
    at = list(
      stimulus_log = stim_log_values,
      Block        = c(massed_test_block, spaced_test_block)
    ),
    type = "response"
  )
  
  emm_resp <- regrid(emm, transform = "response")
  emm_df   <- as.data.frame(emm_resp)
  
  keep <- (emm_df$Training == "massed" & as.character(emm_df$Block) == massed_test_block) |
    (emm_df$Training == "spaced" & as.character(emm_df$Block) == spaced_test_block)
  
  emm_resp[keep]
}


emm_peak_to_per_cell <- function(emm_obj) {
  emm_df <- as.data.frame(confint(emm_obj))
  
  est_col <- get_response_est_col(emm_df)
  lcl_col <- get_lcl_col(emm_df)
  ucl_col <- get_ucl_col(emm_df)
  
  emm_df %>%
    transmute(
      Genotype,
      Training,
      Block,
      stimulus_log,
      stimulus = exp(stimulus_log),
      max_peak = .data[[est_col]],
      SE       = SE,
      lower    = pmax(0, .data[[lcl_col]]),
      upper    = .data[[ucl_col]]
    )
}


proper_spaced_vs_massed_peak_contrast <- function(model, stim_log_values, label) {
  
  emm_sub <- get_test_emm_cells_peak(model, stim_log_values)
  emm_df  <- as.data.frame(emm_sub)
  
  genotypes_local <- levels(droplevels(emm_df$Genotype))
  k <- length(stim_log_values)
  n_cells <- nrow(emm_df)
  
  con_list <- lapply(genotypes_local, function(g) {
    w <- numeric(n_cells)
    
    spaced_idx <- which(emm_df$Genotype == g & emm_df$Training == "spaced")
    massed_idx <- which(emm_df$Genotype == g & emm_df$Training == "massed")
    
    if (length(spaced_idx) != k || length(massed_idx) != k) {
      stop(sprintf(
        "Genotype %s: expected %d spaced and %d massed cells, got %d and %d.",
        g, k, k, length(spaced_idx), length(massed_idx)
      ))
    }
    
    w[spaced_idx] <-  1 / k
    w[massed_idx] <- -1 / k
    w
  })
  
  names(con_list) <- genotypes_local
  
  con <- contrast(emm_sub, method = con_list)
  
  con_ci <- as.data.frame(confint(con))
  con_p  <- as.data.frame(con)
  
  est_col <- c("estimate", "Estimate")[c("estimate", "Estimate") %in% names(con_ci)][1]
  lo_col  <- c("asymp.LCL", "lower.CL", "LCL")[c("asymp.LCL", "lower.CL", "LCL") %in% names(con_ci)][1]
  hi_col  <- c("asymp.UCL", "upper.CL", "UCL")[c("asymp.UCL", "upper.CL", "UCL") %in% names(con_ci)][1]
  z_col   <- c("z.ratio", "t.ratio")[c("z.ratio", "t.ratio") %in% names(con_p)][1]
  p_col   <- c("p.value", "pvalue")[c("p.value", "pvalue") %in% names(con_p)][1]
  
  out <- tibble(
    Genotype  = factor(con_ci$contrast, levels = genotypes_local),
    diff_est  = con_ci[[est_col]],
    diff_se   = con_ci$SE,
    diff_low  = con_ci[[lo_col]],
    diff_high = con_ci[[hi_col]],
    z         = con_p[[z_col]],
    p         = con_p[[p_col]],
    contrast  = label
  )
  
  list(
    contrast = out,
    per_cell = emm_peak_to_per_cell(emm_sub)
  )
}


aggregate_per_cell_peak_for_plot <- function(per_cell_df) {
  per_cell_df %>%
    group_by(Genotype, Training, Block) %>%
    summarise(
      max_peak       = mean(max_peak),
      SE_descriptive = sqrt(sum(SE^2)) / n(),
      lower          = pmax(0, max_peak - 1.96 * SE_descriptive),
      upper          = max_peak + 1.96 * SE_descriptive,
      .groups        = "drop"
    )
}


finite_diff_gradient <- function(f, x, eps = 1e-6) {
  grad <- numeric(length(x))
  
  for (i in seq_along(x)) {
    x_hi <- x
    x_lo <- x
    
    step <- eps * max(1, abs(x[i]))
    
    x_hi[i] <- x_hi[i] + step
    x_lo[i] <- x_lo[i] - step
    
    grad[i] <- (f(x_hi) - f(x_lo)) / (2 * step)
  }
  
  grad
}


# ==============================================================================
# 8. Memory contrasts A, B, C with proper joint-covariance SEs
# ==============================================================================

# A: first test stimulus
res_A <- proper_spaced_vs_massed_peak_contrast(
  m_joint_peak,
  stim_log_values = 0,
  label = "A_stim1"
)

contrast_A <- res_A$contrast
per_cell_A <- res_A$per_cell


# B: mean over test stimuli 1:3
res_B <- proper_spaced_vs_massed_peak_contrast(
  m_joint_peak,
  stim_log_values = log(1:3),
  label = "B_mean_stim1_to_3"
)

contrast_B     <- res_B$contrast
per_cell_B_raw <- res_B$per_cell
per_cell_B     <- aggregate_per_cell_peak_for_plot(per_cell_B_raw)


# C: mean over all shared test stimuli 1:8
test_stims_for_C <- 1:8

res_C <- proper_spaced_vs_massed_peak_contrast(
  m_joint_peak,
  stim_log_values = log(test_stims_for_C),
  label = "C_mean_all_test_stim"
)

contrast_C     <- res_C$contrast
per_cell_C_raw <- res_C$per_cell
per_cell_C     <- aggregate_per_cell_peak_for_plot(per_cell_C_raw)


write.csv(
  per_cell_A,
  file.path(save_results_dir, "per_cell_test_max_peak_stim1.csv"),
  row.names = FALSE
)

write.csv(
  per_cell_B_raw,
  file.path(save_results_dir, "per_cell_test_max_peak_stim1to3.csv"),
  row.names = FALSE
)

write.csv(
  per_cell_B,
  file.path(save_results_dir, "per_cell_test_max_peak_stim1to3_aggregated.csv"),
  row.names = FALSE
)

write.csv(
  per_cell_C_raw,
  file.path(save_results_dir, "per_cell_test_max_peak_meanAllStim_raw.csv"),
  row.names = FALSE
)

write.csv(
  per_cell_C,
  file.path(save_results_dir, "per_cell_test_max_peak_meanAllStim.csv"),
  row.names = FALSE
)

write.csv(
  contrast_A,
  file.path(save_results_dir, "contrast_A_test_stim1_spaced_vs_massed_max_peak.csv"),
  row.names = FALSE
)

write.csv(
  contrast_B,
  file.path(save_results_dir, "contrast_B_test_meanStim1to3_spaced_vs_massed_max_peak.csv"),
  row.names = FALSE
)

write.csv(
  contrast_C,
  file.path(save_results_dir, "contrast_C_test_meanAllStim_spaced_vs_massed_max_peak.csv"),
  row.names = FALSE
)

cat("\n--- Contrast A: first test stimulus ---\n")
print(contrast_A)

cat("\n--- Contrast B: mean stim 1-3 ---\n")
print(contrast_B)

cat("\n--- Contrast C: mean stim 1-8 ---\n")
print(contrast_C)


# ==============================================================================
# 9. Plots for contrasts A, B, C
# ==============================================================================

p_per_cell_A <- per_cell_A %>%
  ggplot(aes(x = Training, y = max_peak,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Predicted max_peak at stim 1 of test block",
    title = "(A) Test peak distance at first stimulus: lower = better retention"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastA_per_cell_test_stim1_max_peak.png"),
  p_per_cell_A, width = 12, height = 4.5, dpi = 300, bg = "white"
)


p_contrast_A <- contrast_A %>%
  ggplot(aes(x = Genotype, y = diff_est,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x = NULL,
    y = "max_peak difference at stim 1: spaced - massed\nNegative = spaced has better retention",
    title = "(A) Headline contrast: stim 1 of test block"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastA_diff_spaced_vs_massed_max_peak.png"),
  p_contrast_A, width = 8, height = 5, dpi = 300, bg = "white"
)


p_per_cell_B <- per_cell_B %>%
  ggplot(aes(x = Training, y = max_peak,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Mean predicted max_peak over stim 1-3 of test block",
    title = "(B) Test peak distance averaged over first 3 stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastB_per_cell_test_stim1to3_max_peak.png"),
  p_per_cell_B, width = 12, height = 4.5, dpi = 300, bg = "white"
)


p_contrast_B <- contrast_B %>%
  ggplot(aes(x = Genotype, y = diff_est,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x = NULL,
    y = "Mean max_peak difference over stim 1-3: spaced - massed\nNegative = spaced has better retention",
    title = "(B) Averaged contrast: mean of stim 1-3"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastB_diff_spaced_vs_massed_max_peak.png"),
  p_contrast_B, width = 8, height = 5, dpi = 300, bg = "white"
)


p_per_cell_C <- per_cell_C %>%
  ggplot(aes(x = Training, y = max_peak,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Mean predicted max_peak over all shared test stimuli 1-8",
    title = "(C) Test peak distance averaged over all shared test stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastC_per_cell_test_meanAllStim_max_peak.png"),
  p_per_cell_C, width = 12, height = 4.5, dpi = 300, bg = "white"
)


p_contrast_C <- contrast_C %>%
  ggplot(aes(x = Genotype, y = diff_est,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x = NULL,
    y = "Mean max_peak difference over stim 1-8: spaced - massed\nNegative = spaced has better retention",
    title = "(C) Averaged contrast: mean of all shared test stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastC_diff_spaced_vs_massed_max_peak.png"),
  p_contrast_C, width = 8, height = 5, dpi = 300, bg = "white"
)


# ==============================================================================
# 10. Contrast D: inter-block recovery with joint-covariance SEs
# ==============================================================================

last_train_stim <- df_all %>%
  filter(BlockRole == "training") %>%
  group_by(Training) %>%
  summarise(
    last_block = as.character(max(as.integer(as.character(Block)))),
    last_stim  = max(stimulus, na.rm = TRUE),
    .groups = "drop"
  )

print(last_train_stim)

last_stim_massed  <- last_train_stim$last_stim[last_train_stim$Training == "massed"]
last_stim_spaced  <- last_train_stim$last_stim[last_train_stim$Training == "spaced"]
last_block_massed <- last_train_stim$last_block[last_train_stim$Training == "massed"]
last_block_spaced <- last_train_stim$last_block[last_train_stim$Training == "spaced"]


recovery_grid_small <- bind_rows(
  tibble(
    Training = "massed",
    point    = c("end_train", "start_test"),
    Block    = c(last_block_massed, massed_test_block),
    stimulus = c(last_stim_massed, 1)
  ),
  tibble(
    Training = "spaced",
    point    = c("end_train", "start_test"),
    Block    = c(last_block_spaced, spaced_test_block),
    stimulus = c(last_stim_spaced, 1)
  )
) %>%
  mutate(
    stimulus_log = log(stimulus),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block, levels = levels(df_all$Block))
  )


emm_D <- emmeans(
  m_joint_peak,
  specs = ~ Genotype * Training * Block * stimulus_log,
  at = list(
    stimulus_log = sort(unique(recovery_grid_small$stimulus_log)),
    Block        = unique(as.character(recovery_grid_small$Block))
  ),
  type = "response"
)

emm_D_resp <- regrid(emm_D, transform = "response")
emm_D_df   <- as.data.frame(emm_D_resp)

emm_D_df$point <- NA_character_

for (i in seq_len(nrow(recovery_grid_small))) {
  r <- recovery_grid_small[i, ]
  
  hit <- which(
    emm_D_df$Training == as.character(r$Training) &
      as.character(emm_D_df$Block) == as.character(r$Block) &
      abs(emm_D_df$stimulus_log - r$stimulus_log) < 1e-8
  )
  
  emm_D_df$point[hit] <- r$point
}

keep_D       <- !is.na(emm_D_df$point)
emm_D_sub    <- emm_D_resp[keep_D]
emm_D_df_sub <- emm_D_df[keep_D, ]

est_col_D_raw <- get_response_est_col(emm_D_df_sub)

per_cell_D <- emm_D_df_sub %>%
  as_tibble() %>%
  mutate(
    point    = emm_D_df_sub$point,
    max_peak = .data[[est_col_D_raw]],
    stimulus = exp(stimulus_log),
    lower    = pmax(0, max_peak - 1.96 * SE),
    upper    = max_peak + 1.96 * SE
  ) %>%
  select(Genotype, Training, Block, point, stimulus, max_peak, SE, lower, upper)

print(per_cell_D)

write.csv(
  per_cell_D,
  file.path(save_results_dir, "per_cell_endTrain_startTest_max_peak.csv"),
  row.names = FALSE
)


genotypes_D <- levels(droplevels(emm_D_df_sub$Genotype))
n_cells_D   <- nrow(emm_D_df_sub)

con_list_D_raw <- lapply(genotypes_D, function(g) {
  w <- numeric(n_cells_D)
  
  idx_ss <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "spaced" & emm_D_df_sub$point == "start_test")
  idx_se <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "spaced" & emm_D_df_sub$point == "end_train")
  idx_ms <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "massed" & emm_D_df_sub$point == "start_test")
  idx_me <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "massed" & emm_D_df_sub$point == "end_train")
  
  if (length(idx_ss) != 1 || length(idx_se) != 1 ||
      length(idx_ms) != 1 || length(idx_me) != 1) {
    stop(sprintf("Genotype %s: missing one of the four D reference cells.", g))
  }
  
  w[idx_ss] <-  1
  w[idx_se] <- -1
  w[idx_ms] <- -1
  w[idx_me] <-  1
  
  w
})

names(con_list_D_raw) <- genotypes_D

con_D_raw <- contrast(emm_D_sub, method = con_list_D_raw)

con_D_raw_ci <- as.data.frame(confint(con_D_raw))
con_D_raw_p  <- as.data.frame(con_D_raw)

est_col_D <- c("estimate", "Estimate")[c("estimate", "Estimate") %in% names(con_D_raw_ci)][1]
lo_col_D  <- c("asymp.LCL", "lower.CL", "LCL")[c("asymp.LCL", "lower.CL", "LCL") %in% names(con_D_raw_ci)][1]
hi_col_D  <- c("asymp.UCL", "upper.CL", "UCL")[c("asymp.UCL", "upper.CL", "UCL") %in% names(con_D_raw_ci)][1]
z_col_D   <- c("z.ratio", "t.ratio")[c("z.ratio", "t.ratio") %in% names(con_D_raw_p)][1]
p_col_D   <- c("p.value", "pvalue")[c("p.value", "pvalue") %in% names(con_D_raw_p)][1]

contrast_D_raw <- tibble(
  Genotype  = factor(con_D_raw_ci$contrast, levels = genotypes_D),
  diff_est  = con_D_raw_ci[[est_col_D]],
  diff_se   = con_D_raw_ci$SE,
  diff_low  = con_D_raw_ci[[lo_col_D]],
  diff_high = con_D_raw_ci[[hi_col_D]],
  z         = con_D_raw_p[[z_col_D]],
  p         = con_D_raw_p[[p_col_D]],
  contrast  = "D_raw_recovery"
)

print(contrast_D_raw)

write.csv(
  contrast_D_raw,
  file.path(save_results_dir, "contrast_D_raw_recovery_spaced_vs_massed_max_peak.csv"),
  row.names = FALSE
)


# Normalized recovery:
#   recovery_norm = (max_peak_start_test - max_peak_end_train) / max_peak_end_train
#                 = max_peak_start_test / max_peak_end_train - 1

beta_D <- as.numeric(emm_D_df_sub[[est_col_D_raw]])
V_D    <- vcov(emm_D_sub)

contrast_D_norm <- bind_rows(lapply(genotypes_D, function(g) {
  
  idx_ss <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "spaced" & emm_D_df_sub$point == "start_test")
  idx_se <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "spaced" & emm_D_df_sub$point == "end_train")
  idx_ms <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "massed" & emm_D_df_sub$point == "start_test")
  idx_me <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "massed" & emm_D_df_sub$point == "end_train")
  
  idx <- c(idx_ss, idx_se, idx_ms, idx_me)
  
  f_norm <- function(x) {
    x_ss <- x[1]
    x_se <- x[2]
    x_ms <- x[3]
    x_me <- x[4]
    
    rec_spaced <- (x_ss - x_se) / x_se
    rec_massed <- (x_ms - x_me) / x_me
    
    rec_spaced - rec_massed
  }
  
  x_hat <- beta_D[idx]
  V_hat <- V_D[idx, idx, drop = FALSE]
  
  est  <- f_norm(x_hat)
  grad <- finite_diff_gradient(f_norm, x_hat)
  se   <- sqrt(as.numeric(t(grad) %*% V_hat %*% grad))
  
  tibble(
    Genotype  = factor(g, levels = genotypes_D),
    diff_est  = est,
    diff_se   = se,
    diff_low  = est - 1.96 * se,
    diff_high = est + 1.96 * se,
    z         = est / se,
    p         = 2 * pnorm(-abs(est / se)),
    contrast  = "D_normalized_recovery"
  )
}))

print(contrast_D_norm)

write.csv(
  contrast_D_norm,
  file.path(save_results_dir, "contrast_D_normalized_recovery_spaced_vs_massed_max_peak.csv"),
  row.names = FALSE
)


recovery_per_cell <- per_cell_D %>%
  select(Genotype, Training, point, max_peak, SE) %>%
  pivot_wider(names_from = point, values_from = c(max_peak, SE)) %>%
  mutate(
    recovery       = max_peak_start_test - max_peak_end_train,
    recovery_SE    = sqrt(SE_start_test^2 + SE_end_train^2),
    recovery_low   = recovery - 1.96 * recovery_SE,
    recovery_high  = recovery + 1.96 * recovery_SE,
    
    recovery_norm      = recovery / max_peak_end_train,
    recovery_norm_SE   = recovery_SE / max_peak_end_train,
    recovery_norm_low  = recovery_norm - 1.96 * recovery_norm_SE,
    recovery_norm_high = recovery_norm + 1.96 * recovery_norm_SE
  )

write.csv(
  recovery_per_cell,
  file.path(save_results_dir, "recovery_per_cell_max_peak.csv"),
  row.names = FALSE
)


# ==============================================================================
# 11. Plots for Contrast D
# ==============================================================================

p_per_cell_D <- per_cell_D %>%
  mutate(
    point_label = factor(
      point,
      levels = c("end_train", "start_test"),
      labels = c("End of training", "Start of test")
    )
  ) %>%
  ggplot(aes(x = point_label, y = max_peak,
             ymin = lower, ymax = upper,
             color = Genotype, group = Training)) +
  facet_grid(Genotype ~ Training) +
  geom_pointrange(linewidth = 0.7, size = 0.6) +
  geom_line(linewidth = 0.6, linetype = "dashed", alpha = 0.7) +
  theme_pubr(base_size = 12) +
  labs(
    x = NULL,
    y = "Predicted max_peak",
    title = "End-of-training vs start-of-test peak distance per protocol"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(save_fig_dir, "contrastD_per_cell_endTrain_startTest_max_peak.png"),
  p_per_cell_D, width = 9, height = 11, dpi = 300, bg = "white"
)


p_recovery_per_cell <- recovery_per_cell %>%
  select(Genotype, Training, recovery, recovery_low, recovery_high) %>%
  ggplot(aes(x = Training, y = recovery,
             ymin = recovery_low, ymax = recovery_high,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Raw recovery = max_peak(start_test) - max_peak(end_train)",
    title = "(D) Inter-block max_peak recovery per protocol: lower = better memory"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_raw_recovery_per_cell_max_peak.png"),
  p_recovery_per_cell, width = 10, height = 4.5, dpi = 300, bg = "white"
)


p_recovery_norm_per_cell <- recovery_per_cell %>%
  select(Genotype, Training, recovery_norm, recovery_norm_low, recovery_norm_high) %>%
  ggplot(aes(x = Training, y = recovery_norm,
             ymin = recovery_norm_low, ymax = recovery_norm_high,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Normalized recovery = recovery / max_peak(end_train)",
    title = "(D-norm) Normalized inter-block max_peak recovery per protocol"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_normalized_recovery_per_cell_max_peak.png"),
  p_recovery_norm_per_cell, width = 10, height = 4.5, dpi = 300, bg = "white"
)


p_contrast_D_raw <- contrast_D_raw %>%
  ggplot(aes(x = Genotype, y = diff_est,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x = NULL,
    y = "Raw max_peak recovery difference: spaced - massed\nNegative = spaced recovers less = better retention",
    title = "(D) Recovery contrast: raw max_peak"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_raw_diff_spaced_vs_massed_max_peak.png"),
  p_contrast_D_raw, width = 8, height = 5, dpi = 300, bg = "white"
)


p_contrast_D_norm <- contrast_D_norm %>%
  ggplot(aes(x = Genotype, y = diff_est,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x = NULL,
    y = "Normalized max_peak recovery difference: spaced - massed\nNegative = spaced recovers less = better retention",
    title = "(D-norm) Recovery contrast: normalized max_peak"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_norm_diff_spaced_vs_massed_max_peak.png"),
  p_contrast_D_norm, width = 8, height = 5, dpi = 300, bg = "white"
)


# ==============================================================================
# 12. Combined contrast table and plot
# ==============================================================================

all_contrasts <- bind_rows(
  contrast_A,
  contrast_B,
  contrast_C,
  contrast_D_raw,
  contrast_D_norm
) %>%
  mutate(
    contrast = factor(
      contrast,
      levels = c(
        "A_stim1",
        "B_mean_stim1_to_3",
        "C_mean_all_test_stim",
        "D_raw_recovery",
        "D_normalized_recovery"
      )
    )
  )

write.csv(
  all_contrasts,
  file.path(save_results_dir, "ALL_contrasts_spaced_vs_massed_max_peak.csv"),
  row.names = FALSE
)

p_all_contrasts <- all_contrasts %>%
  ggplot(aes(x = Genotype, y = diff_est,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y") +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Spaced - massed max_peak difference\nNegative = spaced has better retention",
    title = "All max_peak memory contrasts"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "ALL_contrasts_spaced_vs_massed_max_peak.png"),
  p_all_contrasts, width = 9, height = 12, dpi = 300, bg = "white"
)

cat("\n--- All max_peak contrasts ---\n")
print(all_contrasts)


# ==============================================================================
# 13. Supplementary contrasts
# ==============================================================================

# ------------------------------------------------------------------------------
# 13.1 Training sanity check: start vs end of Block 1
# ------------------------------------------------------------------------------

last_b1_stim <- df_all %>%
  filter(Block == "1") %>%
  group_by(Training) %>%
  summarise(last_stim = max(stimulus, na.rm = TRUE), .groups = "drop")

last_b1_massed <- last_b1_stim$last_stim[last_b1_stim$Training == "massed"]
last_b1_spaced <- last_b1_stim$last_stim[last_b1_stim$Training == "spaced"]

training_grid_small <- tibble(
  Training = rep(c("massed", "spaced"), each = 2),
  point    = rep(c("start_b1", "end_b1"), times = 2),
  Block    = "1",
  stimulus = c(1, last_b1_massed, 1, last_b1_spaced)
) %>%
  mutate(
    stimulus_log = log(stimulus),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block, levels = levels(df_all$Block))
  )


emm_S <- emmeans(
  m_joint_peak,
  specs = ~ Genotype * Training * Block * stimulus_log,
  at = list(
    stimulus_log = sort(unique(training_grid_small$stimulus_log)),
    Block = "1"
  ),
  type = "response"
)

emm_S_resp <- regrid(emm_S, transform = "response")
emm_S_df   <- as.data.frame(emm_S_resp)

emm_S_df$point <- NA_character_

for (i in seq_len(nrow(training_grid_small))) {
  r <- training_grid_small[i, ]
  
  hit <- which(
    emm_S_df$Training == as.character(r$Training) &
      as.character(emm_S_df$Block) == as.character(r$Block) &
      abs(emm_S_df$stimulus_log - r$stimulus_log) < 1e-8
  )
  
  emm_S_df$point[hit] <- r$point
}

keep_S       <- !is.na(emm_S_df$point)
emm_S_sub    <- emm_S_resp[keep_S]
emm_S_df_sub <- emm_S_df[keep_S, ]

est_col_S_raw <- get_response_est_col(emm_S_df_sub)

training_per_cell <- emm_S_df_sub %>%
  as_tibble() %>%
  mutate(
    point    = emm_S_df_sub$point,
    max_peak = .data[[est_col_S_raw]],
    stimulus = exp(stimulus_log),
    lower    = pmax(0, max_peak - 1.96 * SE),
    upper    = max_peak + 1.96 * SE
  ) %>%
  select(Genotype, Training, Block, point, stimulus, max_peak, SE, lower, upper)


combos_S <- training_per_cell %>%
  distinct(Genotype, Training) %>%
  arrange(Genotype, Training)

n_cells_S <- nrow(emm_S_df_sub)

con_list_S <- list()

for (i in seq_len(nrow(combos_S))) {
  g <- as.character(combos_S$Genotype[i])
  t <- as.character(combos_S$Training[i])
  
  idx_start <- which(emm_S_df_sub$Genotype == g & emm_S_df_sub$Training == t & emm_S_df_sub$point == "start_b1")
  idx_end   <- which(emm_S_df_sub$Genotype == g & emm_S_df_sub$Training == t & emm_S_df_sub$point == "end_b1")
  
  if (length(idx_start) != 1 || length(idx_end) != 1) {
    stop(sprintf("Genotype %s, Training %s: missing start_b1 or end_b1.", g, t))
  }
  
  w <- numeric(n_cells_S)
  w[idx_end]   <-  1
  w[idx_start] <- -1
  
  con_list_S[[paste(g, t, sep = "|")]] <- w
}

con_S <- contrast(emm_S_sub, method = con_list_S)

con_S_ci <- as.data.frame(confint(con_S))
con_S_p  <- as.data.frame(con_S)

est_col_S <- c("estimate", "Estimate")[c("estimate", "Estimate") %in% names(con_S_ci)][1]
lo_col_S  <- c("asymp.LCL", "lower.CL", "LCL")[c("asymp.LCL", "lower.CL", "LCL") %in% names(con_S_ci)][1]
hi_col_S  <- c("asymp.UCL", "upper.CL", "UCL")[c("asymp.UCL", "upper.CL", "UCL") %in% names(con_S_ci)][1]
z_col_S   <- c("z.ratio", "t.ratio")[c("z.ratio", "t.ratio") %in% names(con_S_p)][1]
p_col_S   <- c("p.value", "pvalue")[c("p.value", "pvalue") %in% names(con_S_p)][1]

training_learning <- con_S_ci %>%
  as_tibble() %>%
  mutate(
    contrast      = as.character(contrast),
    Genotype      = sapply(strsplit(contrast, "\\|"), `[`, 1),
    Training      = sapply(strsplit(contrast, "\\|"), `[`, 2),
    learning_est  = .data[[est_col_S]],
    learning_se   = SE,
    learning_low  = .data[[lo_col_S]],
    learning_high = .data[[hi_col_S]]
  ) %>%
  left_join(
    con_S_p %>%
      as_tibble() %>%
      mutate(contrast = as.character(contrast)) %>%
      select(contrast, z = all_of(z_col_S), p = all_of(p_col_S)),
    by = "contrast"
  ) %>%
  select(Genotype, Training, learning_est, learning_se,
         learning_low, learning_high, z, p)

print(training_learning)

write.csv(
  training_learning,
  file.path(save_results_dir, "training_sanity_check_block1_learning_max_peak.csv"),
  row.names = FALSE
)

p_training_sanity <- training_learning %>%
  ggplot(aes(x = Training, y = learning_est,
             ymin = learning_low, ymax = learning_high,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "max_peak(end B1) - max_peak(start B1)\nNegative = peak movement decreased",
    title = "Training sanity check: within-Block-1 change in max_peak"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "sanity_check_block1_learning_max_peak.png"),
  p_training_sanity, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# ------------------------------------------------------------------------------
# 13.2 Within-training max_peak slope contrast
# ------------------------------------------------------------------------------

slope_emm <- emtrends(
  m_joint_peak,
  specs = ~ Genotype * Training | Block,
  var   = "stimulus_log",
  at    = list(Block = "1")
)

slope_df <- as_tibble(confint(slope_emm)) %>%
  rename(
    slope = stimulus_log.trend,
    SE    = SE,
    lower = asymp.LCL,
    upper = asymp.UCL
  ) %>%
  select(Genotype, Training, Block, slope, SE, lower, upper)

print(slope_df)

write.csv(
  slope_df,
  file.path(save_results_dir, "training_slope_block1_max_peak.csv"),
  row.names = FALSE
)

slope_pairs <- pairs(slope_emm, by = "Genotype", reverse = TRUE)

slope_pairs_df <- as_tibble(confint(slope_pairs)) %>%
  rename(
    diff_est  = estimate,
    diff_se   = SE,
    diff_low  = asymp.LCL,
    diff_high = asymp.UCL
  ) %>%
  left_join(
    as_tibble(slope_pairs) %>%
      select(Genotype, contrast, z.ratio, p.value),
    by = c("Genotype", "contrast")
  ) %>%
  select(Genotype, contrast, diff_est, diff_se, diff_low, diff_high,
         z.ratio, p.value)

print(slope_pairs_df)

write.csv(
  slope_pairs_df,
  file.path(save_results_dir, "training_slope_diff_spaced_vs_massed_max_peak.csv"),
  row.names = FALSE
)

p_slope <- slope_df %>%
  ggplot(aes(x = Training, y = slope,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Habituation slope d(log max_peak)/d(log stim) in Block 1\nMore negative = faster reduction",
    title = "Within-Block-1 max_peak habituation rate per protocol"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "training_slope_block1_per_protocol_max_peak.png"),
  p_slope, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# ------------------------------------------------------------------------------
# 13.3 Block-by-block learning curve, spaced only
# ------------------------------------------------------------------------------

spaced_last_stim_per_block <- df_all %>%
  filter(Training == "spaced", Block %in% c("1", "2", "3", "4")) %>%
  group_by(Block) %>%
  summarise(last_stim = max(stimulus, na.rm = TRUE), .groups = "drop")

print(spaced_last_stim_per_block)

spaced_endblock_grid <- tidyr::expand_grid(
  Genotype = levels(df_all$Genotype),
  Block    = c("1", "2", "3", "4")
) %>%
  left_join(
    spaced_last_stim_per_block %>% mutate(Block = as.character(Block)),
    by = "Block"
  ) %>%
  mutate(
    Training     = "spaced",
    stimulus     = last_stim,
    stimulus_log = log(stimulus),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block, levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  )


emm_endblock <- emmeans(
  m_joint_peak,
  specs = ~ Genotype * Training * Block * stimulus_log,
  at = list(
    Training = "spaced",
    Block = c("1", "2", "3", "4"),
    stimulus_log = sort(unique(spaced_endblock_grid$stimulus_log))
  ),
  type = "response"
)

emm_endblock_resp <- regrid(emm_endblock, transform = "response")
emm_endblock_df <- as.data.frame(emm_endblock_resp)

emm_endblock_df$keep <- FALSE

for (i in seq_len(nrow(spaced_endblock_grid))) {
  r <- spaced_endblock_grid[i, ]
  
  hit <- which(
    emm_endblock_df$Genotype == as.character(r$Genotype) &
      emm_endblock_df$Training == "spaced" &
      as.character(emm_endblock_df$Block) == as.character(r$Block) &
      abs(emm_endblock_df$stimulus_log - r$stimulus_log) < 1e-8
  )
  
  emm_endblock_df$keep[hit] <- TRUE
}

emm_endblock_sub <- emm_endblock_resp[emm_endblock_df$keep]
emm_endblock_df_sub <- emm_endblock_df[emm_endblock_df$keep, ]

est_col_end <- get_response_est_col(emm_endblock_df_sub)

spaced_endblock_result <- emm_endblock_df_sub %>%
  as_tibble() %>%
  mutate(
    max_peak = .data[[est_col_end]],
    stimulus = exp(stimulus_log),
    block_num = as.integer(as.character(Block)),
    lower = pmax(0, max_peak - 1.96 * SE),
    upper = max_peak + 1.96 * SE
  ) %>%
  select(Genotype, Block, block_num, stimulus, max_peak, SE, lower, upper)

print(spaced_endblock_result)

write.csv(
  spaced_endblock_result,
  file.path(save_results_dir, "spaced_endblock_max_peak.csv"),
  row.names = FALSE
)


emm_end_df <- as.data.frame(emm_endblock_sub)
genotypes_end <- levels(droplevels(emm_end_df$Genotype))

con_list_end <- list()

for (g in genotypes_end) {
  for (transition in c("B1->B2", "B2->B3", "B3->B4")) {
    
    b_from <- str_extract(transition, "(?<=B)[1-4]")
    b_to   <- str_extract(transition, "(?<=->B)[1-4]")
    
    idx_from <- which(emm_end_df$Genotype == g & as.character(emm_end_df$Block) == b_from)
    idx_to   <- which(emm_end_df$Genotype == g & as.character(emm_end_df$Block) == b_to)
    
    if (length(idx_from) == 1 && length(idx_to) == 1) {
      w <- numeric(nrow(emm_end_df))
      w[idx_to]   <-  1
      w[idx_from] <- -1
      
      con_list_end[[paste(g, transition, sep = "|")]] <- w
    }
  }
}

con_end <- contrast(emm_endblock_sub, method = con_list_end)

con_end_ci <- as.data.frame(confint(con_end))

spaced_block_diffs <- con_end_ci %>%
  as_tibble() %>%
  mutate(
    label = as.character(contrast),
    Genotype = sapply(strsplit(label, "\\|"), `[`, 1),
    transition = sapply(strsplit(label, "\\|"), `[`, 2),
    delta = estimate,
    delta_se = SE,
    delta_low = asymp.LCL,
    delta_high = asymp.UCL
  ) %>%
  select(Genotype, transition, delta, delta_se, delta_low, delta_high)

print(spaced_block_diffs)

write.csv(
  spaced_block_diffs,
  file.path(save_results_dir, "spaced_block_to_block_endblock_changes_max_peak.csv"),
  row.names = FALSE
)


p_spaced_endblock <- spaced_endblock_result %>%
  ggplot(aes(x = block_num, y = max_peak,
             ymin = lower, ymax = upper,
             color = Genotype, group = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.7, size = 0.6) +
  geom_line(linewidth = 0.6) +
  scale_x_continuous(breaks = 1:4) +
  theme_pubr(base_size = 13) +
  labs(
    x = "Spaced training block",
    y = "Predicted max_peak at end of block",
    title = "Spaced training: end-of-block max_peak across blocks 1-4"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "spaced_endblock_progression_max_peak.png"),
  p_spaced_endblock, width = 10, height = 4.5, dpi = 300, bg = "white"
)


p_spaced_diffs <- spaced_block_diffs %>%
  ggplot(aes(x = transition, y = delta,
             ymin = delta_low, ymax = delta_high,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.7, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x = "Block-to-block transition",
    y = "Change in end-of-block max_peak\nNegative = further reduction",
    title = "Spaced training: incremental max_peak reduction per block"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(save_fig_dir, "spaced_block_to_block_changes_max_peak.png"),
  p_spaced_diffs, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# ==============================================================================
# Done
# ==============================================================================

message("Joint Gamma GLMM max_peak analysis complete.")
message("Results saved under: ", save_results_dir)
message("Figures saved under: ", save_fig_dir)