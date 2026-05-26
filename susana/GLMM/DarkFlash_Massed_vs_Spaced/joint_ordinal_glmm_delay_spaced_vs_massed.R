###############################################################################
# Joint Ordinal GLMM Analysis: Spaced vs Massed Training - Response Delay
# Author: Nils Brehm
# Date: 2026
#
# Description:
#   Joint ordinal mixed model combining spaced and massed training datasets into
#   one model to test the spaced-vs-massed memory contrast directly for:
#
#     delay = response delay category
#
#   delay can only take ordered values: 0, 1, 2, 3, 4
#
#   Model:
#     ordered(delay) ~ Genotype * Training * Block * stimulus_log + (1 | animal)
#
#   Model family:
#     Cumulative link mixed model, logit link
#
#   Main readout:
#     Expected delay:
#       E(delay) = 0*P(delay=0) + 1*P(delay=1) + ... + 4*P(delay=4)
#
#   Memory contrasts per genotype:
#     (A) expected delay at first test stimulus, spaced vs massed
#     (B) mean expected delay over first 3 test stimuli, spaced vs massed
#     (C) mean expected delay over all shared test stimuli, spaced vs massed
#     (D) inter-block change in expected delay, spaced vs massed
#
#   Positive spaced-massed difference => spaced fish have longer/slower delay
#   in the test block than massed fish.
###############################################################################

# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
library(readr)
library(ordinal)      # clmm(): cumulative link mixed model
library(emmeans)
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

save_results_dir <- file.path(base_dir, "results", "ordinal_joint_delay")
save_fig_dir     <- file.path(base_dir, "figs",    "ordinal_joint_delay")

dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. Load and prepare both datasets
# ==============================================================================
# Delay is only meaningful for responses. This follows your old delay analysis,
# which used df_final_sub. The fallback filters to valid delay values.
res_massed <- load_data(file_massed, move_th = 1, drop = c("th2, tyr", "tyr"))
res_spaced <- load_data(file_spaced, move_th = 1, drop = c("th2, tyr", "tyr"))

df_massed <- res_massed$df_final_sub
df_spaced <- res_spaced$df_final_sub

if (is.null(df_massed)) {
  df_massed <- res_massed$df_final
}
if (is.null(df_spaced)) {
  df_spaced <- res_spaced$df_final
}

# Keep valid ordinal delay values only.
# Use numeric first, then ordered factor for the ordinal model.
df_massed <- df_massed %>%
  mutate(delay_num = as.integer(as.character(delay))) %>%
  filter(!is.na(delay_num), delay_num %in% 0:4)

df_spaced <- df_spaced %>%
  mutate(delay_num = as.integer(as.character(delay))) %>%
  filter(!is.na(delay_num), delay_num %in% 0:4)


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
    animal    = factor(paste0(Training, "_", Video, ".", Well)),
    delay_ord = ordered(delay_num, levels = 0:4)
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

cat("\n--- Delay counts by Training x Block ---\n")
print(
  df_all %>%
    count(Training, Block, BlockRole, delay_num) %>%
    arrange(Training, Block, delay_num)
)

cat("\n--- Mean delay by Training x Block ---\n")
print(
  df_all %>%
    group_by(Training, Block, BlockRole) %>%
    summarise(
      n = n(),
      mean_delay   = mean(delay_num, na.rm = TRUE),
      median_delay = median(delay_num, na.rm = TRUE),
      .groups = "drop"
    )
)


# ==============================================================================
# 3. Exploratory distribution plot
# ==============================================================================
message("Plotting delay distributions...")

p_dist <- ggplot(df_all, aes(x = delay_ord)) +
  geom_bar(fill = "skyblue", color = "black") +
  facet_grid(Training ~ Genotype, scales = "free_y") +
  theme_pubr(base_size = 12) +
  labs(
    title = "Response delay distribution: joint spaced + massed dataset",
    x = "Response delay category",
    y = "Count"
  )

print(p_dist)

ggsave(
  file.path(save_fig_dir, "joint_delay_distribution.png"),
  p_dist, width = 12, height = 6, dpi = 300, bg = "white"
)


# ==============================================================================
# 4. Fit the joint ordinal mixed model
# ==============================================================================
message("Fitting joint ordinal mixed model for delay...")

# clmm() fits a cumulative link mixed model.
# Hess = TRUE is important for SEs, summary(), and emmeans().
if (from_file){
  # Load fitted model if available
  m_joint_delay  <- readRDS(
    file.path(save_results_dir, "joint_ordinal_glmm_spaced_vs_massed_delay.rds")
  )
  
} else {
  m_joint_delay <- clmm(
    delay_ord ~ Genotype * Training * Block * stimulus_log + (1 | animal),
    data = df_all,
    link = "logit",
    Hess = TRUE,
    nAGQ = 1
  )
  
  saveRDS(
    m_joint_delay,
    file.path(save_results_dir, "joint_ordinal_glmm_spaced_vs_massed_delay.rds")
  )
}


# ==============================================================================
# Load fitted model if available
# ==============================================================================
m_joint_delay  <- readRDS(
  file.path(save_results_dir, "joint_ordinal_glmm_spaced_vs_massed_delay.rds")
)

capture.output(
  summary(m_joint_delay),
  file = file.path(save_results_dir, "summary_results_delay.txt")
)

print(summary(m_joint_delay))


# ==============================================================================
# 5. Basic model checks
# ==============================================================================
# DHARMa is excellent for many GLMMs, but ordinal::clmm is not always supported
# by custom validation helpers designed for glmmTMB. For this reason, we save:
#   - model summary
#   - convergence information
#   - fitted threshold/fixed-effect summary
#   - observed vs predicted delay summaries
#
# If your validate_model() supports ordinal::clmm, you can uncomment this:
# validate_model(m_joint_delay, df_all)

cat("\n--- clmm convergence info ---\n")
print(m_joint_delay$convergence)
print(m_joint_delay$optRes)

# ==============================================================================
# 6. Helpers: proper joint-covariance contrasts via emmeans
# ==============================================================================
#
# Per-stimulus contrasts are NOT independent (same model, same fish, shared
# fixed-effect parameters). Aggregating them with sqrt(sum(SE^2))/n gives
# conservative (too-wide) CIs. The fix is to build each averaged contrast as
# a single linear combination of emmeans cells, so the SE uses the full
# covariance matrix of the fixed effects via emmeans::contrast().
#
# For ordinal::clmm, emmeans with mode = "mean.class" returns expected class
# (1-indexed: 1..5 for delay levels 0..4). We subtract 1 to map back to the
# 0..4 delay scale. The SE is invariant to this shift.
# ==============================================================================

delay_scores <- 0:4

# Direct prediction of expected delay from a clmm model (used for habituation
# curves in section 7). No SE, just point predictions.
predict_expected_delay <- function(model, newdata) {
  beta  <- model$beta
  alpha <- model$alpha
  
  Terms <- delete.response(terms(model))
  X <- model.matrix(Terms, newdata)
  X <- X[, names(beta), drop = FALSE]
  
  eta <- as.vector(X %*% beta)
  
  cumprob <- sapply(alpha, function(a) plogis(a - eta))
  
  probs <- cbind(
    cumprob[, 1],
    cumprob[, 2] - cumprob[, 1],
    cumprob[, 3] - cumprob[, 2],
    cumprob[, 4] - cumprob[, 3],
    1 - cumprob[, 4]
  )
  
  as.vector(probs %*% delay_scores)
}


# Helper: get emmeans cells for the test block of each protocol at a set of
# test stimuli, restricted to the matching (Training, Block) combination.
get_test_emm_cells <- function(model, stim_log_values) {
  
  emm <- emmeans(
    model,
    specs = ~ Genotype * Training * Block * stimulus_log,
    at = list(
      stimulus_log = stim_log_values,
      Block        = c(massed_test_block, spaced_test_block)
    ),
    mode = "mean.class"
  )
  
  emm_df <- as.data.frame(emm)
  keep <- (emm_df$Training == "massed" & as.character(emm_df$Block) == massed_test_block) |
    (emm_df$Training == "spaced" & as.character(emm_df$Block) == spaced_test_block)
  emm[keep]
}


# Per-cell summary (estimate + CI) on the 0..4 expected-delay scale.
emm_to_per_cell <- function(emm_obj) {
  emm_df <- as.data.frame(confint(emm_obj))
  
  estimate_col <- c("emmean", "mean.class", "estimate")[
    c("emmean", "mean.class", "estimate") %in% names(emm_df)
  ][1]
  lower_col <- c("asymp.LCL", "lower.CL", "LCL")[
    c("asymp.LCL", "lower.CL", "LCL") %in% names(emm_df)
  ][1]
  upper_col <- c("asymp.UCL", "upper.CL", "UCL")[
    c("asymp.UCL", "upper.CL", "UCL") %in% names(emm_df)
  ][1]
  
  emm_df %>%
    transmute(
      Genotype,
      Training,
      Block,
      stimulus_log,
      stimulus       = exp(stimulus_log),
      expected_delay = .data[[estimate_col]] - 1,
      SE             = SE,
      lower          = .data[[lower_col]] - 1,
      upper          = .data[[upper_col]] - 1
    )
}


# Core: proper joint-covariance "spaced - massed averaged over stimuli" per Genotype.
proper_spaced_vs_massed_contrast <- function(model, stim_log_values, label) {
  
  emm_sub <- get_test_emm_cells(model, stim_log_values)
  emm_df  <- as.data.frame(emm_sub)
  n_cells <- nrow(emm_df)
  
  genotypes <- levels(droplevels(emm_df$Genotype))
  k <- length(stim_log_values)
  
  con_list <- lapply(genotypes, function(g) {
    w <- numeric(n_cells)
    spaced_idx <- which(emm_df$Genotype == g & emm_df$Training == "spaced")
    massed_idx <- which(emm_df$Genotype == g & emm_df$Training == "massed")
    
    if (length(spaced_idx) != k || length(massed_idx) != k) {
      stop(sprintf(
        "Genotype %s: expected %d spaced and %d massed cells, got %d and %d. Check Block filtering.",
        g, k, k, length(spaced_idx), length(massed_idx)
      ))
    }
    
    w[spaced_idx] <-  1 / k
    w[massed_idx] <- -1 / k
    w
  })
  names(con_list) <- genotypes
  
  con <- contrast(emm_sub, method = con_list)
  con_ci <- as.data.frame(confint(con))
  con_p  <- as.data.frame(con)
  
  est_col <- c("estimate", "Estimate")[c("estimate", "Estimate") %in% names(con_ci)][1]
  lo_col  <- c("asymp.LCL", "lower.CL", "LCL")[c("asymp.LCL", "lower.CL", "LCL") %in% names(con_ci)][1]
  hi_col  <- c("asymp.UCL", "upper.CL", "UCL")[c("asymp.UCL", "upper.CL", "UCL") %in% names(con_ci)][1]
  p_col   <- c("p.value", "pvalue")[c("p.value", "pvalue") %in% names(con_p)][1]
  z_col   <- c("z.ratio", "t.ratio")[c("z.ratio", "t.ratio") %in% names(con_p)][1]
  
  out <- tibble(
    Genotype  = factor(con_ci$contrast, levels = genotypes),
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
    per_cell = emm_to_per_cell(emm_sub)
  )
}


# Descriptive per-cell aggregation across stimuli (for plotting only).
aggregate_per_cell_for_plot <- function(per_cell_df) {
  per_cell_df %>%
    group_by(Genotype, Training, Block) %>%
    summarise(
      expected_delay = mean(expected_delay),
      SE_descriptive = sqrt(sum(SE^2)) / n(),
      lower          = expected_delay - 1.96 * SE_descriptive,
      upper          = expected_delay + 1.96 * SE_descriptive,
      .groups = "drop"
    )
}


# ==============================================================================
# 7. Plot habituation curves: expected delay across stimuli
# ==============================================================================
message("Plotting joint expected-delay habituation curves...")

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

new_data_joint <- new_data_joint %>%
  mutate(fit_expected_delay = predict_expected_delay(m_joint_delay, new_data_joint))

raw_summary_joint <- df_all %>%
  group_by(Training, Block, Genotype, stimulus) %>%
  summarise(
    mean_delay = mean(delay_num, na.rm = TRUE),
    .groups = "drop"
  )


p_massed_curves <- ggplot(
  new_data_joint %>% filter(Training == "massed"),
  aes(x = stimulus, color = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_joint %>% filter(Training == "massed"),
    aes(x = stimulus, y = mean_delay, color = Genotype),
    inherit.aes = FALSE, alpha = 0.25, size = 0.6
  ) +
  geom_line(aes(y = fit_expected_delay), linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Expected response delay",
    title = "Joint ordinal GLMM: Massed training"
  ) +
  theme(legend.position = "none")

p_spaced_curves <- ggplot(
  new_data_joint %>% filter(Training == "spaced"),
  aes(x = stimulus, color = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_joint %>% filter(Training == "spaced"),
    aes(x = stimulus, y = mean_delay, color = Genotype),
    inherit.aes = FALSE, alpha = 0.25, size = 0.6
  ) +
  geom_line(aes(y = fit_expected_delay), linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Expected response delay",
    title = "Joint ordinal GLMM: Spaced training"
  ) +
  theme(legend.position = "top")

print(p_spaced_curves)
print(p_massed_curves)

ggsave(
  file.path(save_fig_dir, "joint_ordinal_delay_curves_massed.png"),
  p_massed_curves, width = 10, height = 12, dpi = 300, bg = "white"
)
ggsave(
  file.path(save_fig_dir, "joint_ordinal_delay_curves_spaced.png"),
  p_spaced_curves, width = 14, height = 12, dpi = 300, bg = "white"
)

# Test-block curves only: massed vs spaced side by side for expected response delay

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

# Pick the raw-summary delay column automatically
delay_raw_col <- c("mean_expected_delay", "mean_delay", "expected_delay", "delay")[
  c("mean_expected_delay", "mean_delay", "expected_delay", "delay") %in% names(test_raw_data)
][1]

if (is.na(delay_raw_col)) {
  stop(
    "Could not find delay column in raw_summary_joint. Expected one of: ",
    "mean_expected_delay, mean_delay, expected_delay, delay"
  )
}

p_test_massed_spaced_side_by_side_delay <- ggplot(
  test_curve_data,
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ TestBlock, scales = "free_x") +
  geom_point(
    data = test_raw_data,
    aes(x = stimulus, y = .data[[delay_raw_col]], color = Genotype),
    inherit.aes = FALSE,
    alpha = 0.25,
    size = 0.7
  ) +
  geom_line(aes(y = fit_expected_delay), linewidth = 1.1) +
  scale_color_manual(values = genotype_colors, drop = FALSE) +
  scale_fill_manual(values = genotype_colors, drop = FALSE) +
  theme_pub(base_size = 13) +
  labs(
    x = "Stimulus number within test block",
    y = "Expected response delay (s)",
    title = "Ordinal mixed model: test-block response-delay curves",
    subtitle = "Massed and spaced test blocks shown side by side"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

print(p_test_massed_spaced_side_by_side_delay)

ggsave(
  file.path(save_fig_dir, "joint_ordinal_test_block_massed_vs_spaced_side_by_side_expected_delay.png"),
  p_test_massed_spaced_side_by_side_delay,
  width = 9,
  height = 12,
  dpi = 300,
  bg = "white"
)

# ==============================================================================
# 8. Memory contrasts A, B, C with proper joint-covariance SEs
# ==============================================================================

# --- A: first test stimulus -------------------------------------------------
res_A <- proper_spaced_vs_massed_contrast(
  m_joint_delay,
  stim_log_values = 0,                # log(1) = 0
  label           = "A_stim1"
)
contrast_A <- res_A$contrast
per_cell_A <- res_A$per_cell

# --- B: mean over test stimuli 1..3 -----------------------------------------
res_B <- proper_spaced_vs_massed_contrast(
  m_joint_delay,
  stim_log_values = log(1:3),
  label           = "B_mean_stim1_to_3"
)
contrast_B     <- res_B$contrast
per_cell_B_raw <- res_B$per_cell
per_cell_B     <- aggregate_per_cell_for_plot(per_cell_B_raw)

# --- C: mean over all shared test stimuli (1..8) ----------------------------
test_stims_for_C <- 1:8

res_C <- proper_spaced_vs_massed_contrast(
  m_joint_delay,
  stim_log_values = log(test_stims_for_C),
  label           = "C_mean_all_test_stim"
)
contrast_C     <- res_C$contrast
per_cell_C_raw <- res_C$per_cell
per_cell_C     <- aggregate_per_cell_for_plot(per_cell_C_raw)


# --- Save per-cell and contrast tables --------------------------------------
write.csv(
  per_cell_A,
  file.path(save_results_dir, "per_cell_test_expected_delay_stim1.csv"),
  row.names = FALSE
)
write.csv(
  per_cell_B_raw,
  file.path(save_results_dir, "per_cell_test_expected_delay_stim1to3.csv"),
  row.names = FALSE
)
write.csv(
  per_cell_B,
  file.path(save_results_dir, "per_cell_test_expected_delay_stim1to3_aggregated.csv"),
  row.names = FALSE
)
write.csv(
  per_cell_C,
  file.path(save_results_dir, "per_cell_test_expected_delay_meanAllStim.csv"),
  row.names = FALSE
)

write.csv(
  contrast_A,
  file.path(save_results_dir, "contrast_A_test_stim1_spaced_vs_massed_expected_delay.csv"),
  row.names = FALSE
)
write.csv(
  contrast_B,
  file.path(save_results_dir, "contrast_B_test_meanStim1to3_spaced_vs_massed_expected_delay.csv"),
  row.names = FALSE
)
write.csv(
  contrast_C,
  file.path(save_results_dir, "contrast_C_test_meanAllStim_spaced_vs_massed_expected_delay.csv"),
  row.names = FALSE
)

cat("\n--- Contrast A: first test stimulus ---\n")
print(contrast_A)
cat("\n--- Contrast B: mean of test stim 1-3 ---\n")
print(contrast_B)
cat("\n--- Contrast C: mean of all shared test stim 1-8 ---\n")
print(contrast_C)


# ==============================================================================
# 9. Plots for contrasts A, B, C
# ==============================================================================

p_per_cell_A <- per_cell_A %>%
  ggplot(aes(x = Training, y = expected_delay,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "Expected delay at stim 1 of test block",
    title = "(A) Test response delay at first stimulus"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastA_per_cell_test_stim1_expected_delay.png"),
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
    x     = NULL,
    y     = "Expected delay difference at stim 1 (spaced - massed)",
    title = "(A) Headline contrast: stim 1 of test block"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastA_diff_spaced_vs_massed_expected_delay.png"),
  p_contrast_A, width = 8, height = 5, dpi = 300, bg = "white"
)


p_per_cell_B <- per_cell_B %>%
  ggplot(aes(x = Training, y = expected_delay,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "Mean expected delay over stim 1-3 of test block",
    title = "(B) Test response delay averaged over first 3 stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastB_per_cell_test_stim1to3_expected_delay.png"),
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
    x     = NULL,
    y     = "Mean expected delay difference over stim 1-3 (spaced - massed)",
    title = "(B) Averaged contrast: mean of stim 1-3 of test block"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastB_diff_spaced_vs_massed_expected_delay.png"),
  p_contrast_B, width = 8, height = 5, dpi = 300, bg = "white"
)


p_per_cell_C <- per_cell_C %>%
  ggplot(aes(x = Training, y = expected_delay,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "Mean expected delay over all test-block stim 1-8",
    title = "(C) Test response delay averaged over all shared test stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastC_per_cell_test_meanAllStim_expected_delay.png"),
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
    x     = NULL,
    y     = "Mean expected delay difference over all test stim 1-8 (spaced - massed)",
    title = "(C) Averaged contrast: mean of all shared test stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastC_diff_spaced_vs_massed_expected_delay.png"),
  p_contrast_C, width = 8, height = 5, dpi = 300, bg = "white"
)


# ==============================================================================
# 10. Contrast D: inter-block change in expected delay, joint-covariance SE
# ==============================================================================
# delay_change_g_t  = E[delay | g, t, start_test] - E[delay | g, t, end_train]
# diff_g            = delay_change_g_spaced - delay_change_g_massed
#
# Build a single linear combination per Genotype that captures the entire
# difference-in-differences, so the SE uses the full vcov.
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
    Block    = factor(Block,    levels = levels(df_all$Block))
  )


# Get all relevant cells in one emmeans call so the covariance is fully captured.
emm_D <- emmeans(
  m_joint_delay,
  specs = ~ Genotype * Training * Block * stimulus_log,
  at = list(
    stimulus_log = sort(unique(recovery_grid_small$stimulus_log)),
    Block        = unique(as.character(recovery_grid_small$Block))
  ),
  mode = "mean.class"
)

emm_D_df <- as.data.frame(emm_D)

# Tag each row with the matching reference point (end_train / start_test), or NA.
emm_D_df$point <- NA_character_
for (i in seq_len(nrow(recovery_grid_small))) {
  r <- recovery_grid_small[i, ]
  hit <- which(
    emm_D_df$Training            == as.character(r$Training) &
      as.character(emm_D_df$Block) == as.character(r$Block) &
      abs(emm_D_df$stimulus_log - r$stimulus_log) < 1e-8
  )
  emm_D_df$point[hit] <- r$point
}

keep_D       <- !is.na(emm_D_df$point)
emm_D_sub    <- emm_D[keep_D]
emm_D_df_sub <- emm_D_df[keep_D, ]


# Per-cell summary (for plotting/reporting; inference uses contrast_D).
est_col_raw <- c("emmean", "mean.class", "estimate")[
  c("emmean", "mean.class", "estimate") %in% names(emm_D_df_sub)
][1]

per_cell_D <- emm_D_df_sub %>%
  as_tibble() %>%
  mutate(
    expected_delay = .data[[est_col_raw]] - 1,
    stimulus       = exp(stimulus_log),
    lower          = expected_delay - 1.96 * SE,
    upper          = expected_delay + 1.96 * SE
  ) %>%
  select(Genotype, Training, Block, point, stimulus, expected_delay, SE, lower, upper)

print(per_cell_D)

write.csv(
  per_cell_D,
  file.path(save_results_dir, "per_cell_endTrain_startTest_expected_delay.csv"),
  row.names = FALSE
)


# Build the difference-in-differences contrast per Genotype.
#   contrast_g = (spaced_start_test - spaced_end_train)
#              - (massed_start_test - massed_end_train)
genotypes_D <- levels(droplevels(emm_D_df_sub$Genotype))
n_cells_D   <- nrow(emm_D_df_sub)

con_list_D <- lapply(genotypes_D, function(g) {
  w <- numeric(n_cells_D)
  
  idx_ss <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "spaced" & emm_D_df_sub$point == "start_test")
  idx_se <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "spaced" & emm_D_df_sub$point == "end_train")
  idx_ms <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "massed" & emm_D_df_sub$point == "start_test")
  idx_me <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "massed" & emm_D_df_sub$point == "end_train")
  
  if (length(idx_ss) != 1 || length(idx_se) != 1 ||
      length(idx_ms) != 1 || length(idx_me) != 1) {
    stop(sprintf("Genotype %s: missing one of the four reference cells.", g))
  }
  
  w[idx_ss] <-  1
  w[idx_se] <- -1
  w[idx_ms] <- -1
  w[idx_me] <-  1
  w
})
names(con_list_D) <- genotypes_D

con_D <- contrast(emm_D_sub, method = con_list_D)

con_D_ci <- as.data.frame(confint(con_D))
con_D_p  <- as.data.frame(con_D)

est_col_D <- c("estimate", "Estimate")[c("estimate", "Estimate") %in% names(con_D_ci)][1]
lo_col_D  <- c("asymp.LCL", "lower.CL", "LCL")[c("asymp.LCL", "lower.CL", "LCL") %in% names(con_D_ci)][1]
hi_col_D  <- c("asymp.UCL", "upper.CL", "UCL")[c("asymp.UCL", "upper.CL", "UCL") %in% names(con_D_ci)][1]
p_col_D   <- c("p.value", "pvalue")[c("p.value", "pvalue") %in% names(con_D_p)][1]
z_col_D   <- c("z.ratio", "t.ratio")[c("z.ratio", "t.ratio") %in% names(con_D_p)][1]

contrast_D <- tibble(
  Genotype  = factor(con_D_ci$contrast, levels = genotypes_D),
  diff_est  = con_D_ci[[est_col_D]],
  diff_se   = con_D_ci$SE,
  diff_low  = con_D_ci[[lo_col_D]],
  diff_high = con_D_ci[[hi_col_D]],
  z         = con_D_p[[z_col_D]],
  p         = con_D_p[[p_col_D]],
  contrast  = "D_delay_change"
)

print(contrast_D)

write.csv(
  contrast_D,
  file.path(save_results_dir, "contrast_D_delay_change_spaced_vs_massed.csv"),
  row.names = FALSE
)


# Descriptive per-cell delay change (for plotting; CI here is naive).
delay_change_per_cell <- per_cell_D %>%
  select(Genotype, Training, point, expected_delay, SE) %>%
  pivot_wider(names_from = point, values_from = c(expected_delay, SE)) %>%
  mutate(
    delay_change      = expected_delay_start_test - expected_delay_end_train,
    delay_change_SE   = sqrt(SE_start_test^2 + SE_end_train^2),
    delay_change_low  = delay_change - 1.96 * delay_change_SE,
    delay_change_high = delay_change + 1.96 * delay_change_SE
  )

write.csv(
  delay_change_per_cell,
  file.path(save_results_dir, "delay_change_per_cell.csv"),
  row.names = FALSE
)


p_per_cell_D <- per_cell_D %>%
  mutate(
    point_label = factor(
      point,
      levels = c("end_train", "start_test"),
      labels = c("End of training", "Start of test")
    )
  ) %>%
  ggplot(aes(x = point_label, y = expected_delay,
             ymin = lower, ymax = upper,
             color = Genotype, group = Training)) +
  facet_grid(Genotype ~ Training) +
  geom_pointrange(linewidth = 0.7, size = 0.6) +
  geom_line(linewidth = 0.6, linetype = "dashed", alpha = 0.7) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 12) +
  labs(
    x = NULL,
    y = "Expected response delay",
    title = "End-of-training vs start-of-test expected delay per protocol"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(save_fig_dir, "contrastD_per_cell_endTrain_startTest_expected_delay.png"),
  p_per_cell_D, width = 9, height = 11, dpi = 300, bg = "white"
)


p_delay_change_per_cell <- delay_change_per_cell %>%
  select(Genotype, Training, delay_change, delay_change_low, delay_change_high) %>%
  ggplot(aes(x = Training, y = delay_change,
             ymin = delay_change_low, ymax = delay_change_high,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Delay change = E(delay start_test) - E(delay end_train)",
    title = "(D) Inter-block expected-delay change per protocol"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_delay_change_per_cell.png"),
  p_delay_change_per_cell, width = 10, height = 4.5, dpi = 300, bg = "white"
)


p_contrast_D <- contrast_D %>%
  ggplot(aes(x = Genotype, y = diff_est,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x = NULL,
    y = "Delay-change difference (spaced - massed)",
    title = "(D) Expected-delay change contrast"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_delay_change_diff_spaced_vs_massed.png"),
  p_contrast_D, width = 8, height = 5, dpi = 300, bg = "white"
)


# ==============================================================================
# Combined contrast table A + B + C + D
# ==============================================================================
all_contrasts <- bind_rows(contrast_A, contrast_B, contrast_C, contrast_D) %>%
  mutate(contrast = factor(
    contrast,
    levels = c("A_stim1", "B_mean_stim1_to_3", "C_mean_all_test_stim", "D_delay_change")
  ))

write.csv(
  all_contrasts,
  file.path(save_results_dir, "ALL_contrasts_spaced_vs_massed_expected_delay.csv"),
  row.names = FALSE
)

p_all_contrasts <- all_contrasts %>%
  ggplot(aes(x = Genotype, y = diff_est,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  facet_wrap(~ contrast, ncol = 1) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "Spaced - massed expected-delay difference",
    title = "All expected-delay memory contrasts side-by-side"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "ALL_contrasts_spaced_vs_massed_expected_delay.png"),
  p_all_contrasts, width = 9, height = 10, dpi = 300, bg = "white"
)

cat("\n--- All contrasts (joint-covariance SEs) ---\n")
print(all_contrasts)


# ==============================================================================
# 11. Supplementary contrasts
# ==============================================================================

# ------------------------------------------------------------------------------
# 11.1 Training sanity check: start vs end of Block 1 (joint-covariance SE)
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
    Block    = factor(Block,    levels = levels(df_all$Block))
  )


emm_S <- emmeans(
  m_joint_delay,
  specs = ~ Genotype * Training * Block * stimulus_log,
  at = list(
    stimulus_log = sort(unique(training_grid_small$stimulus_log)),
    Block        = "1"
  ),
  mode = "mean.class"
)

emm_S_df <- as.data.frame(emm_S)

emm_S_df$point <- NA_character_
for (i in seq_len(nrow(training_grid_small))) {
  r <- training_grid_small[i, ]
  hit <- which(
    emm_S_df$Training            == as.character(r$Training) &
      as.character(emm_S_df$Block) == as.character(r$Block) &
      abs(emm_S_df$stimulus_log - r$stimulus_log) < 1e-8
  )
  emm_S_df$point[hit] <- r$point
}

keep_S       <- !is.na(emm_S_df$point)
emm_S_sub    <- emm_S[keep_S]
emm_S_df_sub <- emm_S_df[keep_S, ]


est_col_raw_S <- c("emmean", "mean.class", "estimate")[
  c("emmean", "mean.class", "estimate") %in% names(emm_S_df_sub)
][1]

training_per_cell <- emm_S_df_sub %>%
  as_tibble() %>%
  mutate(
    expected_delay = .data[[est_col_raw_S]] - 1,
    stimulus       = exp(stimulus_log),
    lower          = expected_delay - 1.96 * SE,
    upper          = expected_delay + 1.96 * SE
  ) %>%
  select(Genotype, Training, Block, point, stimulus, expected_delay, SE, lower, upper)


# Within-Block-1 change per (Genotype, Training): end_b1 - start_b1.
# Build per-(Genotype, Training) contrasts so the SE uses the joint vcov.
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
    stop(sprintf("Genotype %s, Training %s: missing start_b1 or end_b1 cell.", g, t))
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
p_col_S   <- c("p.value", "pvalue")[c("p.value", "pvalue") %in% names(con_S_p)][1]
z_col_S   <- c("z.ratio", "t.ratio")[c("z.ratio", "t.ratio") %in% names(con_S_p)][1]

# Split the "Genotype|Training" label back out.
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
  file.path(save_results_dir, "training_sanity_check_block1_expected_delay.csv"),
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
    y = "E(delay end B1) - E(delay start B1)",
    title = "Training sanity check: within-Block 1 change in expected delay"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "sanity_check_block1_expected_delay.png"),
  p_training_sanity, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# ------------------------------------------------------------------------------
# 11.2 Latent-scale slope contrast in Block 1
# ------------------------------------------------------------------------------
# emtrends already uses the full vcov for the slopes and their pairwise
# differences, so no rewrite is needed here.
slope_emm <- emtrends(
  m_joint_delay,
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
  file.path(save_results_dir, "training_slope_block1_delay_latent_scale.csv"),
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
  file.path(save_results_dir, "training_slope_diff_spaced_vs_massed_delay_latent_scale.csv"),
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
    y = "Latent-scale slope d(cumulative logit)/d(log stim) in Block 1",
    title = "Within-Block-1 delay slope per protocol"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "training_slope_block1_per_protocol_delay_latent_scale.png"),
  p_slope, width = 10, height = 4.5, dpi = 300, bg = "white"
)

# ==============================================================================
# Done
# ==============================================================================
message("Joint ordinal GLMM delay analysis complete.")
message("Results saved under: ", save_results_dir)
message("Figures saved under: ", save_fig_dir)
