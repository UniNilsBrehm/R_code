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
m_joint_peak <- glmmTMB(
  max_peak ~ Genotype * Training * Block * stimulus_log + (1 | animal),
  family = Gamma(link = "log"),
  data   = df_all
)

saveRDS(
  m_joint_peak,
  file.path(save_results_dir, "joint_glmm_spaced_vs_massed_max_peak.rds")
)

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


# ==============================================================================
# 7. Helper functions for max_peak memory contrasts
# ==============================================================================

genotypes <- levels(df_all$Genotype)

# Helper: robustly extract the response-scale estimate and CI names from emmeans.
# For Gamma(log), emmeans(type = "response") usually returns column "response",
# but this helper also works if your installed version returns a different name.
get_emm_response_col <- function(emm_df) {
  candidate_cols <- c("response", "rate", "prob", "emmean")
  found <- candidate_cols[candidate_cols %in% names(emm_df)]
  if (length(found) == 0) {
    stop("Could not find response estimate column in emmeans output. Columns are: ",
         paste(names(emm_df), collapse = ", "))
  }
  found[1]
}

get_lcl_col <- function(emm_df) {
  candidate_cols <- c("asymp.LCL", "lower.CL", "LCL")
  found <- candidate_cols[candidate_cols %in% names(emm_df)]
  if (length(found) == 0) {
    stop("Could not find lower CI column in emmeans output. Columns are: ",
         paste(names(emm_df), collapse = ", "))
  }
  found[1]
}

get_ucl_col <- function(emm_df) {
  candidate_cols <- c("asymp.UCL", "upper.CL", "UCL")
  found <- candidate_cols[candidate_cols %in% names(emm_df)]
  if (length(found) == 0) {
    stop("Could not find upper CI column in emmeans output. Columns are: ",
         paste(names(emm_df), collapse = ", "))
  }
  found[1]
}


# Helper to compute spaced-vs-massed contrast at a fixed stimulus_log value.
# Returns:
#   per_cell: predicted max_peak for each Genotype x Training test-block cell
#   contrasts: spaced - massed difference per Genotype
compute_peak_contrast_at <- function(model, stim_log_value) {

  emm <- emmeans(
    model,
    specs = ~ Genotype * Training * Block,
    at = list(
      stimulus_log = stim_log_value,
      Block        = c(massed_test_block, spaced_test_block)
    ),
    type = "response"
  )

  emm_df <- as_tibble(emm) %>%
    filter(
      (Training == "massed" & Block == massed_test_block) |
        (Training == "spaced" & Block == spaced_test_block)
    )

  resp_col <- get_emm_response_col(emm_df)
  lcl_col  <- get_lcl_col(emm_df)
  ucl_col  <- get_ucl_col(emm_df)

  per_cell <- emm_df %>%
    transmute(
      Genotype,
      Training,
      Block,
      max_peak = .data[[resp_col]],
      SE       = SE,
      lower    = .data[[lcl_col]],
      upper    = .data[[ucl_col]],
      stim_log_value = stim_log_value,
      stimulus       = exp(stim_log_value)
    )

  # Build spaced - massed contrast per Genotype.
  # This uses response-scale differences in mm.
  contrasts_df <- per_cell %>%
    select(Genotype, Training, max_peak, SE) %>%
    pivot_wider(
      names_from  = Training,
      values_from = c(max_peak, SE)
    ) %>%
    mutate(
      diff_est  = max_peak_spaced - max_peak_massed,
      diff_se   = sqrt(SE_spaced^2 + SE_massed^2),
      diff_low  = diff_est - 1.96 * diff_se,
      diff_high = diff_est + 1.96 * diff_se,
      z         = diff_est / diff_se,
      p         = 2 * pnorm(-abs(z)),
      stim_log_value = stim_log_value,
      stimulus       = exp(stim_log_value)
    )

  list(
    per_cell  = per_cell,
    contrasts = contrasts_df
  )
}


aggregate_per_cell_predictions <- function(per_cell_df) {
  per_cell_df %>%
    group_by(Genotype, Training, Block) %>%
    summarise(
      max_peak_mean = mean(max_peak),
      SE_mean       = sqrt(sum(SE^2)) / n(),
      lower         = pmax(0, max_peak_mean - 1.96 * SE_mean),
      upper         = max_peak_mean + 1.96 * SE_mean,
      .groups       = "drop"
    ) %>%
    rename(max_peak = max_peak_mean)
}


# ==============================================================================
# 8. Memory contrasts A, B, C
# ==============================================================================
# (A) at stimulus = 1
# (B) mean over stimuli 1, 2, 3
# (C) mean over shared test stimuli 1..8
# ==============================================================================

# A: test stimulus 1
res_A <- compute_peak_contrast_at(m_joint_peak, stim_log_value = 0)

contrast_A <- res_A$contrasts %>%
  select(Genotype, diff_est, diff_se, diff_low, diff_high, z, p,
         stim_log_value, stimulus) %>%
  mutate(contrast = "A_stim1")


# B: mean over test stimuli 1..3
res_B1 <- compute_peak_contrast_at(m_joint_peak, stim_log_value = log(1))
res_B2 <- compute_peak_contrast_at(m_joint_peak, stim_log_value = log(2))
res_B3 <- compute_peak_contrast_at(m_joint_peak, stim_log_value = log(3))

contrast_B <- bind_rows(res_B1$contrasts, res_B2$contrasts, res_B3$contrasts) %>%
  group_by(Genotype) %>%
  summarise(
    diff_est  = mean(diff_est),
    diff_se   = sqrt(sum(diff_se^2)) / n(),
    diff_low  = diff_est - 1.96 * diff_se,
    diff_high = diff_est + 1.96 * diff_se,
    z         = diff_est / diff_se,
    p         = 2 * pnorm(-abs(z)),
    .groups   = "drop"
  ) %>%
  mutate(
    stim_log_value = NA_real_,
    stimulus       = NA_real_,
    contrast       = "B_mean_stim1_to_3"
  )

per_cell_B <- aggregate_per_cell_predictions(
  bind_rows(res_B1$per_cell, res_B2$per_cell, res_B3$per_cell)
)


# C: mean over all shared test stimuli
# Spaced test block has 8 stimuli and massed has 9 in your response-prob script.
# Use 1..8 for both for a clean apples-to-apples comparison.
test_stims_for_C <- 1:8

res_C_list <- lapply(test_stims_for_C, function(s) {
  compute_peak_contrast_at(m_joint_peak, stim_log_value = log(s))
})

contrast_C <- bind_rows(lapply(res_C_list, `[[`, "contrasts")) %>%
  group_by(Genotype) %>%
  summarise(
    diff_est  = mean(diff_est),
    diff_se   = sqrt(sum(diff_se^2)) / n(),
    diff_low  = diff_est - 1.96 * diff_se,
    diff_high = diff_est + 1.96 * diff_se,
    z         = diff_est / diff_se,
    p         = 2 * pnorm(-abs(z)),
    .groups   = "drop"
  ) %>%
  mutate(
    stim_log_value = NA_real_,
    stimulus       = NA_real_,
    contrast       = "C_mean_all_test_stim"
  )

per_cell_C <- aggregate_per_cell_predictions(
  bind_rows(lapply(res_C_list, `[[`, "per_cell"))
)


# Save per-cell predictions and contrasts --------------------------------------
write.csv(
  res_A$per_cell,
  file.path(save_results_dir, "per_cell_test_max_peak_stim1.csv"),
  row.names = FALSE
)

write.csv(
  bind_rows(res_B1$per_cell, res_B2$per_cell, res_B3$per_cell),
  file.path(save_results_dir, "per_cell_test_max_peak_stim1to3.csv"),
  row.names = FALSE
)

write.csv(
  per_cell_B,
  file.path(save_results_dir, "per_cell_test_max_peak_stim1to3_aggregated.csv"),
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

print(contrast_A)
print(contrast_B)
print(contrast_C)


# ==============================================================================
# 9. Plots for contrasts A, B, C
# ==============================================================================

# A: per-cell max_peak at stim 1 of test block
p_per_cell_A <- res_A$per_cell %>%
  ggplot(aes(x = Training, y = max_peak,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "Predicted max_peak at stim 1 of test block",
    title = "(A) Test peak distance at first stimulus: lower = better retention"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastA_per_cell_test_stim1_max_peak.png"),
  p_per_cell_A, width = 12, height = 4.5, dpi = 300, bg = "white"
)


# A: spaced - massed difference
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
    y     = "max_peak difference at stim 1 (spaced - massed)\nNegative = spaced has better retention",
    title = "(A) Headline contrast: stim 1 of test block"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastA_diff_spaced_vs_massed_max_peak.png"),
  p_contrast_A, width = 8, height = 5, dpi = 300, bg = "white"
)


# B: per-cell mean over stim 1..3
p_per_cell_B <- per_cell_B %>%
  ggplot(aes(x = Training, y = max_peak,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "Mean predicted max_peak over stim 1-3 of test block",
    title = "(B) Test peak distance averaged over first 3 stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastB_per_cell_test_stim1to3_max_peak.png"),
  p_per_cell_B, width = 12, height = 4.5, dpi = 300, bg = "white"
)


# B: spaced - massed difference
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
    y     = "Mean max_peak difference over stim 1-3 (spaced - massed)\nNegative = spaced has better retention",
    title = "(B) Averaged contrast: mean of stim 1-3 of test block"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastB_diff_spaced_vs_massed_max_peak.png"),
  p_contrast_B, width = 8, height = 5, dpi = 300, bg = "white"
)


# C: per-cell mean over shared all-test stimuli
p_per_cell_C <- per_cell_C %>%
  ggplot(aes(x = Training, y = max_peak,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "Mean predicted max_peak over all test-block stim 1-8",
    title = "(C) Test peak distance averaged over all shared test stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastC_per_cell_test_meanAllStim_max_peak.png"),
  p_per_cell_C, width = 12, height = 4.5, dpi = 300, bg = "white"
)


# C: spaced - massed difference
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
    y     = "Mean max_peak difference over all test stim 1-8 (spaced - massed)\nNegative = spaced has better retention",
    title = "(C) Averaged contrast: mean of all shared test stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastC_diff_spaced_vs_massed_max_peak.png"),
  p_contrast_C, width = 8, height = 5, dpi = 300, bg = "white"
)


# Combined A/B/C plot
all_contrasts <- bind_rows(contrast_A, contrast_B, contrast_C) %>%
  mutate(contrast = factor(
    contrast,
    levels = c("A_stim1", "B_mean_stim1_to_3", "C_mean_all_test_stim")
  ))

write.csv(
  all_contrasts,
  file.path(save_results_dir, "ALL_contrasts_spaced_vs_massed_max_peak.csv"),
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
    y     = "Spaced - massed max_peak difference\nNegative = spaced has better retention",
    title = "All max_peak memory contrasts side-by-side"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "ALL_contrasts_spaced_vs_massed_max_peak.png"),
  p_all_contrasts, width = 9, height = 10, dpi = 300, bg = "white"
)


# ==============================================================================
# 10. Contrast D: Inter-block recovery in max_peak
# ==============================================================================
# For each (Genotype, Training), compute:
#   recovery      = max_peak(start_test) - max_peak(end_train)
#   recovery_norm = recovery / max_peak(end_train)
#
# Then contrast spaced - massed within each Genotype.
#
# Negative spaced-massed recovery difference means spaced recovers less from
# the end-of-training state, consistent with better retention.
# ==============================================================================

# Determine actual last training stimulus per Training
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


recovery_grid <- bind_rows(
  tidyr::expand_grid(
    Genotype = levels(df_all$Genotype),
    point    = c("end_train", "start_test")
  ) %>%
    mutate(
      Training = "massed",
      Block    = ifelse(point == "end_train", last_block_massed, massed_test_block),
      stimulus = ifelse(point == "end_train", last_stim_massed, 1)
    ),
  tidyr::expand_grid(
    Genotype = levels(df_all$Genotype),
    point    = c("end_train", "start_test")
  ) %>%
    mutate(
      Training = "spaced",
      Block    = ifelse(point == "end_train", last_block_spaced, spaced_test_block),
      stimulus = ifelse(point == "end_train", last_stim_spaced, 1)
    )
) %>%
  mutate(
    stimulus_log = log(stimulus),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  )

print(recovery_grid)

preds <- predict(
  m_joint_peak,
  newdata = recovery_grid,
  re.form = NA,
  se.fit  = TRUE,
  type    = "link"
)

recovery_grid <- recovery_grid %>%
  mutate(
    fit_link = preds$fit,
    SE_link  = preds$se.fit,
    max_peak = exp(fit_link),
    # Delta method for log link: SE_response = response * SE_link
    SE       = max_peak * SE_link,
    lower    = pmax(0, max_peak - 1.96 * SE),
    upper    = max_peak + 1.96 * SE
  )

per_cell_D <- recovery_grid %>%
  select(Genotype, Training, Block, point, stimulus, max_peak, SE, lower, upper)

print(per_cell_D)

write.csv(
  per_cell_D,
  file.path(save_results_dir, "per_cell_endTrain_startTest_max_peak.csv"),
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


contrast_D_raw <- recovery_per_cell %>%
  select(Genotype, Training, recovery, recovery_SE) %>%
  pivot_wider(names_from = Training, values_from = c(recovery, recovery_SE)) %>%
  mutate(
    diff_est  = recovery_spaced - recovery_massed,
    diff_se   = sqrt(recovery_SE_spaced^2 + recovery_SE_massed^2),
    diff_low  = diff_est - 1.96 * diff_se,
    diff_high = diff_est + 1.96 * diff_se,
    z         = diff_est / diff_se,
    p         = 2 * pnorm(-abs(z)),
    contrast  = "D_raw_recovery"
  ) %>%
  select(Genotype, diff_est, diff_se, diff_low, diff_high, z, p, contrast)


contrast_D_norm <- recovery_per_cell %>%
  select(Genotype, Training, recovery_norm, recovery_norm_SE) %>%
  pivot_wider(names_from = Training,
              values_from = c(recovery_norm, recovery_norm_SE)) %>%
  mutate(
    diff_est  = recovery_norm_spaced - recovery_norm_massed,
    diff_se   = sqrt(recovery_norm_SE_spaced^2 + recovery_norm_SE_massed^2),
    diff_low  = diff_est - 1.96 * diff_se,
    diff_high = diff_est + 1.96 * diff_se,
    z         = diff_est / diff_se,
    p         = 2 * pnorm(-abs(z)),
    contrast  = "D_normalized_recovery"
  ) %>%
  select(Genotype, diff_est, diff_se, diff_low, diff_high, z, p, contrast)

print(contrast_D_raw)
print(contrast_D_norm)

write.csv(
  contrast_D_raw,
  file.path(save_results_dir, "contrast_D_raw_recovery_spaced_vs_massed_max_peak.csv"),
  row.names = FALSE
)

write.csv(
  contrast_D_norm,
  file.path(save_results_dir, "contrast_D_normalized_recovery_spaced_vs_massed_max_peak.csv"),
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
    y = "Raw max_peak recovery difference (spaced - massed)\nNegative = spaced recovers less = better retention",
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
    y = "Normalized max_peak recovery difference (spaced - massed)\nNegative = spaced recovers less = better retention",
    title = "(D-norm) Recovery contrast: normalized max_peak"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_norm_diff_spaced_vs_massed_max_peak.png"),
  p_contrast_D_norm, width = 8, height = 5, dpi = 300, bg = "white"
)


# ==============================================================================
# 12. Supplementary contrasts
# ==============================================================================

# ------------------------------------------------------------------------------
# 12.1 Training sanity check: start vs end of Block 1
# ------------------------------------------------------------------------------
last_b1_stim <- df_all %>%
  filter(Block == "1") %>%
  group_by(Training) %>%
  summarise(last_stim = max(stimulus, na.rm = TRUE), .groups = "drop")

last_b1_massed <- last_b1_stim$last_stim[last_b1_stim$Training == "massed"]
last_b1_spaced <- last_b1_stim$last_stim[last_b1_stim$Training == "spaced"]

training_grid <- tidyr::expand_grid(
  Genotype = levels(df_all$Genotype),
  Training = c("massed", "spaced"),
  point    = c("start_b1", "end_b1")
) %>%
  mutate(
    Block    = "1",
    stimulus = case_when(
      point == "start_b1"                      ~ 1,
      point == "end_b1" & Training == "massed" ~ last_b1_massed,
      point == "end_b1" & Training == "spaced" ~ last_b1_spaced
    ),
    stimulus_log = log(stimulus),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  )

train_preds <- predict(
  m_joint_peak,
  newdata = training_grid,
  re.form = NA,
  se.fit  = TRUE,
  type    = "link"
)

training_grid <- training_grid %>%
  mutate(
    fit_link = train_preds$fit,
    SE_link  = train_preds$se.fit,
    max_peak = exp(fit_link),
    SE       = max_peak * SE_link,
    lower    = pmax(0, max_peak - 1.96 * SE),
    upper    = max_peak + 1.96 * SE
  )

training_learning <- training_grid %>%
  select(Genotype, Training, point, max_peak, SE) %>%
  pivot_wider(names_from = point, values_from = c(max_peak, SE)) %>%
  mutate(
    learning_est  = max_peak_end_b1 - max_peak_start_b1,
    learning_se   = sqrt(SE_end_b1^2 + SE_start_b1^2),
    learning_low  = learning_est - 1.96 * learning_se,
    learning_high = learning_est + 1.96 * learning_se,
    z = learning_est / learning_se,
    p = 2 * pnorm(-abs(z))
  )

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
    title = "Training sanity check: within-Block 1 change in max_peak"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "sanity_check_block1_learning_max_peak.png"),
  p_training_sanity, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# ------------------------------------------------------------------------------
# 12.2 Within-training habituation slope contrast
# ------------------------------------------------------------------------------
# Slope is on the log-response scale because the model uses a log link.
# Negative slope = max_peak decreases as stimulus number increases.
slope_emm <- emtrends(
  m_joint_peak,
  specs = ~ Genotype * Training | Block,
  var   = "stimulus_log",
  at    = list(Block = "1")
)

slope_df <- as_tibble(slope_emm) %>%
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
# 12.3 Block-by-block learning curve, spaced only
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
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  )

endblock_preds <- predict(
  m_joint_peak,
  newdata = spaced_endblock_grid,
  re.form = NA,
  se.fit  = TRUE,
  type    = "link"
)

spaced_endblock_grid <- spaced_endblock_grid %>%
  mutate(
    fit_link = endblock_preds$fit,
    SE_link  = endblock_preds$se.fit,
    max_peak = exp(fit_link),
    SE       = max_peak * SE_link,
    lower    = pmax(0, max_peak - 1.96 * SE),
    upper    = max_peak + 1.96 * SE,
    block_num = as.integer(as.character(Block))
  )

print(
  spaced_endblock_grid %>%
    select(Genotype, Block, block_num, stimulus, max_peak, lower, upper)
)

write.csv(
  spaced_endblock_grid %>%
    select(Genotype, Block, block_num, stimulus, max_peak, SE, lower, upper),
  file.path(save_results_dir, "spaced_endblock_max_peak.csv"),
  row.names = FALSE
)


spaced_block_diffs <- spaced_endblock_grid %>%
  arrange(Genotype, block_num) %>%
  group_by(Genotype) %>%
  mutate(
    max_peak_prev = lag(max_peak),
    SE_prev       = lag(SE),
    delta         = max_peak - max_peak_prev,
    delta_se      = sqrt(SE^2 + SE_prev^2),
    delta_low     = delta - 1.96 * delta_se,
    delta_high    = delta + 1.96 * delta_se,
    transition    = paste0("B", block_num - 1, "->B", block_num)
  ) %>%
  ungroup() %>%
  filter(!is.na(max_peak_prev)) %>%
  select(Genotype, transition, delta, delta_se, delta_low, delta_high)

print(spaced_block_diffs)

write.csv(
  spaced_block_diffs,
  file.path(save_results_dir, "spaced_block_to_block_endblock_changes_max_peak.csv"),
  row.names = FALSE
)


p_spaced_endblock <- spaced_endblock_grid %>%
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
