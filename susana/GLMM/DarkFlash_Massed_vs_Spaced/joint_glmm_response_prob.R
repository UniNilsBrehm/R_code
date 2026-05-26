###############################################################################
# Joint GLMM Analysis: Spaced vs Massed Training - Response Probability
# Author: Nils Brehm
# Date: 11/2025
#
# Description:
#   Joint GLMM combining the spaced and massed training datasets into one
#   model to test the spaced-vs-massed memory contrast directly.
#
#   Model:
#     move ~ Genotype * Training * Block * stimulus_log + (1 | animal)
#
#   Memory contrasts (per genotype):
#     (A) P(move) at first test stimulus,    spaced vs massed
#     (B) Mean P(move) over first 3 test stimuli, spaced vs massed
#
#   Negative spaced-massed difference => spaced fish respond less in the test
#   block => better memory retention.
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
source_base_dir <- "C:/UniFreiburg/Code"
# source_base_dir <- "C:/Users/NilsPC/Desktop/Susana"
source(file.path(source_base_dir,"/R_code/susana/GLMM/DarkFlash_Massed_vs_Spaced/utils.R"))

from_file <- FALSE  # load pre-fitted model from file

# ==============================================================================
# 1. Paths
# ==============================================================================
root_dir <- "D:/WorkingData"
# root_dir <- "C:/Users/NilsPC/Desktop/Susana"
base_dir <- file.path(root_dir, "/Susana/GLMM/DarkFlash_Massed_vs_Spaced")

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

save_results_dir <- file.path(base_dir, "results", "glmm_joint_response_prob")
save_fig_dir     <- file.path(base_dir, "figs",    "glmm_joint_response_prob")

dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. Load and prepare both datasets
# ==============================================================================
# All 5 genotypes -> NO keep filter
res_massed <- load_data(file_massed, move_th = 0, drop = c('th2, tyr', 'tyr'))
res_spaced <- load_data(file_spaced, move_th = 0, drop = c('th2, tyr', 'tyr'))

df_massed <- res_massed$df_final
df_spaced <- res_spaced$df_final


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


# CRITICAL: make animal IDs unique across experiments
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

# Remove all non-responders to stimulus 1 in Block 1
df_filtered <- df_all %>%
  group_by(animal) %>%
  filter(!any(Block == 1 & stimulus == 1 & move == 0)) %>%
  ungroup()

summary_compare <- bind_rows(
  df_all %>%
    distinct(animal, Genotype, Training) %>%
    mutate(dataset = "before"),
  
  df_filtered %>%
    distinct(animal, Genotype, Training) %>%
    mutate(dataset = "after")
) %>%
  count(dataset, Genotype, Training)

summary_compare
df_all <- df_filtered


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


# ==============================================================================
# 3. Fit the joint GLMM
# ==============================================================================
message("Fitting joint GLMM (spaced + massed)...")

# We use the same model structure as your single-experiment GLMMs but add
# Training as an additional factor in the full interaction.
# (1 | animal) replaces (1 | Video/Well) because animal IDs are now globally
# unique across the two experiments.
if (from_file) {
  # Load fitted model if available
  m_joint  <- readRDS(
    file.path(save_results_dir, "joint_glmm_spaced_vs_massed.rds")
  )
  
} else {
  # Fit the model
  m_joint <- glmmTMB(
    move ~ Genotype * Training * Block * stimulus_log + (1 | animal),
    family = binomial(link = "logit"),
    data   = df_all
  )
  
  # Save fitted model
  saveRDS(
    m_joint,
    file.path(save_results_dir, "joint_glmm_spaced_vs_massed.rds")
  )
}

capture.output(
  summary(m_joint),
  file = file.path(save_results_dir, "summary_results.txt")
)

# ==============================================================================
# 4. Validate model
# ==============================================================================
message("Validating joint model residuals...")
validate_model(m_joint, df_all)


# ==============================================================================
# 5. Plot habituation curves (joint model)
# ==============================================================================
message("Plotting joint habituation curves...")

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
  m_joint,
  newdata = new_data_joint,
  re.form = NA,
  se.fit  = TRUE
)

new_data_joint <- new_data_joint %>%
  mutate(
    fit     = plogis(pred_joint$fit),
    CI_low  = plogis(pred_joint$fit - 1.96 * pred_joint$se.fit),
    CI_high = plogis(pred_joint$fit + 1.96 * pred_joint$se.fit)
  )

raw_summary_joint <- df_all %>%
  group_by(Training, Block, Genotype, stimulus) %>%
  summarise(p_move = mean(move, na.rm = TRUE), .groups = "drop")


# Plot massed and spaced separately so the panels are readable
p_massed_curves <- ggplot(
  new_data_joint %>% filter(Training == "massed"),
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_joint %>% filter(Training == "massed"),
    aes(x = stimulus, y = p_move, color = Genotype),
    inherit.aes = FALSE, alpha = 0.25, size = 0.6
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Response probability",
    title = "Joint GLMM: Massed training"
  ) +
  theme(legend.position = "none")

p_spaced_curves <- ggplot(
  new_data_joint %>% filter(Training == "spaced"),
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_joint %>% filter(Training == "spaced"),
    aes(x = stimulus, y = p_move, color = Genotype),
    inherit.aes = FALSE, alpha = 0.25, size = 0.6
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Response probability",
    title = "Joint GLMM: Spaced training"
  ) +
  theme(legend.position = "top")

print(p_spaced_curves)
print(p_massed_curves)

ggsave(
  file.path(save_fig_dir, "joint_glmm_curves_massed.png"),
  p_massed_curves, width = 10, height = 12, dpi = 300, bg = "white"
)
ggsave(
  file.path(save_fig_dir, "joint_glmm_curves_spaced.png"),
  p_spaced_curves, width = 14, height = 12, dpi = 300, bg = "white"
)

# ------------------------------------------------------------------------------
# Test-block curves only: massed vs spaced side by side
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
    aes(x = stimulus, y = p_move, color = Genotype),
    inherit.aes = FALSE,
    alpha = 0.25,
    size = 0.7
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_pubr(base_size = 13) +
  labs(
    x = "Stimulus number within test block",
    y = "Response probability",
    title = "Joint GLMM: test-block habituation curves",
    subtitle = "Massed and spaced test blocks shown side by side"
  ) +
  theme(legend.position = "none")

print(p_test_massed_spaced_side_by_side)

ggsave(
  file.path(save_fig_dir, "joint_glmm_test_block_massed_vs_spaced_side_by_side.png"),
  p_test_massed_spaced_side_by_side,
  width = 9,
  height = 12,
  dpi = 300,
  bg = "white"
)

# ==============================================================================
# 6. Memory contrasts with proper joint-covariance SEs on response-probability scale
# ==============================================================================
#
# Important:
#   emmeans(..., type = "response") gives probabilities, but contrasts may still
#   be handled on the link scale unless we regrid.
#
#   Therefore we use:
#     emm_resp <- regrid(emm, transform = "response")
#
#   This makes probability the working linear scale, so contrast() estimates
#   differences in P(move), with SEs based on the full fixed-effect covariance.
# ==============================================================================

genotypes <- levels(df_all$Genotype)

# ------------------------------------------------------------------------------
# Helper: extract test-block emmeans cells at one or more stimulus_log values
# ------------------------------------------------------------------------------

get_test_emm_cells_response <- function(model, stim_log_values) {
  
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


emm_response_to_per_cell <- function(emm_obj) {
  emm_df <- as.data.frame(confint(emm_obj))
  
  estimate_col <- c("prob", "response", "emmean", "estimate")[
    c("prob", "response", "emmean", "estimate") %in% names(emm_df)
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
      stimulus = exp(stimulus_log),
      prob     = .data[[estimate_col]],
      SE       = SE,
      lower    = pmax(0, .data[[lower_col]]),
      upper    = pmin(1, .data[[upper_col]])
    )
}


# ------------------------------------------------------------------------------
# Core helper: proper spaced - massed contrast, averaged across selected stimuli
# ------------------------------------------------------------------------------

proper_spaced_vs_massed_prob_contrast <- function(model, stim_log_values, label) {
  
  emm_sub <- get_test_emm_cells_response(model, stim_log_values)
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
  
  est_col <- c("estimate", "Estimate")[
    c("estimate", "Estimate") %in% names(con_ci)
  ][1]
  
  lo_col <- c("asymp.LCL", "lower.CL", "LCL")[
    c("asymp.LCL", "lower.CL", "LCL") %in% names(con_ci)
  ][1]
  
  hi_col <- c("asymp.UCL", "upper.CL", "UCL")[
    c("asymp.UCL", "upper.CL", "UCL") %in% names(con_ci)
  ][1]
  
  z_col <- c("z.ratio", "t.ratio")[
    c("z.ratio", "t.ratio") %in% names(con_p)
  ][1]
  
  p_col <- c("p.value", "pvalue")[
    c("p.value", "pvalue") %in% names(con_p)
  ][1]
  
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
    per_cell = emm_response_to_per_cell(emm_sub)
  )
}


aggregate_per_cell_prob_for_plot <- function(per_cell_df) {
  per_cell_df %>%
    group_by(Genotype, Training, Block) %>%
    summarise(
      prob           = mean(prob),
      SE_descriptive = sqrt(sum(SE^2)) / n(),
      lower          = pmax(0, prob - 1.96 * SE_descriptive),
      upper          = pmin(1, prob + 1.96 * SE_descriptive),
      .groups        = "drop"
    )
}


# ------------------------------------------------------------------------------
# A, B, C memory contrasts
# ------------------------------------------------------------------------------

# A: first test stimulus
res_A <- proper_spaced_vs_massed_prob_contrast(
  m_joint,
  stim_log_values = 0,
  label = "A_stim1"
)

contrast_A <- res_A$contrast
per_cell_A <- res_A$per_cell


# B: mean over test stimuli 1:3
res_B <- proper_spaced_vs_massed_prob_contrast(
  m_joint,
  stim_log_values = log(1:3),
  label = "B_mean_stim1_to_3"
)

contrast_B     <- res_B$contrast
per_cell_B_raw <- res_B$per_cell
per_cell_B     <- aggregate_per_cell_prob_for_plot(per_cell_B_raw)


# C: mean over all shared test stimuli 1:8
test_stims_for_C <- 1:8

res_C <- proper_spaced_vs_massed_prob_contrast(
  m_joint,
  stim_log_values = log(test_stims_for_C),
  label = "C_mean_all_test_stim"
)

contrast_C     <- res_C$contrast
per_cell_C_raw <- res_C$per_cell
per_cell_C     <- aggregate_per_cell_prob_for_plot(per_cell_C_raw)


# Save outputs
write.csv(per_cell_A, file.path(save_results_dir, "per_cell_test_response_stim1.csv"), row.names = FALSE)
write.csv(per_cell_B_raw, file.path(save_results_dir, "per_cell_test_response_stim1to3.csv"), row.names = FALSE)
write.csv(per_cell_B, file.path(save_results_dir, "per_cell_test_response_stim1to3_aggregated.csv"), row.names = FALSE)
write.csv(per_cell_C_raw, file.path(save_results_dir, "per_cell_test_response_meanAllStim_raw.csv"), row.names = FALSE)
write.csv(per_cell_C, file.path(save_results_dir, "per_cell_test_response_meanAllStim.csv"), row.names = FALSE)

write.csv(contrast_A, file.path(save_results_dir, "contrast_A_test_stim1_spaced_vs_massed.csv"), row.names = FALSE)
write.csv(contrast_B, file.path(save_results_dir, "contrast_B_test_meanStim1to3_spaced_vs_massed.csv"), row.names = FALSE)
write.csv(contrast_C, file.path(save_results_dir, "contrast_C_test_meanAllStim_spaced_vs_massed.csv"), row.names = FALSE)

cat("\n--- Contrast A: first test stimulus ---\n")
print(contrast_A)

cat("\n--- Contrast B: mean stim 1-3 ---\n")
print(contrast_B)

cat("\n--- Contrast C: mean stim 1-8 ---\n")
print(contrast_C)


# ==============================================================================
# 7. Plots: memory contrasts A, B, C
# ==============================================================================

p_per_cell_A <- per_cell_A %>%
  ggplot(aes(x = Training, y = prob,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "P(move) at stim 1 of test block",
    title = "(A) Test response at first stimulus: lower = better retention"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastA_per_cell_test_stim1.png"),
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
    y = "P(move) difference at stim 1: spaced - massed\nNegative = spaced has better retention",
    title = "(A) Headline contrast: stim 1 of test block"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastA_diff_spaced_vs_massed.png"),
  p_contrast_A, width = 8, height = 5, dpi = 300, bg = "white"
)


p_per_cell_B <- per_cell_B %>%
  ggplot(aes(x = Training, y = prob,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Mean P(move) over stim 1-3 of test block",
    title = "(B) Test response averaged over first 3 stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastB_per_cell_test_stim1to3.png"),
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
    y = "Mean P(move) difference over stim 1-3: spaced - massed\nNegative = spaced has better retention",
    title = "(B) Averaged contrast: mean of stim 1-3"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastB_diff_spaced_vs_massed.png"),
  p_contrast_B, width = 8, height = 5, dpi = 300, bg = "white"
)


p_per_cell_C <- per_cell_C %>%
  ggplot(aes(x = Training, y = prob,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Mean P(move) over all shared test stimuli 1-8",
    title = "(C) Test response averaged over all shared test stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastC_per_cell_test_meanAllStim.png"),
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
    y = "Mean P(move) difference over stim 1-8: spaced - massed\nNegative = spaced has better retention",
    title = "(C) Averaged contrast: mean of all shared test stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastC_diff_spaced_vs_massed.png"),
  p_contrast_C, width = 8, height = 5, dpi = 300, bg = "white"
)


# ==============================================================================
# 8. Contrast D: inter-block recovery with joint-covariance SEs
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


emm_D <- emmeans(
  m_joint,
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


estimate_col_D_raw <- c("prob", "response", "emmean", "estimate")[
  c("prob", "response", "emmean", "estimate") %in% names(emm_D_df_sub)
][1]

per_cell_D <- emm_D_df_sub %>%
  as_tibble() %>%
  mutate(
    point    = emm_D_df_sub$point,
    prob     = .data[[estimate_col_D_raw]],
    stimulus = exp(stimulus_log),
    lower    = pmax(0, prob - 1.96 * SE),
    upper    = pmin(1, prob + 1.96 * SE)
  ) %>%
  select(Genotype, Training, Block, point, stimulus, prob, SE, lower, upper)

print(per_cell_D)

write.csv(
  per_cell_D,
  file.path(save_results_dir, "per_cell_endTrain_startTest.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# D raw recovery:
#   recovery = P(start_test) - P(end_train)
#   contrast = recovery_spaced - recovery_massed
# ------------------------------------------------------------------------------

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

est_col_D <- c("estimate", "Estimate")[
  c("estimate", "Estimate") %in% names(con_D_raw_ci)
][1]

lo_col_D <- c("asymp.LCL", "lower.CL", "LCL")[
  c("asymp.LCL", "lower.CL", "LCL") %in% names(con_D_raw_ci)
][1]

hi_col_D <- c("asymp.UCL", "upper.CL", "UCL")[
  c("asymp.UCL", "upper.CL", "UCL") %in% names(con_D_raw_ci)
][1]

z_col_D <- c("z.ratio", "t.ratio")[
  c("z.ratio", "t.ratio") %in% names(con_D_raw_p)
][1]

p_col_D <- c("p.value", "pvalue")[
  c("p.value", "pvalue") %in% names(con_D_raw_p)
][1]

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
  file.path(save_results_dir, "contrast_D_raw_recovery_spaced_vs_massed.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# D normalized recovery: nonlinear delta method using full covariance
#   recovery_norm = (P(start_test) - P(end_train)) / (1 - P(end_train))
# ------------------------------------------------------------------------------

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


beta_D <- as.numeric(emm_D_df_sub[[estimate_col_D_raw]])
V_D    <- vcov(emm_D_sub)

contrast_D_norm <- bind_rows(lapply(genotypes_D, function(g) {
  
  idx_ss <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "spaced" & emm_D_df_sub$point == "start_test")
  idx_se <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "spaced" & emm_D_df_sub$point == "end_train")
  idx_ms <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "massed" & emm_D_df_sub$point == "start_test")
  idx_me <- which(emm_D_df_sub$Genotype == g & emm_D_df_sub$Training == "massed" & emm_D_df_sub$point == "end_train")
  
  idx <- c(idx_ss, idx_se, idx_ms, idx_me)
  
  f_norm <- function(p) {
    p_ss <- p[1]
    p_se <- p[2]
    p_ms <- p[3]
    p_me <- p[4]
    
    rec_spaced <- (p_ss - p_se) / (1 - p_se)
    rec_massed <- (p_ms - p_me) / (1 - p_me)
    
    rec_spaced - rec_massed
  }
  
  p_hat <- beta_D[idx]
  V_hat <- V_D[idx, idx, drop = FALSE]
  
  est <- f_norm(p_hat)
  grad <- finite_diff_gradient(f_norm, p_hat)
  se <- sqrt(as.numeric(t(grad) %*% V_hat %*% grad))
  
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
  file.path(save_results_dir, "contrast_D_normalized_recovery_spaced_vs_massed.csv"),
  row.names = FALSE
)


# Descriptive per-cell recovery summaries for plotting only
recovery_per_cell <- per_cell_D %>%
  select(Genotype, Training, point, prob, SE) %>%
  pivot_wider(names_from = point, values_from = c(prob, SE)) %>%
  mutate(
    recovery       = prob_start_test - prob_end_train,
    recovery_SE    = sqrt(SE_start_test^2 + SE_end_train^2),
    recovery_low   = recovery - 1.96 * recovery_SE,
    recovery_high  = recovery + 1.96 * recovery_SE,
    
    recovery_norm      = recovery / (1 - prob_end_train),
    recovery_norm_SE   = recovery_SE / (1 - prob_end_train),
    recovery_norm_low  = recovery_norm - 1.96 * recovery_norm_SE,
    recovery_norm_high = recovery_norm + 1.96 * recovery_norm_SE
  )

write.csv(
  recovery_per_cell,
  file.path(save_results_dir, "recovery_per_cell.csv"),
  row.names = FALSE
)


# ==============================================================================
# 9. Plots for contrast D
# ==============================================================================

p_per_cell_D <- per_cell_D %>%
  mutate(
    point_label = factor(
      point,
      levels = c("end_train", "start_test"),
      labels = c("End of training", "Start of test")
    )
  ) %>%
  ggplot(aes(x = point_label, y = prob,
             ymin = lower, ymax = upper,
             color = Genotype, group = Training)) +
  facet_grid(Genotype ~ Training) +
  geom_pointrange(linewidth = 0.7, size = 0.6) +
  geom_line(linewidth = 0.6, linetype = "dashed", alpha = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 12) +
  labs(
    x = NULL,
    y = "P(move)",
    title = "End-of-training vs start-of-test response per protocol"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(save_fig_dir, "contrastD_per_cell_endTrain_startTest.png"),
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
    y = "Raw recovery = P(start_test) - P(end_train)",
    title = "(D) Inter-block recovery per protocol: lower = better memory"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_raw_recovery_per_cell.png"),
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
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Normalized recovery",
    title = "(D-norm) Normalized inter-block recovery per protocol"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_normalized_recovery_per_cell.png"),
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
    y = "Raw recovery difference: spaced - massed\nNegative = spaced recovers less = better retention",
    title = "(D) Recovery contrast: raw"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_raw_diff_spaced_vs_massed.png"),
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
    y = "Normalized recovery difference: spaced - massed\nNegative = spaced recovers less = better retention",
    title = "(D-norm) Recovery contrast: normalized"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_norm_diff_spaced_vs_massed.png"),
  p_contrast_D_norm, width = 8, height = 5, dpi = 300, bg = "white"
)


# ==============================================================================
# 10. Combined memory contrast table and plot
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
  file.path(save_results_dir, "ALL_contrasts_spaced_vs_massed.csv"),
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
    y = "Spaced - massed difference\nNegative = spaced has better retention",
    title = "All response-probability memory contrasts"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "ALL_contrasts_spaced_vs_massed.png"),
  p_all_contrasts, width = 9, height = 12, dpi = 300, bg = "white"
)

cat("\n--- All response-probability contrasts ---\n")
print(all_contrasts)


# ==============================================================================
# 11. Supplementary contrasts
# ==============================================================================

# ------------------------------------------------------------------------------
# 11.1 Training sanity check: start vs end of Block 1
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
  m_joint,
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

estimate_col_S_raw <- c("prob", "response", "emmean", "estimate")[
  c("prob", "response", "emmean", "estimate") %in% names(emm_S_df_sub)
][1]

training_per_cell <- emm_S_df_sub %>%
  as_tibble() %>%
  mutate(
    point    = emm_S_df_sub$point,
    prob     = .data[[estimate_col_S_raw]],
    stimulus = exp(stimulus_log),
    lower    = pmax(0, prob - 1.96 * SE),
    upper    = pmin(1, prob + 1.96 * SE)
  ) %>%
  select(Genotype, Training, Block, point, stimulus, prob, SE, lower, upper)


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

est_col_S <- c("estimate", "Estimate")[
  c("estimate", "Estimate") %in% names(con_S_ci)
][1]

lo_col_S <- c("asymp.LCL", "lower.CL", "LCL")[
  c("asymp.LCL", "lower.CL", "LCL") %in% names(con_S_ci)
][1]

hi_col_S <- c("asymp.UCL", "upper.CL", "UCL")[
  c("asymp.UCL", "upper.CL", "UCL") %in% names(con_S_ci)
][1]

z_col_S <- c("z.ratio", "t.ratio")[
  c("z.ratio", "t.ratio") %in% names(con_S_p)
][1]

p_col_S <- c("p.value", "pvalue")[
  c("p.value", "pvalue") %in% names(con_S_p)
][1]

training_learning <- con_S_ci %>%
  as_tibble() %>%
  mutate(
    contrast     = as.character(contrast),
    Genotype     = sapply(strsplit(contrast, "\\|"), `[`, 1),
    Training     = sapply(strsplit(contrast, "\\|"), `[`, 2),
    learning_est = .data[[est_col_S]],
    learning_se  = SE,
    learning_low = .data[[lo_col_S]],
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
  file.path(save_results_dir, "training_sanity_check_block1_learning.csv"),
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
    y = "P(end B1) - P(start B1)\nNegative = learning occurred",
    title = "Training sanity check: within-Block-1 learning"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "sanity_check_block1_learning.png"),
  p_training_sanity, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# ------------------------------------------------------------------------------
# 11.2 Latent-scale slope contrast in Block 1
# ------------------------------------------------------------------------------

slope_emm <- emtrends(
  m_joint,
  specs = ~ Genotype * Training | Block,
  var = "stimulus_log",
  at = list(Block = "1")
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
  file.path(save_results_dir, "training_slope_block1.csv"),
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
  file.path(save_results_dir, "training_slope_diff_spaced_vs_massed.csv"),
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
    y = "Habituation slope d(logit P)/d(log stim) in Block 1",
    title = "Within-Block-1 habituation rate per protocol"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "training_slope_block1_per_protocol.png"),
  p_slope, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# ------------------------------------------------------------------------------
# 11.3 Spaced-only end-of-block progression
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
  m_joint,
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

# Keep only matching block-specific last stimuli
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

estimate_col_end <- c("prob", "response", "emmean", "estimate")[
  c("prob", "response", "emmean", "estimate") %in% names(emm_endblock_df_sub)
][1]

spaced_endblock_result <- emm_endblock_df_sub %>%
  as_tibble() %>%
  mutate(
    prob = .data[[estimate_col_end]],
    stimulus = exp(stimulus_log),
    block_num = as.integer(as.character(Block)),
    lower = pmax(0, prob - 1.96 * SE),
    upper = pmin(1, prob + 1.96 * SE)
  ) %>%
  select(Genotype, Block, block_num, stimulus, prob, SE, lower, upper)

print(spaced_endblock_result)

write.csv(
  spaced_endblock_result,
  file.path(save_results_dir, "spaced_endblock_response.csv"),
  row.names = FALSE
)


# Block-to-block changes using joint covariance
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
  file.path(save_results_dir, "spaced_block_to_block_endblock_changes.csv"),
  row.names = FALSE
)


p_spaced_endblock <- spaced_endblock_result %>%
  ggplot(aes(x = block_num, y = prob,
             ymin = lower, ymax = upper,
             color = Genotype, group = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.7, size = 0.6) +
  geom_line(linewidth = 0.6) +
  scale_x_continuous(breaks = 1:4) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 13) +
  labs(
    x = "Spaced training block",
    y = "P(move) at end of block",
    title = "Spaced training: end-of-block response across blocks 1-4"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "spaced_endblock_progression.png"),
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
    y = "Change in end-of-block P(move)\nNegative = further habituation",
    title = "Spaced training: incremental habituation per block"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(save_fig_dir, "spaced_block_to_block_changes.png"),
  p_spaced_diffs, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# ==============================================================================
# Done
# ==============================================================================

message("Joint GLMM response-probability analysis complete.")
message("Results saved under: ", save_results_dir)
message("Figures saved under: ", save_fig_dir)