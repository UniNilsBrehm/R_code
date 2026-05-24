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

save_results_dir <- file.path(base_dir, "results", "glmm_joint")
save_fig_dir     <- file.path(base_dir, "figs",    "glmm_joint")

dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. Load and prepare both datasets
# ==============================================================================
# All 5 genotypes -> NO keep filter
res_massed <- load_data(file_massed, move_th = 1, drop = c('th2, tyr', 'tyr'))
res_spaced <- load_data(file_spaced, move_th = 1, drop = c('th2, tyr', 'tyr'))

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


# ==============================================================================
# 6. Memory contrasts
# ==============================================================================
# We extract P(move) from the joint GLMM at the start of the test block,
# for each (Genotype, Training) cell. Two definitions:
#
#   (A) at stimulus = 1                       -> stimulus_log = 0
#   (B) mean over stimulus in {1, 2, 3}       -> stimulus_log = log(1), log(2), log(3)
#
# We use emmeans with by-cell predictions at fixed stimulus_log values, then
# compute spaced - massed differences per genotype.
# ==============================================================================

genotypes <- levels(df_all$Genotype)

# Define the test block per Training (needed to pin Block in emmeans)
test_block_per_training <- tibble(
  Training = c("massed", "spaced"),
  Block    = c(massed_test_block, spaced_test_block)
)

# Helper to compute spaced-vs-massed contrast at a fixed stimulus_log value.
# Returns a tibble: Genotype, stimulus_log, est_diff, lower, upper, p
compute_memory_contrast_at <- function(model, stim_log_value) {
  
  # Get predicted P(move) on the response scale at each
  # (Genotype, Training, Block) cell, with Block restricted to the test blocks.
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
  
  per_cell <- emm_df %>%
    select(Genotype, Training, Block,
           prob = prob, SE,
           lower = asymp.LCL, upper = asymp.UCL) %>%
    mutate(stim_log_value = stim_log_value,
           stimulus       = exp(stim_log_value))
  
  # Build the spaced - massed contrast per Genotype.
  # spaced uses spaced_test_block, massed uses massed_test_block -> disjoint
  # animals, so we can combine SEs assuming independence.
  contrasts_df <- per_cell %>%
    select(Genotype, Training, prob, SE) %>%
    pivot_wider(
      names_from  = Training,
      values_from = c(prob, SE)
    ) %>%
    mutate(
      diff_est  = prob_spaced - prob_massed,
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

# (A) at stimulus = 1   (stimulus_log = 0)
res_A <- compute_memory_contrast_at(m_joint, stim_log_value = 0)

# (B) mean over stimulus in {1, 2, 3}: average of three predictions
res_B1 <- compute_memory_contrast_at(m_joint, stim_log_value = log(1))
res_B2 <- compute_memory_contrast_at(m_joint, stim_log_value = log(2))
res_B3 <- compute_memory_contrast_at(m_joint, stim_log_value = log(3))

# Aggregate (B): average the spaced-massed contrasts across the three stimuli
# per Genotype, propagating SE conservatively as the mean of SEs / sqrt(3).
contrast_B <- bind_rows(res_B1$contrasts, res_B2$contrasts, res_B3$contrasts) %>%
  group_by(Genotype) %>%
  summarise(
    diff_est  = mean(diff_est),
    diff_se   = sqrt(sum(diff_se^2)) / n(),     # conservative SE for the mean
    diff_low  = diff_est - 1.96 * diff_se,
    diff_high = diff_est + 1.96 * diff_se,
    z         = diff_est / diff_se,
    p         = 2 * pnorm(-abs(z)),
    .groups   = "drop"
  ) %>%
  mutate(
    stim_log_value = NA,
    stimulus       = NA,
    contrast       = "B_mean_stim1_to_3"
  )

contrast_A <- res_A$contrasts %>%
  select(Genotype, diff_est, diff_se, diff_low, diff_high, z, p,
         stim_log_value, stimulus) %>%
  mutate(contrast = "A_stim1")


# Save per-cell predictions and contrasts ---------------------------
write.csv(
  res_A$per_cell,
  file.path(save_results_dir, "per_cell_test_response_stim1.csv"),
  row.names = FALSE
)

write.csv(
  bind_rows(res_B1$per_cell, res_B2$per_cell, res_B3$per_cell),
  file.path(save_results_dir, "per_cell_test_response_stim1to3.csv"),
  row.names = FALSE
)

write.csv(
  contrast_A,
  file.path(save_results_dir, "contrast_A_test_stim1_spaced_vs_massed.csv"),
  row.names = FALSE
)

write.csv(
  contrast_B,
  file.path(save_results_dir, "contrast_B_test_meanStim1to3_spaced_vs_massed.csv"),
  row.names = FALSE
)

print(contrast_A)
print(contrast_B)


# ==============================================================================
# 7. Plots: per-cell test response, and spaced - massed contrasts
# ==============================================================================

# (A) per-cell P(move) at stim 1 of test block
p_per_cell_A <- res_A$per_cell %>%
  ggplot(aes(x = Training, y = prob,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "P(move) at stim 1 of test block",
    title = "(A) Test response at first stimulus: lower = better retention"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastA_per_cell_test_stim1.png"),
  p_per_cell_A, width = 12, height = 4.5, dpi = 300, bg = "white"
)


# Contrast (A) plot
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
    y     = "P(move) difference at stim 1 (spaced - massed)\nNegative = spaced has better retention",
    title = "(A) Headline contrast: stim 1 of test block"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastA_diff_spaced_vs_massed.png"),
  p_contrast_A, width = 8, height = 5, dpi = 300, bg = "white"
)


# Contrast (B) plot
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
    y     = "Mean P(move) difference over stim 1-3 (spaced - massed)\nNegative = spaced has better retention",
    title = "(B) Averaged contrast: mean of stim 1-3 of test block"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastB_diff_spaced_vs_massed.png"),
  p_contrast_B, width = 8, height = 5, dpi = 300, bg = "white"
)

# ==============================================================================
# Contrast C: mean P(move) over ALL test-block stimuli
# ==============================================================================
# Spaced test block has 8 stimuli (1-8), massed has 9 (1-9). For a clean
# apples-to-apples comparison we use stim 1..8 for both (drop massed stim 9).
# Alternative: average each protocol over its own full range. The difference
# is tiny in practice; we use stim 1..8 for symmetry.
# ==============================================================================

test_stims_for_C <- 1:8

res_C_list <- lapply(test_stims_for_C, function(s) {
  compute_memory_contrast_at(m_joint, stim_log_value = log(s))
})

# Aggregate Contrast C: mean of spaced-massed differences across stim 1..8
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
    stim_log_value = NA,
    stimulus       = NA,
    contrast       = "C_mean_all_test_stim"
  )

print(contrast_C)

write.csv(
  contrast_C,
  file.path(save_results_dir, "contrast_C_test_meanAllStim_spaced_vs_massed.csv"),
  row.names = FALSE
)


# Per-cell predictions for Contrast C (averaged across stim 1..8 within each cell)
per_cell_C <- bind_rows(lapply(res_C_list, `[[`, "per_cell")) %>%
  group_by(Genotype, Training, Block) %>%
  summarise(
    prob_mean = mean(prob),
    # Conservative SE for the average within a cell
    SE_mean   = sqrt(sum(SE^2)) / n(),
    lower     = pmax(0, prob_mean - 1.96 * SE_mean),
    upper     = pmin(1, prob_mean + 1.96 * SE_mean),
    .groups   = "drop"
  ) %>%
  rename(prob = prob_mean)

write.csv(
  per_cell_C,
  file.path(save_results_dir, "per_cell_test_response_meanAllStim.csv"),
  row.names = FALSE
)


# ==============================================================================
# Per-cell plots for Contrast B and Contrast C
# ==============================================================================

# Per-cell for Contrast B (mean over stim 1-3)
# This requires aggregating res_B1, res_B2, res_B3 the same way
per_cell_B <- bind_rows(res_B1$per_cell, res_B2$per_cell, res_B3$per_cell) %>%
  group_by(Genotype, Training, Block) %>%
  summarise(
    prob_mean = mean(prob),
    SE_mean   = sqrt(sum(SE^2)) / n(),
    lower     = pmax(0, prob_mean - 1.96 * SE_mean),
    upper     = pmin(1, prob_mean + 1.96 * SE_mean),
    .groups   = "drop"
  ) %>%
  rename(prob = prob_mean)

write.csv(
  per_cell_B,
  file.path(save_results_dir, "per_cell_test_response_stim1to3_aggregated.csv"),
  row.names = FALSE
)


# Plot: Contrast B per-cell
p_per_cell_B <- per_cell_B %>%
  ggplot(aes(x = Training, y = prob,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "Mean P(move) over stim 1-3 of test block",
    title = "(B) Test response averaged over first 3 stimuli: lower = better retention"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastB_per_cell_test_stim1to3.png"),
  p_per_cell_B, width = 12, height = 4.5, dpi = 300, bg = "white"
)


# Plot: Contrast C per-cell
p_per_cell_C <- per_cell_C %>%
  ggplot(aes(x = Training, y = prob,
             ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "Mean P(move) over all test-block stim (1-8)",
    title = "(C) Test response averaged over all test stimuli: lower = better retention"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastC_per_cell_test_meanAllStim.png"),
  p_per_cell_C, width = 12, height = 4.5, dpi = 300, bg = "white"
)


# Plot: Contrast C diff spaced - massed
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
    y     = "Mean P(move) difference over all test stim (spaced - massed)\nNegative = spaced has better retention",
    title = "(C) Averaged contrast: mean of all test stimuli"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastC_diff_spaced_vs_massed.png"),
  p_contrast_C, width = 8, height = 5, dpi = 300, bg = "white"
)


# ==============================================================================
# Updated combined comparison (A, B, C side-by-side)
# ==============================================================================
all_contrasts <- bind_rows(contrast_A, contrast_B, contrast_C) %>%
  mutate(contrast = factor(
    contrast,
    levels = c("A_stim1", "B_mean_stim1_to_3", "C_mean_all_test_stim")
  ))

write.csv(
  all_contrasts,
  file.path(save_results_dir, "ALL_contrasts_spaced_vs_massed.csv"),
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
    y     = "Spaced - massed P(move) difference\nNegative = spaced has better retention",
    title = "All memory contrasts side-by-side"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "ALL_contrasts_spaced_vs_massed.png"),
  p_all_contrasts, width = 9, height = 10, dpi = 300, bg = "white"
)

# ==============================================================================
# Contrast D: Inter-block recovery
# ==============================================================================
# For each (Genotype, Training), compute:
#   recovery       = P(move|test_stim1)         - P(move|last_training_stim)
#   recovery_norm  = recovery / (1 - P(move|last_training_stim))
#
# Then contrast spaced - massed within each Genotype.
#
# Stimulus positions used:
#   Massed:  last training stim = stim 478 of Block 1; first test stim = stim 1 of Block 2
#   Spaced:  last training stim = stim 119 of Block 4; first test stim = stim 1 of Block 5
# (Use the maximum stimulus actually observed in each block.)
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

last_stim_massed <- last_train_stim$last_stim[last_train_stim$Training == "massed"]
last_stim_spaced <- last_train_stim$last_stim[last_train_stim$Training == "spaced"]
last_block_massed <- last_train_stim$last_block[last_train_stim$Training == "massed"]
last_block_spaced <- last_train_stim$last_block[last_train_stim$Training == "spaced"]


# Build a single grid for end-of-training and start-of-test predictions
recovery_grid <- bind_rows(
  # massed: end-of-training (Block 1, last stim) and start-of-test (Block 2, stim 1)
  tidyr::expand_grid(
    Genotype = levels(df_all$Genotype),
    point    = c("end_train", "start_test")
  ) %>%
    mutate(
      Training  = "massed",
      Block     = ifelse(point == "end_train", last_block_massed, massed_test_block),
      stimulus  = ifelse(point == "end_train", last_stim_massed, 1)
    ),
  # spaced: end-of-training (Block 4, last stim) and start-of-test (Block 5, stim 1)
  tidyr::expand_grid(
    Genotype = levels(df_all$Genotype),
    point    = c("end_train", "start_test")
  ) %>%
    mutate(
      Training  = "spaced",
      Block     = ifelse(point == "end_train", last_block_spaced, spaced_test_block),
      stimulus  = ifelse(point == "end_train", last_stim_spaced, 1)
    )
) %>%
  mutate(
    stimulus_log = log(stimulus),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  )

print(recovery_grid)

# Get population-level predictions (no animal RE)
preds <- predict(
  m_joint,
  newdata = recovery_grid,
  re.form = NA,
  se.fit  = TRUE,
  type    = "link"        # logit scale; we'll transform after
)

recovery_grid <- recovery_grid %>%
  mutate(
    fit_link  = preds$fit,
    SE_link   = preds$se.fit,
    fit_prob  = plogis(fit_link),
    # delta method: SE on prob scale = SE_link * fit_prob * (1 - fit_prob)
    SE_prob   = SE_link * fit_prob * (1 - fit_prob)
  )


# Per-cell summary
per_cell_D <- recovery_grid %>%
  select(Genotype, Training, Block, point, stimulus, fit_prob, SE_prob) %>%
  rename(prob = fit_prob, SE = SE_prob) %>%
  mutate(
    lower = pmax(0, prob - 1.96 * SE),
    upper = pmin(1, prob + 1.96 * SE)
  )

print(per_cell_D)

write.csv(
  per_cell_D,
  file.path(save_results_dir, "per_cell_endTrain_startTest.csv"),
  row.names = FALSE
)


# Compute raw and normalized recovery per (Genotype, Training)
recovery_per_cell <- per_cell_D %>%
  select(Genotype, Training, point, prob, SE) %>%
  pivot_wider(names_from = point, values_from = c(prob, SE)) %>%
  mutate(
    recovery       = prob_start_test - prob_end_train,
    # Independence assumption between the two predictions; conservative but ok
    recovery_SE    = sqrt(SE_start_test^2 + SE_end_train^2),
    recovery_low   = recovery - 1.96 * recovery_SE,
    recovery_high  = recovery + 1.96 * recovery_SE,
    
    recovery_norm     = recovery / (1 - prob_end_train),
    # Approximate SE for normalized recovery via delta method
    recovery_norm_SE  = recovery_SE / (1 - prob_end_train),
    recovery_norm_low  = recovery_norm - 1.96 * recovery_norm_SE,
    recovery_norm_high = recovery_norm + 1.96 * recovery_norm_SE
  )

write.csv(
  recovery_per_cell,
  file.path(save_results_dir, "recovery_per_cell.csv"),
  row.names = FALSE
)


# Spaced - massed contrast on raw recovery
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


# Spaced - massed contrast on normalized recovery
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
  file.path(save_results_dir, "contrast_D_raw_recovery_spaced_vs_massed.csv"),
  row.names = FALSE
)
write.csv(
  contrast_D_norm,
  file.path(save_results_dir, "contrast_D_normalized_recovery_spaced_vs_massed.csv"),
  row.names = FALSE
)


# ==============================================================================
# Plots for Contrast D
# ==============================================================================

# Per-cell: end-of-training vs start-of-test, side by side
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
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  file.path(save_fig_dir, "contrastD_per_cell_endTrain_startTest.png"),
  p_per_cell_D, width = 9, height = 11, dpi = 300, bg = "white"
)


# Per-cell raw recovery summary
recovery_per_cell_long_raw <- recovery_per_cell %>%
  select(Genotype, Training, recovery, recovery_low, recovery_high)

p_recovery_per_cell <- recovery_per_cell_long_raw %>%
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


# Per-cell normalized recovery summary
recovery_per_cell_long_norm <- recovery_per_cell %>%
  select(Genotype, Training, recovery_norm, recovery_norm_low, recovery_norm_high)

p_recovery_norm_per_cell <- recovery_per_cell_long_norm %>%
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
    y = "Normalized recovery = recovery / (1 - P(end_train))",
    title = "(D-norm) Normalized inter-block recovery per protocol\n0 = no recovery, 1 = full recovery to ceiling"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_normalized_recovery_per_cell.png"),
  p_recovery_norm_per_cell, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# Spaced - massed contrast plots (raw and normalized)
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
    y = "Raw recovery difference (spaced - massed)\nNegative = spaced recovers less = better retention",
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
    y = "Normalized recovery difference (spaced - massed)\nNegative = spaced recovers less = better retention",
    title = "(D-norm) Recovery contrast: normalized"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrastD_norm_diff_spaced_vs_massed.png"),
  p_contrast_D_norm, width = 8, height = 5, dpi = 300, bg = "white"
)

# ==============================================================================
# Supplementary contrasts for joint GLMM
# ------------------------------------------------------------------------------
# 1. Training sanity check: did training produce learning?
#    end_of_block_1 vs stim_1_of_block_1, within each (Genotype, Training)
#
# 2. Inter-genotype contrasts within protocol:
#    P(move | ABTL, test) vs P(move | th-tyr, test) etc., within each protocol
#
# 3. Within-training habituation slope contrast:
#    d(log-odds move)/d(stimulus_log) per (Genotype, Training, Block)
#
# 4. Block-by-block learning curve (spaced only):
#    end-of-block response across blocks 1->4
# ==============================================================================

# All four blocks assume the joint GLMM `m_joint` is loaded, df_all is in scope,
# and these objects from the main script exist:
#   massed_test_block, spaced_test_block, save_results_dir, save_fig_dir,
#   compute_memory_contrast_at()


# ==============================================================================
# 1. TRAINING SANITY CHECK: did training produce learning?
# ==============================================================================
# Compare P(move) at start of Block 1 (stim 1) vs end of Block 1 (last stim),
# within each (Genotype, Training). A large negative difference = strong
# within-training habituation. Small or null = training didn't take.
# ==============================================================================

# Last stimulus of Block 1 for each protocol
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
      point == "start_b1"                            ~ 1,
      point == "end_b1" & Training == "massed"       ~ last_b1_massed,
      point == "end_b1" & Training == "spaced"       ~ last_b1_spaced
    ),
    stimulus_log = log(stimulus),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  )

train_preds <- predict(
  m_joint,
  newdata = training_grid,
  re.form = NA,
  se.fit  = TRUE,
  type    = "link"
)

training_grid <- training_grid %>%
  mutate(
    fit_link = train_preds$fit,
    SE_link  = train_preds$se.fit,
    prob     = plogis(fit_link),
    SE_prob  = SE_link * prob * (1 - prob),
    lower    = pmax(0, prob - 1.96 * SE_prob),
    upper    = pmin(1, prob + 1.96 * SE_prob)
  )


# Within-cell learning: end_b1 - start_b1
training_learning <- training_grid %>%
  select(Genotype, Training, point, prob, SE_prob) %>%
  pivot_wider(names_from = point, values_from = c(prob, SE_prob)) %>%
  mutate(
    learning_est = prob_end_b1 - prob_start_b1,
    learning_se  = sqrt(SE_prob_end_b1^2 + SE_prob_start_b1^2),
    learning_low  = learning_est - 1.96 * learning_se,
    learning_high = learning_est + 1.96 * learning_se,
    z = learning_est / learning_se,
    p = 2 * pnorm(-abs(z))
  )

print(training_learning)

write.csv(
  training_learning,
  file.path(save_results_dir, "training_sanity_check_block1_learning.csv"),
  row.names = FALSE
)


# Plot: did training reduce response? (Negative bars = yes, learning occurred)
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
    y = "P(end of B1) - P(start of B1)\nNegative = learning occurred",
    title = "Training sanity check: within-Block 1 learning"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "sanity_check_block1_learning.png"),
  p_training_sanity, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# ==============================================================================
# 2. INTER-GENOTYPE CONTRASTS WITHIN PROTOCOL
# ==============================================================================
# Within each protocol, contrast each pair of genotypes on test response.
# Two readouts:
#   (a) at stim 1 of test block
#   (b) mean over all 8 test stimuli
# Pairwise: ABTL vs th-th2-tyr, ABTL vs th-tyr, th-th2-tyr vs th-tyr
# ==============================================================================

# Pull per-cell test response predictions from our existing res_A (stim 1) and
# from res_C (all test stim) helpers. Already computed earlier in the script.

genotype_pairs <- tribble(
  ~g1,             ~g2,
  "ABTL",          "th, th2, tyr",
  "ABTL",          "th, tyr",
  "th, th2, tyr",  "th, tyr"
)


compute_genotype_contrast <- function(per_cell_df, pair_table) {
  # per_cell_df has columns: Genotype, Training, prob, SE
  # pair_table has columns g1, g2
  
  out <- pair_table %>%
    tidyr::crossing(Training = unique(per_cell_df$Training)) %>%
    rowwise() %>%
    mutate(
      prob_g1 = per_cell_df$prob[per_cell_df$Genotype == g1 & per_cell_df$Training == Training],
      SE_g1   = per_cell_df$SE[  per_cell_df$Genotype == g1 & per_cell_df$Training == Training],
      prob_g2 = per_cell_df$prob[per_cell_df$Genotype == g2 & per_cell_df$Training == Training],
      SE_g2   = per_cell_df$SE[  per_cell_df$Genotype == g2 & per_cell_df$Training == Training],
      diff_est  = prob_g1 - prob_g2,
      diff_se   = sqrt(SE_g1^2 + SE_g2^2),
      diff_low  = diff_est - 1.96 * diff_se,
      diff_high = diff_est + 1.96 * diff_se,
      z         = diff_est / diff_se,
      p         = 2 * pnorm(-abs(z))
    ) %>%
    ungroup() %>%
    mutate(comparison = paste0(g1, " vs ", g2))
  
  out
}


# (a) Stim 1 of test block
# res_A$per_cell already has prob, SE per (Genotype, Training, Block)
gen_contrast_stim1 <- compute_genotype_contrast(
  res_A$per_cell %>% select(Genotype, Training, prob, SE),
  genotype_pairs
) %>%
  mutate(readout = "stim1_test")


# (b) Mean over all 8 test stimuli
# Reuse per_cell_C built in the main script (mean over stim 1..8)
gen_contrast_meanall <- compute_genotype_contrast(
  per_cell_C %>% select(Genotype, Training, prob, SE = SE_mean),
  genotype_pairs
) %>%
  mutate(readout = "mean_all_test")


genotype_contrasts_all <- bind_rows(gen_contrast_stim1, gen_contrast_meanall)

print(genotype_contrasts_all)

write.csv(
  genotype_contrasts_all,
  file.path(save_results_dir, "inter_genotype_contrasts_within_protocol.csv"),
  row.names = FALSE
)


# Plot: inter-genotype differences within each protocol
p_inter_gen <- genotype_contrasts_all %>%
  ggplot(aes(x = comparison, y = diff_est,
             ymin = diff_low, ymax = diff_high,
             color = Training)) +
  facet_wrap(~ readout, ncol = 1, scales = "free_x") +
  geom_pointrange(linewidth = 0.8, size = 0.7,
                  position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Genotype difference in P(move) within protocol",
    title = "Inter-genotype contrasts within protocol",
    color = "Protocol"
  ) +
  theme(legend.position = "top")

ggsave(
  file.path(save_fig_dir, "inter_genotype_contrasts.png"),
  p_inter_gen, width = 10, height = 8, dpi = 300, bg = "white"
)


# ==============================================================================
# 3. WITHIN-TRAINING HABITUATION SLOPE CONTRAST
# ==============================================================================
# emtrends() gives the slope of stimulus_log within each (Genotype, Training,
# Block) cell. Negative slope = response declines as stimulus increases.
# Steeper (more negative) = faster habituation.
# We focus on Block 1 since both protocols have Block 1.
# ==============================================================================

slope_emm <- emtrends(
  m_joint,
  specs = ~ Genotype * Training | Block,
  var   = "stimulus_log",
  at    = list(Block = "1")
)

slope_df <- as_tibble(slope_emm) %>%
  rename(slope = stimulus_log.trend,
         SE    = SE,
         lower = asymp.LCL,
         upper = asymp.UCL) %>%
  select(Genotype, Training, Block, slope, SE, lower, upper)

print(slope_df)

write.csv(
  slope_df,
  file.path(save_results_dir, "training_slope_block1.csv"),
  row.names = FALSE
)


# Spaced vs massed slope contrast within each genotype
# Spaced vs massed slope contrast within each genotype
slope_pairs <- pairs(slope_emm, by = "Genotype", reverse = TRUE)

# confint() adds asymp.LCL / asymp.UCL columns
slope_pairs_df <- as_tibble(confint(slope_pairs)) %>%
  rename(diff_est  = estimate,
         diff_se   = SE,
         diff_low  = asymp.LCL,
         diff_high = asymp.UCL) %>%
  # Add z and p back in by joining with the original pairs object
  left_join(
    as_tibble(slope_pairs) %>%
      select(Genotype, contrast, z.ratio, p.value),
    by = c("Genotype", "contrast")
  ) %>%
  select(Genotype, contrast, diff_est, diff_se, diff_low, diff_high,
         z.ratio, p.value)

print(slope_pairs_df)
print(slope_pairs_df)

write.csv(
  slope_pairs_df,
  file.path(save_results_dir, "training_slope_diff_spaced_vs_massed.csv"),
  row.names = FALSE
)


# Plot: training slopes
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
    y = "Habituation slope d(logit P)/d(log stim) in Block 1\nMore negative = faster habituation",
    title = "Within-Block-1 habituation rate per protocol"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "training_slope_block1_per_protocol.png"),
  p_slope, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# ==============================================================================
# 4. BLOCK-BY-BLOCK LEARNING CURVE (SPACED ONLY)
# ==============================================================================
# Within spaced training, look at end-of-block response across blocks 1-4.
# If asymptote keeps dropping -> training still ongoing.
# If it plateaus -> training has saturated.
# ==============================================================================

# Determine the last stimulus of each spaced training block (1..4)
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
  m_joint,
  newdata = spaced_endblock_grid,
  re.form = NA,
  se.fit  = TRUE,
  type    = "link"
)

spaced_endblock_grid <- spaced_endblock_grid %>%
  mutate(
    fit_link = endblock_preds$fit,
    SE_link  = endblock_preds$se.fit,
    prob     = plogis(fit_link),
    SE_prob  = SE_link * prob * (1 - prob),
    lower    = pmax(0, prob - 1.96 * SE_prob),
    upper    = pmin(1, prob + 1.96 * SE_prob),
    block_num = as.integer(as.character(Block))
  )

print(spaced_endblock_grid %>%
        select(Genotype, Block, block_num, stimulus, prob, lower, upper))

write.csv(
  spaced_endblock_grid %>%
    select(Genotype, Block, block_num, stimulus, prob, SE_prob, lower, upper),
  file.path(save_results_dir, "spaced_endblock_response.csv"),
  row.names = FALSE
)


# Block-to-block change in end-of-block response
spaced_block_diffs <- spaced_endblock_grid %>%
  arrange(Genotype, block_num) %>%
  group_by(Genotype) %>%
  mutate(
    prob_prev    = lag(prob),
    SE_prev      = lag(SE_prob),
    delta        = prob - prob_prev,
    delta_se     = sqrt(SE_prob^2 + SE_prev^2),
    delta_low    = delta - 1.96 * delta_se,
    delta_high   = delta + 1.96 * delta_se,
    transition   = paste0("B", block_num - 1, "->B", block_num)
  ) %>%
  ungroup() %>%
  filter(!is.na(prob_prev)) %>%
  select(Genotype, transition, delta, delta_se, delta_low, delta_high)

print(spaced_block_diffs)

write.csv(
  spaced_block_diffs,
  file.path(save_results_dir, "spaced_block_to_block_endblock_changes.csv"),
  row.names = FALSE
)


# Plot: end-of-block response across spaced blocks
p_spaced_endblock <- spaced_endblock_grid %>%
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


# Plot: block-to-block changes
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
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  file.path(save_fig_dir, "spaced_block_to_block_changes.png"),
  p_spaced_diffs, width = 10, height = 4.5, dpi = 300, bg = "white"
)


message("Supplementary contrasts 1-4 complete.")

# ==============================================================================
# Done
# ==============================================================================
message("Joint GLMM analysis complete.")
message("Results saved under: ", save_results_dir)
message("Figures saved under: ", save_fig_dir)