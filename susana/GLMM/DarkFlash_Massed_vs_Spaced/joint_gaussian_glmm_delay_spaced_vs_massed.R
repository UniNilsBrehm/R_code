###############################################################################
# Joint Gaussian GLMM Analysis: Spaced vs Massed Training - Response Delay
# Author: Nils Brehm
# Date: 2026
#
# Description:
#   Joint Gaussian mixed model combining spaced and massed training datasets into
#   one model to test the spaced-vs-massed memory contrast directly for:
#
#     delay = response delay category, treated numerically
#
#   delay can only take values: 0, 1, 2, 3, 4
#
#   Model:
#     delay_num ~ Genotype * Training * Block * stimulus_log + (1 | animal)
#
#   Model family:
#     Gaussian identity-link mixed model
#
#   Memory contrasts per genotype:
#     (A) predicted delay at first test stimulus, spaced vs massed
#     (B) mean predicted delay over first 3 test stimuli, spaced vs massed
#     (C) mean predicted delay over all shared test stimuli, spaced vs massed
#     (D) inter-block change in predicted delay, spaced vs massed
#
#   Positive spaced-massed difference => spaced fish have longer/slower delay
#   in the test block than massed fish.
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

save_results_dir <- file.path(base_dir, "results", "gaussian_joint_delay")
save_fig_dir     <- file.path(base_dir, "figs",    "gaussian_joint_delay")

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

# Keep valid delay values only.
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

p_dist <- ggplot(df_all, aes(x = factor(delay_num, levels = 0:4))) +
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
# 4. Fit the joint Gaussian GLMM
# ==============================================================================
message("Fitting joint Gaussian GLMM for delay...")

m_joint_delay <- glmmTMB(
  delay_num ~ Genotype * Training * Block * stimulus_log + (1 | animal),
  family = gaussian(link = "identity"),
  data   = df_all
)

saveRDS(
  m_joint_delay,
  file.path(save_results_dir, "joint_gaussian_glmm_spaced_vs_massed_delay.rds")
)

capture.output(
  summary(m_joint_delay),
  file = file.path(save_results_dir, "summary_results_delay.txt")
)

print(summary(m_joint_delay))


# ==============================================================================
# 5. Validate model
# ==============================================================================
message("Validating joint Gaussian delay model residuals...")
validate_model(m_joint_delay, df_all)


# ==============================================================================
# 6. Plot habituation curves: predicted delay across stimuli
# ==============================================================================
message("Plotting joint Gaussian delay habituation curves...")

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
  m_joint_delay,
  newdata = new_data_joint,
  re.form = NA,
  se.fit  = TRUE,
  type    = "response"
)

new_data_joint <- new_data_joint %>%
  mutate(
    fit     = pred_joint$fit,
    CI_low  = pred_joint$fit - 1.96 * pred_joint$se.fit,
    CI_high = pred_joint$fit + 1.96 * pred_joint$se.fit
  )

raw_summary_joint <- df_all %>%
  group_by(Training, Block, Genotype, stimulus) %>%
  summarise(
    mean_delay = mean(delay_num, na.rm = TRUE),
    .groups = "drop"
  )


p_massed_curves <- ggplot(
  new_data_joint %>% filter(Training == "massed"),
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_joint %>% filter(Training == "massed"),
    aes(x = stimulus, y = mean_delay, color = Genotype),
    inherit.aes = FALSE, alpha = 0.25, size = 0.6
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Predicted response delay",
    title = "Joint Gaussian GLMM: Massed training"
  ) +
  theme(legend.position = "none")

p_spaced_curves <- ggplot(
  new_data_joint %>% filter(Training == "spaced"),
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_joint %>% filter(Training == "spaced"),
    aes(x = stimulus, y = mean_delay, color = Genotype),
    inherit.aes = FALSE, alpha = 0.25, size = 0.6
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Predicted response delay",
    title = "Joint Gaussian GLMM: Spaced training"
  ) +
  theme(legend.position = "top")

print(p_spaced_curves)
print(p_massed_curves)

ggsave(
  file.path(save_fig_dir, "joint_gaussian_delay_curves_massed.png"),
  p_massed_curves, width = 10, height = 12, dpi = 300, bg = "white"
)
ggsave(
  file.path(save_fig_dir, "joint_gaussian_delay_curves_spaced.png"),
  p_spaced_curves, width = 14, height = 12, dpi = 300, bg = "white"
)


# ==============================================================================
# 7. Helper functions for delay memory contrasts
# ==============================================================================

genotypes <- levels(df_all$Genotype)

compute_delay_contrast_at <- function(model, stim_log_value) {
  
  emm <- emmeans(
    model,
    specs = ~ Genotype * Training * Block,
    at = list(
      stimulus_log = stim_log_value,
      Block        = c(massed_test_block, spaced_test_block)
    ),
    type = "response"
  )
  
  emm_df <- as_tibble(confint(emm)) %>%
    filter(
      (Training == "massed" & Block == massed_test_block) |
        (Training == "spaced" & Block == spaced_test_block)
    )
  
  estimate_col <- c("emmean", "response", "estimate")[
    c("emmean", "response", "estimate") %in% names(emm_df)
  ][1]
  
  lower_col <- c("lower.CL", "asymp.LCL", "LCL")[
    c("lower.CL", "asymp.LCL", "LCL") %in% names(emm_df)
  ][1]
  
  upper_col <- c("upper.CL", "asymp.UCL", "UCL")[
    c("upper.CL", "asymp.UCL", "UCL") %in% names(emm_df)
  ][1]
  
  if (is.na(estimate_col) || is.na(lower_col) || is.na(upper_col)) {
    stop(
      "Could not identify estimate/CI columns in emmeans output. Columns are: ",
      paste(names(emm_df), collapse = ", ")
    )
  }
  
  per_cell <- emm_df %>%
    transmute(
      Genotype,
      Training,
      Block,
      delay_pred     = .data[[estimate_col]],
      SE             = SE,
      lower          = .data[[lower_col]],
      upper          = .data[[upper_col]],
      stim_log_value = stim_log_value,
      stimulus       = exp(stim_log_value)
    )
  
  contrasts_df <- per_cell %>%
    select(Genotype, Training, delay_pred, SE) %>%
    pivot_wider(
      names_from  = Training,
      values_from = c(delay_pred, SE)
    ) %>%
    mutate(
      diff_est  = delay_pred_spaced - delay_pred_massed,
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
      delay_pred_mean = mean(delay_pred),
      SE_mean         = sqrt(sum(SE^2)) / n(),
      lower           = delay_pred_mean - 1.96 * SE_mean,
      upper           = delay_pred_mean + 1.96 * SE_mean,
      .groups         = "drop"
    ) %>%
    rename(delay_pred = delay_pred_mean)
}


# ==============================================================================
# 8. Memory contrasts A, B, C
# ==============================================================================

res_A <- compute_delay_contrast_at(m_joint_delay, stim_log_value = 0)

contrast_A <- res_A$contrasts %>%
  select(Genotype, diff_est, diff_se, diff_low, diff_high, z, p,
         stim_log_value, stimulus) %>%
  mutate(contrast = "A_stim1")

res_B1 <- compute_delay_contrast_at(m_joint_delay, stim_log_value = log(1))
res_B2 <- compute_delay_contrast_at(m_joint_delay, stim_log_value = log(2))
res_B3 <- compute_delay_contrast_at(m_joint_delay, stim_log_value = log(3))

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

# Shared test stimuli: 1..8 for both protocols.
test_stims_for_C <- 1:8

res_C_list <- lapply(test_stims_for_C, function(s) {
  compute_delay_contrast_at(m_joint_delay, stim_log_value = log(s))
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

write.csv(res_A$per_cell, file.path(save_results_dir, "per_cell_test_delay_stim1.csv"), row.names = FALSE)
write.csv(bind_rows(res_B1$per_cell, res_B2$per_cell, res_B3$per_cell), file.path(save_results_dir, "per_cell_test_delay_stim1to3.csv"), row.names = FALSE)
write.csv(per_cell_B, file.path(save_results_dir, "per_cell_test_delay_stim1to3_aggregated.csv"), row.names = FALSE)
write.csv(per_cell_C, file.path(save_results_dir, "per_cell_test_delay_meanAllStim.csv"), row.names = FALSE)
write.csv(contrast_A, file.path(save_results_dir, "contrast_A_test_stim1_spaced_vs_massed_delay.csv"), row.names = FALSE)
write.csv(contrast_B, file.path(save_results_dir, "contrast_B_test_meanStim1to3_spaced_vs_massed_delay.csv"), row.names = FALSE)
write.csv(contrast_C, file.path(save_results_dir, "contrast_C_test_meanAllStim_spaced_vs_massed_delay.csv"), row.names = FALSE)

print(contrast_A)
print(contrast_B)
print(contrast_C)


# ==============================================================================
# 9. Plots for contrasts A, B, C
# ==============================================================================

p_per_cell_A <- res_A$per_cell %>%
  ggplot(aes(x = Training, y = delay_pred, ymin = lower, ymax = upper, color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 13) +
  labs(x = NULL, y = "Predicted delay at stim 1 of test block", title = "(A) Test response delay at first stimulus") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "contrastA_per_cell_test_stim1_delay.png"), p_per_cell_A, width = 12, height = 4.5, dpi = 300, bg = "white")

p_contrast_A <- contrast_A %>%
  ggplot(aes(x = Genotype, y = diff_est, ymin = diff_low, ymax = diff_high, color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(x = NULL, y = "Delay difference at stim 1 (spaced - massed)", title = "(A) Headline contrast: stim 1 of test block") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "contrastA_diff_spaced_vs_massed_delay.png"), p_contrast_A, width = 8, height = 5, dpi = 300, bg = "white")

p_per_cell_B <- per_cell_B %>%
  ggplot(aes(x = Training, y = delay_pred, ymin = lower, ymax = upper, color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 13) +
  labs(x = NULL, y = "Mean predicted delay over stim 1-3 of test block", title = "(B) Test response delay averaged over first 3 stimuli") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "contrastB_per_cell_test_stim1to3_delay.png"), p_per_cell_B, width = 12, height = 4.5, dpi = 300, bg = "white")

p_contrast_B <- contrast_B %>%
  ggplot(aes(x = Genotype, y = diff_est, ymin = diff_low, ymax = diff_high, color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(x = NULL, y = "Mean delay difference over stim 1-3 (spaced - massed)", title = "(B) Averaged contrast: mean of stim 1-3 of test block") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "contrastB_diff_spaced_vs_massed_delay.png"), p_contrast_B, width = 8, height = 5, dpi = 300, bg = "white")

p_per_cell_C <- per_cell_C %>%
  ggplot(aes(x = Training, y = delay_pred, ymin = lower, ymax = upper, color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 13) +
  labs(x = NULL, y = "Mean predicted delay over all test-block stim 1-8", title = "(C) Test response delay averaged over all shared test stimuli") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "contrastC_per_cell_test_meanAllStim_delay.png"), p_per_cell_C, width = 12, height = 4.5, dpi = 300, bg = "white")

p_contrast_C <- contrast_C %>%
  ggplot(aes(x = Genotype, y = diff_est, ymin = diff_low, ymax = diff_high, color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(x = NULL, y = "Mean delay difference over all test stim 1-8 (spaced - massed)", title = "(C) Averaged contrast: mean of all shared test stimuli") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "contrastC_diff_spaced_vs_massed_delay.png"), p_contrast_C, width = 8, height = 5, dpi = 300, bg = "white")

all_contrasts <- bind_rows(contrast_A, contrast_B, contrast_C) %>%
  mutate(contrast = factor(contrast, levels = c("A_stim1", "B_mean_stim1_to_3", "C_mean_all_test_stim")))

write.csv(all_contrasts, file.path(save_results_dir, "ALL_contrasts_spaced_vs_massed_delay.csv"), row.names = FALSE)

p_all_contrasts <- all_contrasts %>%
  ggplot(aes(x = Genotype, y = diff_est, ymin = diff_low, ymax = diff_high, color = Genotype)) +
  facet_wrap(~ contrast, ncol = 1) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 13) +
  labs(x = NULL, y = "Spaced - massed predicted-delay difference", title = "All Gaussian delay memory contrasts side-by-side") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "ALL_contrasts_spaced_vs_massed_delay.png"), p_all_contrasts, width = 9, height = 10, dpi = 300, bg = "white")


# ==============================================================================
# 10. Contrast D: Inter-block change/recovery in predicted delay
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

recovery_grid <- bind_rows(
  tidyr::expand_grid(Genotype = levels(df_all$Genotype), point = c("end_train", "start_test")) %>%
    mutate(Training = "massed", Block = ifelse(point == "end_train", last_block_massed, massed_test_block), stimulus = ifelse(point == "end_train", last_stim_massed, 1)),
  tidyr::expand_grid(Genotype = levels(df_all$Genotype), point = c("end_train", "start_test")) %>%
    mutate(Training = "spaced", Block = ifelse(point == "end_train", last_block_spaced, spaced_test_block), stimulus = ifelse(point == "end_train", last_stim_spaced, 1))
) %>%
  mutate(
    stimulus_log = log(stimulus),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  )

print(recovery_grid)

preds <- predict(m_joint_delay, newdata = recovery_grid, re.form = NA, se.fit = TRUE, type = "response")

recovery_grid <- recovery_grid %>%
  mutate(delay_pred = preds$fit, SE = preds$se.fit, lower = delay_pred - 1.96 * SE, upper = delay_pred + 1.96 * SE)

per_cell_D <- recovery_grid %>%
  select(Genotype, Training, Block, point, stimulus, delay_pred, SE, lower, upper)

print(per_cell_D)

write.csv(per_cell_D, file.path(save_results_dir, "per_cell_endTrain_startTest_delay.csv"), row.names = FALSE)

delay_change_per_cell <- per_cell_D %>%
  select(Genotype, Training, point, delay_pred, SE) %>%
  pivot_wider(names_from = point, values_from = c(delay_pred, SE)) %>%
  mutate(
    delay_change      = delay_pred_start_test - delay_pred_end_train,
    delay_change_SE   = sqrt(SE_start_test^2 + SE_end_train^2),
    delay_change_low  = delay_change - 1.96 * delay_change_SE,
    delay_change_high = delay_change + 1.96 * delay_change_SE
  )

write.csv(delay_change_per_cell, file.path(save_results_dir, "delay_change_per_cell.csv"), row.names = FALSE)

contrast_D <- delay_change_per_cell %>%
  select(Genotype, Training, delay_change, delay_change_SE) %>%
  pivot_wider(names_from = Training, values_from = c(delay_change, delay_change_SE)) %>%
  mutate(
    diff_est  = delay_change_spaced - delay_change_massed,
    diff_se   = sqrt(delay_change_SE_spaced^2 + delay_change_SE_massed^2),
    diff_low  = diff_est - 1.96 * diff_se,
    diff_high = diff_est + 1.96 * diff_se,
    z         = diff_est / diff_se,
    p         = 2 * pnorm(-abs(z)),
    contrast  = "D_delay_change"
  ) %>%
  select(Genotype, diff_est, diff_se, diff_low, diff_high, z, p, contrast)

print(contrast_D)
write.csv(contrast_D, file.path(save_results_dir, "contrast_D_delay_change_spaced_vs_massed.csv"), row.names = FALSE)

p_per_cell_D <- per_cell_D %>%
  mutate(point_label = factor(point, levels = c("end_train", "start_test"), labels = c("End of training", "Start of test"))) %>%
  ggplot(aes(x = point_label, y = delay_pred, ymin = lower, ymax = upper, color = Genotype, group = Training)) +
  facet_grid(Genotype ~ Training) +
  geom_pointrange(linewidth = 0.7, size = 0.6) +
  geom_line(linewidth = 0.6, linetype = "dashed", alpha = 0.7) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 12) +
  labs(x = NULL, y = "Predicted response delay", title = "End-of-training vs start-of-test predicted delay per protocol") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(save_fig_dir, "contrastD_per_cell_endTrain_startTest_delay.png"), p_per_cell_D, width = 9, height = 11, dpi = 300, bg = "white")

p_delay_change_per_cell <- delay_change_per_cell %>%
  select(Genotype, Training, delay_change, delay_change_low, delay_change_high) %>%
  ggplot(aes(x = Training, y = delay_change, ymin = delay_change_low, ymax = delay_change_high, color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(x = NULL, y = "Delay change = predicted delay(start_test) - predicted delay(end_train)", title = "(D) Inter-block predicted-delay change per protocol") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "contrastD_delay_change_per_cell.png"), p_delay_change_per_cell, width = 10, height = 4.5, dpi = 300, bg = "white")

p_contrast_D <- contrast_D %>%
  ggplot(aes(x = Genotype, y = diff_est, ymin = diff_low, ymax = diff_high, color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(x = NULL, y = "Delay-change difference (spaced - massed)", title = "(D) Predicted-delay change contrast") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "contrastD_delay_change_diff_spaced_vs_massed.png"), p_contrast_D, width = 8, height = 5, dpi = 300, bg = "white")


# ==============================================================================
# 11. Supplementary contrasts
# ==============================================================================

# 11.1 Training sanity check: start vs end of Block 1
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

train_preds <- predict(m_joint_delay, newdata = training_grid, re.form = NA, se.fit = TRUE, type = "response")

training_grid <- training_grid %>%
  mutate(delay_pred = train_preds$fit, SE = train_preds$se.fit, lower = delay_pred - 1.96 * SE, upper = delay_pred + 1.96 * SE)

training_learning <- training_grid %>%
  select(Genotype, Training, point, delay_pred, SE) %>%
  pivot_wider(names_from = point, values_from = c(delay_pred, SE)) %>%
  mutate(
    learning_est  = delay_pred_end_b1 - delay_pred_start_b1,
    learning_se   = sqrt(SE_end_b1^2 + SE_start_b1^2),
    learning_low  = learning_est - 1.96 * learning_se,
    learning_high = learning_est + 1.96 * learning_se,
    z = learning_est / learning_se,
    p = 2 * pnorm(-abs(z))
  )

print(training_learning)
write.csv(training_learning, file.path(save_results_dir, "training_sanity_check_block1_delay.csv"), row.names = FALSE)

p_training_sanity <- training_learning %>%
  ggplot(aes(x = Training, y = learning_est, ymin = learning_low, ymax = learning_high, color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(x = NULL, y = "Predicted delay(end B1) - predicted delay(start B1)", title = "Training sanity check: within-Block 1 change in predicted delay") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "sanity_check_block1_delay.png"), p_training_sanity, width = 10, height = 4.5, dpi = 300, bg = "white")


# ------------------------------------------------------------------------------
# 11.2 Within-training delay slope contrast
# ------------------------------------------------------------------------------
slope_emm <- emtrends(
  m_joint_delay,
  specs = ~ Genotype * Training | Block,
  var   = "stimulus_log",
  at    = list(Block = "1")
)

# Per-cell slopes ---------------------------------------------------------------
slope_df_raw <- as_tibble(confint(slope_emm))

lower_col <- c("lower.CL", "asymp.LCL", "LCL")[
  c("lower.CL", "asymp.LCL", "LCL") %in% names(slope_df_raw)
][1]

upper_col <- c("upper.CL", "asymp.UCL", "UCL")[
  c("upper.CL", "asymp.UCL", "UCL") %in% names(slope_df_raw)
][1]

slope_df <- slope_df_raw %>%
  transmute(
    Genotype,
    Training,
    Block,
    slope = stimulus_log.trend,
    SE,
    lower = .data[[lower_col]],
    upper = .data[[upper_col]]
  )

print(slope_df)

write.csv(
  slope_df,
  file.path(save_results_dir, "training_slope_block1_delay.csv"),
  row.names = FALSE
)


# Spaced vs massed slope contrasts ---------------------------------------------
slope_pairs <- pairs(slope_emm, by = "Genotype", reverse = TRUE)

slope_pairs_raw <- as_tibble(confint(slope_pairs))

lower_col <- c("lower.CL", "asymp.LCL", "LCL")[
  c("lower.CL", "asymp.LCL", "LCL") %in% names(slope_pairs_raw)
][1]

upper_col <- c("upper.CL", "asymp.UCL", "UCL")[
  c("upper.CL", "asymp.UCL", "UCL") %in% names(slope_pairs_raw)
][1]

slope_pairs_df <- slope_pairs_raw %>%
  transmute(
    Genotype,
    contrast,
    diff_est  = estimate,
    diff_se   = SE,
    diff_low  = .data[[lower_col]],
    diff_high = .data[[upper_col]]
  ) %>%
  left_join(
    as_tibble(slope_pairs) %>%
      select(Genotype, contrast, z.ratio, p.value),
    by = c("Genotype", "contrast")
  ) %>%
  select(
    Genotype, contrast,
    diff_est, diff_se, diff_low, diff_high,
    z.ratio, p.value
  )

print(slope_pairs_df)

write.csv(
  slope_pairs_df,
  file.path(save_results_dir, "training_slope_diff_spaced_vs_massed_delay.csv"),
  row.names = FALSE
)


# Plot slopes ------------------------------------------------------------------
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
    y = "Delay slope d(delay)/d(log stim) in Block 1",
    title = "Within-Block-1 delay slope per protocol"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "training_slope_block1_per_protocol_delay.png"),
  p_slope, width = 10, height = 4.5, dpi = 300, bg = "white"
)
# 11.3 Block-by-block learning curve, spaced only
spaced_last_stim_per_block <- df_all %>%
  filter(Training == "spaced", Block %in% c("1", "2", "3", "4")) %>%
  group_by(Block) %>%
  summarise(last_stim = max(stimulus, na.rm = TRUE), .groups = "drop")

print(spaced_last_stim_per_block)

spaced_endblock_grid <- tidyr::expand_grid(
  Genotype = levels(df_all$Genotype),
  Block    = c("1", "2", "3", "4")
) %>%
  left_join(spaced_last_stim_per_block %>% mutate(Block = as.character(Block)), by = "Block") %>%
  mutate(
    Training     = "spaced",
    stimulus     = last_stim,
    stimulus_log = log(stimulus),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  )

endblock_preds <- predict(m_joint_delay, newdata = spaced_endblock_grid, re.form = NA, se.fit = TRUE, type = "response")

spaced_endblock_grid <- spaced_endblock_grid %>%
  mutate(
    delay_pred = endblock_preds$fit,
    SE         = endblock_preds$se.fit,
    lower      = delay_pred - 1.96 * SE,
    upper      = delay_pred + 1.96 * SE,
    block_num  = as.integer(as.character(Block))
  )

print(spaced_endblock_grid %>% select(Genotype, Block, block_num, stimulus, delay_pred, lower, upper))
write.csv(spaced_endblock_grid %>% select(Genotype, Block, block_num, stimulus, delay_pred, SE, lower, upper), file.path(save_results_dir, "spaced_endblock_delay.csv"), row.names = FALSE)

spaced_block_diffs <- spaced_endblock_grid %>%
  arrange(Genotype, block_num) %>%
  group_by(Genotype) %>%
  mutate(
    delay_pred_prev = lag(delay_pred),
    SE_prev         = lag(SE),
    delta           = delay_pred - delay_pred_prev,
    delta_se        = sqrt(SE^2 + SE_prev^2),
    delta_low       = delta - 1.96 * delta_se,
    delta_high      = delta + 1.96 * delta_se,
    transition      = paste0("B", block_num - 1, "->B", block_num)
  ) %>%
  ungroup() %>%
  filter(!is.na(delay_pred_prev)) %>%
  select(Genotype, transition, delta, delta_se, delta_low, delta_high)

print(spaced_block_diffs)
write.csv(spaced_block_diffs, file.path(save_results_dir, "spaced_block_to_block_endblock_changes_delay.csv"), row.names = FALSE)

p_spaced_endblock <- spaced_endblock_grid %>%
  ggplot(aes(x = block_num, y = delay_pred, ymin = lower, ymax = upper, color = Genotype, group = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.7, size = 0.6) +
  geom_line(linewidth = 0.6) +
  scale_x_continuous(breaks = 1:4) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = 0:4) +
  theme_pubr(base_size = 13) +
  labs(x = "Spaced training block", y = "Predicted delay at end of block", title = "Spaced training: end-of-block predicted delay across blocks 1-4") +
  theme(legend.position = "none")

ggsave(file.path(save_fig_dir, "spaced_endblock_progression_delay.png"), p_spaced_endblock, width = 10, height = 4.5, dpi = 300, bg = "white")

p_spaced_diffs <- spaced_block_diffs %>%
  ggplot(aes(x = transition, y = delta, ymin = delta_low, ymax = delta_high, color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.7, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(x = "Block-to-block transition", y = "Change in end-of-block predicted delay", title = "Spaced training: incremental delay change per block") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(save_fig_dir, "spaced_block_to_block_changes_delay.png"), p_spaced_diffs, width = 10, height = 4.5, dpi = 300, bg = "white")


# ==============================================================================
# Done
# ==============================================================================
message("Joint Gaussian GLMM delay analysis complete.")
message("Results saved under: ", save_results_dir)
message("Figures saved under: ", save_fig_dir)
