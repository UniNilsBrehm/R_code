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
# 6. Helper functions: expected delay from ordinal model
# ==============================================================================

delay_scores <- 0:4

# Robust extraction of probabilities from predict.clmm(type = "prob").
# For newdata rows, this should return one probability per delay class.
predict_expected_delay <- function(model, newdata) {
  beta  <- model$beta
  alpha <- model$alpha
  
  Terms <- delete.response(terms(model))
  X <- model.matrix(Terms, newdata)
  
  # keep only fixed-effect columns that exist in model$beta
  X <- X[, names(beta), drop = FALSE]
  
  eta <- as.vector(X %*% beta)
  
  # cumulative probabilities: P(Y <= category)
  cumprob <- sapply(alpha, function(a) plogis(a - eta))
  
  # category probabilities for delay 0,1,2,3,4
  probs <- cbind(
    cumprob[, 1],
    cumprob[, 2] - cumprob[, 1],
    cumprob[, 3] - cumprob[, 2],
    cumprob[, 4] - cumprob[, 3],
    1 - cumprob[, 4]
  )
  
  as.vector(probs %*% 0:4)
}


# emmeans helper for expected delay.
# emmeans mode = "mean.class" returns expected class number for ordinal models.
# With ordered levels 0:4, this is usually class index 1:5, so we subtract 1
# to report expected delay on the original 0:4 scale.
get_expected_delay_emm <- function(model, stim_log_value) {

  emm <- emmeans(
    model,
    specs = ~ Genotype * Training * Block,
    at = list(
      stimulus_log = stim_log_value,
      Block        = c(massed_test_block, spaced_test_block)
    ),
    mode = "mean.class"
  )

  emm_df <- as_tibble(confint(emm)) %>%
    filter(
      (Training == "massed" & Block == massed_test_block) |
        (Training == "spaced" & Block == spaced_test_block)
    )

  # Find estimate and CI columns robustly.
  estimate_col <- c("emmean", "mean.class", "estimate")[c("emmean", "mean.class", "estimate") %in% names(emm_df)][1]
  lower_col    <- c("asymp.LCL", "lower.CL", "LCL")[c("asymp.LCL", "lower.CL", "LCL") %in% names(emm_df)][1]
  upper_col    <- c("asymp.UCL", "upper.CL", "UCL")[c("asymp.UCL", "upper.CL", "UCL") %in% names(emm_df)][1]

  if (is.na(estimate_col) || is.na(lower_col) || is.na(upper_col)) {
    stop("Could not identify estimate/CI columns in emmeans output. Columns are: ",
         paste(names(emm_df), collapse = ", "))
  }

  # Convert from expected class number 1:5 to expected delay 0:4.
  # If your emmeans version already reports 0:4, remove the '- 1' below.
  per_cell <- emm_df %>%
    transmute(
      Genotype,
      Training,
      Block,
      expected_delay = .data[[estimate_col]] - 1,
      SE             = SE,
      lower          = .data[[lower_col]] - 1,
      upper          = .data[[upper_col]] - 1,
      stim_log_value = stim_log_value,
      stimulus       = exp(stim_log_value)
    )

  contrasts_df <- per_cell %>%
    select(Genotype, Training, expected_delay, SE) %>%
    pivot_wider(
      names_from  = Training,
      values_from = c(expected_delay, SE)
    ) %>%
    mutate(
      diff_est  = expected_delay_spaced - expected_delay_massed,
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
      expected_delay_mean = mean(expected_delay),
      SE_mean             = sqrt(sum(SE^2)) / n(),
      lower               = expected_delay_mean - 1.96 * SE_mean,
      upper               = expected_delay_mean + 1.96 * SE_mean,
      .groups             = "drop"
    ) %>%
    rename(expected_delay = expected_delay_mean)
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

# Point predictions only. CI ribbons for ordinal expected values require
# additional delta-method or bootstrap code and are omitted here.
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


# ==============================================================================
# 8. Memory contrasts A, B, C on expected-delay scale
# ==============================================================================

# A: test stimulus 1
res_A <- get_expected_delay_emm(m_joint_delay, stim_log_value = 0)

contrast_A <- res_A$contrasts %>%
  select(Genotype, diff_est, diff_se, diff_low, diff_high, z, p,
         stim_log_value, stimulus) %>%
  mutate(contrast = "A_stim1")


# B: mean over test stimuli 1..3
res_B1 <- get_expected_delay_emm(m_joint_delay, stim_log_value = log(1))
res_B2 <- get_expected_delay_emm(m_joint_delay, stim_log_value = log(2))
res_B3 <- get_expected_delay_emm(m_joint_delay, stim_log_value = log(3))

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
  get_expected_delay_emm(m_joint_delay, stim_log_value = log(s))
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
  file.path(save_results_dir, "per_cell_test_expected_delay_stim1.csv"),
  row.names = FALSE
)

write.csv(
  bind_rows(res_B1$per_cell, res_B2$per_cell, res_B3$per_cell),
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

print(contrast_A)
print(contrast_B)
print(contrast_C)


# ==============================================================================
# 9. Plots for contrasts A, B, C
# ==============================================================================

p_per_cell_A <- res_A$per_cell %>%
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


all_contrasts <- bind_rows(contrast_A, contrast_B, contrast_C) %>%
  mutate(contrast = factor(
    contrast,
    levels = c("A_stim1", "B_mean_stim1_to_3", "C_mean_all_test_stim")
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


# ==============================================================================
# 10. Contrast D: Inter-block change/recovery in expected delay
# ==============================================================================
# For each (Genotype, Training), compute:
#   delay_change = expected_delay(start_test) - expected_delay(end_train)
#
# Then contrast spaced - massed within each Genotype.
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


# Use emmeans helper at two points per protocol.
# We predict expected delay and approximate SE using emmeans at mean.class scale.
get_expected_delay_at_custom_points <- function(model, grid_df) {
  out <- lapply(seq_len(nrow(grid_df)), function(i) {
    row_i <- grid_df[i, ]
    emm <- emmeans(
      model,
      specs = ~ Genotype,
      at = list(
        Training     = as.character(row_i$Training),
        Block        = as.character(row_i$Block),
        stimulus_log = row_i$stimulus_log
      ),
      mode = "mean.class"
    )
    emm_df <- as_tibble(confint(emm))

    estimate_col <- c("emmean", "mean.class", "estimate")[c("emmean", "mean.class", "estimate") %in% names(emm_df)][1]
    lower_col    <- c("asymp.LCL", "lower.CL", "LCL")[c("asymp.LCL", "lower.CL", "LCL") %in% names(emm_df)][1]
    upper_col    <- c("asymp.UCL", "upper.CL", "UCL")[c("asymp.UCL", "upper.CL", "UCL") %in% names(emm_df)][1]

    emm_df %>%
      transmute(
        Genotype,
        Training = row_i$Training,
        Block    = row_i$Block,
        point    = row_i$point,
        stimulus = row_i$stimulus,
        expected_delay = .data[[estimate_col]] - 1,
        SE             = SE,
        lower          = .data[[lower_col]] - 1,
        upper          = .data[[upper_col]] - 1
      )
  })

  bind_rows(out)
}


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

per_cell_D <- get_expected_delay_at_custom_points(m_joint_delay, recovery_grid_small)

print(per_cell_D)

write.csv(
  per_cell_D,
  file.path(save_results_dir, "per_cell_endTrain_startTest_expected_delay.csv"),
  row.names = FALSE
)


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

write.csv(
  contrast_D,
  file.path(save_results_dir, "contrast_D_delay_change_spaced_vs_massed.csv"),
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
    Block    = factor(Block,    levels = levels(df_all$Block))
  )

training_grid <- get_expected_delay_at_custom_points(m_joint_delay, training_grid_small)

training_learning <- training_grid %>%
  select(Genotype, Training, point, expected_delay, SE) %>%
  pivot_wider(names_from = point, values_from = c(expected_delay, SE)) %>%
  mutate(
    learning_est  = expected_delay_end_b1 - expected_delay_start_b1,
    learning_se   = sqrt(SE_end_b1^2 + SE_start_b1^2),
    learning_low  = learning_est - 1.96 * learning_se,
    learning_high = learning_est + 1.96 * learning_se,
    z = learning_est / learning_se,
    p = 2 * pnorm(-abs(z))
  )

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
# 11.2 Latent-scale slope contrast
# ------------------------------------------------------------------------------
# For ordinal models, emtrends gives the slope on the latent cumulative-logit
# scale, not directly on the expected-delay scale.
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
