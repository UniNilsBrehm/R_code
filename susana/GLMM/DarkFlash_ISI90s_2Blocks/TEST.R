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

source("C:/UniFreiburg/Code/R_code/susana/GLMM/DarkFlash_ISI90s_2Blocks/utils.R")
base_dir <- "D:/WorkingData/Susana/GLMM/DarkFlash_ISI90s_2Blocks"

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
res <- load_data_darkflash_60s(file_dir, move_th = 0.2, , take_peak = 0)

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


# base_dir <- "D:/WorkingData/Susana/GLMM/DarkFlash_ISI90s_2Blocks/TEST"
save_fig_dir = file.path(base_dir, "figs")
save_results_dir = file.path(base_dir, "results")
models_dir = file.path(base_dir, "models")

# Create directories if they do not exist
dir.create(save_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Phase Split: Early (1-3) vs Late (4+) stimuli
# ==============================================================================

# Add phase variable to both dataframes
df_final <- df_final %>%
  mutate(
    phase = factor(
      ifelse(as.numeric(stimulus) <= 10, "early", "late"),
      levels = c("early", "late")
    )
  )

df_final_sub <- df_final_sub %>%
  mutate(
    phase = factor(
      ifelse(as.numeric(stimulus) <= 10, "early", "late"),
      levels = c("early", "late")
    )
  )

# Split dataframes
df_early     <- df_final     %>% filter(phase == "early")
df_late      <- df_final     %>% filter(phase == "late")
df_sub_early <- df_final_sub %>% filter(phase == "early")
df_sub_late  <- df_final_sub %>% filter(phase == "late")

# Re-center stimulus for late phase so it starts at 0
# (important for interpretable intercepts and NLS starting values)
df_late <- df_late %>%
  mutate(stimulus_log = log(as.numeric(stimulus) - 3))

df_sub_late <- df_sub_late %>%
  mutate(stimulus_log = log(as.numeric(stimulus) - 3))


# ==============================================================================
# Fit Models
# ==============================================================================
# Early phase: stimulus as categorical (only 3 levels)
df_early$stimulus_fac     <- factor(df_early$stimulus)
df_sub_early$stimulus_fac <- factor(df_sub_early$stimulus)

# Response probability - early
m_prob_early <- glmmTMB(
  move ~ Genotype * Block + Genotype * stimulus_fac + Block * stimulus_fac + 
    (1 | Video/Well),
  family = binomial(link = "logit"),
  data = df_early
)

# Peak distance - early
m_peak_early <- glmmTMB(
  max_peak ~ Genotype * Block * stimulus_fac + (1 | Video/Well),
  family = Gamma(link = "log"),
  data = df_sub_early
)

# Delay - early
m_delay_early <- clmm(
  delay_ord ~ Genotype * Block * stimulus_fac +
    (1 | Video) + (1 | Video:Well),
  data = df_sub_early,
  link = "logit"
)

# Response probability - late
m_prob_late <- glmmTMB(
  move ~ Genotype * Block * stimulus_log + (stimulus_log | Video/Well),
  family = binomial(link = "logit"),
  data = df_late
)

# Peak distance - late
m_peak_late <- glmmTMB(
  max_peak ~ Genotype * Block * stimulus_log + (stimulus_log | Video/Well),
  family = Gamma(link = "log"),
  data = df_sub_late
)

# Delay - late
m_delay_late_ordinal <- clmm(
  delay_ord ~ Genotype * Block * stimulus_log +
    (1 | Video) + (1 | Video:Well),
  data = df_sub_late,
  link = "logit"
)



# --- Early phase predictions ---
# Use the actual stimulus levels (1, 2, 3)
early_grid <- expand.grid(
  Genotype     = levels(df_early$Genotype),
  Block        = levels(df_early$Block),
  stimulus_fac = levels(df_early$stimulus_fac)
) %>%
  mutate(stimulus = as.numeric(as.character(stimulus_fac)))

early_grid$pred_logit <- predict(m_prob_early, 
                                 newdata = early_grid, 
                                 re.form = NA,  # population-level only
                                 type = "link")

early_grid <- early_grid %>%
  mutate(
    fit   = plogis(pred_logit),
    phase = "early"
  )

# --- Late phase predictions ---
# Stimulus sequence 4:60, recentered as log(stimulus - 3)
late_grid <- expand.grid(
  Genotype = levels(df_late$Genotype),
  Block    = levels(df_late$Block),
  stimulus = 4:60
) %>%
  mutate(stimulus_log = log(stimulus - 3))

late_grid$pred_logit <- predict(m_prob_late,
                                newdata = late_grid,
                                re.form = NA,
                                type = "link")

late_grid <- late_grid %>%
  mutate(
    fit   = plogis(pred_logit),
    phase = "late"
  )


# Bind both phases
pred_combined <- bind_rows(
  early_grid %>% select(Genotype, Block, stimulus, fit, phase),
  late_grid  %>% select(Genotype, Block, stimulus, fit, phase)
) %>%
  arrange(Genotype, Block, stimulus)

# SE from early model
early_ci <- early_grid %>%
  mutate(
    se  = predict(m_prob_early, newdata = early_grid, 
                  re.form = NA, type = "link", se.fit = TRUE)$se.fit,
    lwr = plogis(pred_logit - 1.96 * se),
    upr = plogis(pred_logit + 1.96 * se)
  ) %>%
  select(Genotype, Block, stimulus, fit, lwr, upr, phase)

# SE from late model
late_ci <- late_grid %>%
  mutate(
    se  = predict(m_prob_late, newdata = late_grid,
                  re.form = NA, type = "link", se.fit = TRUE)$se.fit,
    lwr = plogis(pred_logit - 1.96 * se),
    upr = plogis(pred_logit + 1.96 * se)
  ) %>%
  select(Genotype, Block, stimulus, fit, lwr, upr, phase)

pred_combined <- bind_rows(early_ci, late_ci) %>%
  arrange(Genotype, Block, stimulus)

df_prob_agg <- df_final %>%
  mutate(stimulus = as.numeric(stimulus)) %>%
  group_by(Genotype, Block, stimulus) %>%
  summarise(
    response_prob = mean(as.numeric(move), na.rm = TRUE),
    .groups = "drop"
  )

# Genotype colour palette — adjust to match your existing plots
geno_colors <- c(
  "ABTL"        = "#E41A1C",
  "th, th2, tyr"= "#8B8B00",
  "th, tyr"     = "#2E8B57",
  "th2, tyr"    = "#00BFFF",
  "tyr"         = "#FF69B4"
)

p_split <- ggplot() +
  # Raw aggregated points
  geom_point(
    data  = df_prob_agg,
    aes(x = stimulus, y = response_prob, color = Genotype),
    alpha = 0.4, size = 1
  ) +
  # Confidence ribbon
  geom_ribbon(
    data = pred_combined,
    aes(x = stimulus, ymin = lwr, ymax = upr, fill = Genotype),
    alpha = 0.2
  ) +
  # Fitted line — dashed for early, solid for late
  geom_line(
    data = pred_combined,
    aes(x = stimulus, y = fit, color = Genotype, 
        linetype = phase),
    linewidth = 1.2
  ) +
  scale_linetype_manual(
    values = c("early" = "dashed", "late" = "solid"),
    guide  = "none"   # hide from legend, or keep if you want it
  ) +
  scale_color_manual(values = geno_colors) +
  scale_fill_manual(values  = geno_colors) +
  facet_grid(Block ~ Genotype) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x     = "Stimulus number within block",
    y     = "Response probability",
    color = "Genotype",
    fill  = "Genotype",
    title = "Habituation curves — split early/late model"
  ) +
  theme_pubr(base_size = 14) +
  theme(legend.position = "top")

print(p_split)


