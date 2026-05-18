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
library(scales)

# Load helper functions (validation, plotting, EMM utilities, reporting)
source("C:/UniFreiburg/Code/R_code/tapping_habituation/utils_tapping.R")

# Base directory for saving results
base_dir_halo <- "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/Haloperidol_96wp/"
base_dir_buta <- "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/MinusButaclamol_96wp/"
base_dir_sch <- "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/SCH_96wp/"

# Load data
res_halo <- load_data(file.path(base_dir_halo, "haloperidol_96_well_plates_data.csv"))
res_buta <- load_data(file.path(base_dir_buta, "minus_butaclamol_96_well_plates_data.csv"))
res_sch <- load_data(file.path(base_dir_sch, "SCH29330_96_well_plates_data.csv"))

df_halo <- res_halo$df
df_buta <- res_buta$df
df_sch <- res_sch$df

df_halo$auc <- df_halo$auc + 0.00001
df_buta$auc <- df_buta$auc + 0.00001
df_sch$auc <- df_sch$auc + 0.00001

# Number of animals
n_haloperidol <- df_halo %>%
  distinct(fish_id, drug_condition, plate_id) %>%   # keep each fish only once per group
  count(drug_condition, name = "n_fish")
  # count(drug_condition, plate_id, name = "n_fish")
n_haloperidol

n_buta <- df_buta %>%
  distinct(fish_id, drug_condition, plate_id) %>%   # keep each fish only once per group
  count(drug_condition, name = "n_fish")
n_buta

n_sch <- df_sch %>%
  distinct(fish_id, drug_condition, plate_id) %>%   # keep each fish only once per group
  count(drug_condition, name = "n_fish")
n_sch

# ==============================================================================
# Response Probability
# ==============================================================================
# GLM Response Prob.
m_prob_halo <- glmmTMB(
  prob ~ Block * drug_condition * stimulus_log
  +(1 | plate_id/fish_id), # random intercept for plate an fish nested in plate
  #+(1 + stimulus_log | plate_id/fish_id), # random slope and intercept for plate an fish nested in plate
  data = df_halo,
  family = binomial
)

m_prob_buta <- glmmTMB(
  prob ~ Block * drug_condition * stimulus_log
  +(1 | plate_id/fish_id), # random intercept for plate an fish nested in plate
  # +(1 + stimulus_log | plate_id/fish_id), # random slope and intercept for plate an fish nested in plate
  data = df_buta,
  family = binomial
)


m_prob_sch <- glmmTMB(
  prob ~ Block * drug_condition * stimulus_log
  +(1 | plate_id/fish_id), # random intercept for plate an fish nested in plate
  #+(1 + stimulus_log | plate_id/fish_id), # random slope and intercept for plate an fish nested in plate
  data = df_sch,
  family = binomial
)

# Validation
validate_model(m_prob_halo, df_halo)
validate_model(m_prob_buta, df_buta)
validate_model(m_prob_sch, df_sch)

# Plots
g_prob_halo <- plot_habituation(df_halo, m_prob_halo, label = "Response prob.", response_var = "prob", Ymin = 0, Ymax = 1, transform = "plogis")
g_prob_halo
save_plot(g_prob_halo, file.path(base_dir_halo, "response_prob", "response_prob_habituation_curves"))

g_prob_buta <- plot_habituation(df_buta, m_prob_buta, label = "Response prob.", response_var = "prob", Ymin = 0, Ymax = 1, transform = "plogis")
g_prob_buta
save_plot(g_prob_buta, file.path(base_dir_buta, "response_prob", "response_prob_habituation_curves"))

g_prob_sch <- plot_habituation(df_sch, m_prob_sch, label = "Response prob.", response_var = "prob", Ymin = 0, Ymax = 1, transform = "plogis")
g_prob_sch
save_plot(g_prob_sch, file.path(base_dir_sch, "response_prob", "response_prob_habituation_curves"))


# ==============================================================================
# Estimated Marginal Means (EMMs)
# ==============================================================================
# Haloperidol
compute_emm_response_prob(m_prob_halo, df_halo, base_dir_halo)

# (-)Butaclalmol
compute_emm_response_prob(m_prob_buta, df_buta, base_dir_buta)

# SCH-29330
compute_emm_response_prob(m_prob_sch, df_sch, base_dir_sch)


# ==============================================================================
# Summed Distance
# ==============================================================================
# GLMs

m_distance_halo <- glmmTMB(
  auc ~ Block * drug_condition * stimulus_log 
  +(1 | plate_id/fish_id),
  data = df_halo,
  family = Gamma(link = "log")
)

m_distance_buta <- glmmTMB(
  auc ~ Block * drug_condition * stimulus_log 
  +(1 | plate_id/fish_id),
  data = df_buta,
  family = Gamma(link = "log")
)

m_distance_sch <- glmmTMB(
  auc ~ Block * drug_condition * stimulus_log 
  +(1 | plate_id/fish_id),
  data = df_sch,
  family = Gamma(link = "log")
)

# Validation
# Simulate residuals
res <- simulateResiduals(m_distance_halo)
plot(res)

res <- simulateResiduals(m_distance_buta)
plot(res)

res <- simulateResiduals(m_distance_sch)
plot(res)

# Plots
g_distance_halo <- plot_habituation(df_halo, m_distance_halo, label = "Distance moved (mm)", response_var = "auc", Ymin = 10, Ymax = 60, transform = "exp")
g_distance_halo
save_plot(g_distance_halo, file.path(base_dir_halo, "distance", "distance_habituation_curves"))

g_distance_buta <- plot_habituation(df_buta, m_distance_buta, label = "Distance moved (mm)", response_var = "auc", Ymin = 10, Ymax = 80, transform = "exp")
g_distance_buta
save_plot(g_distance_buta, file.path(base_dir_buta, "distance", "distance_habituation_curves"))

g_distance_sch <- plot_habituation(df_sch, m_distance_sch, label = "Distance moved (mm)", response_var = "auc", Ymin = 10, Ymax = 60, transform = "exp")
g_distance_sch
save_plot(g_distance_sch, file.path(base_dir_sch, "distance", "distance_habituation_curves"))


# ==============================================================================
# Estimated Marginal Means (EMMs)
# ==============================================================================
# Haloperidol
compute_emm_distance(m_distance_halo, df_halo, base_dir_halo)

# (-)Butaclalmol
compute_emm_distance(m_distance_buta, df_buta, base_dir_buta)

# SCH-29330
compute_emm_distance(m_distance_sch, df_sch, base_dir_sch)

# ==============================================================================
# ==============================================================================
# Light/Dark
library(lme4)
library(lmerTest)

# Haloperidol
# ==============================================================================

halo_ld <- read_csv(file.path(base_dir_halo, "haloperidol_96_well_light_dark.csv"))
# halo_ld <- read_csv(file.path(base_dir_buta, "MinusButaclamol_96_well_light_dark.csv"))
# halo_ld <- read_csv(file.path(base_dir_sch, "SCH_96_well_light_dark.csv"))

halo_ld_df <- halo_ld %>%
  select(-`...1`) %>%
  mutate(
    treatment = factor(treatment, levels = c("Control","Acute","DK")),
    phase     = factor(phase, levels = c("Light","Dark")),
    block     = factor(block, levels = c("L1","D1","L2","D2")),
    plate     = factor(plate),
    larva_id  = factor(larva_id),
    # good idea for time modeling:
    time_c    = as.numeric(scale(time_sec, center = TRUE, scale = TRUE))
  ) %>%
  filter(!is.na(block))   # drop anything outside the 40 min schedule

#PLOTS
sum_df <- halo_ld_df %>%
  group_by(treatment, time_bin, time_sec, phase, block) %>%
  summarise(
    n = n_distinct(larva_id),
    mean = mean(distance, na.rm = TRUE),
    sd   = sd(distance, na.rm = TRUE),
    se   = sd / sqrt(n),
    ci95 = 1.96 * se,
    .groups = "drop"
  )

# Shade dark periods (10–20 min and 30–40 min)
dark_rects <- data.frame(
  xmin = c(10*60, 30*60),
  xmax = c(20*60, 40*60)
)

# Figure 1
halo_ld_plot01 <- ggplot(sum_df, aes(x = time_sec, y = mean, color = treatment)) +
  geom_rect(
    data = dark_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.10
  ) +
  geom_ribbon(
    aes(ymin = mean - ci95, ymax = mean + ci95, fill = treatment),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(
    name = "Time (min)",
    breaks = seq(0, 40, by = 5) * 60,
    labels = seq(0, 40, by = 5)
  ) +
  ylab("Distance moved per bin") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
halo_ld_plot01
save_plot(halo_ld_plot01, file.path(base_dir_halo, "lightdark", "light_dark_curves_01"))


# Figure 2
halo_ld_plot02 <- ggplot(sum_df, aes(x = time_sec, y = mean)) +
  geom_rect(
    data = dark_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.10
  ) +
  geom_ribbon(aes(ymin = mean - ci95, ymax = mean + ci95), alpha = 0.18) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ treatment, ncol = 1) +
  scale_x_continuous(
    name = "Time (min)",
    breaks = seq(0, 40, by = 5) * 60,
    labels = seq(0, 40, by = 5)
  ) +
  ylab("Distance moved per bin") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
halo_ld_plot02
save_plot(halo_ld_plot02, file.path(base_dir_halo, "lightdark", "light_dark_curves_02"))


# Minus Butaclamol
# ==============================================================================
buta_ld <- read_csv(file.path(base_dir_buta, "MinusButaclamol_96_well_light_dark.csv"))

buta_ld_df <- buta_ld %>%
  select(-`...1`) %>%
  mutate(
    treatment = factor(treatment, levels = c("Control","Acute","DK")),
    phase     = factor(phase, levels = c("Light","Dark")),
    block     = factor(block, levels = c("L1","D1","L2","D2")),
    plate     = factor(plate),
    larva_id  = factor(larva_id),
    # good idea for time modeling:
    time_c    = as.numeric(scale(time_sec, center = TRUE, scale = TRUE))
  ) %>%
  filter(!is.na(block))   # drop anything outside the 40 min schedule

#PLOTS
sum_df <- buta_ld_df %>%
  group_by(treatment, time_bin, time_sec, phase, block) %>%
  summarise(
    n = n_distinct(larva_id),
    mean = mean(distance, na.rm = TRUE),
    sd   = sd(distance, na.rm = TRUE),
    se   = sd / sqrt(n),
    ci95 = 1.96 * se,
    .groups = "drop"
  )

# Shade dark periods (10–20 min and 30–40 min)
dark_rects <- data.frame(
  xmin = c(10*60, 30*60),
  xmax = c(20*60, 40*60)
)

# Figure 1
buta_ld_plot01 <- ggplot(sum_df, aes(x = time_sec, y = mean, color = treatment)) +
  geom_rect(
    data = dark_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.10
  ) +
  geom_ribbon(
    aes(ymin = mean - ci95, ymax = mean + ci95, fill = treatment),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(
    name = "Time (min)",
    breaks = seq(0, 40, by = 5) * 60,
    labels = seq(0, 40, by = 5)
  ) +
  ylab("Distance moved per bin") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
buta_ld_plot01
save_plot(buta_ld_plot01, file.path(base_dir_buta, "lightdark", "light_dark_curves_01"))


# Figure 2
buta_ld_plot02 <- ggplot(sum_df, aes(x = time_sec, y = mean)) +
  geom_rect(
    data = dark_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.10
  ) +
  geom_ribbon(aes(ymin = mean - ci95, ymax = mean + ci95), alpha = 0.18) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ treatment, ncol = 1) +
  scale_x_continuous(
    name = "Time (min)",
    breaks = seq(0, 40, by = 5) * 60,
    labels = seq(0, 40, by = 5)
  ) +
  ylab("Distance moved per bin") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
buta_ld_plot02
save_plot(buta_ld_plot02, file.path(base_dir_buta, "lightdark", "light_dark_curves_02"))

# SCH-29330
# ==============================================================================
sch_ld <- read_csv(file.path(base_dir_sch, "SCH_96_well_light_dark.csv"))

sch_ld_df <- sch_ld %>%
  select(-`...1`) %>%
  mutate(
    treatment = factor(treatment, levels = c("Control","Acute","DK")),
    phase     = factor(phase, levels = c("Light","Dark")),
    block     = factor(block, levels = c("L1","D1","L2","D2")),
    plate     = factor(plate),
    larva_id  = factor(larva_id),
    # good idea for time modeling:
    time_c    = as.numeric(scale(time_sec, center = TRUE, scale = TRUE))
  ) %>%
  filter(!is.na(block))   # drop anything outside the 40 min schedule

#PLOTS
sum_df <- sch_ld_df %>%
  group_by(treatment, time_bin, time_sec, phase, block) %>%
  summarise(
    n = n_distinct(larva_id),
    mean = mean(distance, na.rm = TRUE),
    sd   = sd(distance, na.rm = TRUE),
    se   = sd / sqrt(n),
    ci95 = 1.96 * se,
    .groups = "drop"
  )

# Shade dark periods (10–20 min and 30–40 min)
dark_rects <- data.frame(
  xmin = c(10*60, 30*60),
  xmax = c(20*60, 40*60)
)

# Figure 1
sch_ld_plot01 <- ggplot(sum_df, aes(x = time_sec, y = mean, color = treatment)) +
  geom_rect(
    data = dark_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.10
  ) +
  geom_ribbon(
    aes(ymin = mean - ci95, ymax = mean + ci95, fill = treatment),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(
    name = "Time (min)",
    breaks = seq(0, 40, by = 5) * 60,
    labels = seq(0, 40, by = 5)
  ) +
  ylab("Distance moved per bin") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
sch_ld_plot01
save_plot(sch_ld_plot01, file.path(base_dir_sch, "lightdark", "light_dark_curves_01"))


# Figure 2
sch_ld_plot02 <- ggplot(sum_df, aes(x = time_sec, y = mean)) +
  geom_rect(
    data = dark_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.10
  ) +
  geom_ribbon(aes(ymin = mean - ci95, ymax = mean + ci95), alpha = 0.18) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ treatment, ncol = 1) +
  scale_x_continuous(
    name = "Time (min)",
    breaks = seq(0, 40, by = 5) * 60,
    labels = seq(0, 40, by = 5)
  ) +
  ylab("Distance moved per bin") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
sch_ld_plot02
save_plot(sch_ld_plot02, file.path(base_dir_sch, "lightdark", "light_dark_curves_02"))



# ==============================================================================
# ==============================================================================
# Linear Mixed Model
m_halo_ld <- lmer(
  distance ~ treatment * phase +
    (1 | plate) +
    (1 + phase | larva_id),
  REML = FALSE,
  data = halo_ld_df
)

# 1) Light vs Dark within each treatment: Does each treatment show an LD response?
# Tests average activity
emm_phase_within_trt <- emmeans(m_halo_ld, ~ phase | treatment)
pairs(emm_phase_within_trt, adjust = "holm")

# 2) Treatment differences within Light and within Dark: 
# Do treatments differ during Light? Do treatments differ during Dark?
# Tests average activity.
emm_trt_within_phase <- emmeans(m_halo_ld, ~ treatment | phase)
pairs(emm_trt_within_phase, adjust = "holm")

# 3) Tests early vs late mean level.
# 3.1) Same tests but for initial response
initial_df <- halo_ld_df %>%
  group_by(larva_id, block) %>%
  mutate(time_within = time_sec - min(time_sec)) %>%
  ungroup() %>%
  filter(time_within < 5*60)   # first half only

m_initial <- lmer(
  distance ~ treatment * phase +
    (1 + phase | larva_id),
  data = initial_df
)

emmeans(m_initial, ~ treatment | phase) |> pairs(adjust = "holm")

# 3.2) Same tests but for steady state only
steady_df <- halo_ld_df %>%
  group_by(larva_id, block) %>%
  mutate(time_within = time_sec - min(time_sec)) %>%
  ungroup() %>%
  filter(time_within >= 5*60)   # second half only

m_steady <- lmer(
  distance ~ treatment * phase +
    (1 + phase | larva_id),
  data = steady_df
)

# emmeans(m_steady, ~ phase | treatment) |> pairs(adjust = "holm")
emmeans(m_steady, ~ treatment | phase) |> pairs(adjust = "holm")

# 4) Test slopes within each phase for each treatment
# Tests minute-scale adaptation dynamics.
# treatment × phase × time: Does the within-phase slope differ between treatments AND between light vs dark?
# Each larva gets: its own baseline activity (intercept) and its own Light–Dark sensitivity and its own habituation slope over time
m_slope <- lmer(
  distance ~ treatment * phase * time_c +
    (1 + phase + time_c | larva_id),
  data = halo_ld_df,
  control = lmerControl(optimizer = "bobyqa")
)

# emtrends: So instead of comparing means, you compare rates of change.
# slope = change in distance per 1 unit of time_c
# Interpretation: 
# - negative slope -> activity decreases (habituation)
# - positive slope -> activity increases
# - zero slope -> stable steady state

emtrends(m_slope, ~ treatment * phase, var = "time_c")

# Are the slopes different between groups?
pairs(emtrends(m_slope, ~ treatment * phase, var = "time_c"),
      adjust = "holm")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Peri-transition model (on and off responses)
# Does treatment change sensory-evoked onset magnitude and/or post-transition 
# dynamics, and does this depend on flash direction?
# Tests seconds-scale sensory bursts.
W <- 60  # seconds before/after transition

# FALSE = before transition
# TRUE = after transition
peri_df <- halo_ld_df %>%
  mutate(
    within_block = time_sec %% 600
  ) %>%
  filter(within_block <= W | within_block >= (600 - W)) %>%
  mutate(
    # relative time to transition
    rel_time = ifelse(within_block <= W, within_block, within_block - 600),
    
    # indicator: before vs after transition
    after = rel_time >= 0
  )

peri_df <- peri_df %>%
  mutate(
    transition = case_when(
      block %in% c("L1", "L2") & rel_time < 0  ~ "LightToDark",
      block %in% c("D1", "D2") & rel_time >= 0 ~ "LightToDark",
      block %in% c("D1") & rel_time < 0        ~ "DarkToLight",
      block %in% c("L2") & rel_time >= 0       ~ "DarkToLight",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(transition)) %>%
  mutate(
    transition = factor(transition,
                        levels = c("LightToDark", "DarkToLight"))
  )

peri_df$rel_time_z <- scale(peri_df$rel_time)

m_peri <- lmer(
  distance ~ treatment * transition * after +
    treatment * transition * rel_time_z +
    (1 + after | larva_id),
  data = peri_df,
  control = lmerControl(optimizer = "bobyqa")
)

anova(m_peri)
summary(m_peri)

# Jump size per treatment & transition
# How big is the difference between before vs after?
emm <- emmeans(m_peri, ~ treatment * transition * after,
               at = list(rel_time_z = 0))

# Compute jump (after - before), then compare across treatments
jump_contrasts <- contrast(emm,
                           method = "revpairwise",
                           by = c("treatment", "transition"))

# Now compare those jumps between treatments
pairs(jump_contrasts, by = "transition", adjust = "holm")

# Recovery slopes
emtrends(m_peri, ~ treatment * transition, var = "rel_time_z") |>
  pairs(adjust = "holm")

# Visualization (PLOT)
peri_sum <- peri_df %>%
  group_by(transition, treatment, rel_time) %>%
  summarise(
    mean = mean(distance),
    se = sd(distance)/sqrt(n_distinct(larva_id)),
    .groups = "drop"
  )

ggplot(peri_sum, aes(rel_time, mean, color = treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - 1.96*se,
                  ymax = mean + 1.96*se,
                  fill = treatment),
              alpha = 0.15, color = NA) +
  facet_wrap(~ transition, scales = "free_y") +
  theme_bw() +
  labs(x = "Time relative to transition (s)",
       y = "Distance moved")

# ==============================================================================
# FINISHED =====================================================================
# ******************************************************************************
# ******************************************************************************

# ==============================================================================
# Compute response prob from raw data (ignoring mixed structure)
# Just for checking the exp decay (the log(stimulus) in the GLMM)
# ==============================================================================
# Response Prob.
df_summary <- df_sch %>%
  group_by(Block, drug_condition, stimulus) %>%  # include drug_condition
  summarise(
    response_prob = mean(prob, na.rm = TRUE),
    n = n(),
    se = sqrt(response_prob * (1 - response_prob) / n),
    .groups = "drop"
  )


ggplot(df_summary, aes(x = stimulus, y = response_prob,
                       color = drug_condition, group = drug_condition)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = response_prob - se,
                    ymax = response_prob + se),
                width = 0.15, alpha = 0.4) +
  facet_wrap(~Block, scales = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(
    title = "Response Probability per Stimulus",
    subtitle = "Split by Block and Drug Condition",
    x = "Stimulus number (within block)",
    y = "Response probability",
    color = "Drug condition"
  ) +
  theme_pubr(base_size = 14)

################################################################################
# Distance moved
df_summary <- df %>%
  group_by(Block, drug_condition, stimulus) %>%  # include drug_condition
  summarise(
    v_mean = mean(auc, na.rm = TRUE),
    n = n(),
    se = sd(auc, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(df_summary, aes(x = stimulus, y = v_mean,
                       color = drug_condition, group = drug_condition)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = v_mean - se,
                    ymax = v_mean + se),
                width = 0.15, alpha = 0.4) +
  facet_wrap(~Block, scales = "free_x") +
  labs(
    title = "Habituation Curves",
    subtitle = "Split by Block and Drug Condition",
    x = "Stimulus number (within block)",
    y = "Metric",
    color = "Drug condition"
  ) +
  theme_pubr(base_size = 14)
