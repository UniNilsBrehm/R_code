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
library(data.table)

# Load helper functions (validation, plotting, EMM utilities, reporting)
source("C:/UniFreiburg/Code/R_code/tapping_habituation/utils_tapping.R")

# Base directory for saving results
base_dir <- "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/48wp"

# Load data

# Load data from csv file
res <- load_48_well_plates(file.path(base_dir, "all_treatments_48wellplates.csv"))
df <- res$df
df$auc <- df$auc + 0.0001
df_sub <- res$df_sub

df_5s <- subset(df, Block %in% c("ISI_5s_block"))
df_90s <- subset(df, Block %in% c("ISI_90s"))

fish_counts <- df %>%
  distinct(fish_id, treatment, drug_condition) %>%   # keep each fish only once per group
  count(treatment, drug_condition, name = "n_fish")
print(fish_counts)

# ==============================================================================
# Response Probability
# ==============================================================================
# GLM Response Prob.
m_prob <- glmmTMB(
  response ~ drug_condition * stimulus_log * treatment
  +(1 | plate_uid/fish_id), # random intercept for plate an fish nested in plate
  #+ (1 + stimulus_log | plate_uid/fish_id), # random slope and intercept for plate an fish nested in plate
  data = df_5s,
  family = binomial
)

# Validation
validate_model(m_prob, df_5s)
# validate_model(m_prob_all, df)
# compare_performance(m_prob_test, m_prob)
# anova(m_prob_test, m_prob)

# ==============================================================================
# Prepare for Plotting
# ==============================================================================
# Build a data grid for all combinations that exist in your dataset
valid_combos <- df_5s %>%
  distinct(drug_condition, treatment)

new_data <- expand.grid(
  stimulus = 1:20
) %>%
  mutate(stimulus_log = log(stimulus)) %>%
  tidyr::crossing(valid_combos)   # only valid drug_condition × treatment pairs

# Predictions from Model
pred <- predict(
  m_prob,
  newdata = new_data,
  re.form = NA,   # fixed effects only (population-level)
  se.fit = TRUE
)

# Back-transform from the logit scale
new_data$fit <- plogis(pred$fit)
new_data$CI_low <- plogis(pred$fit - 1.96 * pred$se.fit)
new_data$CI_high <- plogis(pred$fit + 1.96 * pred$se.fit)

# Compute observed mean responses
p_dt_5s <- data.table(df_5s)
p_dt_5s <- p_dt_5s[, .(
  response = mean(response == 1)
), by = .(treatment, stimulus, drug_condition)]

# Plot one treatment per panel
ggplot(new_data, aes(x = stimulus, color = drug_condition)) +
  geom_line(aes(y = fit), size = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = drug_condition),
              alpha = 0.2, color = NA) +
  geom_point(data = p_dt_5s,
             aes(x = stimulus, y = response, color = drug_condition),
             alpha = 0.6, size = 1.5) +
  facet_wrap(~ treatment, scales = "fixed") +
  # Set consistent tick marks for both axes
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20),
    limits = c(0, 20)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1)
  ) +
  labs(
    x = "Tap number (within block)",
    y = "Response probability",
    color = "Drug condition",
    fill  = "Drug condition",
    # title = "Habituation curves by treatment",
    # subtitle = "Predicted from GLMM (binomial logit link)"
  ) +
  theme_pubr(base_size = 14)

ggsave(
  filename = file.path(base_dir, "Habituation_GLMM_response_prob.jpg"),
  plot = last_plot(),       # saves the most recent ggplot
  width = 10,
  height = 7,
  units = "in",
  dpi=600
)

# ==============================================================================
# Estimated Marginal Means (EMMs)
# ==============================================================================
# Habituation slope
slopes <- emtrends(
  m_prob, 
  var = "stimulus_log",
  ~ drug_condition * treatment
)
pairs(slopes, by = "treatment")

# First 3 stimuli
emm_first_3 <- emmeans(
  m_prob,
  ~ drug_condition * treatment,
  at = list(stimulus_log = log(1:3)),
  cov.reduce = mean
)
get_effects_logit(emm_first_3, p_th=0.1)


# Middle stimuli
emm_center <- emmeans(
  m_prob,
  ~ drug_condition * treatment,
  at = list(stimulus_log = log(8:12)),
  cov.reduce = mean
)
get_effects_logit(emm_center, p_th=0.1)

# Last 3 stimuli
emm_last_3 <- emmeans(
  m_prob,
  ~ drug_condition * treatment,
  at = list(stimulus_log = log(17:20)),
  cov.reduce = mean
)
get_effects_logit(emm_last_3, p_th=0.1)

# Over all stimuli (average, y-shift)
emm_avg <- emmeans(m_prob, ~ drug_condition | treatment)
get_effects_logit(emm_avg, p_th=0.1)

# Some Notes
# odds = p / (1-p)
# odds ratio = odds A / odds B
# Same mathematical meaning:
# "On average, acute Apomorphine increases the relative likelihood of responding versus not responding by about 62%."
# "On average, acute Apomorphine increases the response probability by about 9 percentage points (from 22% to 31%)."

# ==============================================================================
# Summed distance moved
# ==============================================================================
m_distance <- glmmTMB(
  auc ~ drug_condition * stimulus_log * treatment
  +(1 | plate_uid/fish_id), # random intercept for plate an fish nested in plate
  # +(1 + stimulus_log | plate_uid/fish_id), # random slope and intercept for plate an fish nested in plate
  data = df_5s,
  family = Gamma(link = "log")
  )

# Validation
validate_model(m_distance, df_5s)

hist(df_5s$auc, breaks=100)

res <- residuals(m_distance)
ggplot(data.frame(res), aes(x = res)) +
  geom_histogram(bins = 30, color = "black", fill = "grey") +
  labs(title = "Residuals Histogram", x = "Residuals", y = "Count")

# ==============================================================================
# Estimated Marginal Means (EMMs)
# ==============================================================================
# Habituation slope
slopes <- emtrends(
  m_distance, 
  var = "stimulus_log",
  ~ drug_condition * treatment
)
get_slope_effects_log(slopes, p_th=0.1)

########################################
# First 3 stimuli
emm_first_3 <- emmeans(
  m_distance,
  ~ drug_condition * treatment,
  at = list(stimulus_log = log(1:3)),
  cov.reduce = mean
)
get_effects_log(emm_first_3, p_th=0.1)

# Middle stimuli
emm_center <- emmeans(
  m_distance,
  ~ drug_condition * treatment,
  at = list(stimulus_log = log(8:12)),
  cov.reduce = mean
)
get_effects_log(emm_center, p_th=0.1)

# Last 3 stimuli
emm_last_3 <- emmeans(
  m_distance,
  ~ drug_condition * treatment,
  at = list(stimulus_log = log(17:20)),
  cov.reduce = mean
)
get_effects_log(emm_last_3, p_th=0.1)

# Over all stimuli (average, y-shift)
emm_avg <- emmeans(m_distance, ~ drug_condition | treatment)
get_effects_log(emm_avg, p_th=0.1)

# ==============================================================================
# Prepare for Plotting
# ==============================================================================
# Build a data grid for all combinations that exist in your dataset
valid_combos <- df_5s %>%
  distinct(drug_condition, treatment)

new_data <- expand.grid(
  stimulus = 1:20
) %>%
  mutate(stimulus_log = log(stimulus)) %>%
  tidyr::crossing(valid_combos)   # only valid drug_condition × treatment pairs

# Predictions from Model
pred <- predict(
  m_distance,
  newdata = new_data,
  re.form = NA,   # fixed effects only (population-level)
  se.fit = TRUE
)

# Back-transform from the logit scale
new_data$fit <- exp(pred$fit)
new_data$CI_low <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$CI_high <- exp(pred$fit + 1.96 * pred$se.fit)

# Compute observed mean responses
p_dt_5s <- data.table(df_5s)
p_dt_5s <- p_dt_5s[, .(
  response = mean(auc)
), by = .(treatment, stimulus, drug_condition)]

# Plot one treatment per panel
ggplot(new_data, aes(x = stimulus, color = drug_condition)) +
  geom_line(aes(y = fit), size = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = drug_condition),
              alpha = 0.2, color = NA) +
  geom_point(data = p_dt_5s,
             aes(x = stimulus, y = response, color = drug_condition),
             alpha = 0.6, size = 1.5) +
  facet_wrap(~ treatment, scales = "fixed") +
  # Set consistent tick marks for both axes
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20),
    limits = c(0, 20)
  ) +
  labs(
    x = "Tap number (within block)",
    y = "Summed distance moved",
    color = "Drug condition",
    fill  = "Drug condition",
    #title = "Habituation curves by treatment",
    #subtitle = "Predicted from GLMM (Gamma, log)"
  ) +
  theme_pubr(base_size = 14)

ggsave(
  filename = file.path(base_dir, "Habituation_GLMM_distance_moved.jpg"),
  plot = last_plot(),       # saves the most recent ggplot
  width = 10,
  height = 7,
  units = "in",
  dpi=600
)

##
##
##
# ==============================================================================
# ISI 90s
# ==============================================================================
# Response Prob.
m_prob_90 <- glmmTMB(
  response ~ drug_condition * treatment + stimulus_log
  +(1 | plate_id/fish_id), # random intercept for plate an fish nested in plate
  data = df_90s,
  family = binomial
)

# Validation
validate_model(m_prob_90, df_90s)

emm_90 <- emmeans(m_prob_90, ~ drug_condition| treatment, type='response')
pairs(emm_90)

# Build a data grid for all combinations that exist in your dataset
valid_combos <- df_90s %>%
  distinct(drug_condition, treatment)

new_data <- expand.grid(
  stimulus = 1:3
) %>%
  mutate(stimulus_log = log(stimulus)) %>%
  tidyr::crossing(valid_combos)   # only valid drug_condition × treatment pairs

# Predictions from Model
pred <- predict(
  m_prob_90,
  newdata = new_data,
  re.form = NA,   # fixed effects only (population-level)
  se.fit = TRUE
)

# Back-transform from the logit scale
new_data$fit <- plogis(pred$fit)
new_data$CI_low <- plogis(pred$fit - 1.96 * pred$se.fit)
new_data$CI_high <- plogis(pred$fit + 1.96 * pred$se.fit)

# Compute observed mean responses
p_dt_5s <- data.table(df_90s)
p_dt_5s <- p_dt_5s[, .(
  response = mean(response == 1)
), by = .(treatment, stimulus, drug_condition)]

# Plot one treatment per panel
ggplot(new_data, aes(x = stimulus, color = drug_condition)) +
  geom_line(aes(y = fit), size = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = drug_condition),
              alpha = 0.2, color = NA) +
  geom_point(data = p_dt_5s,
             aes(x = stimulus, y = response, color = drug_condition),
             alpha = 0.6, size = 1.5) +
  facet_wrap(~ treatment, scales = "fixed") +
  # Set consistent tick marks for both axes
  scale_x_continuous(
    breaks = c(1, 2, 3),
    limits = c(1, 3)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0.5, 1)
  ) +
  labs(
    x = "Tap number (within block)",
    y = "Response probability",
    color = "Drug condition",
    fill  = "Drug condition",
    # title = "Habituation curves by treatment",
    # subtitle = "Predicted from GLMM (binomial logit link)"
  ) +
  theme_pubr(base_size = 14)

ggsave(
  filename = file.path(base_dir, "Habituation_GLMM_reponse_prob_ISI90s.jpg"),
  plot = last_plot(),       # saves the most recent ggplot
  width = 10,
  height = 7,
  units = "in",
  dpi=600
)

# ==============================================================================
# Distance moved
#
# ==============================================================================
df_90s$auc <- df_90s$auc + 0.0001
m_distance_90 <- glmmTMB(
  auc ~ drug_condition * stimulus_log * treatment
  +(1 | plate_id/fish_id), # random intercept for plate an fish nested in plate
  data = df_90s,
  family = Gamma(link = "log")
)

emm_90 <- emmeans(m_prob_90, ~ drug_condition| treatment, type='response')
pairs(emm_90)

# Build a data grid for all combinations that exist in your dataset
valid_combos <- df_90s %>%
  distinct(drug_condition, treatment)

new_data <- expand.grid(
  stimulus = 1:3
) %>%
  mutate(stimulus_log = log(stimulus)) %>%
  tidyr::crossing(valid_combos)   # only valid drug_condition × treatment pairs

# Predictions from Model
pred <- predict(
  m_distance_90,
  newdata = new_data,
  re.form = NA,   # fixed effects only (population-level)
  se.fit = TRUE
)

# Back-transform from the log scale
new_data$fit <- exp(pred$fit)
new_data$CI_low <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$CI_high <- exp(pred$fit + 1.96 * pred$se.fit)

# Compute observed mean responses
p_dt_5s <- data.table(df_90s)
p_dt_5s <- p_dt_5s[, .(
  response = mean(auc)
), by = .(treatment, stimulus, drug_condition)]

# Plot one treatment per panel
ggplot(new_data, aes(x = stimulus, color = drug_condition)) +
  geom_line(aes(y = fit), size = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = drug_condition),
              alpha = 0.2, color = NA) +
  geom_point(data = p_dt_5s,
             aes(x = stimulus, y = response, color = drug_condition),
             alpha = 0.6, size = 1.5) +
  facet_wrap(~ treatment, scales = "fixed") +
  # Set consistent tick marks for both axes
  scale_x_continuous(
    breaks = c(1, 2, 3),
    limits = c(1, 3)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 110)
  ) +
  labs(
    x = "Tap number (within block)",
    y = "Distance moved (mm)",
    color = "Drug condition",
    fill  = "Drug condition",
    # title = "Habituation curves by treatment",
    # subtitle = "Predicted from GLMM (binomial logit link)"
  ) +
  theme_pubr(base_size = 14)

ggsave(
  filename = file.path(base_dir, "Habituation_GLMM_distance_ISI90s.jpg"),
  plot = last_plot(),       # saves the most recent ggplot
  width = 10,
  height = 7,
  units = "in",
  dpi=600
)

##
##
##
# ==============================================================================
# Light/Dark
# ==============================================================================
df_ld <- read_csv(file.path(base_dir, "light_dark_all_treatments_48wellplates.csv"))

ld_df <- df_ld %>%
  select(-`...1`) %>%
  mutate(
    treatment     = factor(treatment, levels = c("Control","Acute","DK")),
    phase         = factor(phase, levels = c("Light","Dark")),
    block         = factor(block, levels = c("L1","D1","L2","D2")),
    plate         = factor(plate),
    larva_id      = factor(larva_id),
    drug          = factor(drug),
    concentration = factor(concentration),
    time_c        = as.numeric(scale(time_sec, center = TRUE, scale = TRUE))
  ) %>%
  filter(!is.na(block))

# summary WITH drug + concentration
sum_df <- ld_df %>%
  group_by(drug, concentration, treatment, time_bin, time_sec, phase, block) %>%
  summarise(
    n    = n_distinct(larva_id),
    mean = mean(distance, na.rm = TRUE),
    sd   = sd(distance, na.rm = TRUE),
    se   = sd / sqrt(n),
    ci95 = 1.96 * se,
    .groups = "drop"
  ) %>%
  # one-line facet strip label like: "apomorphine (20uM)"
  mutate(drug_conc = paste0(as.character(drug), " (", as.character(concentration), ")"))

# Shade dark periods
dark_rects <- data.frame(
  xmin = c(20*60, 60*60),
  xmax = c(40*60, 80*60)
)

plot_ld <- ggplot(sum_df, aes(x = time_sec, y = mean, color = treatment)) +
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
    name   = "Time (min)",
    breaks = seq(0, 80, by = 10) * 60,
    labels = seq(0, 80, by = 10)
  ) +
  ylab("Distance moved per 30 seconds") +
  facet_wrap(~ drug_conc, scales = "fixed") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA)
  )



plot_ld
save_plot(plot_ld, file.path(base_dir, "light_dark_curves"))
