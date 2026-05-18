library(readr)
library(lme4)
library(performance)
library(sjPlot)
library(emmeans)
library(lmerTest)
library(dplyr)
library(glmmTMB)
library(ggplot2)
library(scales)
library(DHARMa)
library(data.table)
library(broom)

# Load data from csv file
# df <- read_csv("D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/apomorphine_20uM_48wellplates.csv")
# df <- read_csv("D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/SCH_48wellplates.csv")
df <- read_csv("D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/all_treatments_48wellplates.csv")

# Make sure response is numeric (0/1)
df$response <- as.numeric(df$response)       # 0/1
df$fish_id <- as.factor(df$fish_id)
df$plate_id <- as.factor(df$plate_id)
df$drug_condition <- as.factor(df$drug_condition)
df$ISI_block <- as.factor(df$ISI_block)

# Define tap number within block
df <- df %>%
  group_by(ISI_block, fish_id, plate_id) %>%
  mutate(tap_in_block = tap_num - min(tap_num) + 1) %>%
  ungroup()

# ======================================================================================================================================
# Only the ISI 5s Blocks
df_5s <- subset(df, ISI_block %in% c("ISI_5s_block"))

# Set References
df_5s$drug_condition <- factor(df_5s$drug_condition, levels = c("control", "acute", "dk"))
df_5s$ISI_block <- factor(df_5s$ISI_block, levels = c("ISI_5s_block"))

# Log transform tap number (for exponential decay)
df_5s$tap_in_block_log <- log(df_5s$tap_in_block)

# Create unique FISH ID
df_5s <- df_5s %>%
  dplyr::group_by(treatment, plate_id, fish_id) %>%
  dplyr::mutate(fish_uid_num = dplyr::cur_group_id()) %>%
  dplyr::ungroup()

check <- df_5s %>%
  dplyr::count(treatment, plate_id, fish_id, fish_uid_num)
print(check, n=1000)

# Center and Scale taps in block
# df_5s$tap_in_block <- scale(df_5s$tap_num, center = TRUE, scale = TRUE)

# GLMM:
# fish_id is nested within plate_id and treatment.
# (1 | treatment/plate_id/fish_id) expands to:
# (1 | treatment) + (1 | treatment:plate_id) + (1 | treatment:plate_id:fish_id)

# Treat treatment as a fixed factor (if you want to compare drugs)
glmm_block <- glmmTMB(
  response ~ drug_condition * tap_in_block_log * treatment
  # +(1 + tap_in_block_log | fish_id) + (1 | treatment/plate_id),
  +(1 | plate_id/fish_id),
  #+ (1 | plate_id) + (1 | fish_uid_num),
  data = df_5s,
  family = binomial
)

# # Model 2: Delay (only for responding fish)
# # Keep only responders
# df_5s_resp <- subset(df_5s, response == TRUE) %>% droplevels()
# 
# # Log-transform tap number
# df_5s_resp$tap_in_block_log <- log(df_5s_resp$tap_in_block)
# 
# # Remove values >= 0
# df_5s_resp$delay[df_5s_resp$delay <= 0] <- 0.00001
# df_5s_resp$delay_ms <- df_5s_resp$delay * 1000
# df_5s_resp <- df_5s_resp %>%
#   filter(delay_ms >= 30)
# 
# df_5s_resp_trim <- df_5s_resp %>% filter(delay_ms < quantile(delay_ms, 0.95, na.rm = TRUE))
# # Model for delay (Gaussian)
# glmm_delay <- glmmTMB(
#   log(delay_ms) ~ drug_condition * tap_in_block_log * treatment + (1 | plate_id/fish_id),
#   data = df_5s_resp,
#   family = Gamma(link = "log")
# )

# Compare Models
# anova(glmm_block_simpler, glmm_block)

# Test Model
# Simulate Residuals (QQ Plot)
# simulationOutput = simulateResiduals(glmm_block)
# plot(simulationOutput)
# plotResiduals(simulationOutput, df_5s$drug_condition)
# testDispersion(simulationOutput)     # Over/underdispersion check
# testZeroInflation(simulationOutput)  # Extra zeros?

# # Empirical Dist. of response variable
# ggplot(df_5s_resp, aes(x = delay)) +
#   geom_histogram(bins = 30, fill = "gray70", color = "black") +
#   #scale_x_log10() +
#   theme_minimal() +
#   labs(title = "Distribution of delay (ms)")

# Get Summary of GLMM
summary(glmm_block)
# VarCorr(glmm_block)

# Visualize Random Intercepts
# sjPlot::plot_model(glmm_block, type = "re")

# Compare to Models
# anova(glmm_block_no, glmm_block)

# ==============================================================================
# Prepare for Plotting
# Build a data grid for all combinations that exist in your dataset
valid_combos <- df_5s %>%
  distinct(drug_condition, treatment)

new_data <- expand.grid(
  tap_in_block = 1:20
) %>%
  mutate(tap_in_block_log = log(tap_in_block)) %>%
  tidyr::crossing(valid_combos)   # only valid drug_condition × treatment pairs

# Predictions from Model
pred <- predict(
  glmm_block,
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
), by = .(treatment, tap_in_block, drug_condition)]

# Plot one treatment per panel
ggplot(new_data, aes(x = tap_in_block, color = drug_condition)) +
  geom_line(aes(y = fit), size = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = drug_condition),
              alpha = 0.2, color = NA) +
  geom_point(data = p_dt_5s,
             aes(x = tap_in_block, y = response, color = drug_condition),
             alpha = 0.6, size = 1.5) +
  facet_wrap(~ treatment, scales = "free_y") +
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
    y = "Predicted response probability",
    color = "Drug condition",
    fill  = "Drug condition",
    title = "Habituation curves by treatment",
    subtitle = "Predicted from GLMM (binomial logit link)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )
ggsave(
  filename = "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/Habituation_GLMM_response_prob.pdf",
  plot = last_plot(),       # saves the most recent ggplot
  width = 10,
  height = 7,
  units = "in"
)
ggsave(
  filename = "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/Habituation_GLMM_response_prob.jpg",
  plot = last_plot(),
  width = 10,
  height = 7,
  units = "in",
  dpi = 600
)
# ==============================================================================
# Pair-wise Comparisons
# ------------------------------------------------------------------------------
# Habituation rate (slope)

# ------------------------------------------------------------------------------
# Overall responsiveness (vertical offset)


# # Pairs Tests
# pairs(emm_first_tap)
# pairs(emm_last_taps)
# pairs(emm_overall)


# ------------------------------------------------------------------------------
# Combine results and back-transform to probability scale
# Collect results



# View(results_all)
# ------------------------------------------------------------------------------
# Store to Disk
readr::write_csv(results_all,
                 "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/GLMM_response_prob_contrasts.csv")
model_tidy <- broom.mixed::tidy(glmm_block, effects = "fixed", conf.int = TRUE)
readr::write_csv(model_tidy,
                 "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/GLMM_response_prob_summary.csv")
ranefs <- broom.mixed::tidy(glmm_block, effects = "ran_vals")
readr::write_csv(ranefs,
                 "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/GLMM_response_prob_random_effects.csv")
sink("D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/GLMM_response_prob_summary.txt")
summary(glmm_block)
sink()

sink("D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/GLMM_response_prob_contrasts.txt")
kable(results_all, digits = 3, format = "simple")
sink()

# ==============================================================================
# ==============================================================================
#  Exponential fits to GLMM-predicted habituation curves (with error bars)
# ==============================================================================
# ---------------------------------------------------------------
# 1 Predict GLMM means and 95% CIs for each tap & drug
# ---------------------------------------------------------------
# 
# new_data <- expand.grid(
#   tap_in_block = 1:20,
#   drug_condition = c("control", "acute", "dk")
# )
# new_data$tap_in_block_log <- log(new_data$tap_in_block)
# 
# # Predict on link (logit) scale so we can get SEs
# pred <- predict(glmm_block, newdata = new_data, re.form = NA,
#                 se.fit = TRUE, type = "link")
# 
# new_data <- new_data %>%
#   mutate(
#     fit_link = pred$fit,
#     se_link  = pred$se.fit,
#     fit      = plogis(fit_link),
#     CI_low   = plogis(fit_link - 1.96 * se_link),
#     CI_high  = plogis(fit_link + 1.96 * se_link)
#   )
# 
# # ---------------------------------------------------------------
# # 2 Fit exponential decay curves per condition
# #       p = A * exp(-t/tau) + C
# # ---------------------------------------------------------------
# 
# fits_glmm <- new_data %>%
#   group_by(drug_condition) %>%
#   do({
#     fit <- nls(fit ~ A * exp(-tap_in_block / tau) + C,
#                data = .,
#                start = list(A = 0.8, tau = 3, C = 0.1))
#     df <- broom::augment(fit)
#     df$tau  <- coef(fit)[["tau"]]
#     df$se_tau <- summary(fit)$coefficients["tau", "Std. Error"]
#     df
#   })
# 
# # ---------------------------------------------------------------
# # 3 Summarize τ ± SE per drug condition (convert to seconds)
# # ---------------------------------------------------------------
# 
# ISI_seconds <- 5   # inter-stimulus interval in seconds
# 
# tau_table <- fits_glmm %>%
#   group_by(drug_condition) %>%
#   summarise(
#     tau  = mean(tau),
#     se_tau = mean(se_tau)
#   ) %>%
#   mutate(
#     tau_seconds = tau * ISI_seconds,
#     se_tau_seconds = se_tau * ISI_seconds,
#     label = sprintf(
#       "%s (τ = %.2f ± %.2f taps = %.1f ± %.1f s)",
#       drug_condition, tau, se_tau, tau_seconds, se_tau_seconds
#     )
#   )
# 
# # Merge labels for plotting legend
# fits_glmm <- left_join(fits_glmm, tau_table, by = "drug_condition")
# 
# # ---------------------------------------------------------------
# # 4 Plot: GLMM points + error bars + exponential fits + τ labels
# # ---------------------------------------------------------------
# 
# ggplot() +
#   # 95 % CI error bars for GLMM predictions
#   geom_errorbar(
#     data = new_data,
#     aes(x = tap_in_block, ymin = CI_low, ymax = CI_high, color = drug_condition),
#     width = 0.3, size = 0.8, alpha = 0.7
#   ) +
#   
#   # GLMM predicted probabilities (points)
#   geom_point(
#     data = new_data,
#     aes(x = tap_in_block, y = fit, color = drug_condition),
#     size = 2, alpha = 0.7
#   ) +
#   
#   # Exponential fits (lines)
#   geom_line(
#     data = fits_glmm,
#     aes(x = tap_in_block, y = .fitted, color = drug_condition),
#     linewidth = 1.2
#   ) +
#   
#   # Custom color palette and legend with τ ± SE
#   scale_color_manual(
#     name = "Drug condition (τ in taps and seconds)",
#     values = c("control" = "#1B9E77",
#                "acute"   = "#D95F02",
#                "dk"      = "#7570B3"),
#     labels = tau_table$label
#   ) +
#   
#   scale_y_continuous(labels = percent_format(accuracy = 1)) +
#   scale_x_continuous(breaks = seq(1, 20, 5)) +
#   labs(
#     x = "Tap number (within block)",
#     y = "Predicted response probability",
#     title = "Exponential fits to GLMM-predicted habituation curves"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "top",
#     panel.grid.minor = element_blank(),
#     strip.text = element_text(face = "bold")
#   )

# SOME NOTES:
# Biological intuition:
# Raw tap number: implies the fish loses a fixed amount of responsiveness each tap 
# (same drop between tap 2→3 as between tap 18→19).
# log(tap): implies the largest change happens in the first few taps, 
# and the rate of change slows down — exactly what we expect for neural or behavioral habituation processes.
