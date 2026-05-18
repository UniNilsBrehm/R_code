library(readr)
library(lme4)
library(performance)
library(sjPlot)
library(emmeans)
library(lmerTest)   # optional, adds p-values for fixed effects
library(dplyr)
library(glmmTMB)
library(ggplot2)
library(scales)
library(DHARMa)
library(data.table)
library(tibble)


# Load data from csv file
df <- read_csv("D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/Haloperidol_20uM_96wellplates.csv")

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
df_5s <- subset(df, ISI_block %in% c("ISI_5s_block1", "ISI_5s_block2"))

# Set References
df_5s$drug_condition <- factor(df_5s$drug_condition, levels = c("control", "acute", "dk"))
df_5s$ISI_block <- factor(df_5s$ISI_block, levels = c("ISI_5s_block1", "ISI_5s_block2"))

# Log transform tap number (for exponential decay)
df_5s$tap_in_block_log <- log(df_5s$tap_in_block)

# GLMM:
# each fish belongs to only one plate: (1 | plate_id/fish_id)
glmm_block <- glmmTMB(
  response ~ ISI_block * drug_condition * tap_in_block_log 
  #+(1 + tap_in_block_log | fish_id) + (1 | plate_id),
  +(1 | plate_id/fish_id),
  #+(1 | fish_id) + (1 | plate_id),
  #+(1 | plate_id/fish_id),
  #,
  data = df_5s,
  family = binomial
)

# glmm_block_simpler <- glmmTMB(
#   response ~ ISI_block * drug_condition * tap_in_block_log
#   +(1 | plate_id/fish_id),
#   data = df_5s,
#   family = binomial
# )
# 
# # Compare Models
# AIC(glmm_block, glmm_block_simpler)
# BIC(glmm_block, glmm_block_simpler)
# anova(glmm_block_simpler, glmm_block)

# Test Model
# Simulate Residuals (QQ Plot)
simulationOutput = simulateResiduals(glmm_block)
plot(simulationOutput)
plotResiduals(simulationOutput, df_5s$drug_condition)
testDispersion(simulationOutput)     # Over/underdispersion check
testZeroInflation(simulationOutput)  # Extra zeros?

# Get Summary of GLMM
summary(glmm_block)
# VarCorr(glmm_block)

# Visualize Random Intercepts
sjPlot::plot_model(glmm_block, type = "re")

# # Extract random slope effects for fish
# re_fish <- ranef(glmm_block)$cond$fish_id %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column("fish_id") %>%
#   rename(
#     intercept = `(Intercept)`,
#     slope = tap_in_block_log
#   )

# table(df_5s$ISI_block)
# table(df_5s$drug_condition)
# Arange Data for Plotting
new_data <- expand.grid(
  ISI_block = c("ISI_5s_block1", "ISI_5s_block2"),
  tap_in_block = 1:20,
  drug_condition = c("control", "acute", "dk")
)

# Match the model predictor: add the log-transformed version
new_data$tap_in_block_log <- log(new_data$tap_in_block)


# Get Predictions from the Model
pred = predict(glmm_block, newdata = new_data, re.form = NA, se.fit = T)
# str(pred)

# Back-Transform results
new_data$fit = plogis(pred$fit)
new_data$CI_low = plogis(pred$fit - 1.95*pred$se.fit)
new_data$CI_high = plogis(pred$fit + 1.95*pred$se.fit)

# Plotting
p_dt_5s <- data.table(df_5s)
p_dt_5s <- p_dt_5s[,.(response = sum(response == 1)/.N), by = .(tap_in_block, ISI_block, drug_condition)]

ggplot(new_data, aes(x = tap_in_block, color = drug_condition)) +
  facet_wrap(~ISI_block)+
  geom_line(aes(y = fit), size = 2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = drug_condition),
              alpha = 0.2, color = NA) +
  geom_point(data = p_dt_5s, aes(x = tap_in_block, y = response), alpha = 0.6) +
  labs(
    x = "Tap number (within block)",
    y = "Predicted response probability",
    color = "Drug condition",
    fill  = "Drug condition"
  ) +
  theme_minimal(base_size = 14)
  

# Plot random intercept-slope relationship for each fish
# ggplot(re_fish, aes(x = slope, y = intercept)) +
#   geom_point(size = 2, alpha = 0.8, color = "steelblue") +
#   geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
#   labs(
#     x = "Random slope for tap_in_block_log (habituation rate deviation)",
#     y = "Random intercept (baseline response deviation)",
#     title = "Random intercept–slope relationship per fish"
#   ) +
#   theme_minimal(base_size = 14)

# ==============================================================================
# Pair-wise Comparisons
# ------------------------------------------------------------------------------
# EMM: Estimated Marginal Means of the linear predictors (log-odds of response)

emm_by_tap <- emmeans(
  glmm_block,
  ~ drug_condition | ISI_block * tap_in_block_log,   # include tap_in_block_log as a by-factor
  at = list(tap_in_block_log = log(taps)),
  type = "response",
  cov.reduce = FALSE
)

emm_by_tap_df <- as.data.frame(emm_by_tap) %>%
  mutate(tap_in_block = exp(tap_in_block_log))   # back-transform log(tap) → tap numbe
# First Tap, Intercept of Habituation Curve
emm_first_tap <- emmeans(glmm_block,
                       ~ drug_condition | ISI_block,
                       at = list(tap_in_block_log = 0),
                       type = "response")

emm_last_taps <- emmeans(glmm_block,
                         ~ drug_condition | ISI_block,
                         at = list(tap_in_block_log = log(15:20)),
                         type = "response")

# All taps (mean y shift of curve): mean response probability across the entire block
emm_overall <- emmeans(
  glmm_block,
  ~ drug_condition | ISI_block,
  at = list(tap_in_block_log = log(1:20)),   # all taps within block
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)

# Habituation slopes
emm_slopes <- emtrends(glmm_block,
                       ~ drug_condition | ISI_block,
                       var = "tap_in_block_log")
# emtrends() output:
# Estimated slope of the habituation curve (change in log-odds of response per 1-unit increase in log(tap number))

# BLOCK 1 vs 2
# Overall responsiveness between blocks: First Tap
emm_block_first_tap <- emmeans(glmm_block,
                          ~ ISI_block | drug_condition,
                          at = list(tap_in_block_log = 0),
                          type = "response")

# Overall responsiveness between blocks: All taps (mean y shift of curve): mean response probability across the entire block
emm_block_overall <- emmeans(glmm_block,
                               ~ ISI_block | drug_condition,
                               at = list(tap_in_block_log = log(1:20)),   # all taps within block
                               cov.reduce = mean,                         # average over those predictions
                               type = "response")


# Habituation slope between blocks
emm_block_slope <- emtrends(glmm_block,
                            ~ ISI_block | drug_condition,
                            var = "tap_in_block_log")

pairs(emm_first_tap)
pairs(emm_last_taps)
pairs(emm_overall)
pairs(emm_slopes)
pairs(emm_block_first_tap)
pairs(emm_block_overall)
pairs(emm_block_slope)

pairs_by_tap <- contrast(emm_by_tap, method = "pairwise", type = "response") %>%
  as.data.frame() %>%
  mutate(
    tap_in_block_log = as.numeric(as.character(tap_in_block_log)),
    tap_in_block = exp(tap_in_block_log),
    ISI_block = as.factor(ISI_block),
    sig_clean = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "n.s."
    )
  )

# Collect results
res_overall_within <- pairs(emm_overall) %>%
  as.data.frame() %>%
  mutate(contrast_type = "DrugCondition_withinBlock",
         test_type = "Overall_Responsiveness")

res_slope_within <- pairs(emm_slopes) %>%
  as.data.frame() %>%
  mutate(contrast_type = "DrugCondition_withinBlock",
         test_type = "Habituation_Slope")

res_overall_between <- pairs(emm_block_resp) %>%
  as.data.frame() %>%
  mutate(contrast_type = "Block_withinCondition",
         test_type = "Overall_Responsiveness")

res_slope_between <- pairs(emm_block_slope) %>%
  as.data.frame() %>%
  mutate(contrast_type = "Block_withinCondition",
         test_type = "Habituation_Slope")

results_all <- bind_rows(
  res_overall_within,
  res_slope_within,
  res_overall_between,
  res_slope_between
) %>%
  # Add transformed & readable values
  mutate(
    odds_ratio = exp(estimate),
    OR_low  = exp(estimate - 1.96 * SE),
    OR_high = exp(estimate + 1.96 * SE),
    prob_est   = plogis(estimate),
    prob_low   = plogis(estimate - 1.96 * SE),
    prob_high  = plogis(estimate + 1.96 * SE),
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  mutate(across(
    c(estimate, SE, odds_ratio, OR_low, OR_high,
      prob_est, prob_low, prob_high, z.ratio, p.value),
    ~ round(., 4)
  )) %>%
  # Arrange nicely
  arrange(contrast_type, test_type) %>%
  # Select key columns
  select(
    contrast_type, test_type,
    contrast, ISI_block, drug_condition,
    estimate, SE, z.ratio, p.value, sig,
    odds_ratio, OR_low, OR_high,
    prob_est, prob_low, prob_high
  )

# Convert Habituation rates to half-life time constants (tap where it reaches 50%)
# Assume emm_slopes already exists (from emtrends)
emm_slopes_df <- as.data.frame(emm_slopes)

# Compute time constant (tau)
emm_slopes_df <- emm_slopes_df %>%
  mutate(
    tau_taps = 2^(1 / abs(tap_in_block_log.trend)),   # half-life in taps
    slope_sign = ifelse(tap_in_block_log.trend < 0, "decay", "growth")
  )

# Round for readability
emm_slopes_df <- emm_slopes_df %>%
  mutate(across(c(tap_in_block_log.trend, SE, tau_taps), ~ round(., 3)))

print(emm_slopes_df)

# Convert Habituation rate to tau time constants (1/e)
emm_slopes_df <- as.data.frame(emm_slopes)

emm_slopes_df <- emm_slopes_df %>%
  mutate(
    tau_exp = exp(1) / abs(tap_in_block_log.trend),   # 1/e time constant (≈2.718/|slope|)
    tau_half = log(2) / abs(tap_in_block_log.trend)   # 50% half-life for comparison
  )

emm_slopes_df %>%
  mutate(across(c(tap_in_block_log.trend, tau_exp, tau_half), round, 3))

# Compute CIs for taus
emm_slopes_df <- as.data.frame(emm_slopes)

emm_slopes_df <- emm_slopes_df %>%
  mutate(
    # --- Time constants ---
    tau_exp  = exp(1) / abs(tap_in_block_log.trend),   # 1/e time constant
    tau_half = log(2) / abs(tap_in_block_log.trend),   # 50% half-life
    
    # --- Confidence intervals (from slope CI bounds) ---
    tau_exp_low  = exp(1) / abs(asymp.UCL),  # smaller |slope| → larger tau
    tau_exp_high = exp(1) / abs(asymp.LCL),
    tau_half_low  = log(2) / abs(asymp.UCL),
    tau_half_high = log(2) / abs(asymp.LCL),
    
    slope_sign = ifelse(tap_in_block_log.trend < 0, "decay", "growth")
  ) %>%
  mutate(across(
    c(tap_in_block_log.trend, SE, tau_exp, tau_half,
      tau_exp_low, tau_exp_high, tau_half_low, tau_half_high),
    ~ round(., 3)
  ))

view(emm_slopes_df)

# Plot

ggplot(emm_slopes_df,
       aes(x = drug_condition, y = tau_exp,
           color = ISI_block, group = ISI_block)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = tau_exp_low, ymax = tau_exp_high),
                width = 0.2, position = position_dodge(width = 0.3)) +
  labs(
    x = "Drug condition",
    y = expression(tau~"(1/e time constant, taps)"),
    color = "ISI Block",
    title = "Habituation time constants with 95% CIs"
  ) +
  theme_minimal(base_size = 14)
