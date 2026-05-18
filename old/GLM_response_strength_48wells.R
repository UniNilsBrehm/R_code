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
library(broom.mixed)

# Load data from csv file
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
df_5s$auc <- df_5s$auc + 1e-6  # avoid zeros
glmm_block <- glmmTMB(
  auc ~ drug_condition * tap_in_block_log * treatment
  +(1 | plate_id/fish_id),
  data = df_5s,
  family = Gamma(link = "log")
)

# Compare Models
# AIC(glmm_block, glmm_tweedie)
# BIC(glmm_block, glmm_tweedie)
# anova(glmm_tweedie, glmm_block)

# Validate Model
# Simulate Residuals (QQ Plot)
# simulationOutput = simulateResiduals(glmm_block)
# plot(simulationOutput)
# plotResiduals(simulationOutput, df_5s$drug_condition)
# 
# testDispersion(simulationOutput)     # Over/underdispersion check
# testZeroInflation(simulationOutput)  # Extra zeros?
# testResiduals(simulationOutput)

# Get Summary of GLMM
summary(glmm_block)

# Visualize Random Intercepts
# sjPlot::plot_model(glmm_block, type = "re")


# ==============================================================================
# Prepare for Plotting
# Build valid combinations of drug_condition × treatment present in the data
valid_combos <- df_5s %>%
  distinct(drug_condition, treatment)

# Build prediction grid
new_data <- expand.grid(
  tap_in_block = 1:20
) %>%
  mutate(tap_in_block_log = log(tap_in_block)) %>%
  tidyr::crossing(valid_combos)   # only valid combos

# Predictions from the Gamma model
pred <- predict(
  glmm_block,
  newdata = new_data,
  re.form = NA,  # fixed effects only
  se.fit = TRUE
)

# Back-transform from log link → expected mean AUC
new_data$fit <- exp(pred$fit)
new_data$CI_low <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$CI_high <- exp(pred$fit + 1.96 * pred$se.fit)

# Relative Scale
# new_data <- new_data %>%
#   group_by(drug_condition, treatment) %>%
#   mutate(fit_rel = fit / first(fit))

# Compute observed mean AUCs for each combination
p_dt_5s <- data.table(df_5s)
p_dt_5s <- p_dt_5s[, .(
  auc_mean = mean(auc, na.rm = TRUE)
), by = .(treatment, tap_in_block, drug_condition)]

# Plot one treatment per panel
ggplot(new_data, aes(x = tap_in_block, color = drug_condition)) +
  geom_line(aes(y = fit), size = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = drug_condition),
              alpha = 0.2, color = NA) +
  geom_point(data = p_dt_5s,
             aes(x = tap_in_block, y = auc_mean, color = drug_condition),
             alpha = 0.6, size = 1.5) +
  facet_wrap(~ treatment, scales = "free_y") +
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20),
    limits = c(0, 20)
  ) +
  # Optional: adapt Y scale to your observed range
  scale_y_continuous(
    limits = c(0, 80),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Tap number (within block)",
    y = "Predicted AUC (Gamma GLMM, log link)",
    color = "Drug condition",
    fill  = "Drug condition",
    title = "Habituation of escape response (AUC) by treatment",
    subtitle = "Predicted from GLMM (Gamma, log link)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )
ggsave(
  filename = "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/Habituation_GLMM_AUC.pdf",
  plot = last_plot(),       # saves the most recent ggplot
  width = 10,
  height = 7,
  units = "in"
)
ggsave(
  filename = "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/Habituation_GLMM_AUC.jpg",
  plot = last_plot(),
  width = 10,
  height = 7,
  units = "in",
  dpi = 600
)
# ==============================================================================
# Pair-wise Comparisons

# ------------------------------------------------------------------------------
# Habituation Rate (slope)

# Estimate marginal slopes for each combination
emtr <- emtrends(
  glmm_block,
  var = "tap_in_block_log",
  specs = c("drug_condition", "treatment")
)

# Look at the estimated slopes
# summary(emtr)

# Test pairwise differences of slopes within each treatment
slope_tests <- pairs(emtr, by = "treatment") %>%
  as.data.frame() %>%
  mutate(test_type = "habituation")

# Negative slopes: Faster habituation
# Positive Difference: "A" has less negative slope -> slower habituation
# Negative Difference: "A" has steeper decline -> faster habituation 

# ------------------------------------------------------------------------------
# Overall responsiveness (vertical offset)
# Estimated mean AUC at average tap_in_block_log (overall responsiveness)
emm <- emmeans(
  glmm_block,
  ~ drug_condition | treatment,
  at = list(tap_in_block_log = mean(df_5s$tap_in_block_log))
)

# Compare dk vs control within each treatment
level_tests <- pairs(emm, by = "treatment") %>%
  as.data.frame() %>%
  mutate(test_type = "overall")

# ------------------------------------------------------------------------------
# Combine slope and level test results

results_all <- bind_rows(slope_tests, level_tests) %>%
  # (No filtering — keep all contrasts)
  # Add fold-change and percent-change columns
  mutate(
    fold_change = exp(estimate),                         # back-transform from log scale
    percent_change = (exp(estimate) - 1) * 100           # convert to %
  ) %>%
  # Round numeric values
  mutate(across(c(estimate, fold_change, percent_change), round, 3)) %>%
  # Add significance markers
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  # Arrange neatly
  arrange(test_type, treatment, contrast) %>%
  # Select relevant columns
  select(treatment, test_type, contrast, estimate, SE, z.ratio, p.value,
         fold_change, percent_change, sig)

# View(results_all)

# Store to Disk
readr::write_csv(results_all,
                 "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/GLMM_response_strength_contrasts.csv")
model_tidy <- broom.mixed::tidy(glmm_block, effects = "fixed", conf.int = TRUE)
readr::write_csv(model_tidy,
                 "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/GLMM_response_strength_summary.csv")
ranefs <- broom.mixed::tidy(glmm_block, effects = "ran_vals")
readr::write_csv(ranefs,
                 "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/GLMM_response_strength_random_effects.csv")
sink("D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/GLMM_response_strength_summary.txt")
summary(glmm_block)
sink()

sink("D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/results/GLMM_response_strength_contrasts.txt")
kable(results_all, digits = 3, format = "simple")
sink()

