library(readr)
library(lme4)
library(lmerTest)  # for p-values
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(scales)
library(data.table)
library(tibble)
library(tidyr)
library(purrr)
library(furrr)


# Load data from csv file
df <- read_csv("D:/WorkingData/PrTecDA_Data/PrDA_somas_Ca_imaging/R_data/stim_motor_data.csv")

df$stim <- as.factor(df$stim)
df$motor <- as.factor(df$motor)
df$fish <- as.factor(df$fish)
df$roi <- as.factor(df$roi)
df$condition <- as.factor(df$condition)
df$stim_name <- as.factor(df$stim_name)

df <- df %>%
  mutate(fish = recode(fish,
                       "1" = "1",
                       "3" = "2",
                       "4" = "3",
                       "5" = "4",
                       "6" = "5",
                       "9" = "6"))
df <- df %>%
  mutate(fish = factor(fish, levels = c("1", "2", "3", "4", "5", "6")))


# Remove baseline condition (pseudo)
df <- df %>%
  filter(condition != "baseline")

# ==============================================================================
# Mixed Models
# =============================================================================
# ==============================================================================
# Categorical GLMM
df$condition <- factor(df$condition, 
                       levels = c("stim_only", "motor_only", "stim+motor"))

df$score <- df$score + 1e-6  # avoid zeros
df$score_norm99 <- ifelse(df$score_norm99 < 0, 0, df$score_norm99)
df$score_norm99 <- df$score_norm99 + 1e-6
df$log1p_score <- log1p(df$score)
df$condition <- factor(df$condition, 
                       levels = c("stim_only", "motor_only", "stim+motor"))

glmm_cat <- glmmTMB(
  score ~ condition + (1 | fish) + (1 | roi),
  data = df,
  family = Gamma(link = "log")
)

summary(glmm_cat)

# ------------------------------------------------------------------------------
# Validation
simres <- simulateResiduals(glmm_cat)
plot(simres)

# ------------------------------------------------------------------------------
# Pairwise comparisons
emm_glmm <- emmeans(glmm_cat, ~ condition, type = "response")
pairs(emm_glmm)

# Synergy Contrast
contrast(
  emm_glmm,
  list(synergy = c(-0.5, -0.5, 1)),  # weights correspond to (stim_only, motor_only, stim+motor)
  adjust = "none"  # no multiple comparison correction
)

emm_glmm_link <- emmeans(glmm_cat, ~ condition, type = "link")
contrast(emm_glmm_link, list(synergy = c(-0.5, -0.5, 1)), adjust = "none")

# Expected additive mean (for illustration):
pred_df <- as.data.frame(emm_glmm)

# use additive mean so that the sum of weihgts = 0, since emmeans expects normalized weights
pred_df$additive_mean <- (pred_df$response[pred_df$condition == "stim_only"] +
                            pred_df$response[pred_df$condition == "motor_only"]) / 2
pred_df$additive_sum <- (
  pred_df$response[pred_df$condition == "stim_only"] +
    pred_df$response[pred_df$condition == "motor_only"]
)

pred_df$response[pred_df$condition == "stim+motor"]
# ------------------------------------------------------------------------------
# Plots

# Extract numeric additive mean for plotting reference
add_mean <- unique(pred_df$additive_mean)
add_sum <- unique(pred_df$additive_sum)

ggplot(pred_df, aes(x = condition, y = response)) +
  geom_point(size = 4, color = "steelblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.1, color = "steelblue") +
  geom_hline(yintercept = add_mean, linetype = "dashed", color = "firebrick", linewidth = 1) +
  geom_hline(yintercept = add_sum, linetype = "dashed", color = "firebrick", linewidth = 1) +
  annotate("text", x = 2.5, y = add_mean * 1.08, label = "Additive mean", color = "firebrick", hjust = 0.5) +
  annotate("text", x = 2.5, y = add_sum * 1.05, label = "Additive sum", color = "firebrick", hjust = 0.5) +
  scale_y_continuous("Predicted score (Gamma GLMM, response scale)", expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, title = "Synergistic enhancement of neuronal responses") +
  theme_classic(base_size = 14)



