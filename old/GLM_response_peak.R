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
# Load data from csv file
df <- read_csv("D:/WorkingData/NoldusBehavior/Tapping_Habituation/tapping_peaks.csv")

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

# GLMM:
# each fish belongs to only one plate: (1 | plate_id/fish_id)
glmm_block <- glmmTMB(
  response ~ ISI_block * drug_condition * tap_in_block 
  +(1 | plate_id/fish_id),
  data = df_5s,
  family = binomial
)

# Simulate Residuals (QQ Plot)
res = simulateResiduals(glmm_block)
plot(res)

# Get Summary of GLMM
summary(glmm_block)