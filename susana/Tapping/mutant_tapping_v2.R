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
library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)

source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/Tapping/utils.R")
# ==============================================================================
# Load Data
# ==============================================================================
# Base directory for saving results
base_dir<- "C:/Users/NilsPC/Desktop/Susana/Susana/Tapping"
file_dir <- file.path(base_dir, "data_files", 'SPZ_ISI5_Nils2Blocks_Tap.csv')
res <- load_data_tapping(file_dir, move_th = 1, take_peak = 0, ref = "tyr")

df_final <- res$df_final
df_final_sub <- res$df_final_sub

# Model
# Response Probability
responses <- responses %>%
  mutate(response_yes_no_num = as.numeric(response_yes_no))
responses_sub <- responses %>%
  filter(response_yes_no == TRUE)

m_response_prob <- glmmTMB(
  response_yes_no_num ~ genotype * phase * log(block_tap) + (1 | Well),
  data = responses %>% filter(phase != "isolated"),
  family = binomial()
)
summary(m_response_prob)

# Peak distance moved
m_peak <- glmmTMB(
  peak_distance  ~ genotype * phase * log(block_tap) + (1 | Well),
  data = responses_sub %>% filter(phase != "isolated"),
  family = Gamma(link="log")
)
summary(m_peak)

# Summed distance moved
m_sum <- glmmTMB(
  summed_distance  ~ genotype * phase * log(block_tap) + (1 | Well),
  data = responses_sub %>% filter(phase != "isolated"),
  family = Gamma(link="log")
)
summary(m_sum)

# delay
m_delay <- glmmTMB(
  response_delay  ~ genotype * phase * log(block_tap) + (1 | Well),
  data = responses_sub %>% filter(phase != "isolated"),
  family = gaussian(link='identity')
)
summary(m_delay)


# Check Residuals
library(DHARMa)
sim <- simulateResiduals(m_response_prob)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

# Model Plots
plot_model_response(m_peak, responses_sub, "peak_distance", "Peak distance moved")

plot_model_response(m_sum, responses_sub, "summed_distance", "Summed distance moved")

plot_model_response(m_delay, responses_sub, "response_delay", "Response delay")

plot_model_response(m_response_prob, responses, "response_yes_no_num", "Response probability")
