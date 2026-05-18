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

# Load helper functions (validation, plotting, EMM utilities, reporting)
source("utils.R")

# Base directory for saving results
base_dir <- "D:/WorkingData/Susana/results/"


# ==============================================================================
# 1. Load and Prepare Data
# ==============================================================================
res_massed <- load_data(
  "D:/WorkingData/Susana/SPZ_Massed_Training_7Nov2025.csv",
  # keep = c("ABTL", "th, th2, tyr", "th, tyr")
)
res_spaced <- load_data(
  "D:/WorkingData/Susana/SPZ_Spaced_Training_Nov2025.csv",
  # keep = c("ABTL", "th, th2, tyr", "th, tyr")
)

df_massed <- res_massed$df_final
df_massed_sub <- res_massed$df_final_sub
df_spaced <- res_spaced$df_final
df_spaced_sub <- res_spaced$df_final_sub


# Number of animals per genotype
n_massed <- df_massed %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)
n_spaced <- df_spaced %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

# --- Filter by Block ---
df_massed_b2 <- df_massed %>%
  filter(Block == "2") %>%
  mutate(Training = "Massed")

df_spaced_b5 <- df_spaced %>%
  filter(Block == "5") %>%
  mutate(Training = "Spaced")

# --- Combine ---
df_combined <- bind_rows(df_massed_b2, df_spaced_b5)
df_combined <- droplevels(df_combined)

colSums(is.na(df_combined))
df_combined <- df_combined %>%
  mutate(Block = recode(Block,
                        "2" = "Massed",
                        "5" = "Spaced"))

# ==============================================================================
# Fit GLMM Model
# ==============================================================================
# Response Probability
# ==============================================================================

m_prob <- glmmTMB(
  move ~ Genotype * Block * stimulus_log + (1 | Video/Well),
  family = binomial(link = "logit"),
  data = df_combined
)

validate_model(m_prob, df_combined)

g_prob <- plot_habituation(df_combined, m_prob, label = "Response prob.", transform = "plogis")
g_prob


# Difference of Genotypes in massed and spaced
emm <- emmeans(
  m_prob,
  ~ Genotype | Block,
  type = "response"
)
pairs(emm)

# Differences of Genotypes between massed and spaced
emm_between <- emmeans(
  m_prob,
  ~ Block | Genotype,
  type = "response"
)

pairs(emm_between)

emm_between_slopes <- emtrends(m_prob, ~ Block | Genotype, var = "stimulus_log")
pairs(emm_between_slopes)

# Peak Distance
# ==============================================================================
df_combined_sub <- df_combined[df_combined$move > 0, ]

m_peak <- glmmTMB(
  max_peak ~ Genotype * stimulus_log * Block + (1 | Video/Well),
  family = Gamma(link = "log"),
  data = df_combined_sub
)

validate_model(m_peak, df_combined_sub)

g_peak <- plot_habituation(df_combined_sub, m_peak, label = "Peak distance moved (mm)", transform = "exp")
g_peak


# Difference of Genotypes in massed and spaced
emm_peak <- emmeans(
  m_peak,
  ~ Genotype | Block,
  type = "response"
)
pairs(emm_peak)

# Differences of Genotypes between massed and spaced
emm_peak_between <- emmeans(
  m_peak,
  ~ Block | Genotype,
  #at = list(stimulus_log = log(1:8)),
  type = "response"
)

pairs(emm_peak_between)

emm_peak_between_slopes <- emtrends(m_peak, ~ Block | Genotype, var = "stimulus_log")
pairs(emm_peak_between_slopes)

# Summed Distance
# ==============================================================================
m_sum <- glmmTMB(
  max_cumsum ~ Genotype * stimulus_log * Block + (1 | Video/Well),
  family = Gamma(link = "log"),
  data = df_combined_sub
)

validate_model(m_sum, df_combined_sub)

g_sum <- plot_habituation(df_combined_sub, m_sum, label = "Summed distance moved (mm)", transform = "exp")
g_sum


# Difference of Genotypes in massed and spaced
emm_sum <- emmeans(
  m_sum,
  ~ Genotype | Block,
  type = "response"
)
pairs(emm_sum)

# Differences of Genotypes between massed and spaced
emm_sum_between <- emmeans(
  m_sum,
  ~ Block | Genotype,
  #at = list(stimulus_log = log(1:8)),
  type = "response"
)

pairs(emm_sum_between)

emm_sum_between_slopes <- emtrends(m_sum, ~ Block | Genotype, var = "stimulus_log")
pairs(emm_sum_between_slopes)


# Response Delay
# ==============================================================================
# --- Model 3: Response Delay (Gaussian LMM) ----------------------------------
m_delay_gaussian <- glmmTMB(
  delay ~ Genotype * stimulus_log * Block + (1 | Video/Well),
  family = gaussian(link = "identity"),
  data = df_combined_sub
)

# --- Model 4: Response Delay (Poisson LMM) ----------------------------------
m_delay_poisson <- glmmTMB(
  delay ~ Genotype * stimulus_log * Block + (1 | Video/Well),
  family = poisson(link = "log"),
  data = df_combined_sub
)

validate_model(m_delay_gaussian, df_combined_sub)
validate_model(m_delay_poisson, df_combined_sub)

g_delay <- plot_habituation(df_combined_sub, m_delay_poisson, label = "Response delay (s)", transform = "exp")
g_delay


# Difference of Genotypes in massed and spaced
emm_delay <- emmeans(
  m_delay_poisson,
  ~ Genotype | Block,
  type = "response"
)
pairs(emm_delay)

# Differences of Genotypes between massed and spaced
emm_delay_between <- emmeans(
  m_delay_poisson,
  ~ Block | Genotype,
  #at = list(stimulus_log = log(1:8)),
  type = "response"
)

pairs(emm_delay_between)

emm_delay_between_slopes <- emtrends(m_delay_poisson, ~ Block | Genotype, var = "stimulus_log")
pairs(emm_delay_between_slopes)

