###############################################################################
# GLMM Analysis: Dark Flash 60 s – Response Probability
# Author: Nils Brehm
# Date: 02/2026
#
# Description:

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
library(here)

# Load helper functions (validation, plotting, EMM utilities, reporting)
source("C:/UniFreiburg/Code/R_code/susana/utils.R")

# Base directory for saving results
base_dir <- "D:/WorkingData/Susana/results/darkflash_60s"


# ==============================================================================
# 1. Load and Prepare Data
# ==============================================================================
# Optionally subset to certain genotypes
df_raw <- read_csv("D:/WorkingData/Susana/SPZ_ISI60_DF.csv")
res <- load_data_ISI60_DF("D:/WorkingData/Susana/SPZ_ISI60_DF.csv", move_th = 0)

df_final <- res$df_final
df_final_sub <- res$df_final_sub

# ----------
df_final <- df_final %>%
  group_by(Block, Well, Video) %>%
  mutate(animal_id = cur_group_id()) %>%
  ungroup()


df_one_animal <- df_final %>%
  filter(animal_id == 200)


# Number of animals per genotype
n_per_genotype <- df_final %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

print(n_per_genotype)


# ==============================================================================
# 2. Exploratory Data Visualization
# ==============================================================================
# Histogram of response distribution by genotype
h1 <- ggplot(df_final, aes(x = move)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Distribution of Binary Movement Responses",
    x = "Response (0 = No, 1 = Yes)",
    y = "Count"
  )

# Display and save plot
h1
ggsave(
  filename = file.path(base_dir, "response_prob", "spaced_response_binary_dist.pdf"),
  plot = h1,
  width = 6, height = 4
)


# ==============================================================================
# 3. Fit GLMM Model
# ==============================================================================
message("Fitting GLMM model for response probability...")

m_prob <- glmmTMB(
  move ~ Genotype * Block * stimulus_log + (1 | Video/Well),
  family = binomial(link = "logit"),
  data = df_final
)

# --- Model 1: Peak Movement (Gamma GLMM) -------------------------------------
m_peak <- glmmTMB(
  max_peak ~ Genotype * stimulus_log * Block + (1 | Video/Well),
  family = Gamma(link = "log"),
  data = df_final_sub
)

# --- Model 2: Summed Distance (Gamma GLMM) -----------------------------------
m_sum <- glmmTMB(
  max_cumsum ~ Genotype * stimulus_log * Block + (1 | Video/Well),
  family = Gamma(link = "log"),
  data = df_final_sub)

# ==============================================================================
# 4. Model Validation
# ==============================================================================
message("Validating model residuals...")

validate_model(m_prob, df_final)


# ==============================================================================
# 5. Plot Habituation Curves
# ==============================================================================
message("Plotting habituation curves...")

g_prob <- plot_habituation(df_final, m_prob, label='Response prob.', transform = "plogis")
g_prob

g_peak <- plot_habituation(df_final_sub, m_peak, label='peak distance', transform = "exp")
g_peak

g_sum <- plot_habituation(df_final_sub, m_sum, label='summed distance', transform = "exp")
g_sum

ggsave(
  filename = file.path(base_dir, "response_prob", "spaced_response_prob_habituation_curves.pdf"),
  plot = g_prob,
  width = 20, height = 8
)


# ==============================================================================
# 6. Estimated Marginal Means (EMMs)
# ==============================================================================
message("Calculating EMMs for response probability...")

# --- Block-level EMMs (dynamic across all blocks) ------------------------------
emm_prob <- get_emm_blocks(m_prob, df_final)

# --- Habituation slopes (within blocks) --------------------------------------
emm_prob_slopes <- emtrends(m_prob, ~ Genotype | Block, var = "stimulus_log")
emm_prob_slopes_pairs <- pairs(emm_prob_slopes)

# --- Between-block comparisons (1 vs 2, 2 vs 3, etc.) -------------------------
emm_prob_between <- get_emm_consecutive_blocks(m_prob, df_final, n_stim = 9)

# --- Between-block slopes -----------------------------------------------------
emm_prob_between_blocks_slopes <- emtrends(m_prob, ~ Block | Genotype, var = "stimulus_log")
emm_prob_between_blocks_slopes_pairs <- pairs(emm_prob_between_blocks_slopes)


# ==============================================================================
# 7. Export EMM Results
# ==============================================================================
write_emm_report(
  model = m_prob,
  emm_blocks = emm_prob,
  emm_slopes = emm_prob_slopes,
  emm_slopes_pairs = emm_prob_slopes_pairs,
  emm_between = emm_prob_between,
  emm_between_blocks_slopes = emm_prob_between_blocks_slopes,
  emm_between_blocks_slopes_pairs = emm_prob_between_blocks_slopes_pairs,
  outfile = file.path(base_dir, "response_prob", "spaced_glmm_response_prob_comparisons.txt"),
  model_name = "GLMM Response Probability Model"
)

# =============================================================================
# 8. Completion Message
# ==============================================================================
message("All analyses completed successfully — results saved to disk.")

