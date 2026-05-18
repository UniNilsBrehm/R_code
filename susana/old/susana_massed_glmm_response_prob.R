###############################################################################
# GLMM Analysis: Massed Training – Response Probability
# Author: Nils Brehm
# Date: 11/2025
#
# Description:
#   This script analyzes habituation behavior data from massed training
#   experiments. It fits a GLMM model predicting the probability of movement
#   (response probability) across stimulus blocks, validates the model,
#   visualizes habituation curves, and computes estimated marginal means (EMMs)
#   and contrasts between genotypes and blocks.
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

# Load helper functions (validation, plotting, EMM utilities, reporting)
source("C:/UniFreiburg/Code/R_code/susana/utils.R")

# Base directory for saving results
base_dir <- "D:/WorkingData/Susana/results/massed"


# ==============================================================================
# 1. Load and Prepare Data
# ==============================================================================
res <- load_data(
  "D:/WorkingData/Susana/SPZ_Massed_Training_7Nov2025.csv",
  # keep = c("ABTL", "th, th2, tyr", "th, tyr"),
  keep = c("ABTL", "th, th2, tyr")
)

df_final <- res$df_final
df_final_sub <- res$df_final_sub

# Number of animals per genotype
n_per_genotype <- df_final %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

print(n_per_genotype)


# ==============================================================================
# 2. Exploratory Data Visualization
# ==============================================================================
h1 <- ggplot(df_final, aes(x = move)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Distribution of Binary Movement Responses (Massed)",
    x = "Response (0 = No, 1 = Yes)",
    y = "Count"
  )

h1
ggsave(
  filename = file.path(base_dir, "response_prob", "massed_response_binary_dist.pdf"),
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
save_plot(g_prob, file.path(base_dir, "response_prob", "massed_response_prob_habituation_curves_v2"), width=10, height=5, dpi=600)


# ==============================================================================
# 6. Estimated Marginal Means (EMMs)
# ==============================================================================
message("Calculating EMMs for response probability...")

# --- Block-level EMMs ---------------------------------------------------------
emm_prob <- get_emm_blocks(m_prob, df_final)

# --- Habituation slopes (within blocks) --------------------------------------
emm_prob_slopes <- emtrends(m_prob, ~ Genotype | Block, var = "stimulus_log")
emm_prob_slopes_pairs <- pairs(emm_prob_slopes)

# --- Between-block comparisons (1 vs 2, 2 vs 3, etc.) -------------------------
# Compare the n first an n last of blocks
emm_prob_between <- get_emm_consecutive_blocks(m_prob, df_final, n_stim = 9)

# Compare the n first an n first of blocks
emm_prob_between_first <- get_emm_consecutive_blocks_first(m_prob, df_final, n_stim = 9)

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
  emm_between_first = emm_prob_between_first,
  emm_between_blocks_slopes = emm_prob_between_blocks_slopes,
  emm_between_blocks_slopes_pairs = emm_prob_between_blocks_slopes_pairs,
  outfile = file.path(base_dir, "response_prob", "massed_glmm_response_prob_comparisons.txt"),
  model_name = "GLMM Response Probability Model (Massed Training)"
)

# ==============================================================================
# 8. Completion Message
# ==============================================================================
message("All analyses completed successfully — results saved to disk.")
