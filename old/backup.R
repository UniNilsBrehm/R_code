################################################################################
# Script: spaced_habituation_glmm.R
# Author: Nils Brehm
# Date:   10.11.2025
#
# Description:
#   This script analyzes habituation behavior data from spaced training experiments.
#   It fits a GLMM model to predict movement probability across stimulus blocks,
#   validates the model, generates diagnostic plots, and computes estimated marginal
#   means (EMMs) and contrasts between genotypes and blocks.
################################################################################

# ==============================================================================
# 1. Load Required Packages -----------------------------------------------------
# ==============================================================================
library(readr)        # Reading CSV files
library(DHARMa)       # Residual diagnostics for (G)LMMs
library(emmeans)      # Estimated marginal means (EMMs) and contrasts
library(glmmTMB)      # Generalized linear mixed models
library(ggplot2)      # Visualization
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(performance)  # Model diagnostics (AIC, R², etc.)
library(ggpubr)       # Publication-ready themes (theme_pubr)

# ==============================================================================
# 2. Load Helper Functions -----------------------------------------------------
# ==============================================================================
source("utils.R")  # Custom utility functions (validate_model, plot_habituation, etc.)

# Base directory for output
base_dir <- "D:/WorkingData/Susana/results/spaced"

# ==============================================================================
# 3. Load and Prepare Data -----------------------------------------------------
# ==============================================================================
# Load data from CSV
res <- load_data("D:/WorkingData/Susana/SPZ_Spaced_Training_Nov2025.csv")
df_final <- res$df_final
df_final_sub <- res$df_final_sub

# Number of animals per genotype
n_per_genotype <- df_final %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

# ==============================================================================
# 4. Exploratory Data Visualization --------------------------------------------
# ==============================================================================
# Histogram of response distribution by genotype
h1 <- ggplot(df_final, aes(x = move)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Distribution of Responses",
    x = "Response (0 = No, 1 = Yes)",
    y = "Count"
  )

# Display and save histogram
h1
ggsave(
  filename = file.path(base_dir, "response_prob", "spaced_response_binary_dist.pdf"),
  plot = h1,
  width = 6, height = 4
)

# ==============================================================================
# 5. Fit GLMM Model ------------------------------------------------------------
# ==============================================================================
# Response: move (binary)
# Fixed effects: Genotype * Block * stimulus_log
# Random effects: nested random intercepts for Video/Well

m1 <- glmmTMB(
  move ~ Genotype * Block * stimulus_log + (1 | Video/Well),
  family = binomial(link = "logit"),
  data = df_final
)

# ==============================================================================
# 6. Model Validation ----------------------------------------------------------
# ==============================================================================
validate_model(m1, df_final)

# ==============================================================================
# 7. Habituation Curves --------------------------------------------------------
# ==============================================================================
# Generate predicted response probabilities and plot habituation curves
g <- plot_habituation(df_final, m1, transform = "plogis")

# Display and save the plot
g
ggsave(
  filename = file.path(base_dir, "response_prob", "spaced_response_prob_habituation_curves.pdf"),
  plot = g,
  width = 20, height = 8
)

# ==============================================================================
# 8. Estimated Marginal Means (EMMs) -------------------------------------------
# ==============================================================================
# Compute and compare EMMs for responsiveness across blocks and stimuli subsets

# --- Block-level responsiveness ----------------------------------------------
emm_block1 <- emmeans(m1, ~ Genotype,
                      at = list(stimulus_log = log(1:477), Block = "1"),
                      cov.reduce = mean, type = "response"
)

emm_block2 <- emmeans(m1, ~ Genotype,
                      at = list(stimulus_log = log(1:9), Block = "2"),
                      cov.reduce = mean, type = "response"
)

# --- Within-block subsets -----------------------------------------------------
# First, middle, and last stimuli within Block 1
emm_block1_first_stim <- emmeans(m1, ~ Genotype,
                                 at = list(stimulus_log = log(1:10), Block = "1"),
                                 cov.reduce = mean, type = "response"
)

emm_block1_middle_stim <- emmeans(m1, ~ Genotype,
                                  at = list(stimulus_log = log(230:240), Block = "1"),
                                  cov.reduce = mean, type = "response"
)

emm_block1_last_stim <- emmeans(m1, ~ Genotype,
                                at = list(stimulus_log = log(466:477), Block = "1"),
                                cov.reduce = mean, type = "response"
)

# --- Habituation slopes -------------------------------------------------------
emm_slopes <- emtrends(m1, ~ Genotype | Block, var = "stimulus_log")

# ==============================================================================
# 9. Between-Block Comparisons -------------------------------------------------
# ==============================================================================
# Compare responsiveness and slopes between Block 1 and Block 2

# --- First stimulus in each block --------------------------------------------
emm_between_blocks_first_stimulus <- emmeans(
  m1, ~ Block | Genotype,
  at = list(stimulus_log = log(1)),
  type = "response"
)

# --- Recovery test: last 9 stim in Block 1 vs. first 9 stim in Block 2 --------
emm_block1_recovery <- emmeans(
  m1, ~ Block | Genotype,
  at = list(stimulus_log = log(466:477), Block = "1"),
  type = "response"
)

emm_block2_recovery <- emmeans(
  m1, ~ Block | Genotype,
  at = list(stimulus_log = log(1:9), Block = "2"),
  type = "response"
)

# Combine and contrast
combined_emms <- rbind(emm_block1_recovery, emm_block2_recovery)
block_comparisons <- contrast(
  combined_emms,
  method = "revpairwise",
  by = "Genotype",
  adjust = "Tukey"
)

# --- Last stimulus of Block 1 vs first stimulus of Block 2 --------------------
emm_block1_last_recovery <- emmeans(
  m1, ~ Block | Genotype,
  at = list(stimulus_log = log(477), Block = "1"),
  type = "response"
)

emm_block2_first_recovery <- emmeans(
  m1, ~ Block | Genotype,
  at = list(stimulus_log = log(1), Block = "2"),
  type = "response"
)

combined_last_first_emms <- rbind(emm_block1_last_recovery, emm_block2_first_recovery)
block_last_first_comparisons <- contrast(
  combined_last_first_emms,
  method = "revpairwise",
  by = "Genotype",
  adjust = "Tukey"
)

# --- Slope comparison between blocks -----------------------------------------
emm_between_blocks_slopes <- emtrends(
  m1, ~ Block | Genotype,
  var = "stimulus_log"
)

# ==============================================================================
# 10. Save Results to Text File ------------------------------------------------
# ==============================================================================
# Convert EMMs to readable data frames using custom helper functions
results <- list(
  responsiveness_block1                = pretty_pairs(emm_block1),
  responsiveness_block2                = pretty_pairs(emm_block2),
  responsiveness_block1_first          = pretty_pairs(emm_block1_first_stim),
  responsiveness_block1_middle         = pretty_pairs(emm_block1_middle_stim),
  responsiveness_block1_last           = pretty_pairs(emm_block1_last_stim),
  habituation_slope_block1_and_block2  = pretty_pairs(emm_slopes),
  between_blocks_responsiveness_first  = pretty_pairs(emm_between_blocks_first_stimulus),
  between_blocks_habituation_slope     = pretty_pairs(emm_between_blocks_slopes),
  between_block_10_stimuli             = pretty_blocks(block_comparisons),
  between_block_last_first             = pretty_blocks(block_last_first_comparisons)
)

# Descriptive section headers for the report
descriptions <- list(
  responsiveness_block1 = "Pairwise contrasts of responsiveness within Block 1.",
  responsiveness_block2 = "Pairwise contrasts of responsiveness within Block 2.",
  responsiveness_block1_first = "Responsiveness contrasts within Block 1 (first stimulus subset).",
  responsiveness_block1_middle = "Responsiveness contrasts within Block 1 (middle stimuli subset).",
  responsiveness_block1_last = "Responsiveness contrasts within Block 1 (last stimulus subset).",
  habituation_slope_block1_and_block2 = "Comparisons of habituation slopes between Block 1 and Block 2.",
  between_blocks_responsiveness_first = "Between-block comparison for the first stimulus.",
  between_blocks_habituation_slope = "Between-block comparison for habituation slopes.",
  between_block_10_stimuli = "Block-wise comparison: last 9 stimuli of Block 1 vs. first 9 of Block 2.",
  between_block_last_first = "Comparison between last stimulus of Block 1 and first stimulus of Block 2."
)

# Info block for report header
info <- "
Interpretation for estimate (odds ratio):
---------------------------------------
A ratio of 1.00 → no difference between genotypes.
A ratio < 1.00 → the first genotype (left of “/”) has a lower predicted response.
A ratio > 1.00 → the first genotype has a higher predicted response.

Habituation Rate (Slope):
Negative slope → faster habituation.
"

# Output file path
outfile <- file.path(base_dir, "spaced_glmm_response_prob_comparisons.txt")

# Write report to text file
sink(outfile)

cat("### GLMM Response Probability (Spaced) ###\n\n")
print(summary(m1), row.names = FALSE)
cat("\n")
print(n_per_genotype)
cat("\n\n")

cat("### Estimated Marginal Means Comparisons ###\n\n")
cat(info)

for (nm in names(results)) {
  cat("----", nm, "----\n")
  if (!is.null(descriptions[[nm]])) {
    cat(descriptions[[nm]], "\n\n")
  }
  print(results[[nm]], row.names = FALSE)
  cat("\n\n")
}

sink()
cat("All results written to:", outfile, "\n")
