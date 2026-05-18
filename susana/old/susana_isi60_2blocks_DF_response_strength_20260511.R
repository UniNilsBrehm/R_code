###############################################################################
# GLMM Analysis: Dark Flash Blocks – Response Strength
# Author: Nils Brehm
# Date: 05/2026
#
# Description:
#   This script analyzes habituation behavior of larval zebrafish to dark flash 
#   experiments. It fits a GLMM model predicting the response strength across 
#   stimulus blocks, validates the model, visualizes habituation curves, and 
#   computes estimated marginal means (EMMs) and contrasts between genotypes 
#   and blocks.
#
# Experimental Design:
#   In each block dark flash (DF: brief period of darkness) is presented every
#   60 seconds. There are two blocks with a inter-block pause of 1 hour. Each
#   block gas 60 DF stimuli. The analysis is based on the "distance moved" in
#   response to each DF.
#
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

# Custom utility functions: validation, plotting, EMM helpers
source("C:/UniFreiburg/Code/R_code/susana/utils.R")

# Base directory for results
base_dir <- "D:/WorkingData/Susana/results/darkflash_60s/"


# ==============================================================================
# 1. Load Data
# ==============================================================================
message("Loading data...")
# Optionally subset to certain genotypes
# Load data from csv file, ignore warning
file_dir <-  "D:/WorkingData/Susana/SPZ_ISI60_ALL_GENOTYPES_Block1_Block2_baseline_subtracted_metrics.csv"
# df <- read_csv(file_dir)

res <- load_data_darkflash_60s(file_dir, move_th = 0.2, , take_peak = 0)

df_final <- res$df_final
df_final_sub <- res$df_final_sub

# Number of animals per genotype
n_per_genotype <- df_final %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

print(n_per_genotype)


# ==============================================================================
# 2. Explore Distributions (Histograms)
# ==============================================================================

# Peak Distance Distribution
h1 <- ggplot(df_final_sub, aes(x = max_peak)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Peak Distance Distribution",
    x = "Peak Distance Moved (mm)",
    y = "Count"
  )
h1
ggsave(
  filename = file.path(base_dir, "distance_moved", "peak_distance_dist.pdf"),
  plot = h1, width = 6, height = 4
)

# Summed Distance Distribution
h2 <- ggplot(df_final_sub, aes(x = max_cumsum)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Summed Distance Distribution",
    x = "Summed Distance Moved (mm)",
    y = "Count"
  )
h2
ggsave(
  filename = file.path(base_dir, "distance_moved", "summed_distance_dist.pdf"),
  plot = h2, width = 6, height = 4
)

# Summed Distance Distribution
h3 <- ggplot(df_final_sub, aes(x = delay)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Delay Distribution",
    x = "Delay (s)",
    y = "Count"
  )
h3
ggsave(
  filename = file.path(base_dir, "distance_moved", "delay_dist.pdf"),
  plot = h3, width = 6, height = 4
)


# ==============================================================================
# 3. Fit Models
# ==============================================================================
message("Fitting GLMM models...")
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
  data = df_final_sub
)

# --- Model 3: Response Delay (Gaussian GLMM) ----------------------------------
m_delay <- glmmTMB(
  delay ~ Genotype * stimulus_log * Block + (1 | Video/Well),
  family = gaussian(link = "identity"),
  data = df_final_sub
)

# ==============================================================================
# 4. Model Validation
# ==============================================================================
message("Validating models...")

# Residual check
model_residuals_check(m_peak, df_final_sub)
model_residuals_check(m_sum, df_final_sub)
model_residuals_check(m_delay, df_final_sub)

# Full validation
validate_model(m_peak, df_final_sub)
validate_model(m_sum, df_final_sub)
validate_model(m_delay, df_final_sub)


# ==============================================================================
# 5. Plot Habituation Curves
# ==============================================================================
message("Plotting habituation curves...")

# Peak movement
g_peak <- plot_habituation(df_final_sub, m_peak,
                           label = 'Peak distance moved (mm)',
                           transform = "exp", raw_var = "max_peak")
g_peak
save_plot(g_peak,  file.path(base_dir, "distance_moved", "peak_distance_habituation_curves_v2"), width=10, height=5, dpi=600)

# Summed distance
g_sum <- plot_habituation(df_final_sub, m_sum,
                          label = 'Summed distance moved (mm)',
                          transform = "exp", raw_var = "max_cumsum")
g_sum
save_plot(g_sum,  file.path(base_dir, "distance_moved", "summed_distance_habituation_curves_v2"), width=10, height=5, dpi=600)


# Response delay
g_delay <- plot_habituation(df_final_sub, m_delay,
                            label = 'Response delay (s)',
                            transform = "none", raw_var = "delay")
g_delay
save_plot(g_delay,  file.path(base_dir, "delay", "delay_habituation_curves_v2"), width=10, height=5, dpi=600)


# ==============================================================================
# TEST
# ==============================================================================

# --- Block-level EMMs (dynamic across all blocks) ------------------------------
emm_prob <- get_emm_blocks(m_prob, df_final)

# --- Habituation slopes (within blocks) --------------------------------------
emm_prob_slopes <- emtrends(m_prob, ~ Genotype | Block, var = "stimulus_log")
emm_prob_slopes_pairs <- pairs(emm_prob_slopes)

# --- Between-block comparisons (1 vs 2, 2 vs 3, etc.) -------------------------
emm_prob_between <- get_emm_consecutive_blocks(m_prob, df_final, n_stim = 9)

# Compare the n first an n first of blocks
emm_prob_between_first <- get_emm_consecutive_blocks_first(m_prob, df_final, n_stim = 9)

# --- Between-block slopes -----------------------------------------------------
emm_prob_between_blocks_slopes <- emtrends(m_prob, ~ Block | Genotype, var = "stimulus_log")
emm_prob_between_blocks_slopes_pairs <- pairs(emm_prob_between_blocks_slopes)


write_emm_report(
  model = m_prob,
  emm_blocks = emm_prob,
  emm_slopes = emm_prob_slopes,
  emm_slopes_pairs = emm_prob_slopes_pairs,
  emm_between = emm_prob_between,
  emm_between_first = emm_prob_between_first,
  emm_between_blocks_slopes = emm_prob_between_blocks_slopes,
  emm_between_blocks_slopes_pairs = emm_prob_between_blocks_slopes_pairs,
  outfile = file.path(base_dir, "response_prob", "spaced_glmm_response_prob_comparisons.txt"),
  model_name = "GLMM Response Probability Model"
)


## ---------------------------------------------------------------------------
## 6.1. Peak Distance
## ---------------------------------------------------------------------------
message("Calculating EMMs for peak distance...")

emm_peak <- get_emm_blocks(m_peak, df_final_sub)
emm_peak_slopes <- emtrends(m_peak, ~ Genotype | Block, var = "stimulus_log")
emm_peak_slopes_pairs <- pairs(emm_peak_slopes)

# Between-block comparisons
emm_peak_between <- get_emm_consecutive_blocks(m_peak, df_final_sub, n_stim = 9)
emm_peak_between_blocks_slopes <- emtrends(m_peak, ~ Block | Genotype, var = "stimulus_log")
emm_peak_between_blocks_slopes_pairs <- pairs(emm_peak_between_blocks_slopes)

write_emm_report(
  model = m_peak,
  emm_blocks = emm_peak,
  emm_slopes = emm_peak_slopes,
  emm_slopes_pairs = emm_peak_slopes_pairs,
  emm_between = emm_peak_between,
  emm_between_blocks_slopes = emm_peak_between_blocks_slopes,
  emm_between_blocks_slopes_pairs = emm_peak_between_blocks_slopes_pairs,
  outfile = file.path(base_dir, "distance_moved", "glmm_peak_distance_comparisons.txt"),
  model_name = "GLMM Peak Model"
)


## ---------------------------------------------------------------------------
## 6.2. Summed Distance
## ---------------------------------------------------------------------------
message("Calculating EMMs for summed distance...")

emm_sum <- get_emm_blocks(m_sum, df_final_sub)
emm_sum_slopes <- emtrends(m_sum, ~ Genotype | Block, var = "stimulus_log")
emm_sum_slopes_pairs <- pairs(emm_sum_slopes)

# Between-block comparisons
emm_sum_between <- get_emm_consecutive_blocks(m_sum, df_final_sub, n_stim = 9)
emm_sum_between_blocks_slopes <- emtrends(m_sum, ~ Block | Genotype, var = "stimulus_log")
emm_sum_between_blocks_slopes_pairs <- pairs(emm_sum_between_blocks_slopes)

write_emm_report(
  model = m_sum,
  emm_blocks = emm_sum,
  emm_slopes = emm_sum_slopes,
  emm_slopes_pairs = emm_sum_slopes_pairs,
  emm_between = emm_sum_between,
  emm_between_blocks_slopes = emm_sum_between_blocks_slopes,
  emm_between_blocks_slopes_pairs = emm_sum_between_blocks_slopes_pairs,
  outfile = file.path(base_dir, "distance_moved", "glmm_sum_distance_comparisons.txt"),
  model_name = "GLMM Sum Model"
)


## ---------------------------------------------------------------------------
## 6.3. Response Delay
## ---------------------------------------------------------------------------
message("Calculating EMMs for response delay...")

emm_delay <- get_emm_blocks(m_delay_poisson, df_final_sub)
emm_delay_slopes <- emtrends(m_delay_poisson, ~ Genotype | Block, var = "stimulus_log")
emm_delay_slopes_pairs <- pairs(emm_delay_slopes)

# Between-block comparisons
emm_delay_between <- get_emm_consecutive_blocks(m_delay_poisson, df_final_sub, n_stim = 9)
emm_delay_between_blocks_slopes <- emtrends(m_delay_poisson, ~ Block | Genotype, var = "stimulus_log")
emm_delay_between_blocks_slopes_pairs <- pairs(emm_delay_between_blocks_slopes)

write_emm_report(
  model = m_delay_poisson,
  emm_blocks = emm_delay,
  emm_slopes = emm_delay_slopes,
  emm_slopes_pairs = emm_delay_slopes_pairs,
  emm_between = emm_delay_between,
  emm_between_blocks_slopes = emm_delay_between_blocks_slopes,
  emm_between_blocks_slopes_pairs = emm_delay_between_blocks_slopes_pairs,
  outfile = file.path(base_dir, "distance_moved", "glmm_delay_distance_comparisons.txt"),
  model_name = "GLMM Delay Model"
)

# ------------------------------------------------------------------------------
# Delays Violin Plot (Dynamic by Block)
# ------------------------------------------------------------------------------

# Create dynamic facet labels for each Block
block_labels <- df_final_sub %>%
  group_by(Block) %>%
  summarise(n_stim = n_distinct(stimulus)) %>%
  mutate(label = paste0("Block ", Block, ": ", n_stim, " flashes")) %>%
  { setNames(.$label, .$Block) }  # Convert to named vector for labeller()

# Plot
response_delay_violin_plot <- ggplot(df_final_sub, aes(x = as.factor(delay), y = stimulus, fill = as.factor(delay))) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "grey30") +
  facet_wrap(
    ~Block, ncol = 2, scales = "free_y",
    labeller = as_labeller(block_labels)
  ) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1, name = "Delay (s)") +
  labs(
    x = "Delay category (s)",
    y = "Stimulus number",
    title = "Stimulus distribution per delay category and block"
  ) +
  theme_pubr(base_size = 14)
response_delay_violin_plot
ggsave(
  filename = file.path(base_dir, "distance_moved", "response_delay_violin_plot.pdf"),
  plot = response_delay_violin_plot, width = 15, height = 15
)

# ==============================================================================
# Done!
# ==============================================================================
message("All models fitted, validated, and EMM results exported successfully.")

