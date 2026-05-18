###############################################################################
# GLMM Analysis: Massed Training – Distance Moved and Response Delay
# Author: Nils Brehm
# Date: 11/2025
#
# Description:
#   This script fits Generalized Linear Mixed Models (GLMMs) to measure
#   zebrafish behavior across repeated stimuli in massed training experiments.
#   It computes:
#     - Peak movement distance
#     - Summed movement distance
#     - Response delay
#   For each metric, it performs model validation, visualizes habituation
#   curves, extracts estimated marginal means (EMMs), computes slopes,
#   and writes results to disk.
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
base_dir <- "D:/WorkingData/Susana/results/massed"


# ==============================================================================
# 1. Load and Prepare Data
# ==============================================================================
message("Loading and preparing data...")

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
# 2. Exploratory Distributions
# ==============================================================================

# Peak Distance Distribution
h1 <- ggplot(df_final_sub, aes(x = max_peak)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Peak Distance Distribution (Massed)",
    x = "Peak Distance Moved (mm)",
    y = "Count"
  )
h1
ggsave(
  filename = file.path(base_dir, "distance_moved", "massed_peak_distance_dist.pdf"),
  plot = h1, width = 6, height = 4
)

# Summed Distance Distribution
h2 <- ggplot(df_final_sub, aes(x = max_cumsum)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Summed Distance Distribution (Massed)",
    x = "Summed Distance Moved (mm)",
    y = "Count"
  )
h2
ggsave(
  filename = file.path(base_dir, "distance_moved", "massed_summed_distance_dist.pdf"),
  plot = h2, width = 6, height = 4
)


# ==============================================================================
# 3. Fit GLMM Models
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

# --- Model 3: Response Delay (Gaussian LMM) ----------------------------------
m_delay_gaussian <- glmmTMB(
  delay ~ Genotype * stimulus_log * Block + (1 | Video/Well),
  family = gaussian(link = "identity"),
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

validate_model(m_peak, df_final_sub)
validate_model(m_sum, df_final_sub)
validate_model(m_delay, df_final_sub)

# ==============================================================================
# 5. Plot Habituation Curves
# ==============================================================================
message("Plotting habituation curves...")

# --- Peak Movement ------------------------------------------------------------
g_peak <- plot_habituation(df_final_sub, m_peak,
                           label = 'Peak distance moved (mm)',
                           transform = "exp", raw_var = "max_peak")
g_peak
save_plot(g_peak,  file.path(base_dir, "distance_moved", "massed_peak_distance_habituation_curves_v2"), width=10, height=5, dpi=600)



# --- Summed Distance ----------------------------------------------------------
g_sum <- plot_habituation(df_final_sub, m_sum,
                          label = 'Summed distance moved (mm)',
                          transform = "exp", raw_var = "max_cumsum")
g_sum
save_plot(g_sum,  file.path(base_dir, "distance_moved", "massed_summed_distance_habituation_curves_v2"), width=10, height=5, dpi=600)


# --- Response Delay -----------------------------------------------------------
g_delay <- plot_habituation(df_final_sub, m_delay,
                            label = 'Response delay (s)',
                            transform = "none", raw_var = "delay")
g_delay
save_plot(g_delay,  file.path(base_dir, "delay", "massed_response_delay_habituation_curves_v2"), width=10, height=5, dpi=600)

# ==============================================================================
# 6. Estimated Marginal Means (EMMs)
# ==============================================================================

## ---------------------------------------------------------------------------
## 6.1 Peak Distance
## ---------------------------------------------------------------------------
message("Calculating EMMs for peak distance...")

emm_peak <- get_emm_blocks(m_peak, df_final_sub)
emm_peak_slopes <- emtrends(m_peak, ~ Genotype | Block, var = "stimulus_log")
emm_peak_slopes_pairs <- pairs(emm_peak_slopes)

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
  outfile = file.path(base_dir, "distance_moved", "massed_glmm_peak_distance_comparisons.txt"),
  model_name = "GLMM Peak Model (Massed)"
)


## ---------------------------------------------------------------------------
## 6.2 Summed Distance
## ---------------------------------------------------------------------------
message("Calculating EMMs for summed distance...")

emm_sum <- get_emm_blocks(m_sum, df_final_sub)
emm_sum_slopes <- emtrends(m_sum, ~ Genotype | Block, var = "stimulus_log")
emm_sum_slopes_pairs <- pairs(emm_sum_slopes)

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
  outfile = file.path(base_dir, "distance_moved", "massed_glmm_sum_distance_comparisons.txt"),
  model_name = "GLMM Summed Distance Model (Massed)"
)


## ---------------------------------------------------------------------------
## 6.3 Response Delay
## ---------------------------------------------------------------------------
message("Calculating EMMs for response delay...")

emm_delay <- get_emm_blocks(m_delay_gaussian, df_final_sub)
emm_delay_slopes <- emtrends(m_delay_gaussian, ~ Genotype | Block, var = "stimulus_log")
emm_delay_slopes_pairs <- pairs(emm_delay_slopes)

emm_delay_between <- get_emm_consecutive_blocks(m_delay_gaussian, df_final_sub, n_stim = 9)
emm_delay_between_blocks_slopes <- emtrends(m_delay_gaussian, ~ Block | Genotype, var = "stimulus_log")
emm_delay_between_blocks_slopes_pairs <- pairs(emm_delay_between_blocks_slopes)

write_emm_report(
  model = m_delay_gaussian,
  emm_blocks = emm_delay,
  emm_slopes = emm_delay_slopes,
  emm_slopes_pairs = emm_delay_slopes_pairs,
  emm_between = emm_delay_between,
  emm_between_blocks_slopes = emm_delay_between_blocks_slopes,
  emm_between_blocks_slopes_pairs = emm_delay_between_blocks_slopes_pairs,
  outfile = file.path(base_dir, "distance_moved", "massed_glmm_delay_distance_comparisons.txt"),
  model_name = "GLMM Delay Model (Massed)"
)


# ==============================================================================
# 7. Additional Visualization – Delay Violin Plot
# ==============================================================================
message("Generating violin plots for delay distributions...")

block_labels <- df_final_sub %>%
  group_by(Block) %>%
  summarise(n_stim = n_distinct(stimulus)) %>%
  mutate(label = paste0("Block ", Block, ": ", n_stim, " flashes")) %>%
  { setNames(.$label, .$Block) }

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
    title = "Stimulus distribution per delay category and block (Massed)"
  ) +
  theme_pubr(base_size = 14)

response_delay_violin_plot
ggsave(
  filename = file.path(base_dir, "distance_moved", "massed_response_delay_violin_plot.pdf"),
  plot = response_delay_violin_plot, width = 15, height = 15
)


# ==============================================================================
# 8. Completion Message
# ==============================================================================
message("All models fitted, validated, and EMM results exported successfully.")

