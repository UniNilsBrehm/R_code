# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
library(readr)        # Reading CSV files
library(ggplot2)      # Visualization
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(ggpubr)       # Publication-ready plots
library(purrr)

source("C:/UniFreiburg/Code/R_code/susana/Aggregated//DarkFlash_ISI90s_2Blocks/utils.R")
base_dir <- "D:/WorkingData/Susana/Aggregated/DarkFlash_ISI90s_2Blocks"

file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
  # "SPZ_ISI60_ALL_GENOTYPES_Block1_Block2_baseline_subtracted_metrics.csv"
)
# ==============================================================================
# Load Data
# ==============================================================================
message("Loading data...")
res <- load_data_darkflash_60s(file_dir, move_th = 0, , take_peak = 0)

df_final <- res$df_final
df_final_sub <- res$df_final_sub
df_final_sub$delay_ord <- factor(
  df_final_sub$delay,
  ordered = TRUE
)

# Number of animals per genotype
n_per_genotype <- df_final %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

print(n_per_genotype)


save_fig_dir = file.path(base_dir, "figs")
save_results_dir = file.path(base_dir, "results")
models_dir = file.path(base_dir, "models")

# Create directories if they do not exist
dir.create(save_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

geno_colors <- c(
  "ABTL" = "#440154",
  "tyr" = "#3B528B",
  "th2, tyr" = "#21908C",
  "th, tyr" = "#5DC863",
  "th, th2, tyr" = "#FDE725"
)

# ==============================================================================
# Compute Exponential fits on aggreated data variables
# Response probability (binary move)
res_prob <- fit_aggregate_exp(
  df_final, 
  outcome = "response_prob", 
  genotype_order = names(geno_colors),
  colors         = geno_colors,
  y_limits  = c(0, 1),
  y_break   = 0.2
)
res_prob$plot
res_prob$params
ggsave(
  filename = file.path(save_fig_dir,"Aggregated_response_prob_exp_fit.png"),
  plot   = res_prob$plot,
  width  = 14,
  height = 7,
  dpi    = 300,
  bg     = "white"
)


# Peak distance
res_peak <- fit_aggregate_exp(
  df_final_sub, 
  outcome = "distance",
  genotype_order = names(geno_colors),
  colors         = geno_colors,
  y_limits  = c(0, 10),
  y_break   = 2
)
res_peak$plot
res_peak$params
ggsave(
  filename = file.path(save_fig_dir,"Aggregated_peak_exp_fit.png"),
  plot   = res_peak$plot,
  width  = 14,
  height = 7,
  dpi    = 300,
  bg     = "white"
)

# Cumulative distance
res_summed <- fit_aggregate_exp(
  df_final_sub, 
  outcome = "cumsum_distance",
  genotype_order = names(geno_colors),
  colors         = geno_colors,
  y_limits  = c(0, 20),
  y_break   = 5
)
res_summed$plot
res_summed$params
ggsave(
  filename = file.path(save_fig_dir,"Aggregated_summed_exp_fit.png"),
  plot   = res_summed$plot,
  width  = 14,
  height = 7,
  dpi    = 300,
  bg     = "white"
)

# Delay (ordinal 0–4)
res_delay <- fit_aggregate_exp(
  df_final_sub, 
  outcome = "delay",
  genotype_order = names(geno_colors),
  colors         = geno_colors,
  y_limits  = c(0, 3),
  y_break   = 1
)
res_delay$plot
res_delay$params
ggsave(
  filename = file.path(save_fig_dir,"Aggregated_delay_exp_fit.png"),
  plot   = res_delay$plot,
  width  = 14,
  height = 7,
  dpi    = 300,
  bg     = "white"
)

message("========== FINISHED ==========")

# res <- bootstrap_exp_ci(df_final, outcome = "response_prob", n_boot = 1000)
# res$plot
# res$params
# res$boot_ci  # raw CI bounds per stimulus × Genotype × Block
