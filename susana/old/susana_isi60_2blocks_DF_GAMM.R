###############################################################################
# NLME Analysis: Dark Flash Blocks – Response Probability
# Author: Nils Brehm
# Date: 05/2026
#
# Description:
#   This script analyzes habituation behavior of larval zebrafish to dark flash
#   experiments. It fits a Bayesian nonlinear mixed model predicting the
#   probability of movement (response probability) across stimulus blocks,
#   validates the model, visualizes habituation curves, and computes posterior
#   contrasts between genotypes and blocks.
#
# Experimental Design:
#   In each block a dark flash (DF: brief period of darkness) is presented every
#   60 seconds. There are two blocks with an inter-block pause of 1 hour. Each
#   block has 60 DF stimuli. The analysis is based on the "distance moved" in
#   response to each DF.
#
# Model:
#   P(respond) = alpha + (1 - alpha) * exp(-exp(loglambda) * t)
#   where:
#     alpha     = asymptote (floor response probability)
#     loglambda = log of habituation rate
#     t         = stimulus number - 1 (so t=0 at stimulus 1, guaranteeing
#                 P=1.0 at first stimulus by construction)
#
###############################################################################


# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
library(readr)        # Reading CSV files
library(brms)         # Bayesian nonlinear mixed models
library(ggplot2)      # Visualization
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(ggpubr)       # Publication-ready plots
library(mgcv)      # GAMMs — already battle-tested, no installation needed
library(gratia)    # Modern ggplot2-based visualization for mgcv models

# Load helper functions
source("C:/UniFreiburg/Code/R_code/susana/utils.R")

# Base directory for saving results
base_dir <- "D:/WorkingData/Susana/results/darkflash_60s/"

# ==============================================================================
# 1. Load and Prepare Data
# ==============================================================================
file_dir <- "D:/WorkingData/Susana/SPZ_ISI60_removed_non_responders_2stimuli.csv"

res <- load_data_darkflash_60s(file_dir, move_th = 0, take_peak = 0)

df_final     <- res$df_final
df_final_sub <- res$df_final_sub

# Number of animals per genotype
n_per_genotype <- df_final %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

print(n_per_genotype)

# Prepare model data:
# t = Stimulus_New - 1 so that t=0 at stimulus 1
# This guarantees P(respond) = 1.0 at t=0 by construction,
# consistent with the selection criterion (only fish that responded
# to stimulus 1 are included)
df_model <- df_final %>%
  mutate(t = Stimulus_New - 1)

