###############################################################################
# GLMM Analysis: Dark Flash Blocks
# for response probability, response strength and response delay
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

source("C:/UniFreiburg/Code/R_code/susana/utils.R")

base_dir <- "D:/WorkingData/Susana/DarkFlash_ISI90s_2Blocks"
file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

# ==============================================================================
# Load Data
# ==============================================================================
message("Loading data...")
res <- load_data_darkflash_60s(file_dir, move_th = 0.2, , take_peak = 0)

df_final <- res$df_final
df_final_sub <- res$df_final_sub

# Number of animals per genotype
n_per_genotype <- df_final %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

print(n_per_genotype)


# ==============================================================================
# Explore Distributions (Histograms)
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
  filename = file.path(base_dir, "figs", "peak_distance_dist.png"),
  plot = h1,
  width = 6,
  height = 4,
  dpi = 300
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
  filename = file.path(base_dir, "figs", "summed_distance_dist.pdf"),
  plot = h2,
  width = 6,
  height = 4,
  dpi = 300
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
  filename = file.path(base_dir, "figs", "delay_dist.pdf"),
  plot = h3,
  width = 6,
  height = 4,
  dpi = 300
)


# ==============================================================================
# Fit Models
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
df_final_sub$delay_non_zero <- df_final_sub$delay + 0.001
m_delay <- glmmTMB(
  delay_non_zero ~ Genotype * stimulus_log * Block + (1 | Video/Well),
  family = gaussian(link = "identity"),
  data = df_final_sub
)

# --- Model 4: Response Prob (Binomial GLMM) -----------------------------------
m_prob <- glmmTMB(
  move ~ Genotype * Block * stimulus_log + (1 | Video/Well),
  family = binomial(link = "logit"),
  data = df_final
)

# ==============================================================================
# Model Validation
# ==============================================================================
message("Validating models...")

# Residual check
model_residuals_check(m_peak, df_final_sub)
model_residuals_check(m_sum, df_final_sub)
model_residuals_check(m_delay, df_final_sub)
model_residuals_check(m_prob, df_final)

# Full validation
validate_model(m_peak, df_final_sub)
validate_model(m_sum, df_final_sub)
validate_model(m_delay, df_final_sub)
validate_model(m_prob, df_final)

# ==============================================================================
# COMPARISON TESTS
# ==============================================================================
## ---------------------------------------------------------------------------
## Response Prob.
## ---------------------------------------------------------------------------
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
## Peak Distance
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
## Summed Distance
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
## Response Delay
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

# ==============================================================================
# Plot Habituation Curves
# ==============================================================================
message("Plotting habituation curves...")

# Peak movement
g_peak <- plot_habituation(df_final_sub, m_peak,
                           label = 'Peak distance moved (mm)',
                           transform = "exp", raw_var = "max_peak")
g_peak
save_plot(g_peak,  file.path(base_dir, "figs", "GLMM_peak_distance_habituation_curves"), width=10, height=5, dpi=600)

# Summed distance
g_sum <- plot_habituation(df_final_sub, m_sum,
                          label = 'Summed distance moved (mm)',
                          transform = "exp", raw_var = "max_cumsum")
g_sum
save_plot(g_sum,  file.path(base_dir, "figs", "GLMM_summed_distance_habituation_curves"), width=10, height=5, dpi=600)


# Response delay
g_delay <- plot_habituation(df_final_sub, m_delay,
                            label = 'Response delay (s)',
                            transform = "none", raw_var = "delay")
g_delay
save_plot(g_delay,  file.path(base_dir, "figs", "GLMM_delay_habituation_curves"), width=10, height=5, dpi=600)

# Response prob
g_prob <- plot_habituation(df_final, m_prob, label='Response prob.', transform = "plogis")
g_prob
save_plot(g_prob, file.path(base_dir, "figs", "GLMM_response_prob_habituation_curves"), width=10, height=5, dpi=600)

# ==============================================================================
# Plot Habituation Curves Separated
# ==============================================================================
# Peak movement
p_peak_by_genotype_block <- plot_habituation_glmm_by_genotype_block(
  df_final = df_final_sub,
  model = m_peak,
  label = "Peak distance moved (mm)",
  raw_var = "max_peak",
  transform = "exp"
)

print(p_peak_by_genotype_block)

ggsave(
  filename = file.path(base_dir, "figs", "GLMM_peak_distance_habituation_each_genotype_each_block.png"),
  plot = p_peak_by_genotype_block,
  width = 14,
  height = 7,
  dpi = 300
)

# Summed movement
p_summed_by_genotype_block <- plot_habituation_glmm_by_genotype_block(
  df_final = df_final_sub,
  model = m_sum,
  label = "Summed distance moved (mm)",
  raw_var = "max_cumsum",
  transform = "exp"
)

print(p_summed_by_genotype_block)

ggsave(
  filename = file.path(base_dir, "figs", "GLMM_summed_distance_habituation_each_genotype_each_block.png"),
  plot = p_summed_by_genotype_block,
  width = 14,
  height = 7,
  dpi = 300
)

# Response delay
p_delay_by_genotype_block <- plot_habituation_glmm_by_genotype_block(
  df_final = df_final_sub,
  model = m_delay,
  label = "Response delay (s)",
  raw_var = "delay",
  transform = "none"
)

print(p_delay_by_genotype_block)

ggsave(
  filename = file.path(base_dir, "figs", "GLMM_delay_habituation_each_genotype_each_block.png"),
  plot = p_delay_by_genotype_block,
  width = 14,
  height = 7,
  dpi = 300
)

# Response prob.
p_response_prob_by_genotype_block <- plot_habituation_glmm_by_genotype_block(
  df_final = df_final,
  model = m_prob,
  label = "response prob",
  raw_var = "move",
  transform = "plogis"
)

print(p_response_prob_by_genotype_block)

ggsave(
  filename = file.path(base_dir, "figs", "GLMM_response_prob_habituation_each_genotype_each_block.png"),
  plot = p_response_prob_by_genotype_block,
  width = 14,
  height = 7,
  dpi = 300
)


# ==============================================================================
# Simple aggregate response-probability exponential fits
# ==============================================================================

library(purrr)
library(broom)

# ------------------------------------------------------------------------------
# 1. Compute response probability per Genotype × Block × stimulus
# ------------------------------------------------------------------------------

df_prob_agg <- df_final %>%
  mutate(
    stimulus = as.numeric(stimulus),
    stimulus0 = stimulus - 1,
    Genotype = factor(Genotype),
    Block = factor(Block),
    move = as.integer(move)
  ) %>%
  group_by(Genotype, Block, stimulus, stimulus0) %>%
  summarise(
    n_animals = n(),
    n_response = sum(move, na.rm = TRUE),
    response_prob = mean(move, na.rm = TRUE),
    .groups = "drop"
  )

print(df_prob_agg)

# ------------------------------------------------------------------------------
# 2. Fit simple exponential per Genotype × Block
# response_prob = Asym + (R0 - Asym) * exp(-k * stimulus0)
# ------------------------------------------------------------------------------

fit_simple_exp <- function(dat) {
  tryCatch(
    {
      nls(
        response_prob ~ Asym + (R0 - Asym) * exp(-k * stimulus0),
        data = dat,
        weights = n_animals,
        start = list(
          Asym = max(0.01, min(dat$response_prob, na.rm = TRUE)),
          R0   = min(0.99, max(dat$response_prob, na.rm = TRUE)),
          k    = 0.1
        ),
        algorithm = "port",
        lower = c(
          Asym = 0.001,
          R0   = 0.001,
          k    = 0.0001
        ),
        upper = c(
          Asym = 0.999,
          R0   = 0.999,
          k    = 10
        ),
        control = nls.control(
          maxiter = 500,
          warnOnly = TRUE
        )
      )
    },
    error = function(e) {
      message(
        "NLS failed for ",
        unique(dat$Genotype), " ",
        unique(dat$Block), ": ",
        e$message
      )
      return(NULL)
    }
  )
}

simple_fits <- df_prob_agg %>%
  group_by(Genotype, Block) %>%
  group_split() %>%
  map(fit_simple_exp)

fit_keys <- df_prob_agg %>%
  distinct(Genotype, Block) %>%
  arrange(Genotype, Block)

# ------------------------------------------------------------------------------
# 3. Extract fitted parameters
# ------------------------------------------------------------------------------

simple_exp_params <- map2_dfr(
  simple_fits,
  seq_along(simple_fits),
  function(mod, i) {
    if (is.null(mod)) {
      tibble(
        Genotype = fit_keys$Genotype[i],
        Block = fit_keys$Block[i],
        Asym = NA_real_,
        R0 = NA_real_,
        k = NA_real_
      )
    } else {
      cc <- coef(mod)
      
      tibble(
        Genotype = fit_keys$Genotype[i],
        Block = fit_keys$Block[i],
        Asym = unname(cc["Asym"]),
        R0 = unname(cc["R0"]),
        k = unname(cc["k"])
      )
    }
  }
) %>%
  mutate(
    half_life_stimuli = log(2) / k
  )

print(simple_exp_params)

write.csv(
  simple_exp_params,
  file = file.path(base_dir, "results", "simple_aggregate_exp_response_prob_parameters.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Generate smooth fitted curves
# ------------------------------------------------------------------------------

simple_exp_pred <- df_prob_agg %>%
  group_by(Genotype, Block) %>%
  summarise(
    stim_min = min(stimulus, na.rm = TRUE),
    stim_max = max(stimulus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Genotype, Block) %>%
  reframe(
    stimulus = seq(stim_min, stim_max, length.out = 100)
  ) %>%
  mutate(
    stimulus0 = stimulus - 1
  ) %>%
  left_join(simple_exp_params, by = c("Genotype", "Block")) %>%
  mutate(
    fit = Asym + (R0 - Asym) * exp(-k * stimulus0)
  )

# ------------------------------------------------------------------------------
# 5. Plot simple aggregate exponential curves
# ------------------------------------------------------------------------------

p_simple_exp <- ggplot(simple_exp_pred, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  
  geom_point(
    data = df_prob_agg,
    aes(x = stimulus, y = response_prob, color = Genotype),
    alpha = 0.45,
    size = 1,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    aes(y = fit),
    linewidth = 1.3
  ) +
  
  theme_pubr(base_size = 14) +
  labs(
    x = "Stimulus number within block",
    y = "Response probability",
    color = "Genotype",
    fill = "Genotype",
    title = "Simple aggregate exponential fits"
  ) +
  theme(
    legend.position = "top",
    panel.spacing = unit(1.2, "lines")
  )

print(p_simple_exp)

ggsave(
  filename = file.path(base_dir, "figs", "Simple_aggregate_exp_response_prob_curves.png"),
  plot = p_simple_exp,
  width = 14,
  height = 7,
  dpi = 300
)

