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
library(tidybayes)    # Posterior draws and summaries
library(ggplot2)      # Visualization
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(ggpubr)       # Publication-ready plots

# Load helper functions
source("C:/UniFreiburg/Code/R_code/susana/utils.R")

# Base directory for saving results
base_dir <- "D:/WorkingData/Susana/results/darkflash_60s/"
dir.create(file.path(base_dir, "nlme"), recursive = TRUE, showWarnings = FALSE)

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


# ==============================================================================
# 2. Define Priors
# ==============================================================================
# Visualize what different parameter values imply before committing to priors
t_seq <- 0:59
expand.grid(
  t      = t_seq,
  alpha  = c(0.2, 0.4, 0.6),
  lambda = c(0.02, 0.08, 0.20)
) %>%
  mutate(prob = alpha + (1 - alpha) * exp(-lambda * t)) %>%
  ggplot(aes(t + 1, prob,
             color = factor(lambda),
             linetype = factor(alpha))) +
  geom_line() +
  labs(
    title  = "Prior predictive check: implied curves",
    x      = "Stimulus number",
    y      = "Response probability",
    color  = "Rate (lambda)",
    linetype = "Asymptote (alpha)"
  ) +
  theme_minimal()

# Priors:
#   alpha:     logit scale, logit(0.3) ≈ -0.85, logit(0.5) ≈ 0.00
#   loglambda: log scale, log(0.08) ≈ -2.5
priors <- c(
  prior(normal(-0.5, 1.0), nlpar = "alpha"),
  prior(normal(-2.5, 0.8), nlpar = "loglambda"),
  prior(exponential(1),    class = "sd", nlpar = "alpha"),
  prior(exponential(1),    class = "sd", nlpar = "loglambda")
)


# ==============================================================================
# 3. Define Model Formula
# ==============================================================================
# P(respond) = alpha + (1 - alpha) * exp(-exp(loglambda) * t)
#
# alpha is on the logit scale (brms applies logit link automatically
# via the bernoulli family) — constrains asymptote to [0, 1]
#
# loglambda uses exp() in the equation to keep rate positive,
# avoiding the need for explicit constraints

nlform <- bf(
  move ~ alpha + (1 - alpha) * exp(-exp(loglambda) * t),
  alpha     ~ Genotype * Block + (1 | Video/Well),
  loglambda ~ Genotype * Block + (1 | Video/Well),
  nl = TRUE
)


# ==============================================================================
# 4. Full Model Fit
# ==============================================================================
# Only run after test fit passes all checks above

message("Fitting full model...")

m_nlme <- brm(
  formula = nlform,
  data    = df_model,
  family  = bernoulli(link = "logit"),
  prior   = priors,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 2000,
  control = list(
    adapt_delta   = 0.95,
    max_treedepth = 12
  ),
  seed = 42,
  file = file.path(base_dir, "nlme", "m_nlme_habituation")
)

# Repeat convergence checks on full model
summary(m_nlme)

nuts_params(m_nlme) %>%
  filter(Parameter == "divergent__") %>%
  summarise(divergences = sum(Value))

pp_check(m_nlme, ndraws = 100)


# ==============================================================================
# 5. Extract Posterior Parameters
# ==============================================================================
# First inspect actual parameter names from the model
# before writing spread_draws calls
# parnames(m_nlme)
variables(m_nlme)

# --- Asymptote (alpha) per genotype x block -----------------------------------
# Back-transform from logit scale with plogis()
alpha_draws <- m_nlme %>%
  spread_draws(
    b_alpha_Intercept,
    `b_alpha_GenotypeXth.th2.tyr`,   # replace with actual names from parnames()
    `b_alpha_GenotypeXth.tyr`,
    `b_alpha_GenotypeXth2.tyr`,
    `b_alpha_GenotypeXtyr`,
    b_alpha_BlockBlock2
  ) %>%
  mutate(
    alpha_ABTL_B1       = plogis(b_alpha_Intercept),
    alpha_ABTL_B2       = plogis(b_alpha_Intercept + b_alpha_BlockBlock2)
    # Add other genotype combinations as needed after checking parnames()
  )

# --- Rate (loglambda) per genotype x block ------------------------------------
# Back-transform from log scale with exp()
lambda_draws <- m_nlme %>%
  spread_draws(
    b_loglambda_Intercept,
    b_loglambda_BlockBlock2
    # Add genotype terms after checking parnames()
  ) %>%
  mutate(
    lambda_ABTL_B1 = exp(b_loglambda_Intercept),
    lambda_ABTL_B2 = exp(b_loglambda_Intercept + b_loglambda_BlockBlock2)
  )

# --- Summary table of posterior medians and 95% credible intervals ------------
m_nlme %>%
  gather_draws(b_alpha_Intercept,
               b_loglambda_Intercept) %>%
  median_qi(.width = 0.95)


# ==============================================================================
# 6. Posterior Contrasts Between Genotypes
# ==============================================================================
# Compare genotypes on habituation rate (loglambda) within each block
# Positive contrast = first genotype has higher log-rate = habituates faster

# Extract all loglambda population-level draws and compute contrasts manually
lambda_genotype_draws <- m_nlme %>%
  spread_draws(
    b_loglambda_Intercept,
    b_loglambda_BlockBlock2,
    `b_loglambda_GenotypeXth.th2.tyr`,  # replace with actual names
    `b_loglambda_GenotypeXth.tyr`,
    `b_loglambda_GenotypeXth2.tyr`,
    `b_loglambda_GenotypeXtyr`
  ) %>%
  mutate(
    # Block 1 log-rates per genotype
    loglambda_ABTL_B1        = b_loglambda_Intercept,
    loglambda_th.th2.tyr_B1  = b_loglambda_Intercept + `b_loglambda_GenotypeXth.th2.tyr`,
    loglambda_th.tyr_B1      = b_loglambda_Intercept + `b_loglambda_GenotypeXth.tyr`,
    loglambda_th2.tyr_B1     = b_loglambda_Intercept + `b_loglambda_GenotypeXth2.tyr`,
    loglambda_tyr_B1         = b_loglambda_Intercept + `b_loglambda_GenotypeXtyr`,
    
    # Block 2 log-rates per genotype (add Block2 main effect +
    # genotype x block interaction if present — check parnames())
    loglambda_ABTL_B2        = b_loglambda_Intercept + b_loglambda_BlockBlock2,
    
    # Contrasts: ABTL vs each genotype in Block 1
    # Negative = mutant habituates faster than ABTL
    contrast_ABTL_vs_th.th2.tyr_B1 = loglambda_ABTL_B1 - loglambda_th.th2.tyr_B1,
    contrast_ABTL_vs_th.tyr_B1     = loglambda_ABTL_B1 - loglambda_th.tyr_B1,
    contrast_ABTL_vs_th2.tyr_B1    = loglambda_ABTL_B1 - loglambda_th2.tyr_B1,
    contrast_ABTL_vs_tyr_B1        = loglambda_ABTL_B1 - loglambda_tyr_B1
  )

# Summarise contrasts with posterior median and 95% credible interval
lambda_genotype_draws %>%
  select(starts_with("contrast_")) %>%
  pivot_longer(everything(), names_to = "contrast", values_to = "estimate") %>%
  group_by(contrast) %>%
  median_qi(estimate, .width = 0.95) %>%
  # Flag contrasts where 95% CI excludes zero
  mutate(significant = (.lower > 0 | .upper < 0))


# ==============================================================================
# 7. Visualize Fitted Curves
# ==============================================================================
# Population-level predictions (re_formula = NA ignores random effects)
newdata <- expand.grid(
  t        = 0:59,                        # t=0 → stimulus 1, t=59 → stimulus 60
  Genotype = unique(df_model$Genotype),
  Block    = unique(df_model$Block)
)

fitted_curves <- newdata %>%
  add_epred_draws(m_nlme, ndraws = 200, re_formula = NA)

# Observed proportions per stimulus for overlay
obs_props <- df_model %>%
  group_by(Genotype, Block, Stimulus_New) %>%
  summarise(obs_prob = mean(move), .groups = "drop")

# Plot
g_nlme <- fitted_curves %>%
  ggplot(aes(x = t + 1, y = .epred,
             color = Genotype, fill = Genotype)) +
  stat_lineribbon(.width = 0.95, alpha = 0.2) +
  geom_point(
    data    = obs_props,
    mapping = aes(x = Stimulus_New, y = obs_prob, color = Genotype),
    size    = 0.8,
    alpha   = 0.5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Block) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x     = "Stimulus number",
    y     = "Response probability",
    title = "Habituation curves — Bayesian nonlinear mixed model",
    color = "Genotype",
    fill  = "Genotype"
  ) +
  theme_minimal()

g_nlme

ggsave(
  filename = file.path(base_dir, "nlme", "nlme_habituation_curves.pdf"),
  plot     = g_nlme,
  width    = 10,
  height   = 5,
  dpi      = 600
)


# ==============================================================================
# 8. Completion Message
# ==============================================================================
message("All analyses completed successfully — results saved to disk.")