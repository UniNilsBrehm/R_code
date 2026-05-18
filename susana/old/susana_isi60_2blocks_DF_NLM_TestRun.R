###############################################################################
# NLME Analysis: Dark Flash Blocks – Response Probability (TEST VERSION)
# Author: Nils Brehm
# Date: 05/2026
#
# Description:
#   Fast test version for model validation only.
#   Simplified formula (no interactions, no random slopes) and
#   subsampled data to keep runtime under 5 minutes.
#   Do NOT use for inference — run full script once test passes.
###############################################################################


# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
library(readr)
library(brms)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

source("C:/UniFreiburg/Code/R_code/susana/utils.R")

base_dir <- "D:/WorkingData/Susana/results/darkflash_60s/"


# ==============================================================================
# 1. Load and Prepare Data
# ==============================================================================
file_dir <- "D:/WorkingData/Susana/SPZ_ISI60_removed_non_responders_2stimuli.csv"

res      <- load_data_darkflash_60s(file_dir, move_th = 0, take_peak = 0)
df_final <- res$df_final

df_model <- df_final %>%
  mutate(t = Stimulus_New - 1)


# ==============================================================================
# 2. Priors (same as full model)
# ==============================================================================
# Priors for TEST model only
# sd for loglambda removed because loglambda has no random effect
# in the simplified test formula
priors_test <- c(
  prior(normal(-0.5, 1.0), nlpar = "alpha"),
  prior(normal(-2.5, 0.8), nlpar = "loglambda"),
  prior(exponential(1),    class = "sd", nlpar = "alpha")
  # no sd prior for loglambda — no (1 | Video/Well) on loglambda in nlform_test
)

# ==============================================================================
# 3. Simplified Formula for Testing
# ==============================================================================
# Differences from full model:
#   - No Genotype:Block interaction (main effects only)
#   - Random intercept on alpha only — loglambda random effect dropped
#   - This is the minimal structure that still validates the nonlinear core

nlform_test <- bf(
  move ~ alpha + (1 - alpha) * exp(-exp(loglambda) * t),
  alpha     ~ Genotype + Block + (1 | Video/Well),  # main effects only
  loglambda ~ Genotype + Block,                     # no random effects
  nl = TRUE
)


# ==============================================================================
# 4. Test Fit
# ==============================================================================
message("Fitting TEST model — expected runtime: 2-5 minutes...")
message("Started at: ", Sys.time())

m_test <- brm(
  formula = nlform_test,
  data    = df_model %>% slice_sample(n = 1000),
  family  = bernoulli(link = "logit"),
  prior   = priors_test,   # <-- changed
  chains  = 2,
  cores   = 2,
  iter    = 800,
  warmup  = 400,
  control = list(
    adapt_delta   = 0.90,
    max_treedepth = 10
  ),
  seed = 42
)

message("Completed at: ", Sys.time())


# ==============================================================================
# 5. Quick Diagnostics
# ==============================================================================

# --- Convergence --------------------------------------------------------------
summary(m_test)
# Look for: all Rhat < 1.01
# ESS will be low — expected given few iterations, not a concern here

# --- Divergences --------------------------------------------------------------
n_div <- nuts_params(m_test) %>%
  filter(Parameter == "divergent__") %>%
  summarise(divergences = sum(Value)) %>%
  pull(divergences)

message("Divergent transitions: ", n_div)

if (n_div > 10) {
  message("WARNING: too many divergences — increase adapt_delta to 0.95 in full model")
} else {
  message("OK: divergences acceptable")
}

# --- Posterior predictive check -----------------------------------------------
pp_check(m_test, ndraws = 50)
# Simulated lines should roughly follow the observed distribution


# ==============================================================================
# 6. Quick Curve Plot
# ==============================================================================
# Population-level predictions to visually check curve shape

newdata_test <- expand.grid(
  t        = 0:59,
  Genotype = unique(df_model$Genotype),
  Block    = unique(df_model$Block)
)

fitted_test <- newdata_test %>%
  add_epred_draws(m_test, ndraws = 50, re_formula = NA)

fitted_test %>%
  ggplot(aes(x = t + 1, y = .epred,
             color = Genotype, fill = Genotype)) +
  stat_lineribbon(.width = 0.95, alpha = 0.2) +
  facet_wrap(~ Block) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x        = "Stimulus number",
    y        = "Response probability",
    title    = "TEST FIT — visual check only, not for inference",
    subtitle = "Simplified formula, subsampled data"
  ) +
  theme_minimal()


# ==============================================================================
# 7. Decision Gate
# ==============================================================================
message("
=== TEST FIT CHECKLIST ===
Before running the full model, verify:

[ ] All Rhat < 1.01 in summary()
[ ] Divergent transitions < 10
[ ] pp_check() looks reasonable
[ ] Curves in plot start near 1.0 and decay plausibly
[ ] Genotype ordering matches biological expectation

If all pass -> run the full script (nlme_full.R)
If Rhat > 1.01 -> check which parameters, may need stronger priors
If divergences > 10 -> increase adapt_delta to 0.95 or 0.99
If curves look wrong -> revisit prior visualization in section 2
=========================
")
