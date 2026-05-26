# ==============================================================================
# Joint Bayesian NLME: spaced vs massed training, memory test
# ------------------------------------------------------------------------------
# Combines the massed dataset (1 long training block + 1 test block) and the
# spaced dataset (4 short training blocks + 1 test block) into one model so
# that the spaced-vs-massed memory contrast can be estimated as a single
# posterior contrast with proper uncertainty propagation.
#
# Model (per animal, per block):
#   p(move) = inv_logit(A) + (inv_logit(R0) - inv_logit(A)) * exp(-exp(logk)*stim0)
#
# Fixed-effect structure on each nonlinear parameter:
#   ~ Genotype * Training * BlockRole + Genotype * Training * Block + (1|animal)
#
# where:
#   Training   in {massed, spaced}
#   Block      = literal block index within each experiment
#   BlockRole  in {training, test}    -- used to define the memory contrast
#
# Key derived quantity (computed from posterior draws after fitting):
#   recovery        = inv_logit(R0_test) - inv_logit(A_lastTrainingBlock)
#   recovery_diff   = recovery_spaced - recovery_massed
#
# A negative recovery_diff means spaced training produced LESS recovery
# in the test block, i.e. BETTER retention than massed.
# ==============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(brms)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(loo)
library(DHARMa)
library(bayesplot)
library(purrr)
library(stringr)
library(patchwork)


# ==============================================================================
# Paths
# ==============================================================================
# source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_Massed_vs_Spaced//utils.R")
source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/NLME/DarkFlash_Massed_vs_Spaced//utils.R")

# base_dir <- "D:/WorkingData/Susana/NLME/DarkFlash_Joint_SpacedVsMassed/"
base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/NLME/DarkFlash_Joint_SpacedVsMassed/"

file_massed <- file.path(
  base_dir,
  "data_files",
  "SPZ_Massed_Training_7Nov2025.csv"
)

file_spaced <- file.path(
  base_dir,
  "data_files",
  "SPZ_Spaced_Training_Nov2025.csv"    
)

var_name <- "response_prob"
col_name <- "move"

save_fig_dir     <- file.path(base_dir, "figs",    "nlme_joint_response_prob", var_name)
save_results_dir <- file.path(base_dir, "results", "nlme_joint_response_prob", var_name)
save_model_dir   <- file.path(base_dir, "models")

dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_model_dir,   recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# Load and prepare both datasets
# ==============================================================================
res_massed <- load_data(file_massed, move_th = 0, drop = c('th2, tyr', 'tyr'))
res_spaced <- load_data(file_spaced, move_th = 0, drop = c('th2, tyr', 'tyr'))

# Use df_final (NOT df_final_sub) -- df_final_sub is responders-only.
df_massed <- res_massed$df_final
df_spaced <- res_spaced$df_final


# ------------------------------------------------------------------------------
# Tag each row with its Training condition and define BlockRole
# ------------------------------------------------------------------------------
# In each experiment, the LAST block is the memory test. All preceding blocks
# are training blocks.
massed_blocks_train <- "1"            # massed: block 1 = training
massed_block_test   <- "2"            # massed: block 2 = test

spaced_blocks_train <- c("1","2","3","4")   # spaced: blocks 1-4 = training
spaced_block_test   <- "5"                  # spaced: block 5 = test


df_massed_tagged <- df_massed %>%
  mutate(
    Training  = "massed",
    BlockRole = ifelse(as.character(Block) == massed_block_test, "test", "training")
  )

df_spaced_tagged <- df_spaced %>%
  mutate(
    Training  = "spaced",
    BlockRole = ifelse(as.character(Block) == spaced_block_test, "test", "training")
  )


# ------------------------------------------------------------------------------
# CRITICAL: animal IDs must be unique across experiments.
# Same Video x Well combo from different experiments would otherwise collide.
# We prefix the animal label with the Training condition.
# ------------------------------------------------------------------------------
df_all <- bind_rows(df_massed_tagged, df_spaced_tagged) %>%
  mutate(
    stimulus    = as.numeric(stimulus),
    stimulus0   = stimulus - 1,
    Training    = factor(Training,  levels = c("massed", "spaced")),
    Block       = factor(Block),
    BlockRole   = factor(BlockRole, levels = c("training", "test")),
    Genotype    = factor(Genotype),
    Video       = factor(Video),
    Well        = factor(Well),
    animal      = factor(paste0(Training, "_", Video, ".", Well))
  )

# Remove all non-responders to stimulus 1 in Block 1
# df_filtered <- df_all %>%
#   group_by(animal) %>%
#   filter(!any(Block == 1 & stimulus == 1 & move == 0)) %>%
#   ungroup()
# 
# summary_compare <- bind_rows(
#   df_all %>%
#     distinct(animal, Genotype, Training) %>%
#     mutate(dataset = "before"),
#   
#   df_filtered %>%
#     distinct(animal, Genotype, Training) %>%
#     mutate(dataset = "after")
# ) %>%
#   count(dataset, Genotype, Training)
# 
# summary_compare
# df_all <- df_filtered

write.csv(df_all, file.path(base_dir, 'data_files', 'final_data.csv'), row.names = FALSE)

# Sanity checks ----------------------------------------------------------------
cat("\n--- Rows per Training x Block ---\n")
print(df_all %>% count(Training, Block, BlockRole))

cat("\n--- Animals per Genotype x Training ---\n")
print(
  df_all %>%
    distinct(animal, Genotype, Training) %>%
    count(Genotype, Training)
)

cat("\n--- Marginal response rate per Training x BlockRole ---\n")
print(
  df_all %>%
    group_by(Training, BlockRole) %>%
    summarise(p_move = mean(.data[[col_name]]), n = n(), .groups = "drop")
)

cat("\n--- Stimulus range per Training x Block ---\n")
print(
  df_all %>%
    group_by(Training, Block) %>%
    summarise(min_stim = min(stimulus0), max_stim = max(stimulus0), .groups = "drop")
)


# ==============================================================================
# The Model
# ==============================================================================
# Three-way interaction Genotype * Training * Block on every nonlinear parameter.
# This is the most flexible structure and is what you want for the spaced-vs-
# massed memory contrast: each (Genotype, Training, Block) cell gets its own
# R0, A, and logk, with partial pooling across animals within cells.
#
# Note: we use Block (not BlockRole) as the fixed effect. BlockRole is only
# used afterwards when extracting the recovery contrast from posterior draws,
# because "the test block" maps to different literal Block labels in the two
# experiments (Block 2 in massed, Block 5 in spaced).
model_joint <- bf(
  # move ~ inv_logit(A) + (inv_logit(R0) - inv_logit(A)) * exp(-exp(logk) * stimulus0)
  # move ~ inv_logit(A) + (inv_logit(R0) - inv_logit(A)) * (stimulus0 + 1)^(-exp(logalpha)),  # power law
  move ~ inv_logit(A) + (inv_logit(R0) - inv_logit(A)) * exp(-exp(logk) * stimulus0^exp(logbeta)),  # strechted exp
  
  A    ~ 1 + Genotype * Training * Block + (1 | animal),
  R0   ~ 1 + Genotype * Training * Block + (1 | animal),
  logk ~ 1 + Genotype * Training * Block + (1 | animal),
  # logalpha ~ 1 + Genotype * Training * Block + (1 | animal),
  logbeta ~ 1 + Genotype * Training + (1 | animal),
  
  
  nl = TRUE
)


# ==============================================================================
# Priors
# ==============================================================================
# Massed Block 1 has ~478 stimuli; spaced training blocks have ~120 stimuli;
# test blocks have ~9. So plausible logk spans roughly:
#   k in [0.003, 0.1]  -> logk in [-5.8, -2.3]
# Center logk at -4 with sd=1.5 to comfortably cover both extremes.
#
# A and R0 priors are on the logit probability scale and unchanged from the
# single-experiment scripts.
priors_joint <- c(
  prior(normal(0,    1.0), class = "b", coef = "Intercept", nlpar = "A"),
  prior(normal(1.5,  1.0), class = "b", coef = "Intercept", nlpar = "R0"),
  prior(normal(-1.5, 1.5), class = "b", coef = "Intercept", nlpar = "logk"),
  # prior(normal(-0.7, 0.6), class = "b", coef = "Intercept", nlpar = "logalpha"),
  prior(normal(-0.3, 0.5), class = "b", coef = "Intercept", nlpar = "logbeta"),
  
  prior(normal(0, 0.75), class = "b", nlpar = "A"),
  prior(normal(0, 0.75), class = "b", nlpar = "R0"),
  prior(normal(0, 0.75), class = "b", nlpar = "logk"),
  # prior(normal(0, 0.4), class = "b", nlpar = "logalpha"),
  prior(normal(0, 0.3), class = "b", nlpar = "logbeta"),
  
  prior(exponential(4), class = "sd", nlpar = "A"),
  prior(exponential(4), class = "sd", nlpar = "R0"),
  prior(exponential(4), class = "sd", nlpar = "logk"),
  # prior(exponential(6), class = "sd", nlpar = "logalpha")
  prior(exponential(8), class = "sd", nlpar = "logbeta")
)


# ==============================================================================
# Prior predictive check
# ==============================================================================
fit_prior <- brm(
  formula      = model_joint,
  data         = df_all,
  family       = bernoulli(link = "identity"),
  prior        = priors_joint,
  sample_prior = "only",
  backend      = "cmdstanr",
  chains       = 4,
  cores        = 4,
  iter         = 2000,
  warmup       = 500,
  seed         = 42
)

# ------------------------------------------------------------------------------
# Extract posterior draws (= prior draws here)
# ------------------------------------------------------------------------------
draws_df <- as_draws_df(fit_prior)

names(draws_df)

# Plot Priors
intercept_draws <- draws_df %>%
  select(
    starts_with("b_A_Intercept"),
    starts_with("b_R0_Intercept"),
    starts_with("b_logk_Intercept")
  ) %>%
  pivot_longer(everything(),
               names_to = "parameter",
               values_to = "value")

ggplot(intercept_draws, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  facet_wrap(~parameter, scales = "free") +
  theme_pubr(base_size = 12) +
  labs(
    title = "Prior distributions: intercepts",
    x = "Parameter value",
    y = "Density"
  )

prob_draws <- draws_df %>%
  transmute(
    A_prob  = plogis(b_A_Intercept),
    R0_prob = plogis(b_R0_Intercept)
  ) %>%
  pivot_longer(everything(),
               names_to = "parameter",
               values_to = "probability")

ggplot(prob_draws, aes(x = probability)) +
  geom_density(fill = "purple", alpha = 0.4) +
  facet_wrap(~parameter, scales = "free") +
  theme_pubr(base_size = 12) +
  labs(
    title = "Intercept priors on probability scale",
    x = "Probability",
    y = "Density"
  )

sd_draws <- draws_df %>%
  select(starts_with("sd_animal")) %>%
  pivot_longer(everything(),
               names_to = "parameter",
               values_to = "value")

ggplot(sd_draws, aes(x = value)) +
  geom_density(fill = "forestgreen", alpha = 0.4) +
  facet_wrap(~parameter, scales = "free") +
  theme_pubr(base_size = 12) +
  labs(
    title = "Prior distributions: random-effect SDs",
    x = "SD value",
    y = "Density"
  )

x <- seq(-10, 10, length.out = 2000)

priors_df <- bind_rows(
  
  tibble(
    x = x,
    density = dnorm(x, 0, 1),
    prior = "A intercept ~ Normal(0,1)"
  ),
  
  tibble(
    x = x,
    density = dnorm(x, 1.5, 1),
    prior = "R0 intercept ~ Normal(1.5,1)"
  ),
  
  tibble(
    x = x,
    density = dnorm(x, -4, 1.5),
    prior = "logk intercept ~ Normal(-4,1.5)"
  ),
  
  tibble(
    x = x,
    density = dnorm(x, 0, 0.75),
    prior = "Fixed effects A/R0 ~ Normal(0,0.75)"
  ),
  
  tibble(
    x = x,
    density = dnorm(x, 0, 1),
    prior = "Fixed effects logk ~ Normal(0,1)"
  )
)

ggplot(priors_df, aes(x = x, y = density)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~prior, scales = "free") +
  theme_pubr(base_size = 12)

x_sd <- seq(0, 3, length.out = 1000)

sd_df <- tibble(
  x = x_sd,
  density = dexp(x_sd, rate = 4)
)

ggplot(sd_df, aes(x = x, y = density)) +
  geom_line(linewidth = 1.2) +
  theme_pubr(base_size = 12) +
  labs(
    title = "SD prior: Exponential(4)",
    x = "SD",
    y = "Density"
  )


# Generate prior predictive curves
curve_grid <- tibble(
  stimulus0 = seq(0, 500, length.out = 300),
  Genotype = levels(df_all$Genotype)[1],
  Training = levels(df_all$Training)[1],
  Block = levels(df_all$Block)[1],
  animal = NA
)

# curve_grid <- tibble(
#   stimulus0 = seq(0, 10, length.out = 100),
#   Genotype = levels(df_all$Genotype)[1],
#   Training = levels(df_all$Training)[1],
#   Block = levels(df_all$Block)[1],
#   animal = NA
# )

ep_prior <- posterior_epred(
  fit_prior,
  newdata = curve_grid,
  re_formula = NA,
  ndraws = 200
)

curve_df <- apply(ep_prior, 1, function(y) {
  tibble(
    stimulus0 = curve_grid$stimulus0,
    response = y
  )
}, simplify = FALSE) %>%
  bind_rows(.id = "draw")

ggplot(curve_df,
       aes(x = stimulus0, y = response, group = draw)) +
  geom_line(alpha = 0.08, color = "steelblue") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 12) +
  labs(
    title = "Prior predictive habituation curves",
    x = "Stimulus number",
    y = "Predicted response probability"
  )

ep <- posterior_epred(fit_prior, ndraws = 500)
cat("\nPrior predictive quantiles (expected p(move)):\n")
print(quantile(as.vector(ep), c(.001, .01, .05, .5, .95, .99, .999), na.rm = TRUE))
cat("\nObserved marginal P(move):", round(mean(df_all$move), 3), "\n")


# ==============================================================================
# Fast VI test fit (sanity check before NUTS)
# ==============================================================================
fit_vi_joint <- brm(
  formula   = model_joint,
  data      = df_all,
  family    = bernoulli(link = "identity"),
  prior     = priors_joint,
  backend   = "cmdstanr",
  algorithm = "meanfield",
  iter      = 10000,
  init      = 0
)

fit_joint <- fit_vi_joint

# ==============================================================================
# Full NUTS fit
# ==============================================================================
# With ~25k (massed) + ~25k (spaced) rows and a 3-way interaction on 3 nlpars,
# expect several hours of sampling. Use file= so you can resume / rerun without
# refitting on unchanged formula.
fit_joint <- brm(
  formula    = model_joint,
  data       = df_all,
  family     = bernoulli(link = "identity"),
  prior      = priors_joint,
  save_pars  = save_pars(all = TRUE),
  backend    = "cmdstanr",
  chains     = 4,
  cores      = 4,
  threads    = threading(6),
  iter       = 6000,
  warmup     = 3000,
  seed       = 42,
  control    = list(adapt_delta = 0.99, max_treedepth = 15),
  init       = 0,
)

saveRDS(
  fit_joint,
  file = file.path(save_model_dir, "bayesian_nlme_joint_respond_prob.rds")
)


# fit_joint <- readRDS(file.path(save_model_dir, "bayesian_nlme_joint_NUTS.rds"))


# ==============================================================================
# Summary and diagnostics
# ==============================================================================
sink(file.path(save_results_dir, "joint_model_response_prob_summary.txt"))
print(summary(fit_joint))
sink()

cat("\nN divergent transitions:\n")
print(sum(subset(nuts_params(fit_joint), Parameter == "divergent__")$Value))

pp_check(fit_joint, type = "dens_overlay", ndraws = 100)


# ==============================================================================
# Plot habituation curves -- per Training, per Block, per Genotype
# ==============================================================================
# Build prediction grid: per Training x Block x Genotype
new_data_joint <- df_all %>%
  group_by(Training, Block, Genotype) %>%
  summarise(
    stim_min = min(stimulus0, na.rm = TRUE),
    stim_max = max(stimulus0, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    grid = list(tibble(stimulus0 = seq(stim_min, stim_max, length.out = 200)))
  ) %>%
  unnest(grid) %>%
  select(-stim_min, -stim_max) %>%
  ungroup() %>%
  mutate(
    stimulus = stimulus0,
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype)),
    animal   = NA
  )

# Posterior expected response probabilities
ep_joint <- posterior_epred(
  fit_joint,
  newdata = new_data_joint,
  re_formula = NA,
  allow_new_levels = TRUE
)

pred_summary <- t(apply(
  ep_joint,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975),
  na.rm = TRUE
))

new_data_joint <- new_data_joint %>%
  mutate(
    CI_low  = pred_summary[, 1],
    fit     = pred_summary[, 2],
    CI_high = pred_summary[, 3]
  )

# Raw observed response probability per stimulus
raw_summary_joint <- df_all %>%
  group_by(Training, Block, Genotype, stimulus0) %>%
  summarise(
    p_move = mean(move, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(stimulus = stimulus0)

# Massed plot
p_massed_curves <- ggplot(
  new_data_joint %>% filter(Training == "massed"),
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_joint %>% filter(Training == "massed"),
    aes(x = stimulus, y = p_move, color = Genotype),
    inherit.aes = FALSE,
    alpha = 0.25,
    size = 0.6
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Response probability",
    title = "Joint Bayesian NLME: Massed training"
  ) +
  theme(legend.position = "none")

# Spaced plot
p_spaced_curves <- ggplot(
  new_data_joint %>% filter(Training == "spaced"),
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_joint %>% filter(Training == "spaced"),
    aes(x = stimulus, y = p_move, color = Genotype),
    inherit.aes = FALSE,
    alpha = 0.25,
    size = 0.6
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Response probability",
    title = "Joint Bayesian NLME: Spaced training"
  ) +
  theme(legend.position = "top")

print(p_spaced_curves)
print(p_massed_curves)

# ==============================================================================
# Contrasts on the NLME (computed from posterior draws)
# ==============================================================================
# All contrasts are computed on the response probability scale by transforming
# inv_logit(R0) and inv_logit(A) from posterior draws. Credible intervals come
# from posterior quantiles -- no SE aggregation or joint-covariance issues.
#
# Sign conventions (response probability scale):
#   training_drop > 0      : habituation during the block
#   memory_dropR0 < 0      : memory retained (test start lower than naive start)
#   recovery > 0           : bounce-back from trained state
#   spaced_vs_massed < 0   : spaced retains more / habituates more / etc.
# ==============================================================================

# ------------------------------------------------------------------------------
# Define which Block label is which role per Training
# ------------------------------------------------------------------------------
training_blocks_per_protocol <- list(
  massed = c("1"),                  # massed: block 1 = training
  spaced = c("1", "2", "3", "4")    # spaced: blocks 1-4 = training
)
test_block_per_protocol <- list(
  massed = "2",
  spaced = "5"
)
last_train_block_per_protocol <- list(
  massed = "1",
  spaced = "4"
)
# "naive" baseline = the first stimulus of Block 1 in either protocol
baseline_block <- "1"


genotypes_levels <- levels(df_all$Genotype)
trainings_levels <- levels(df_all$Training)


# ------------------------------------------------------------------------------
# Helper: posterior draws of inv_logit(R0) and inv_logit(A) per
# (Genotype, Training, Block) on the response-probability scale.
# Returns one row per draw x cell.
# ------------------------------------------------------------------------------
get_param_draws <- function(model, param = c("R0", "A"), genotypes, trainings, blocks) {
  param <- match.arg(param)
  
  grid <- expand.grid(
    Genotype = genotypes,
    Training = trainings,
    Block    = blocks,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Training  = factor(Training, levels = levels(df_all$Training)),
      Block     = factor(Block,    levels = levels(df_all$Block)),
      Genotype  = factor(Genotype, levels = levels(df_all$Genotype)),
      stimulus0 = if (param == "R0") 0 else 1e6,   # 0 = R0, infinity ~ A
      animal    = NA,
      cell_id   = paste(Genotype, Training, Block, sep = "|")
    )
  
  ep <- posterior_epred(
    model,
    newdata           = grid,
    re_formula        = NA,
    allow_new_levels  = TRUE
  )
  # ep is [draws x rows]; assign cell_id column names
  colnames(ep) <- grid$cell_id
  
  draws_df <- as_tibble(ep) %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(
      cols = -.draw,
      names_to  = "cell_id",
      values_to = "prob"
    ) %>%
    separate(cell_id, into = c("Genotype", "Training", "Block"), sep = "\\|") %>%
    mutate(
      param = param,
      Genotype = factor(Genotype, levels = levels(df_all$Genotype)),
      Training = factor(Training, levels = levels(df_all$Training)),
      Block    = factor(Block,    levels = levels(df_all$Block))
    )
  
  draws_df
}


# Helper: summarise posterior draws of a contrast vector
summarise_contrast <- function(draws_vec, label = NA) {
  tibble(
    estimate = median(draws_vec),
    lower    = quantile(draws_vec, 0.025),
    upper    = quantile(draws_vec, 0.975),
    p_dir    = min(mean(draws_vec > 0), mean(draws_vec < 0)) * 2,  # two-sided "p"
    contrast = label
  )
}


# ------------------------------------------------------------------------------
# 0. Pull all the parameter draws we need (R0 and A at every relevant block)
# ------------------------------------------------------------------------------
# All blocks present in the data
all_blocks <- levels(df_all$Block)

R0_draws <- get_param_draws(fit_joint, "R0", genotypes_levels, trainings_levels, all_blocks)
A_draws  <- get_param_draws(fit_joint, "A",  genotypes_levels, trainings_levels, all_blocks)

# Wide format: one column per cell, one row per draw -- easier for arithmetic
R0_wide <- R0_draws %>%
  mutate(cell = paste(Genotype, Training, Block, sep = "|")) %>%
  select(.draw, cell, prob) %>%
  pivot_wider(names_from = cell, values_from = prob)

A_wide <- A_draws %>%
  mutate(cell = paste(Genotype, Training, Block, sep = "|")) %>%
  select(.draw, cell, prob) %>%
  pivot_wider(names_from = cell, values_from = prob)

# Helper to pull a column by (Genotype, Training, Block)
pull_cell <- function(wide_df, g, t, b) {
  col <- paste(g, t, b, sep = "|")
  if (!col %in% names(wide_df)) {
    stop(sprintf("Cell not found in draws: %s", col))
  }
  wide_df[[col]]
}


# ==============================================================================
# (L) Within-block learning: did the response drop during each training block?
# ==============================================================================
# learning_{g,t,b} = inv_logit(R0_{g,t,b}) - inv_logit(A_{g,t,b})
# Positive = habituation occurred in that block.
# Computed only for training blocks.

learning_rows <- list()
for (g in genotypes_levels) {
  for (t in trainings_levels) {
    train_blocks <- training_blocks_per_protocol[[t]]
    for (b in train_blocks) {
      r0_d <- pull_cell(R0_wide, g, t, b)
      a_d  <- pull_cell(A_wide,  g, t, b)
      drop_d <- r0_d - a_d
      s <- summarise_contrast(
        drop_d,
        label = sprintf("learning_drop[%s|%s|B%s]", g, t, b)
      )
      learning_rows[[length(learning_rows) + 1]] <- s %>%
        mutate(Genotype = g, Training = t, Block = b)
    }
  }
}
contrast_L <- bind_rows(learning_rows) %>%
  select(Genotype, Training, Block, estimate, lower, upper, p_dir, contrast)

cat("\n--- (L) Within-block learning: R0 - A per training block ---\n")
print(contrast_L)

write.csv(
  contrast_L,
  file.path(save_results_dir, "contrast_L_within_block_learning.csv"),
  row.names = FALSE
)


# ==============================================================================
# (M_R0) Within-protocol memory at test start vs naive Block 1 start
# ==============================================================================
# memory_{g,t} = inv_logit(R0_{test_block}) - inv_logit(R0_{block 1})
# Negative = test start lower than naive start = memory retained.
# Zero = no memory.

memory_R0_rows <- list()
for (g in genotypes_levels) {
  for (t in trainings_levels) {
    b_test <- test_block_per_protocol[[t]]
    r0_test <- pull_cell(R0_wide, g, t, b_test)
    r0_b1   <- pull_cell(R0_wide, g, t, baseline_block)
    mem_d   <- r0_test - r0_b1
    s <- summarise_contrast(
      mem_d,
      label = sprintf("memory_R0_test_vs_baseline[%s|%s]", g, t)
    )
    memory_R0_rows[[length(memory_R0_rows) + 1]] <- s %>%
      mutate(Genotype = g, Training = t)
  }
}
contrast_M_R0 <- bind_rows(memory_R0_rows) %>%
  select(Genotype, Training, estimate, lower, upper, p_dir, contrast)

cat("\n--- (M_R0) Memory: R0(test) - R0(naive Block 1) per (Genotype, Training) ---\n")
print(contrast_M_R0)

write.csv(
  contrast_M_R0,
  file.path(save_results_dir, "contrast_M_R0_test_vs_baseline.csv"),
  row.names = FALSE
)


# ==============================================================================
# (M_recovery) Recovery from trained state to test start
# ==============================================================================
# recovery_{g,t} = inv_logit(R0_{test_block}) - inv_logit(A_{last_train_block})
# Positive = response bounced back from trained state. Smaller = better retention.
# Note this is the same shape as your old contrast D, but defined on the
# nonlinear parameters directly.

recovery_rows <- list()
for (g in genotypes_levels) {
  for (t in trainings_levels) {
    b_test  <- test_block_per_protocol[[t]]
    b_train <- last_train_block_per_protocol[[t]]
    r0_test <- pull_cell(R0_wide, g, t, b_test)
    a_train <- pull_cell(A_wide,  g, t, b_train)
    rec_d   <- r0_test - a_train
    s <- summarise_contrast(
      rec_d,
      label = sprintf("recovery_R0test_minus_AlastTrain[%s|%s]", g, t)
    )
    recovery_rows[[length(recovery_rows) + 1]] <- s %>%
      mutate(Genotype = g, Training = t)
  }
}
contrast_M_recovery <- bind_rows(recovery_rows) %>%
  select(Genotype, Training, estimate, lower, upper, p_dir, contrast)

cat("\n--- (M_recovery) Recovery: R0(test) - A(last training block) ---\n")
print(contrast_M_recovery)

write.csv(
  contrast_M_recovery,
  file.path(save_results_dir, "contrast_M_recovery_R0test_vs_AlastTrain.csv"),
  row.names = FALSE
)


# ==============================================================================
# (SvM_memory) Spaced vs Massed: between-protocol memory contrast on R0
# ==============================================================================
# Difference of differences per Genotype:
#   [R0_test_spaced  - R0_baseline_spaced]   (= memory in spaced cohort)
# - [R0_test_massed  - R0_baseline_massed]   (= memory in massed cohort)
# Negative = spaced retains MORE memory than massed.
# This controls for any baseline differences between the two fish cohorts.

svm_memory_rows <- list()
for (g in genotypes_levels) {
  r0_test_s <- pull_cell(R0_wide, g, "spaced", test_block_per_protocol[["spaced"]])
  r0_b1_s   <- pull_cell(R0_wide, g, "spaced", baseline_block)
  r0_test_m <- pull_cell(R0_wide, g, "massed", test_block_per_protocol[["massed"]])
  r0_b1_m   <- pull_cell(R0_wide, g, "massed", baseline_block)
  
  diff_diff <- (r0_test_s - r0_b1_s) - (r0_test_m - r0_b1_m)
  s <- summarise_contrast(
    diff_diff,
    label = sprintf("SvM_memory[%s]", g)
  )
  svm_memory_rows[[length(svm_memory_rows) + 1]] <- s %>%
    mutate(Genotype = g)
}
contrast_SvM_memory <- bind_rows(svm_memory_rows) %>%
  select(Genotype, estimate, lower, upper, p_dir, contrast)

cat("\n--- (SvM_memory) Spaced vs Massed memory: (R0_test - R0_baseline) diff ---\n")
print(contrast_SvM_memory)

write.csv(
  contrast_SvM_memory,
  file.path(save_results_dir, "contrast_SvM_memory_diff_of_diff.csv"),
  row.names = FALSE
)


# ==============================================================================
# (SvM_recovery) Spaced vs Massed: between-protocol recovery contrast
# ==============================================================================
# Per Genotype:
#   recovery_spaced - recovery_massed
# Negative = spaced recovers less = better retention.
# This is exactly your old contrast D, on the NLME parameters.

svm_recovery_rows <- list()
for (g in genotypes_levels) {
  r0_test_s <- pull_cell(R0_wide, g, "spaced", test_block_per_protocol[["spaced"]])
  a_last_s  <- pull_cell(A_wide,  g, "spaced", last_train_block_per_protocol[["spaced"]])
  r0_test_m <- pull_cell(R0_wide, g, "massed", test_block_per_protocol[["massed"]])
  a_last_m  <- pull_cell(A_wide,  g, "massed", last_train_block_per_protocol[["massed"]])
  
  rec_diff <- (r0_test_s - a_last_s) - (r0_test_m - a_last_m)
  s <- summarise_contrast(
    rec_diff,
    label = sprintf("SvM_recovery[%s]", g)
  )
  svm_recovery_rows[[length(svm_recovery_rows) + 1]] <- s %>%
    mutate(Genotype = g)
}
contrast_SvM_recovery <- bind_rows(svm_recovery_rows) %>%
  select(Genotype, estimate, lower, upper, p_dir, contrast)

cat("\n--- (SvM_recovery) Spaced vs Massed recovery difference ---\n")
print(contrast_SvM_recovery)

write.csv(
  contrast_SvM_recovery,
  file.path(save_results_dir, "contrast_SvM_recovery_diff.csv"),
  row.names = FALSE
)


# ==============================================================================
# (SvM_trained_state) Spaced vs Massed: depth of training (end-state confound check)
# ==============================================================================
# Per Genotype:
#   inv_logit(A_{last_train_spaced}) - inv_logit(A_{last_train_massed})
# Negative = spaced ended training at a deeper habituated state.
# This is the "end-state confound" diagnostic: if SvM_memory looks good but
# SvM_trained_state is strongly negative, then spaced fish "remember" mainly
# because they ended training in a more habituated state, not because they
# resisted forgetting.

svm_trained_rows <- list()
for (g in genotypes_levels) {
  a_last_s <- pull_cell(A_wide, g, "spaced", last_train_block_per_protocol[["spaced"]])
  a_last_m <- pull_cell(A_wide, g, "massed", last_train_block_per_protocol[["massed"]])
  diff_d   <- a_last_s - a_last_m
  s <- summarise_contrast(
    diff_d,
    label = sprintf("SvM_trained_state[%s]", g)
  )
  svm_trained_rows[[length(svm_trained_rows) + 1]] <- s %>%
    mutate(Genotype = g)
}
contrast_SvM_trained_state <- bind_rows(svm_trained_rows) %>%
  select(Genotype, estimate, lower, upper, p_dir, contrast)

cat("\n--- (SvM_trained_state) Spaced vs Massed: A at last training block ---\n")
print(contrast_SvM_trained_state)

write.csv(
  contrast_SvM_trained_state,
  file.path(save_results_dir, "contrast_SvM_trained_state.csv"),
  row.names = FALSE
)


# ==============================================================================
# (Test_C) Test block average: mean response over test block 1..8 stimuli
# ==============================================================================
# This is the closest analog to your old Contrast C, but here we compute the
# expected p(move) by averaging the NLME prediction across stim 1..8 of the
# test block for each (Genotype, Training).
# Then take spaced - massed per Genotype.

test_stims <- 1:8

# Build the prediction grid
test_grid <- expand.grid(
  Genotype  = genotypes_levels,
  Training  = trainings_levels,
  stimulus0 = test_stims,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Block = ifelse(Training == "massed",
                   test_block_per_protocol[["massed"]],
                   test_block_per_protocol[["spaced"]]),
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype)),
    animal   = NA,
    cell_id  = paste(Genotype, Training, stimulus0, sep = "|")
  )

ep_test <- posterior_epred(
  fit_joint,
  newdata          = test_grid,
  re_formula       = NA,
  allow_new_levels = TRUE
)
colnames(ep_test) <- test_grid$cell_id

# Average across the 8 stim per (Genotype, Training) draw-by-draw
test_long <- as_tibble(ep_test) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(cols = -.draw, names_to = "cell_id", values_to = "prob") %>%
  separate(cell_id, into = c("Genotype", "Training", "stimulus0"),
           sep = "\\|", convert = TRUE) %>%
  group_by(.draw, Genotype, Training) %>%
  summarise(prob_mean = mean(prob), .groups = "drop")

# Spaced - massed per Genotype
test_C_rows <- list()
for (g in genotypes_levels) {
  s_d <- test_long$prob_mean[test_long$Genotype == g & test_long$Training == "spaced"]
  m_d <- test_long$prob_mean[test_long$Genotype == g & test_long$Training == "massed"]
  diff_d <- s_d - m_d
  s <- summarise_contrast(diff_d, label = sprintf("Test_C[%s]", g))
  test_C_rows[[length(test_C_rows) + 1]] <- s %>% mutate(Genotype = g)
}
contrast_Test_C <- bind_rows(test_C_rows) %>%
  select(Genotype, estimate, lower, upper, p_dir, contrast)

cat("\n--- (Test_C) Mean P(move) over test stim 1-8, spaced - massed ---\n")
print(contrast_Test_C)

write.csv(
  contrast_Test_C,
  file.path(save_results_dir, "contrast_TestC_meanStim_spaced_vs_massed.csv"),
  row.names = FALSE
)


# Also keep per-cell test-block means for plotting
test_per_cell <- test_long %>%
  group_by(Genotype, Training) %>%
  summarise(
    prob_mean = median(prob_mean),
    lower     = quantile(prob_mean, 0.025),
    upper     = quantile(prob_mean, 0.975),
    .groups   = "drop"
  )
# (Note: the within-summarise quantile uses the just-rewritten column.
#  If your R/dplyr is strict about this, replace with the line below:)
# test_per_cell <- test_long %>%
#   group_by(Genotype, Training) %>%
#   summarise(
#     prob_median = median(prob_mean),
#     lower       = quantile(prob_mean, 0.025),
#     upper       = quantile(prob_mean, 0.975),
#     .groups     = "drop"
#   ) %>% rename(prob_mean = prob_median)

write.csv(
  test_per_cell,
  file.path(save_results_dir, "per_cell_test_block_mean.csv"),
  row.names = FALSE
)


# ==============================================================================
# Plots
# ==============================================================================

# (L) Within-block learning
p_L <- contrast_L %>%
  ggplot(aes(x = Block, y = estimate, ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_grid(Genotype ~ Training, scales = "free_x") +
  geom_pointrange(linewidth = 0.7, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_cartesian(ylim = c(-0.2, 1)) +
  theme_pubr(base_size = 13) +
  labs(
    x = "Training block",
    y = "P(R0) - P(A): how much the response dropped\nPositive = habituation occurred",
    title = "(L) Within-block learning"
  ) +
  theme(legend.position = "none")
ggsave(
  file.path(save_fig_dir, "contrastL_within_block_learning.png"),
  p_L, width = 9, height = 9, dpi = 300, bg = "white"
)


# (M_R0) Memory: R0 test vs baseline
p_M_R0 <- contrast_M_R0 %>%
  ggplot(aes(x = Training, y = estimate, ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "P(R0 test) - P(R0 naive Block 1)\nNegative = memory retained, 0 = no memory",
    title = "(M_R0) Test start vs naive baseline"
  ) +
  theme(legend.position = "none")
ggsave(
  file.path(save_fig_dir, "contrastM_R0_test_vs_baseline.png"),
  p_M_R0, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# (M_recovery) Recovery: R0_test - A_lastTrain
p_M_recovery <- contrast_M_recovery %>%
  ggplot(aes(x = Training, y = estimate, ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 3) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "P(R0 test) - P(A last training block)\nSmaller = better retention",
    title = "(M_recovery) Recovery from trained state to test"
  ) +
  theme(legend.position = "none")
ggsave(
  file.path(save_fig_dir, "contrastM_recovery_R0test_vs_AlastTrain.png"),
  p_M_recovery, width = 10, height = 4.5, dpi = 300, bg = "white"
)


# (SvM) Spaced vs Massed summary
svm_panel <- bind_rows(
  contrast_SvM_memory       %>% mutate(metric = "Memory\n(R0_test - R0_baseline)"),
  contrast_SvM_recovery     %>% mutate(metric = "Recovery\n(R0_test - A_lastTrain)"),
  contrast_SvM_trained_state %>% mutate(metric = "Trained state\n(A_lastTrain)")
)

p_SvM <- svm_panel %>%
  ggplot(aes(x = Genotype, y = estimate, ymin = lower, ymax = upper,
             color = Genotype)) +
  facet_wrap(~ metric, ncol = 1, scales = "free_x") +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 13) +
  labs(
    x = NULL,
    y = "Spaced - Massed difference (probability scale)",
    title = "Spaced vs Massed across NLME parameter contrasts"
  ) +
  theme(legend.position = "none")
ggsave(
  file.path(save_fig_dir, "contrast_SvM_summary.png"),
  p_SvM, width = 9, height = 9, dpi = 300, bg = "white"
)


# (Test_C)
p_Test_C <- contrast_Test_C %>%
  ggplot(aes(x = Genotype, y = estimate, ymin = lower, ymax = upper,
             color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x = NULL,
    y = "Spaced - Massed mean P(move) over test stim 1-8\nNegative = spaced has better retention",
    title = "(Test_C) Mean test-block response, spaced vs massed"
  ) +
  theme(legend.position = "none")
ggsave(
  file.path(save_fig_dir, "contrast_TestC_meanStim_spaced_vs_massed.png"),
  p_Test_C, width = 8, height = 5, dpi = 300, bg = "white"
)


message("All NLME contrasts complete.")
message("Results saved under: ", save_results_dir)
message("Figures saved under: ", save_fig_dir)
