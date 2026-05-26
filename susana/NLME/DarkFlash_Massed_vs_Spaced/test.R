# ==============================================================================
# Bayesian double-exponential habituation model for zebrafish DF training
# ==============================================================================

library(dplyr)
library(brms)
library(emmeans)
library(ggplot2)
library(readr)
library(tidyr)
library(ggplot2)
library(ggpubr)
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
# ------------------------------------------------------------------------------
# 1. Prepare training data
# ------------------------------------------------------------------------------

df_train <- df_all %>%
  filter(BlockRole == "training")


# Optional but useful check
df_train %>%
  distinct(animal, Training, Genotype) %>%
  count(Training, Genotype)

# ------------------------------------------------------------------------------
# 2. Nonlinear double-exponential habituation model
# ------------------------------------------------------------------------------
bf_hab <- bf(
  move ~ p,
  
  nlf(
    p ~ c + amp *
      (
        w * exp(-stimulus0 / taufast) +
          (1 - w) * exp(-stimulus0 / tauslow)
      )
  ),
  
  nlf(c ~ inv_logit(logitc)),
  nlf(amp ~ inv_logit(logitamp) * (1 - c)),
  nlf(w ~ inv_logit(logitw)),
  nlf(taufast ~ exp(logtaufast)),
  nlf(tauslow ~ taufast + exp(logtaudelta)),
  
  logitc ~ Genotype * Training + (1 | animal),
  logitamp ~ Genotype * Training + (1 | animal),
  logitw ~ Genotype * Training,
  logtaufast ~ Genotype * Training,
  logtaudelta ~ Genotype * Training,
  
  nl = TRUE
)

priors_hab <- c(
  # Asymptotic response probability c
  prior(normal(-2, 1), nlpar = "logitc", class = "b"),
  
  # Habituation amplitude
  prior(normal(1, 1), nlpar = "logitamp", class = "b"),
  
  # Fast-component weight
  prior(normal(0, 1), nlpar = "logitw", class = "b"),
  
  # Fast time constant
  prior(normal(log(20), 0.7), nlpar = "logtaufast", class = "b"),
  
  # Additional slow-time component
  prior(normal(log(150), 0.8), nlpar = "logtaudelta", class = "b"),
  
  # Animal-level variation
  prior(exponential(2), nlpar = "logitc", class = "sd"),
  prior(exponential(2), nlpar = "logitamp", class = "sd")
)

# ==============================================================================
# Fast VI test fit (sanity check before NUTS)
# ==============================================================================
fit_vi_joint <- brm(
  formula   = bf_hab,
  data      = df_train,
  family    = bernoulli(link = "identity"),
  prior     = priors_hab,
  backend   = "cmdstanr",
  algorithm = "meanfield",
  iter      = 10000,
  init      = 0
)


# ------------------------------------------------------------------------------
# 4. Fit model
# ------------------------------------------------------------------------------
fit_hab <- brm(
  formula = bf_hab,
  data = df_train,
  family = bernoulli(link = "identity"),
  prior = priors_hab,
  chains = 2,
  cores = 2,
  threads = threading(6),
  iter = 1000,
  warmup = 500,
  seed = 42,
  control = list(
    adapt_delta = 0.90,
    max_treedepth = 10
  ),
  backend = "cmdstanr",
  file = file.path(save_model_dir, "brms_double_exp_habituation_training")
)

# ------------------------------------------------------------------------------
# 5. Diagnostics
# ------------------------------------------------------------------------------

summary(fit_hab)

plot(fit_hab)

pp_check(fit_hab, ndraws = 100)

bayes_R2(fit_hab)

# Check sampler diagnostics
nuts_params <- nuts_params(fit_hab)

table(nuts_params$Parameter, nuts_params$Value)

# More useful:
np <- nuts_params(fit_hab)

np %>%
  filter(Parameter == "divergent__") %>%
  summarise(n_divergent = sum(Value))

np %>%
  filter(Parameter == "treedepth__") %>%
  summarise(max_treedepth_hit = sum(Value >= 15))


# ==============================================================================
# Plot double-exponential habituation curves
# per Training, per Block, per Genotype
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(brms)
library(ggpubr)

# ------------------------------------------------------------------------------
# Choose model to plot
# ------------------------------------------------------------------------------
fit_use <- fit_vi_joint
fit_use <- fit_hab

# ------------------------------------------------------------------------------
# Use training data only
# ------------------------------------------------------------------------------

df_plot_base <- df_train %>%
  mutate(
    stimulus0 = as.numeric(stimulus0),
    Training = factor(Training, levels = levels(df_train$Training)),
    Genotype = factor(Genotype, levels = levels(df_train$Genotype)),
    Block = factor(Block)
  )

# ------------------------------------------------------------------------------
# Build prediction grid: Training x Block x Genotype
# ------------------------------------------------------------------------------

new_data_hab <- df_plot_base %>%
  group_by(Training, Block, Genotype) %>%
  summarise(
    stim_min = min(stimulus0, na.rm = TRUE),
    stim_max = max(stimulus0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    grid = list(
      tibble(
        stimulus0 = seq(stim_min, stim_max, length.out = 200)
      )
    )
  ) %>%
  unnest(grid) %>%
  select(-stim_min, -stim_max) %>%
  ungroup() %>%
  mutate(
    Training = factor(Training, levels = levels(df_plot_base$Training)),
    Genotype = factor(Genotype, levels = levels(df_plot_base$Genotype)),
    Block = factor(Block, levels = levels(df_plot_base$Block)),
    animal = NA
  )

# ------------------------------------------------------------------------------
# Important for binned binomial models:
# posterior_epred() returns expected successes.
# Setting ntrials = 1 makes it return response probability.
# This is harmless for Bernoulli models if the variable is unused.
# ------------------------------------------------------------------------------

new_data_hab <- new_data_hab %>%
  mutate(
    ntrials = 1,
    stimulus0bin = stimulus0
  )

# ------------------------------------------------------------------------------
# Posterior expected response probabilities
# ------------------------------------------------------------------------------

ep_hab <- posterior_epred(
  fit_use,
  newdata = new_data_hab,
  re_formula = NA,
  allow_new_levels = TRUE
)

pred_summary <- t(apply(
  ep_hab,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975),
  na.rm = TRUE
))

new_data_hab <- new_data_hab %>%
  mutate(
    CI_low = pred_summary[, 1],
    fit = pred_summary[, 2],
    CI_high = pred_summary[, 3]
  )

# ------------------------------------------------------------------------------
# Raw observed response probability per stimulus
# ------------------------------------------------------------------------------

raw_summary_hab <- df_plot_base %>%
  group_by(Training, Block, Genotype, stimulus0) %>%
  summarise(
    p_move = mean(move, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Optional: binned raw points for cleaner plotting
# ------------------------------------------------------------------------------

bin_width_plot <- 10

raw_summary_hab_binned <- df_plot_base %>%
  mutate(
    stimbin = floor(stimulus0 / bin_width_plot),
    stimulus0mid = stimbin * bin_width_plot + bin_width_plot / 2
  ) %>%
  group_by(Training, Block, Genotype, stimulus0mid) %>%
  summarise(
    p_move = mean(move, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_massed_curves <- ggplot(
  new_data_hab %>% filter(Training == "massed"),
  aes(x = stimulus0, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_hab_binned %>% filter(Training == "massed"),
    aes(x = stimulus0mid, y = p_move, color = Genotype),
    inherit.aes = FALSE,
    alpha = 0.35,
    size = 0.8
  ) +
  geom_ribbon(
    aes(ymin = CI_low, ymax = CI_high),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(
    aes(y = fit),
    linewidth = 1.1
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Response probability",
    title = "Double-exponential habituation model: Massed training"
  ) +
  theme(
    legend.position = "none"
  )

print(p_massed_curves)

p_spaced_curves <- ggplot(
  new_data_hab %>% filter(Training == "spaced"),
  aes(x = stimulus0, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_hab_binned %>% filter(Training == "spaced"),
    aes(x = stimulus0mid, y = p_move, color = Genotype),
    inherit.aes = FALSE,
    alpha = 0.35,
    size = 0.8
  ) +
  geom_ribbon(
    aes(ymin = CI_low, ymax = CI_high),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(
    aes(y = fit),
    linewidth = 1.1
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Stimulus number within block",
    y = "Response probability",
    title = "Double-exponential habituation model: Spaced training"
  ) +
  theme(
    legend.position = "top"
  )

print(p_spaced_curves)

