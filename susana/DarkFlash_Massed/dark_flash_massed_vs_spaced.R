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


# ==============================================================================
# Paths
# ==============================================================================
source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/utils.R")
source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/plot_utils.R")

base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_Joint_SpacedVsMassed/"

file_massed <- file.path(
  "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_Massed",
  "data_files",
  "SPZ_Massed_Training_7Nov2025.csv"
)

file_spaced <- file.path(
  "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_Spaced",
  "data_files",
  "SPZ_Spaced_Training_Nov2025.csv"    
)

var_name <- "response_prob"
col_name <- "move"

save_fig_dir     <- file.path(base_dir, "figs",    "nlme_joint", var_name)
save_results_dir <- file.path(base_dir, "results", "nlme_joint", var_name)
save_model_dir   <- file.path(base_dir, "models")

dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_model_dir,   recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# Load and prepare both datasets
# ==============================================================================
res_massed <- load_data_darkflash_60s(file_massed, move_th = 1, take_peak = 1)
res_spaced <- load_data_darkflash_60s(file_spaced, move_th = 1, take_peak = 1)

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
    stimulus0   = stimulus,                          # already 0-indexed per block
    Training    = factor(Training,  levels = c("massed", "spaced")),
    Block       = factor(Block),
    BlockRole   = factor(BlockRole, levels = c("training", "test")),
    Genotype    = factor(Genotype),
    Video       = factor(Video),
    Well        = factor(Well),
    animal      = factor(paste0(Training, "_", Video, ".", Well))
  )

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
  move ~ inv_logit(A) + (inv_logit(R0) - inv_logit(A)) * exp(-exp(logk) * stimulus0),
  
  A    ~ 1 + Genotype * Training * Block + (1 | animal),
  R0   ~ 1 + Genotype * Training * Block + (1 | animal),
  logk ~ 1 + Genotype * Training * Block + (1 | animal),
  
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
  prior(normal(0,    1),     class = "b", coef = "Intercept", nlpar = "A"),
  prior(normal(1.5,  1),     class = "b", coef = "Intercept", nlpar = "R0"),
  prior(normal(-4,   1.5),   class = "b", coef = "Intercept", nlpar = "logk"),
  
  prior(normal(0, 0.75), class = "b", nlpar = "A"),
  prior(normal(0, 0.75), class = "b", nlpar = "R0"),
  prior(normal(0, 1.0),  class = "b", nlpar = "logk"),  # slightly wider to cover both scales
  
  prior(exponential(4), class = "sd", nlpar = "A"),
  prior(exponential(4), class = "sd", nlpar = "R0"),
  prior(exponential(4), class = "sd", nlpar = "logk")
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
  file = file.path(save_model_dir, "bayesian_nlme_joint.rds")
)


# fit_joint <- readRDS(file.path(save_model_dir, "bayesian_nlme_joint_NUTS.rds"))


# ==============================================================================
# Summary and diagnostics
# ==============================================================================
sink(file.path(save_results_dir, "joint_model_summary.txt"))
print(summary(fit_joint))
sink()

cat("\nN divergent transitions:\n")
print(sum(subset(nuts_params(fit_joint), Parameter == "divergent__")$Value))

pp_check(fit_joint, type = "dens_overlay", ndraws = 100)


# ==============================================================================
# Plot habituation curves -- per Training, per Block, per Genotype
# ==============================================================================
library(patchwork)

p_massed <- pred_grid %>% filter(Training == "massed") %>%
  ggplot(aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_prob_all %>% filter(Training == "massed"),
    aes(x = stimulus, y = response_prob, color = Genotype),
    inherit.aes = FALSE, alpha = 0.15, size = 0.6
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 12) +
  labs(title = "Massed training", x = NULL, y = "Response probability") +
  theme(legend.position = "none")

p_spaced <- pred_grid %>% filter(Training == "spaced") %>%
  ggplot(aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_prob_all %>% filter(Training == "spaced"),
    aes(x = stimulus, y = response_prob, color = Genotype),
    inherit.aes = FALSE, alpha = 0.15, size = 0.6
  ) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 12) +
  labs(title = "Spaced training", x = "Stimulus within block", y = "Response probability") +
  theme(legend.position = "top")

p_combined <- p_massed / p_spaced + plot_layout(heights = c(1, 1))
print(p_combined)


# ==============================================================================
# MEMORY CONTRASTS  (revised - asymptote-free)
# ------------------------------------------------------------------------------
# Replaces the old "recovery vs end-of-training asymptote" block.
#
# The asymptote-based contrast assumed both protocols reach a comparable
# "trained state" before the test block. The data don't support this:
#   - Massed reaches A ~ 0.23 in 1 long block (478 stimuli)
#   - Spaced reaches A ~ 0.14 across 4 short blocks (480 stimuli total)
#   - Within spaced training, fish recover between blocks
# So "recovery" measured against the last block's asymptote conflates "more
# learning" with "better retention". We avoid this with three contrasts:
#
#   (1) ABSOLUTE TEST RESPONSE
#       Just R0_test per (Genotype, Training). No asymptote involved.
#       Lower test response = better retention.
#       Headline contrast: diff = R0_test_spaced - R0_test_massed.
#
#   (5) BETWEEN-BLOCK RECOVERY  (within-protocol, asymptote-free)
#       Massed:   R0_block2 - R0_block1    (test vs training start)
#       Spaced:   R0_block5 - R0_block4    (test vs last training block start)
#       Then compare across protocols.
#       This asks "did the inter-block pause cause more recovery in one
#       protocol than the other".
#       Bonus: for spaced we also compute R0_block_(i+1) - R0_block_i for
#       all i = 1..4 to characterize within-training recovery rhythm.
#
#   (4) RE-HABITUATION SAVINGS
#       k_test - k_block1 per protocol.
#       Savings = trained fish habituates faster than naive fish.
#       Compare across protocols. Test blocks have only ~9 stimuli so k_test
#       is poorly identified -- report this as supporting, not primary.
#
# All three contrasts are computed from posterior_epred draws on tiny custom
# grids, exactly the same machinery as the old recovery block.
# ==============================================================================

LARGE <- 1e6   # for sampling asymptote points (kept here for completeness)

last_train_block <- list(massed = "1", spaced = "4")
test_block       <- list(massed = "2", spaced = "5")
first_train_block <- list(massed = "1", spaced = "1")

genotypes <- levels(df_all$Genotype)


# ==============================================================================
# CONTRAST (1): ABSOLUTE TEST RESPONSE   --- HEADLINE
# ==============================================================================
# P(move) at the very start of the test block = inv_logit(R0_test).
# Smaller = better retention (fish stays habituated despite the pause).
# ==============================================================================

abs_test_grid <- tidyr::expand_grid(
  Genotype = genotypes,
  Training = c("massed", "spaced")
) %>%
  rowwise() %>%
  mutate(
    Block     = test_block[[Training]],
    stimulus0 = 0,
    stimulus  = 0
  ) %>%
  ungroup() %>%
  mutate(
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  )

abs_test_raw <- posterior_epred(
  fit_joint,
  newdata    = abs_test_grid,
  re_formula = NA,
  ndraws     = 1000
)

abs_test_long <- abs_test_grid %>%
  mutate(row_id = row_number()) %>%
  select(row_id, Genotype, Training)

abs_test_draws <- as_tibble(t(abs_test_raw)) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(-row_id, names_to = ".draw", values_to = "test_response") %>%
  mutate(.draw = as.integer(str_remove(.draw, "V"))) %>%
  left_join(abs_test_long, by = "row_id") %>%
  select(Genotype, Training, .draw, test_response)


# Per-cell summary
abs_test_summary <- abs_test_draws %>%
  group_by(Genotype, Training) %>%
  summarise(
    median = median(test_response),
    low    = quantile(test_response, 0.025),
    high   = quantile(test_response, 0.975),
    .groups = "drop"
  )

print(abs_test_summary)

write.csv(
  abs_test_summary,
  file.path(save_results_dir, "absolute_test_response_per_genotype.csv"),
  row.names = FALSE
)


# Spaced vs massed contrast on absolute test response
abs_test_diff <- abs_test_draws %>%
  pivot_wider(names_from = Training, values_from = test_response) %>%
  mutate(diff = spaced - massed)

abs_test_diff_summary <- abs_test_diff %>%
  group_by(Genotype) %>%
  summarise(
    diff_median        = median(diff),
    diff_low           = quantile(diff, 0.025),
    diff_high          = quantile(diff, 0.975),
    p_spaced_lt_massed = mean(diff < 0),   # P(spaced responds LESS = better retention)
    p_spaced_gt_massed = mean(diff > 0),
    .groups = "drop"
  )

print(abs_test_diff_summary)

write.csv(
  abs_test_diff_summary,
  file.path(save_results_dir, "absolute_test_response_diff_spaced_vs_massed.csv"),
  row.names = FALSE
)


# Plot: absolute test response per genotype x training
p_abs_test <- abs_test_summary %>%
  ggplot(aes(x = Training, y = median,
             ymin = low, ymax = high,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "P(move) at first test stimulus",
    title = "Absolute test response: lower = better memory retention"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrast1_absolute_test_response.png"),
  p_abs_test, width = 12, height = 4.5, dpi = 300, bg = "white"
)

# Plot: spaced - massed contrast
p_abs_test_diff <- abs_test_diff_summary %>%
  ggplot(aes(x = Genotype, y = diff_median,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x     = NULL,
    y     = "P(move) test difference (spaced - massed)\nNegative = spaced has better retention",
    title = "Headline contrast: spaced vs massed test response"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrast1_absolute_test_response_diff.png"),
  p_abs_test_diff, width = 8, height = 5, dpi = 300, bg = "white"
)




# ==============================================================================
# CONTRAST (5): BETWEEN-BLOCK RECOVERY  (within-protocol, asymptote-free)
# ==============================================================================
# For each protocol, look at how much R0 jumps between the start of the LAST
# training block and the start of the TEST block:
#
#   massed:  jump_test = R0_block2 - R0_block1
#   spaced:  jump_test = R0_block5 - R0_block4
#
# This is "how much did the fish recover during the inter-block pause
# (= the memory test pause)" expressed entirely in terms of R0 values, with
# no asymptote estimation needed.
#
# For SPACED we additionally compute the within-training jumps
#   spaced: jump_i = R0_block(i+1) - R0_block(i)   for i = 1..3
# to characterize the rhythm of recovery between consecutive training blocks.
# This contextualizes the test-pause jump: is it bigger than between-training
# jumps (suggests more forgetting from the longer pause), or comparable
# (suggests test pause is similar to training pauses)?
# ==============================================================================

# Build grid of (Training, Block) cells at stimulus0 = 0  -> gives inv_logit(R0)
r0_grid_massed <- tidyr::expand_grid(
  Genotype = genotypes,
  Block    = c("1", "2")
) %>%
  mutate(Training = "massed")

r0_grid_spaced <- tidyr::expand_grid(
  Genotype = genotypes,
  Block    = c("1", "2", "3", "4", "5")
) %>%
  mutate(Training = "spaced")

r0_grid <- bind_rows(r0_grid_massed, r0_grid_spaced) %>%
  mutate(
    stimulus0 = 0,
    stimulus  = 0,
    Training  = factor(Training, levels = levels(df_all$Training)),
    Block     = factor(Block,    levels = levels(df_all$Block)),
    Genotype  = factor(Genotype, levels = levels(df_all$Genotype))
  )

r0_raw <- posterior_epred(
  fit_joint,
  newdata    = r0_grid,
  re_formula = NA,
  ndraws     = 1000
)

r0_long_meta <- r0_grid %>%
  mutate(row_id = row_number()) %>%
  select(row_id, Genotype, Training, Block)

r0_draws <- as_tibble(t(r0_raw)) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(-row_id, names_to = ".draw", values_to = "R0") %>%
  mutate(.draw = as.integer(str_remove(.draw, "V"))) %>%
  left_join(r0_long_meta, by = "row_id") %>%
  select(Genotype, Training, Block, .draw, R0)


# ------------------------------------------------------------------------------
# Test-pause jump per protocol
# ------------------------------------------------------------------------------
massed_jump <- r0_draws %>%
  filter(Training == "massed") %>%
  select(Genotype, .draw, Block, R0) %>%
  pivot_wider(names_from = Block, values_from = R0, names_prefix = "R0_block") %>%
  mutate(
    Training  = "massed",
    jump_test = R0_block2 - R0_block1
  ) %>%
  select(Genotype, Training, .draw, jump_test)

spaced_jump <- r0_draws %>%
  filter(Training == "spaced") %>%
  select(Genotype, .draw, Block, R0) %>%
  pivot_wider(names_from = Block, values_from = R0, names_prefix = "R0_block") %>%
  mutate(
    Training      = "spaced",
    jump_test     = R0_block5 - R0_block4,
    jump_1_to_2   = R0_block2 - R0_block1,
    jump_2_to_3   = R0_block3 - R0_block2,
    jump_3_to_4   = R0_block4 - R0_block3
  )

# Per-protocol test-pause jump summaries
test_jump_draws <- bind_rows(
  massed_jump %>% select(Genotype, Training, .draw, jump_test),
  spaced_jump %>% select(Genotype, Training, .draw, jump_test) %>% mutate(Training = "spaced")
)

test_jump_summary <- test_jump_draws %>%
  group_by(Genotype, Training) %>%
  summarise(
    median = median(jump_test),
    low    = quantile(jump_test, 0.025),
    high   = quantile(jump_test, 0.975),
    .groups = "drop"
  )

print(test_jump_summary)

write.csv(
  test_jump_summary,
  file.path(save_results_dir, "between_block_jump_test_pause.csv"),
  row.names = FALSE
)


# Spaced vs massed contrast on the test-pause jump
test_jump_diff <- test_jump_draws %>%
  pivot_wider(names_from = Training, values_from = jump_test) %>%
  mutate(diff = spaced - massed)

test_jump_diff_summary <- test_jump_diff %>%
  group_by(Genotype) %>%
  summarise(
    diff_median        = median(diff),
    diff_low           = quantile(diff, 0.025),
    diff_high          = quantile(diff, 0.975),
    p_spaced_lt_massed = mean(diff < 0),
    p_spaced_gt_massed = mean(diff > 0),
    .groups = "drop"
  )

print(test_jump_diff_summary)

write.csv(
  test_jump_diff_summary,
  file.path(save_results_dir, "between_block_jump_diff_spaced_vs_massed.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Within-spaced-training jumps (context for the test-pause jump)
# ------------------------------------------------------------------------------
spaced_all_jumps <- spaced_jump %>%
  select(Genotype, .draw, jump_1_to_2, jump_2_to_3, jump_3_to_4, jump_test) %>%
  pivot_longer(
    cols = starts_with("jump_"),
    names_to = "transition",
    values_to = "jump"
  ) %>%
  mutate(
    transition = factor(
      transition,
      levels = c("jump_1_to_2", "jump_2_to_3", "jump_3_to_4", "jump_test"),
      labels = c("B1->B2", "B2->B3", "B3->B4", "B4->Test")
    )
  )

spaced_all_jumps_summary <- spaced_all_jumps %>%
  group_by(Genotype, transition) %>%
  summarise(
    median = median(jump),
    low    = quantile(jump, 0.025),
    high   = quantile(jump, 0.975),
    .groups = "drop"
  )

print(spaced_all_jumps_summary)

write.csv(
  spaced_all_jumps_summary,
  file.path(save_results_dir, "spaced_within_training_jumps.csv"),
  row.names = FALSE
)


# Plots --------------------------------------------------------------------
p_test_jump <- test_jump_summary %>%
  ggplot(aes(x = Training, y = median,
             ymin = low, ymax = high,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "R0(test) - R0(last training block)\non probability scale",
    title = "Inter-block jump at the memory test pause"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrast5_test_pause_jump.png"),
  p_test_jump, width = 12, height = 4.5, dpi = 300, bg = "white"
)


p_test_jump_diff <- test_jump_diff_summary %>%
  ggplot(aes(x = Genotype, y = diff_median,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x     = NULL,
    y     = "Test-pause jump difference (spaced - massed)\nNegative = spaced has smaller jump = better retention",
    title = "Test-pause recovery: spaced vs massed"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrast5_test_pause_jump_diff.png"),
  p_test_jump_diff, width = 8, height = 5, dpi = 300, bg = "white"
)


p_spaced_jumps <- spaced_all_jumps_summary %>%
  ggplot(aes(x = transition, y = median, ymin = low, ymax = high,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.7) +
  geom_line(aes(group = Genotype), linewidth = 0.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 12) +
  labs(
    x     = "Inter-block transition (spaced)",
    y     = "R0(later block) - R0(earlier block)",
    title = "Spaced: between-block recovery across all training transitions"
  ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(save_fig_dir, "contrast5_spaced_all_jumps.png"),
  p_spaced_jumps, width = 14, height = 5, dpi = 300, bg = "white"
)


# ==============================================================================
# CONTRAST (4): RE-HABITUATION SAVINGS  (k_test vs k_block1)
# ==============================================================================
# Savings = trained fish habituates faster than a naive one.
# Per protocol: k_test - k_block1. Then compare across protocols.
# CAVEAT: test blocks have only ~8-9 stimuli, so k_test is weakly identified.
# Treat as supporting analysis, not headline.
#
# We extract logk by predicting at two carefully chosen stimulus points within
# a (Genotype, Training, Block) cell. The model is:
#   p(stim) = inv_logit(A) + (inv_logit(R0) - inv_logit(A)) * exp(-k * stim)
# At stim = 0:           p_0 = inv_logit(R0)
# At stim = LARGE:       p_inf = inv_logit(A)
# At any intermediate s: p_s   = inv_logit(A) + (inv_logit(R0)-inv_logit(A)) * exp(-k*s)
#
# Solve for k from p_s:
#   exp(-k*s) = (p_s - p_inf) / (p_0 - p_inf)
#   k = -log((p_s - p_inf) / (p_0 - p_inf)) / s
#
# We pick s = 5 (or block-specific) so the prediction is firmly in the decay
# region and the algebra is numerically stable.
# ==============================================================================

# Choose s per block based on block length
s_per_block <- function(block_label, training) {
  if (training == "massed" && block_label == "1") return(50)
  if (training == "massed" && block_label == "2") return(4)
  if (training == "spaced" && block_label == "1") return(20)
  if (training == "spaced" && block_label == "5") return(4)
  return(20)
}

# We only need k_block1 and k_test per protocol
k_cells <- bind_rows(
  tibble(Training = "massed", Block = "1", role = "k_train"),
  tibble(Training = "massed", Block = "2", role = "k_test"),
  tibble(Training = "spaced", Block = "1", role = "k_train"),
  tibble(Training = "spaced", Block = "5", role = "k_test")
)

# Build a grid with three rows per (Genotype, Training, Block):
# stim=0, stim=s, stim=LARGE
k_grid <- tidyr::expand_grid(
  Genotype = genotypes,
  cell_id  = seq_len(nrow(k_cells))
) %>%
  left_join(k_cells %>% mutate(cell_id = row_number()), by = "cell_id") %>%
  rowwise() %>%
  mutate(
    s = s_per_block(Block, Training)
  ) %>%
  ungroup() %>%
  tidyr::crossing(point = c("p0", "ps", "pinf")) %>%
  mutate(
    stimulus0 = case_when(
      point == "p0"   ~ 0,
      point == "ps"   ~ s,
      point == "pinf" ~ LARGE
    ),
    stimulus = stimulus0,
    Training = factor(Training, levels = levels(df_all$Training)),
    Block    = factor(Block,    levels = levels(df_all$Block)),
    Genotype = factor(Genotype, levels = levels(df_all$Genotype))
  )

k_raw <- posterior_epred(
  fit_joint,
  newdata    = k_grid,
  re_formula = NA,
  ndraws     = 1000
)

k_long_meta <- k_grid %>%
  mutate(row_id = row_number()) %>%
  select(row_id, Genotype, Training, Block, role, point, s)

k_draws_raw <- as_tibble(t(k_raw)) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(-row_id, names_to = ".draw", values_to = "p") %>%
  mutate(.draw = as.integer(str_remove(.draw, "V"))) %>%
  left_join(k_long_meta, by = "row_id")

# Pivot so each (Genotype, Training, Block, role, .draw) has p0, ps, pinf
# and compute k.
k_solved <- k_draws_raw %>%
  select(Genotype, Training, Block, role, s, .draw, point, p) %>%
  pivot_wider(names_from = point, values_from = p) %>%
  mutate(
    # Numerical guards: avoid log of negatives or zero
    num   = ps - pinf,
    denom = p0 - pinf,
    valid = denom > 1e-4 & num > 1e-4,
    k     = ifelse(valid, -log(num / denom) / s, NA_real_)
  )

# Per-cell summary
k_cell_summary <- k_solved %>%
  filter(!is.na(k)) %>%
  group_by(Genotype, Training, role) %>%
  summarise(
    median = median(k),
    low    = quantile(k, 0.025),
    high   = quantile(k, 0.975),
    n_valid = n(),
    .groups = "drop"
  )

print(k_cell_summary)

write.csv(
  k_cell_summary,
  file.path(save_results_dir, "k_per_cell_train_and_test.csv"),
  row.names = FALSE
)


# Savings per protocol: k_test - k_train, per (Genotype, Training)
savings_draws <- k_solved %>%
  select(Genotype, Training, role, .draw, k) %>%
  pivot_wider(names_from = role, values_from = k) %>%
  filter(!is.na(k_train), !is.na(k_test)) %>%
  mutate(savings = k_test - k_train)

savings_summary <- savings_draws %>%
  group_by(Genotype, Training) %>%
  summarise(
    median = median(savings),
    low    = quantile(savings, 0.025),
    high   = quantile(savings, 0.975),
    .groups = "drop"
  )

print(savings_summary)

write.csv(
  savings_summary,
  file.path(save_results_dir, "rehabituation_savings_per_protocol.csv"),
  row.names = FALSE
)


# Spaced vs massed savings contrast
savings_diff <- savings_draws %>%
  select(Genotype, Training, .draw, savings) %>%
  pivot_wider(names_from = Training, values_from = savings) %>%
  filter(!is.na(spaced), !is.na(massed)) %>%
  mutate(diff = spaced - massed)

savings_diff_summary <- savings_diff %>%
  group_by(Genotype) %>%
  summarise(
    diff_median        = median(diff),
    diff_low           = quantile(diff, 0.025),
    diff_high          = quantile(diff, 0.975),
    p_spaced_gt_massed = mean(diff > 0),  # P(spaced has MORE savings = stronger memory at rate level)
    p_spaced_lt_massed = mean(diff < 0),
    .groups = "drop"
  )

print(savings_diff_summary)

write.csv(
  savings_diff_summary,
  file.path(save_results_dir, "rehabituation_savings_diff_spaced_vs_massed.csv"),
  row.names = FALSE
)



# Plots --------------------------------------------------------------------
p_savings <- savings_summary %>%
  ggplot(aes(x = Training, y = median,
             ymin = low, ymax = high,
             color = Genotype)) +
  facet_wrap(~ Genotype, ncol = 5) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "k_test - k_block1\nPositive = re-habituates faster than naive",
    title = "Re-habituation savings (CAVEAT: k_test weakly identified)"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrast4_savings_per_protocol.png"),
  p_savings, width = 12, height = 4.5, dpi = 300, bg = "white"
)


p_savings_diff <- savings_diff_summary %>%
  ggplot(aes(x = Genotype, y = diff_median,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x     = NULL,
    y     = "Savings difference (spaced - massed)\nPositive = spaced has stronger rate memory",
    title = "Re-habituation savings: spaced vs massed"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(save_fig_dir, "contrast4_savings_diff_spaced_vs_massed.png"),
  p_savings_diff, width = 8, height = 5, dpi = 300, bg = "white"
)


# ==============================================================================
# UNIFIED COMPARISON TABLE  --- read all three contrasts side by side
# ==============================================================================
# Per genotype, gather the three "spaced - massed" contrasts so you can see
# whether they agree.
# ==============================================================================

all_contrasts <- bind_rows(
  abs_test_diff_summary    %>% mutate(contrast = "1_absolute_test_response"),
  test_jump_diff_summary   %>% mutate(contrast = "5_test_pause_jump"),
  savings_diff_summary     %>% mutate(contrast = "4_rehabituation_savings")
) %>%
  select(contrast, Genotype, diff_median, diff_low, diff_high,
         p_spaced_lt_massed, p_spaced_gt_massed)

print(all_contrasts)

write.csv(
  all_contrasts,
  file.path(save_results_dir, "ALL_contrasts_spaced_vs_massed.csv"),
  row.names = FALSE
)


# Unified figure: rows = contrasts, x = Genotype
p_all_contrasts <- all_contrasts %>%
  ggplot(aes(x = Genotype, y = diff_median,
             ymin = diff_low, ymax = diff_high,
             color = Genotype)) +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y") +
  geom_pointrange(linewidth = 0.8, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_pubr(base_size = 13) +
  labs(
    x     = NULL,
    y     = "Spaced - massed contrast (sign varies by contrast)",
    title = "All three memory contrasts side-by-side"
  ) +
  theme(legend.position = "none", strip.text = element_text(size = 11))

ggsave(
  file.path(save_fig_dir, "ALL_contrasts_spaced_vs_massed.png"),
  p_all_contrasts, width = 10, height = 12, dpi = 300, bg = "white"
)


# ==============================================================================
# Done
# ==============================================================================
message("Contrast analysis complete.")
message("Headline: contrast (1) absolute_test_response_diff_spaced_vs_massed.csv")
message("Supporting: contrast (5) between_block_jump, contrast (4) savings.")
message("All three: ALL_contrasts_spaced_vs_massed.csv")