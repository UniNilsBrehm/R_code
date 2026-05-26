
library(readr)
library(dplyr)
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
library(mgcv)
library(emmeans)
library(gratia)

# ==============================================================================
# Paths
# ==============================================================================
source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_Massed_vs_Spaced//utils.R")

# base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_Joint_SpacedVsMassed/"
base_dir <- "D:/WorkingData/Susana/NLME/DarkFlash_Joint_SpacedVsMassed/"

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

# # Remove all non-responders to stimulus 1 in Block 1
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
# GAM
# ==============================================================================
# df_all is created exactly as in your Bayesian NLME script
df_gam <- df_all %>%
  mutate(
    move      = as.integer(move),
    Training  = factor(Training, levels = c("massed", "spaced")),
    Block     = factor(Block),
    BlockRole = factor(BlockRole, levels = c("training", "test")),
    Genotype  = factor(Genotype),
    animal    = factor(animal),
    cell      = interaction(Genotype, Training, Block, drop = TRUE),
    
    # KEY CHANGE
    stimulus_log = log1p(stimulus0)
  )

# Frequentist GAM
# - separate habituation curve per Genotype x Training x Block
# - parametric Genotype * Training * Block differences
# - random intercept per animal
# gam_joint <- bam(
#   move ~
#     Genotype * Training * Block +
#     s(stimulus_log, by = cell, k = 5) +
#     s(animal, bs = "re"),
#   
#   data = df_gam,
#   family = binomial(link = "logit"),
#   method = "fREML",
#   discrete = TRUE
# )

# Simpler
# gam_joint <- bam(
#   move ~
#     cell +
#     s(stimulus_log, by = cell, k = 5) +
#     s(animal, bs = "re"),
#   
#   data = df_gam,
#   family = binomial(link = "logit"),
#   method = "fREML",
#   discrete = TRUE
# )

# Partial pooling / shared shape structure:
# because habituation almost certainly has:
#   - a shared biological structure
#   - plus genotype/training modifications
gam_joint <- bam(
  move ~
    cell +
    s(stimulus_log, k = 5) +
    s(stimulus_log, by = cell, k = 5) +
    s(animal, bs = "re"),
  
  data = df_gam,
  family = binomial(link = "logit"),
  method = "fREML",
  discrete = TRUE
)

summary(gam_joint)
gam.check(gam_joint)

# ==============================================================================
# Plot
# Build prediction grid: per Training x Block x Genotype
new_data_gam <- df_gam %>%
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
    Training = factor(Training, levels = levels(df_gam$Training)),
    Block    = factor(Block,    levels = levels(df_gam$Block)),
    Genotype = factor(Genotype, levels = levels(df_gam$Genotype)),
    animal   = df_gam$animal[1],
    cell     = interaction(Genotype, Training, Block, drop = TRUE)
  )

new_data_gam <- new_data_gam %>%
  mutate(
    stimulus_log = log1p(stimulus0)
  )

# Frequentist GAM predictions on link scale
pred_gam <- predict(
  gam_joint,
  newdata = new_data_gam,
  type = "link",
  se.fit = TRUE,
  exclude = "s(animal)"
)

new_data_gam <- new_data_gam %>%
  mutate(
    CI_low  = plogis(pred_gam$fit - 1.96 * pred_gam$se.fit),
    fit     = plogis(pred_gam$fit),
    CI_high = plogis(pred_gam$fit + 1.96 * pred_gam$se.fit)
  )

# Raw observed response probability per stimulus
raw_summary_gam <- df_gam %>%
  group_by(Training, Block, Genotype, stimulus0) %>%
  summarise(
    p_move = mean(move, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(stimulus = stimulus0)

# Massed plot
p_massed_curves_gam <- ggplot(
  new_data_gam %>% filter(Training == "massed"),
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_gam %>% filter(Training == "massed"),
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
    title = "Frequentist GAM: Massed training"
  ) +
  theme(legend.position = "none")

# Spaced plot
p_spaced_curves_gam <- ggplot(
  new_data_gam %>% filter(Training == "spaced"),
  aes(x = stimulus, color = Genotype, fill = Genotype)
) +
  facet_grid(Genotype ~ Block, scales = "free_x") +
  geom_point(
    data = raw_summary_gam %>% filter(Training == "spaced"),
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
    title = "Frequentist GAM: Spaced training"
  ) +
  theme(legend.position = "top")

print(p_spaced_curves_gam)
print(p_massed_curves_gam)

