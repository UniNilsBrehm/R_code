# ==============================================================================
# Compare Bayesian hierarchical model vs aggregate NLS
# ==============================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(brms)
library(matrixStats)   # for colLogSumExps
library(readr)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(loo)
library(DHARMa)
library(bayesplot)

source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/nlme_utils.R")
source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/plot_utils.R")

# source("C:/UniFreiburg/Code/R_code/susana/nlme_utils.R")
# source("C:/UniFreiburg/Code/R_code/susana/plot_utils.R")

# source("D:/Behavior_Data/R_code/susana/nlme_utils.R")
# source("D:/Behavior_Data/R_code/susana/plot_utils.R")

base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <- "D:/WorkingData/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <- "D:/Behavior_Data/DarkFlash_ISI90s_2Blocks"

file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

var_name = 'response_prob'
col_name = 'move'

save_fig_dir = file.path(base_dir, "figs", "nlme", var_name)
save_results_dir = file.path(base_dir, "results", "nlme", var_name)

# Create directories if they do not exist
dir.create(save_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Load and prepare data
# ==============================================================================
res <- load_data_darkflash_60s(file_dir, move_th = 0, take_peak = 0)

df_final     <- res$df_final

df_resp <- df_final %>%
  mutate(
    stimulus = as.numeric(stimulus),
    stimulus0 = stimulus - 1,
    Block = factor(Block),
    Genotype = factor(Genotype),
    Video = factor(Video),
    Well = factor(Well),
    animal = interaction(Video, Well, drop = TRUE)
    
  )


# ==============================================================================
# The Model
# ==============================================================================
model_logit <- bf(
  move ~ A + (R0 - A) * exp(-exp(logk) * stimulus0),
  
  A     ~ 0 + Genotype * Block + (1 | animal), # final/asymptotic response level on logit scale
  R0    ~ 0 + Genotype * Block + (1 | animal), # initial response level on logit scale
  logk ~ 0 + Genotype * Block + (1 | animal),  # log decay rate: k = exp(logk), half_life = log(2) / exp(logk)
  
  nl = TRUE
)

# on the prob scale
model_identity <- bf(
  move ~ inv_logit(A) +
    (inv_logit(R0) - inv_logit(A)) *
    exp(-exp(logk) * stimulus0),
  
  A    ~ 0 + Genotype * Block + (1 | animal),
  R0   ~ 0 + Genotype * Block + (1 | animal),
  logk ~ 0 + Genotype * Block + (1 | animal),
  
  nl = TRUE
)

# ==============================================================================
# The Priors
# ==============================================================================
priors <- c(
  prior(normal(0, 1.0), class = "b", nlpar = "A"),
  prior(normal(1.5, 1), class = "b", nlpar = "R0"),
  prior(normal(-3, 1), class = "b", nlpar = "logk"),
  
  prior(exponential(2), class = "sd", nlpar = "A"),
  prior(exponential(2), class = "sd", nlpar = "R0"),
  prior(exponential(2), class = "sd", nlpar = "logk")
)

# ==============================================================================
set.seed(42)
compare_dir <- file.path(base_dir, "results", "model_comparison")
compare_fig_dir <- file.path(base_dir, "figs", "model_comparison")
dir.create(compare_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(compare_fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 0. Helper: fit your NLS on an aggregated dataset (same as your fit_simple_exp)
# ------------------------------------------------------------------------------
fit_simple_exp_safe <- function(dat) {
  tryCatch(
    nls(
      response_prob ~ A + (R0 - A) * exp(-exp(logk) * stimulus0),
      data = dat, weights = n_animals,
      start = list(
        A    = max(0.01, min(dat$response_prob, na.rm = TRUE)),
        R0   = min(0.99, max(dat$response_prob, na.rm = TRUE)),
        logk = log(0.1)
      ),
      algorithm = "port",
      lower = c(A = 0.001, R0 = 0.001, logk = log(1e-4)),
      upper = c(A = 0.999, R0 = 0.999, logk = log(10)),
      control = nls.control(maxiter = 500, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )
}

aggregate_to_probs <- function(df) {
  df %>%
    mutate(stimulus = as.numeric(stimulus),
           stimulus0 = stimulus - 1,
           move = as.integer(move)) %>%
    group_by(Genotype, Block, stimulus, stimulus0) %>%
    summarise(n_animals = n(),
              response_prob = mean(move, na.rm = TRUE),
              .groups = "drop")
}

fit_nls_per_cell <- function(df_agg) {
  keys <- df_agg %>% distinct(Genotype, Block) %>% arrange(Genotype, Block)
  fits <- df_agg %>% group_by(Genotype, Block) %>% group_split() %>%
    map(fit_simple_exp_safe)
  map2_dfr(fits, seq_along(fits), function(mod, i) {
    if (is.null(mod)) {
      tibble(Genotype = keys$Genotype[i], Block = keys$Block[i],
             A = NA_real_, R0 = NA_real_, logk = NA_real_, converged = FALSE)
    } else {
      cc <- coef(mod)
      tibble(Genotype = keys$Genotype[i], Block = keys$Block[i],
             A = unname(cc["A"]), R0 = unname(cc["R0"]),
             logk = unname(cc["logk"]),
             converged = mod$convInfo$isConv %||% TRUE)
    }
  })
}

# ------------------------------------------------------------------------------
# 1. Hold out 20% of animals, stratified by Genotype x Block
# ------------------------------------------------------------------------------
holdout_animals <- df_resp %>%
  distinct(animal, Genotype, Block) %>%
  group_by(Genotype, Block) %>%
  slice_sample(prop = 0.20) %>%
  pull(animal) %>% unique()

df_train <- df_resp %>% filter(!animal %in% holdout_animals)
df_test  <- df_resp %>% filter( animal %in% holdout_animals)

cat("Train animals:", n_distinct(df_train$animal),
    " Test animals:", n_distinct(df_test$animal), "\n")

# ------------------------------------------------------------------------------
# 2. Refit Bayesian model on training data
# ------------------------------------------------------------------------------
fit_bayes_train <- brm(
  formula = model_identity,
  data = df_train,
  family = bernoulli(link = "identity"),
  prior = priors,
  backend = "cmdstanr",
  chains = 2, cores = 2, threads = threading(6),
  iter = 1000, warmup = 500, seed = 42,
  control = list(adapt_delta = 0.90, max_treedepth = 10),
  save_pars = save_pars(all = TRUE)
)

# Per-trial log-likelihood on held-out animals
# allow_new_levels + sample_new_levels integrates over the animal RE prior
ll_bayes <- log_lik(
  fit_bayes_train, newdata = df_test,
  allow_new_levels = TRUE, sample_new_levels = "gaussian"
)
# log mean exp over draws per trial = per-trial lpd
lpd_bayes_per_trial <- colLogSumExps(ll_bayes) - log(nrow(ll_bayes))
elpd_bayes <- sum(lpd_bayes_per_trial)

# ------------------------------------------------------------------------------
# 3. Fit NLS on training set and score on held-out trials
# ------------------------------------------------------------------------------
df_train_agg  <- aggregate_to_probs(df_train)
nls_params_tr <- fit_nls_per_cell(df_train_agg)

df_test_scored <- df_test %>%
  mutate(stimulus0 = stimulus - 1, move = as.integer(move)) %>%
  left_join(nls_params_tr, by = c("Genotype", "Block")) %>%
  mutate(
    p_hat = A + (R0 - A) * exp(-exp(logk) * stimulus0),
    p_hat = pmin(pmax(p_hat, 1e-6), 1 - 1e-6),
    ll    = move * log(p_hat) + (1 - move) * log(1 - p_hat)
  )
lpd_nls_per_trial <- df_test_scored$ll
elpd_nls <- sum(lpd_nls_per_trial, na.rm = TRUE)

# SE of the difference (paired, per trial)
delta <- lpd_bayes_per_trial - lpd_nls_per_trial
elpd_diff    <- sum(delta)
elpd_diff_se <- sqrt(length(delta)) * sd(delta)

elpd_table <- tibble(
  model    = c("Bayesian hierarchical", "Aggregate NLS", "Difference (Bayes - NLS)"),
  elpd     = c(elpd_bayes, elpd_nls, elpd_diff),
  se_diff  = c(NA, NA, elpd_diff_se)
)
print(elpd_table)
write.csv(elpd_table, file.path(compare_dir, "elpd_heldout_animals.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Calibration on held-out trials
# ------------------------------------------------------------------------------
# Bayesian per-trial mean predicted probability (marginalized over draws)
p_bayes_test <- posterior_epred(
  fit_bayes_train, newdata = df_test,
  allow_new_levels = TRUE, sample_new_levels = "gaussian"
) %>% colMeans()

calib_df <- bind_rows(
  tibble(method = "Bayesian", p_hat = p_bayes_test, y = as.integer(df_test$move)),
  tibble(method = "NLS",      p_hat = df_test_scored$p_hat, y = df_test_scored$move)
) %>%
  mutate(bin = ntile(p_hat, 10)) %>%
  group_by(method, bin) %>%
  summarise(p_pred = mean(p_hat), p_obs = mean(y), n = n(), .groups = "drop")

p_calib <- ggplot(calib_df, aes(p_pred, p_obs, color = method, size = n)) +
  geom_abline(linetype = "dashed") +
  geom_point() + geom_line(aes(group = method), linewidth = 0.7) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  scale_size_continuous(range = c(2, 6)) +
  theme_pubr(base_size = 14) +
  labs(x = "Predicted probability", y = "Observed proportion",
       title = "Calibration on held-out animals")
ggsave(file.path(compare_fig_dir, "calibration_heldout.png"),
       p_calib, width = 7, height = 6, dpi = 300, bg = "white")

# ------------------------------------------------------------------------------
# 5. Cluster bootstrap of NLS over animals -> honest CIs for k, half_life
#    Use the FULL data (not just training) for the parameter comparison
# ------------------------------------------------------------------------------
animals_by_cell <- df_resp %>%
  distinct(animal, Genotype, Block) %>%
  group_by(Genotype, Block) %>%
  summarise(animals = list(unique(animal)), .groups = "drop")

n_boot <- 500   # bump to 1000+ for final figures

boot_one <- function(b) {
  resampled <- animals_by_cell %>%
    rowwise() %>%
    mutate(sampled = list(sample(animals, length(animals), replace = TRUE))) %>%
    ungroup() %>%
    select(Genotype, Block, sampled) %>%
    unnest(sampled) %>%
    rename(animal = sampled) %>%
    left_join(df_resp, by = c("Genotype", "Block", "animal"),
              relationship = "many-to-many")
  agg  <- aggregate_to_probs(resampled)
  pars <- fit_nls_per_cell(agg)
  pars %>% mutate(boot = b)
}

boot_params <- map_dfr(seq_len(n_boot), boot_one)

nls_boot_summary <- boot_params %>%
  mutate(k = exp(logk), half_life = log(2) / k) %>%
  group_by(Genotype, Block) %>%
  summarise(
    n_fail      = sum(is.na(logk)),
    pct_fail    = mean(is.na(logk)) * 100,
    k_median    = median(k, na.rm = TRUE),
    k_low       = quantile(k, 0.025, na.rm = TRUE),
    k_high      = quantile(k, 0.975, na.rm = TRUE),
    hl_median   = median(half_life, na.rm = TRUE),
    hl_low      = quantile(half_life, 0.025, na.rm = TRUE),
    hl_high     = quantile(half_life, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(method = "NLS (cluster bootstrap)")

# Bayesian summary from your full-data fit_model
bayes_k_summary <- k_summary %>%
  transmute(Genotype, Block,
            k_median, k_low, k_high,
            hl_median = half_life_median,
            hl_low    = half_life_low,
            hl_high   = half_life_high,
            method    = "Bayesian hierarchical",
            n_fail = 0, pct_fail = 0)

k_compare_tbl <- bind_rows(bayes_k_summary, nls_boot_summary)
print(k_compare_tbl)
write.csv(k_compare_tbl, file.path(compare_dir, "k_halflife_method_comparison.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. Side-by-side intervals for k and half_life
# ------------------------------------------------------------------------------
p_k_methods <- ggplot(k_compare_tbl,
                      aes(x = interaction(Genotype, Block, sep = " / "),
                          y = k_median, ymin = k_low, ymax = k_high,
                          color = method)) +
  geom_pointrange(position = position_dodge(width = 0.5), linewidth = 0.8) +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(x = NULL, y = "Habituation rate k",
       title = "k: Bayesian CrI vs NLS cluster-bootstrap CI")
ggsave(file.path(compare_fig_dir, "k_intervals_methods.png"),
       p_k_methods, width = 9, height = 6, dpi = 300, bg = "white")

p_hl_methods <- ggplot(k_compare_tbl,
                       aes(x = interaction(Genotype, Block, sep = " / "),
                           y = hl_median, ymin = hl_low, ymax = hl_high,
                           color = method)) +
  geom_pointrange(position = position_dodge(width = 0.5), linewidth = 0.8) +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(x = NULL, y = "Half-life (stimuli)",
       title = "Half-life: Bayesian CrI vs NLS cluster-bootstrap CI")
ggsave(file.path(compare_fig_dir, "halflife_intervals_methods.png"),
       p_hl_methods, width = 9, height = 6, dpi = 300, bg = "white")

# ------------------------------------------------------------------------------
# 7. Curve overlay: Bayesian band vs NLS bootstrap band vs aggregate points
# ------------------------------------------------------------------------------
new_grid <- df_resp %>%
  group_by(Genotype, Block) %>%
  summarise(stim_min = min(stimulus), stim_max = max(stimulus), .groups = "drop") %>%
  reframe(stimulus = seq(stim_min, stim_max, length.out = 100),
          .by = c(Genotype, Block)) %>%
  mutate(stimulus0 = stimulus - 1,
         Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
         Block    = factor(Block,    levels = levels(df_resp$Block)))

bayes_curves <- fitted(fit_model, newdata = new_grid, re_formula = NA,
                       probs = c(0.025, 0.975)) %>%
  as.data.frame() %>%
  bind_cols(new_grid) %>%
  rename(fit = Estimate, lo = Q2.5, hi = Q97.5) %>%
  mutate(method = "Bayesian hierarchical")

# NLS bootstrap curves: for each draw, evaluate the curve on the grid
nls_boot_curves <- boot_params %>%
  filter(!is.na(logk)) %>%
  inner_join(new_grid, by = c("Genotype", "Block"),
             relationship = "many-to-many") %>%
  mutate(fit = A + (R0 - A) * exp(-exp(logk) * stimulus0)) %>%
  group_by(Genotype, Block, stimulus) %>%
  summarise(lo  = quantile(fit, 0.025),
            hi  = quantile(fit, 0.975),
            fit = median(fit), .groups = "drop") %>%
  mutate(method = "NLS (bootstrap)")

curves <- bind_rows(bayes_curves, nls_boot_curves)

p_curves <- ggplot(curves, aes(stimulus, fit,
                               color = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(data = df_prob_agg,
             aes(stimulus, response_prob), inherit.aes = FALSE,
             alpha = 0.4, size = 1) +
  facet_grid(Block ~ Genotype) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_pubr(base_size = 14) +
  labs(x = "Stimulus", y = "Response probability",
       title = "Fitted curves with 95% intervals")
ggsave(file.path(compare_fig_dir, "curves_methods_overlay.png"),
       p_curves, width = 14, height = 8, dpi = 300, bg = "white")