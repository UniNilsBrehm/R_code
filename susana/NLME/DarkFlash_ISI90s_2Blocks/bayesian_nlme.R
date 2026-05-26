# ==============================================================================
# Bayesian nonlinear mixed model for zebrafish dark-flash habituation
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

source("C:/Users/NilsPC/Desktop/Susana/NLME/utils.R")

base_dir <- "C:/Users/NilsPC/Desktop/Susana/NLME/"
file_dir <- "C:/Users/NilsPC/Desktop/Susana/NLME/SPZ_ISI60_removed_non_responders_2stimuli.csv"

# ==============================================================================
# 1. Load and prepare data
# ==============================================================================

res <- load_data_darkflash_60s(file_dir, move_th = 0, take_peak = 0)

df_final     <- res$df_final
df_final_sub <- res$df_final_sub

df_brm <- df_final_sub %>%
  mutate(
    animal = interaction(Video, Well, drop = TRUE),
    stimulus = as.numeric(stimulus),
    stimulus0 = stimulus - 1,
    Block = factor(Block),
    Genotype = factor(Genotype),
    y = log(max_peak + 1)
  ) %>%
  filter(
    is.finite(y),
    !is.na(stimulus),
    !is.na(stimulus0),
    !is.na(Block),
    !is.na(Genotype),
    !is.na(animal)
  )

n_per_genotype <- df_brm %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

print(n_per_genotype)

# ==============================================================================
# 2. Bayesian nonlinear model
# ==============================================================================

bf_hab <- bf(
  y ~ Asym + (R0 - Asym) * exp(-exp(lrc) * stimulus0),
  
  Asym ~ Genotype * Block + (1 | animal),
  R0   ~ Genotype * Block + (1 | animal),
  lrc  ~ Genotype * Block + (1 | animal),
  
  nl = TRUE
)

# ==============================================================================
# 3. Starting values for priors
# ==============================================================================

m_start <- nls(
  y ~ SSasymp(stimulus0, Asym, R0, lrc),
  data = df_brm
)

start_vals <- coef(m_start)
print(start_vals)

Asym_start <- unname(start_vals["Asym"])
R0_start   <- unname(start_vals["R0"])
lrc_start  <- unname(start_vals["lrc"])

# ==============================================================================
# 4. Priors
# ==============================================================================

priors_hab <- c(
  set_prior(paste0("normal(", Asym_start, ", 1)"),
            nlpar = "Asym", class = "b", coef = "Intercept"),
  set_prior("normal(0, 0.5)",
            nlpar = "Asym", class = "b"),
  set_prior("exponential(1)",
            nlpar = "Asym", class = "sd"),
  
  set_prior(paste0("normal(", R0_start, ", 1)"),
            nlpar = "R0", class = "b", coef = "Intercept"),
  set_prior("normal(0, 0.5)",
            nlpar = "R0", class = "b"),
  set_prior("exponential(1)",
            nlpar = "R0", class = "sd"),
  
  set_prior(paste0("normal(", lrc_start, ", 1)"),
            nlpar = "lrc", class = "b", coef = "Intercept"),
  set_prior("normal(0, 0.5)",
            nlpar = "lrc", class = "b"),
  set_prior("exponential(1)",
            nlpar = "lrc", class = "sd"),
  
  set_prior("exponential(1)", class = "sigma")
)

get_prior(
  formula = bf_hab,
  data = df_brm,
  family = gaussian()
)

# ==============================================================================
# 5. Test model
# ==============================================================================

options(mc.cores = parallel::detectCores())

fit_test <- brm(
  formula = bf_hab,
  data = df_brm,
  family = gaussian(),
  prior = priors_hab,
  
  backend = "cmdstanr",
  
  chains = 2,
  cores = 2,
  threads = threading(4),
  
  iter = 1000,
  warmup = 500,
  
  seed = 42,
  
  control = list(
    adapt_delta = 0.90,
    max_treedepth = 10
  )
)

summary(fit_test)

# ==============================================================================
# 6. Final model
# ==============================================================================

fit_brm <- brm(
  formula = bf_hab,
  data = df_brm,
  family = gaussian(),
  prior = priors_hab,
  
  backend = "cmdstanr",
  
  chains = 4,
  cores = 4,
  threads = threading(4),
  
  iter = 4000,
  warmup = 1000,
  seed = 42,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  )
)

summary(fit_brm)
plot(fit_brm)
pp_check(fit_brm)

saveRDS(
  fit_brm,
  file = file.path(base_dir, "BRMS_habituation_model.rds")
)

# ==============================================================================
# 7. Plot habituation curves
# ==============================================================================

plot_habituation_brm <- function(df_brm, fit,
                                 label = "Peak distance moved",
                                 raw_var = "max_peak",
                                 save_path = NULL) {
  
  raw_summary <- df_brm %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      y_raw = mean(.data[[raw_var]], na.rm = TRUE),
      .groups = "drop"
    )
  
  new_data <- df_brm %>%
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
      stimulus0 = stimulus - 1,
      Genotype = factor(Genotype, levels = levels(df_brm$Genotype)),
      Block = factor(Block, levels = levels(df_brm$Block))
    )
  
  mu_draws <- posterior_epred(
    fit,
    newdata = new_data,
    re_formula = NA,
    summary = FALSE
  )
  
  draws_df <- as_draws_df(fit)
  sigma_draws <- draws_df$sigma
  
  pred_orig <- sweep(
    mu_draws,
    1,
    sigma_draws^2 / 2,
    FUN = "+"
  )
  
  pred_orig <- exp(pred_orig) - 1
  
  pred_data <- new_data %>%
    mutate(
      fit = apply(pred_orig, 2, median),
      CI_low = apply(pred_orig, 2, quantile, probs = 0.025),
      CI_high = apply(pred_orig, 2, quantile, probs = 0.975)
    )
  
  p <- ggplot(pred_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
    facet_grid(Block ~ Genotype, scales = "fixed") +
    
    geom_point(
      data = raw_summary,
      aes(x = stimulus, y = y_raw, color = Genotype),
      alpha = 0.25,
      size = 1,
      inherit.aes = FALSE
    ) +
    
    geom_ribbon(
      aes(ymin = CI_low, ymax = CI_high),
      alpha = 0.15,
      color = NA
    ) +
    
    geom_line(
      aes(y = fit),
      linewidth = 1.3
    ) +
    
    theme_pubr(base_size = 14) +
    labs(
      x = "Stimulus number within block",
      y = label,
      color = "Genotype",
      fill = "Genotype"
    ) +
    theme(
      panel.spacing = unit(1.2, "lines"),
      legend.position = "top"
    )
  
  print(p)
  
  if (!is.null(save_path)) {
    ggsave(
      filename = save_path,
      plot = p,
      width = 14,
      height = 7,
      dpi = 300
    )
  }
  
  return(p)
}

p_brm <- plot_habituation_brm(
  df_brm = df_brm,
  fit = fit_brm,
  label = "Peak distance moved",
  raw_var = "max_peak",
  save_path = file.path(
    base_dir,
    "BRMS_habituation_curves.png"
  )
)

# ==============================================================================
# 8. Observation-level fitted diagnostic plot
# ==============================================================================

epred_obs <- posterior_epred(
  fit_brm,
  newdata = df_brm,
  re_formula = NULL,
  summary = TRUE
)

df_obs_fit <- bind_cols(df_brm, as.data.frame(epred_obs)) %>%
  mutate(
    fit_obs = exp(Estimate + sigma(fit_brm)^2 / 2) - 1
  )

obs_fit_summary <- df_obs_fit %>%
  group_by(Genotype, Block, stimulus) %>%
  summarise(
    raw = mean(max_peak, na.rm = TRUE),
    fitted = mean(fit_obs, na.rm = TRUE),
    .groups = "drop"
  )

p_obs_fit <- ggplot(obs_fit_summary, aes(x = stimulus)) +
  facet_grid(Block ~ Genotype, scales = "fixed") +
  geom_point(aes(y = raw), alpha = 0.4, size = 1) +
  geom_line(aes(y = fitted), linewidth = 1.2) +
  theme_pubr(base_size = 14) +
  labs(
    x = "Stimulus number within block",
    y = "Peak distance moved",
    title = "Observation-level fitted values"
  )

print(p_obs_fit)

ggsave(
  filename = file.path(base_dir, "BRMS_observation_level_fit_check.png"),
  plot = p_obs_fit,
  width = 14,
  height = 7,
  dpi = 300
)

# ==============================================================================
# 9. Extract biological curve parameters
# ==============================================================================

curve_grid <- expand.grid(
  Genotype = levels(df_brm$Genotype),
  Block = levels(df_brm$Block),
  stimulus = 1
) %>%
  mutate(
    stimulus0 = stimulus - 1,
    Genotype = factor(Genotype, levels = levels(df_brm$Genotype)),
    Block = factor(Block, levels = levels(df_brm$Block))
  )

get_nlpar_draws <- function(fit, newdata, nlpar) {
  posterior_linpred(
    fit,
    newdata = newdata,
    nlpar = nlpar,
    re_formula = NA,
    transform = FALSE
  )
}

Asym_draws <- get_nlpar_draws(fit_brm, curve_grid, "Asym")
R0_draws   <- get_nlpar_draws(fit_brm, curve_grid, "R0")
lrc_draws  <- get_nlpar_draws(fit_brm, curve_grid, "lrc")

param_draws <- bind_rows(lapply(seq_len(nrow(curve_grid)), function(i) {
  tibble(
    Genotype = curve_grid$Genotype[i],
    Block = curve_grid$Block[i],
    Asym_log = Asym_draws[, i],
    R0_log = R0_draws[, i],
    lrc = lrc_draws[, i]
  )
})) %>%
  mutate(
    k = exp(lrc),
    half_life_stimuli = log(2) / k,
    Asym_max_peak = exp(Asym_log + sigma(fit_brm)^2 / 2) - 1,
    R0_max_peak = exp(R0_log + sigma(fit_brm)^2 / 2) - 1
  )

curve_params_bayes <- param_draws %>%
  group_by(Genotype, Block) %>%
  summarise(
    Asym_median = median(Asym_max_peak),
    Asym_low = quantile(Asym_max_peak, 0.025),
    Asym_high = quantile(Asym_max_peak, 0.975),
    
    R0_median = median(R0_max_peak),
    R0_low = quantile(R0_max_peak, 0.025),
    R0_high = quantile(R0_max_peak, 0.975),
    
    k_median = median(k),
    k_low = quantile(k, 0.025),
    k_high = quantile(k, 0.975),
    
    half_life_median = median(half_life_stimuli),
    half_life_low = quantile(half_life_stimuli, 0.025),
    half_life_high = quantile(half_life_stimuli, 0.975),
    
    .groups = "drop"
  )

print(curve_params_bayes)

write.csv(
  curve_params_bayes,
  file = file.path(base_dir, "BRMS_curve_parameters.csv"),
  row.names = FALSE
)