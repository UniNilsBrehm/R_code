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
library(DHARMa)
library(bayesplot)

source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/nlme_utils.R")
source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/plot_utils.R")

# source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/utils.R")
# source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_ISI90s_2Blocks/plot_utils.R")

# source("D:/Behavior_Data/R_code/susana/nlme_utils.R")
# source("D:/Behavior_Data/R_code/susana/plot_utils.R")

# base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/DarkFlash_ISI90s_2Blocks"
# base_dir <- "D:/WorkingData/Susana/NLME/DarkFlash_ISI90s_2Blocks"
base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/NLME//DarkFlash_ISI90s_2Blocks"

file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

var_name = 'peak_distance'
col_name = 'max_peak'

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
df_final_sub <- res$df_final_sub

df_resp <- df_final_sub %>%
  mutate(
    stimulus = as.numeric(stimulus),
    stimulus0 = stimulus - 1,
    Block = factor(Block),
    Genotype = factor(Genotype),
    Video = factor(Video),
    Well = factor(Well),
    animal = interaction(Video, Well, drop = TRUE)
    
  )

# Check if all values are larger than zero (for Gamma)
summary(df_resp[[col_name]])
any(df_resp[[col_name]] <= 0, na.rm = TRUE)

# Histogram
ggplot(df_resp, aes(x = .data[[col_name]])) +
  geom_histogram(
    bins = 100,
    fill = "lightblue",
    color = "white"
  ) +
  labs(
    title = paste("Histogram of", col_name),
    x = col_name,
    y = "Count"
  ) +
  theme_minimal()

# ==============================================================================
# The Model
# ==============================================================================
model <- bf(
  max_peak ~ exp(A) + (exp(R0) - exp(A)) * exp(-exp(logk) * stimulus0),
  A    ~ 1 + Genotype * Block + (1 | animal),
  R0   ~ 1 + Genotype * Block + (1 | animal),
  logk ~ 1 + Genotype * Block + (1 | animal),
  nl = TRUE
)

# ==============================================================================
# The Priors
# ==============================================================================
priors <- c(
  # Reference cell baseline
  prior(normal(2.0, 0.3),  nlpar = "R0",   class = "b", coef = "Intercept"),
  prior(normal(0.5, 0.3),  nlpar = "A",    class = "b", coef = "Intercept"),
  prior(normal(-2.3, 0.4), nlpar = "logk", class = "b", coef = "Intercept"),
  # Deviations from reference (all genotype contrasts, block, interactions)
  prior(normal(0, 0.2), nlpar = "R0",   class = "b"),
  prior(normal(0, 0.2), nlpar = "A",    class = "b"),
  prior(normal(0, 0.2), nlpar = "logk", class = "b"),
  # Animal REs
  prior(student_t(3, 0, 0.2), nlpar = "R0",   class = "sd"),
  prior(student_t(3, 0, 0.2), nlpar = "A",    class = "sd"),
  prior(student_t(3, 0, 0.2), nlpar = "logk", class = "sd"),
  prior(gamma(2, 0.1), class = "shape")
)

# ==============================================================================
# Validate Priors
# ==============================================================================
fit_prior_only <- brm(
  formula = model,
  data = df_resp,
  family = Gamma(link = "identity"),
  prior = priors,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  seed = 42,
  sample_prior = "only"
)

ep <- posterior_epred(fit_prior_only, ndraws = 500)

quantile(as.vector(ep), c(.001,.01,.05,.5,.95,.99,.999), na.rm = TRUE)
quantile(df_resp[[col_name]], c(.001,.01,.05,.5,.95,.99,.999), na.rm = TRUE)

prior_draws <- as_draws_df(fit_prior_only)
summary(prior_draws)


# ==============================================================================
# Fit Fast Test Model
# ==============================================================================
# Fit using Meanfield Variational Inference (VI)
fit_vi_test <- brm(
  formula = model,
  data = df_resp,               
  family  = Gamma(link = "identity"),
  prior = priors,
  backend = "cmdstanr",
  algorithm = "meanfield",          # <--- Reliable, fast VI method
  # algorithm = "fullrank",
  iter = 10000                      # VI likes higher iterations (it's still lightning fast)
)

fit_nuts_test <- brm(
  formula = model,
  data = df_resp,
  family  = Gamma(link = "identity"),
  prior = priors,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr",
  chains = 2,
  cores = 2,
  threads = threading(6),
  iter = 1000,
  warmup = 500,
  seed = 42,
  control = list(adapt_delta = 0.90, max_treedepth = 10),
  init    = 0,                        # important for nonlinear models
)

fit_model <- fit_vi_test
fit_model <- fit_nuts_test

# ==============================================================================
# Fit the model
# ==============================================================================
fit_model <- brm(
  model,
  data    = df_resp,
  family  = Gamma(link = "identity"),
  prior   = priors,
  save_pars = save_pars(all = TRUE),
  chains  = 4, 
  cores = 4,
  threads = threading(6),
  iter    = 4000, 
  warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  init    = 0,
  backend = "cmdstanr"
)

# ==============================================================================
# Save fitted model to HDD
# ==============================================================================
saveRDS(
  fit_model,
  file = file.path(base_dir, "models", paste0("bayesian_nlme_", var_name,"_results.rds"))
)

# ==============================================================================
# Load fitted model if available
# ==============================================================================
fit_model <- readRDS(
  file.path(base_dir, "models", paste0("bayesian_nlme_", var_name,"_results.rds"))
)


# ==============================================================================
# Get summary and diagnostics
# ==============================================================================

diag_dir <- file.path(base_dir, "models", "diagnostics", var_name)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Get summary
# ==============================================================================
# --- Summary table: R-hat, ESS ---
sink(file.path(save_results_dir, paste0(var_name, '_model_summary.txt')))
fit_summary <- summary(fit_model)
print(fit_summary)
sink()

# Divergences
sum(subset(nuts_params(fit_model), Parameter == "divergent__")$Value)
# 0 ideal, single digits tolerable, more is concerning

# Energy diagnostics
energy_plot <- bayesplot::mcmc_nuts_energy(nuts_params(fit_model))
ggsave(
  filename = file.path(diag_dir, "nuts_energy.png"),
  plot     = energy_plot,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# Posterior predictive check
ppc_marginal <- pp_check(fit_model, type = "dens_overlay", ndraws = 200)
ppc_decay <- pp_check(fit_model, type = "stat_grouped", group = "stimulus0",
                      stat = "median", ndraws = 200)
ppc_genotype <- pp_check(fit_model, type = "stat_grouped", group = "Genotype",
                         stat = "mean", ndraws = 200)
ppc_list <- list(
  ppc_marginal  = ppc_marginal,
  ppc_decay     = ppc_decay,
  ppc_genotype  = ppc_genotype
)
purrr::iwalk(ppc_list, ~ ggsave(
  filename = file.path(diag_dir, paste0(.y, ".png")),
  plot     = .x,
  width    = 8,
  height   = 6,
  dpi      = 300
))


mcmc_pairs_plot <- bayesplot::mcmc_pairs(
  fit_model, np = nuts_params(fit_model),
  pars = c("b_logk_Intercept", "sd_animal__logk_Intercept"),
  off_diag_args = list(size = 0.5)
)
ggsave(
  filename = file.path(diag_dir, "mcmc_pairs.png"),
  plot     = mcmc_pairs_plot,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# Residuals
yrep <- posterior_predict(fit_model, ndraws = 1000)

sim_res <- createDHARMa(
  simulatedResponse = t(yrep),
  observedResponse = df_resp$max_peak,
  fittedPredictedResponse = fitted(fit_model)[, "Estimate"]
)

print(plot(sim_res))
plotResiduals(sim_res, df_resp$stimulus)
plotResiduals(sim_res, df_resp$stimulus0)
plotResiduals(sim_res, df_resp$Genotype)
plotResiduals(sim_res, df_resp$Block)
plotResiduals(sim_res, df_resp$animal)

save_plot_as_png(
  paste0("nlme_", var_name, "_DHARMa_residuals.png"),
  quote(plot(sim_res))
)


# Compute leave-one-out cross-validation:
loo_var <- loo(fit_model)

print(loo_var)
plot(loo_var)
save_plot_as_png(
  paste0("nlme_", var_name, "_loo_plot.png"),
  quote(plot(loo_var))
)
# Identify which observations they are
plot(loo_var, diagnostic = "k")

# Or extract directly
problematic <- which(loo_var$diagnostics$pareto_k > 0.7)
df_resp[problematic, ]

# loo_compare(loo1, loo2)  # compare models


# ==============================================================================
# Plot Habituation curves
# ==============================================================================
p_hab_ind <- plot_habituation_response(
  df_resp, fit_model,
  response_var = "max_peak",
  y_label = "Peak distance moved",
  y_limits = c(0, 10),
  raw_data = "individual",
  ndraws        = 200, 
  save_fig_dir = save_fig_dir
)

p_hab_ind

p_hab <- plot_habituation_response(
  df_resp, fit_model,
  response_var = "max_peak",
  y_label = "Peak distance moved",
  y_limits = c(0, 10),
  raw_data = "trials",
  save_fig_dir = save_fig_dir
)

p_hab

p_hab_agg <- plot_habituation_response(
  df_resp, fit_model,
  response_var = "max_peak",
  y_label = "Peak distance moved",
  y_limits = c(0, 10),
  raw_data = "aggregate",
  save_fig_dir = save_fig_dir
)

p_hab_agg


p_hab_averaged <- plot_habituation_response_animal_averaged(
  df_resp, fit_model,
  response_var = "max_peak",
  ndraws = 500,
  y_label = "Peak distance moved",
  y_limits = c(0, 10),
  save_fig_dir = save_fig_dir
)

p_hab_averaged


# ==============================================================================
# Compare nonlinear model parameters - peak distance
# ==============================================================================

# ==============================================================================
# 1. Habituation rate k = exp(logk)
# ==============================================================================
k_draws_df <- make_nlpar_draws(
  fit_model     = fit_model,
  df_resp       = df_resp,
  nlpar         = "logk",
  transform_fun = exp
) %>%
  rename(k = value, logk = value_raw)

k_summary <- k_draws_df %>%
  group_by(Genotype, Block) %>%
  summarise(
    k_median         = median(k),
    k_low            = quantile(k, 0.025),
    k_high           = quantile(k, 0.975),
    half_life_median = median(log(2) / k),
    half_life_low    = quantile(log(2) / k, 0.025),
    half_life_high   = quantile(log(2) / k, 0.975),
    .groups = "drop"
  )

print(k_summary)

write.csv(
  k_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_k_summary.csv")),
  row.names = FALSE
)

k_comparison <- compare_nlpar(
  k_draws_df %>% rename(value = k),
  value_name = "k",
  ratio      = TRUE
)

k_comparison_summary <- k_comparison$summary
print(k_comparison_summary)

write.csv(
  k_comparison_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_k_all_pairwise_comparisons.csv")),
  row.names = FALSE
)

# ==============================================================================
# 2. Asymptote A_inf = exp(A)  (distance moved at infinite stimuli)
# ==============================================================================
Ainf_draws_df <- make_nlpar_draws(
  fit_model     = fit_model,
  df_resp       = df_resp,
  nlpar         = "A",
  transform_fun = exp
) %>%
  rename(Ainf = value, A_log = value_raw)

Ainf_summary <- summarise_nlpar(
  Ainf_draws_df %>% rename(value = Ainf),
  "Ainf"
)

print(Ainf_summary)

write.csv(
  Ainf_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_Ainf_summary.csv")),
  row.names = FALSE
)

Ainf_comparison <- compare_nlpar(
  Ainf_draws_df %>% rename(value = Ainf),
  value_name = "Ainf",
  ratio      = TRUE   # multiplicative scale makes ratios more natural here
)

Ainf_comparison_summary <- Ainf_comparison$summary
print(Ainf_comparison_summary)

write.csv(
  Ainf_comparison_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_Ainf_all_pairwise_comparisons.csv")),
  row.names = FALSE
)

# ==============================================================================
# 3. Initial response R0 = exp(R0)  (distance moved at stimulus 1)
# ==============================================================================
R0_draws_df <- make_nlpar_draws(
  fit_model     = fit_model,
  df_resp       = df_resp,
  nlpar         = "R0",
  transform_fun = exp
) %>%
  rename(R0 = value, R0_log = value_raw)

R0_summary <- summarise_nlpar(
  R0_draws_df %>% rename(value = R0),
  "R0"
)

print(R0_summary)

write.csv(
  R0_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_R0_summary.csv")),
  row.names = FALSE
)

R0_comparison <- compare_nlpar(
  R0_draws_df %>% rename(value = R0),
  value_name = "R0",
  ratio      = TRUE
)

R0_comparison_summary <- R0_comparison$summary
print(R0_comparison_summary)

write.csv(
  R0_comparison_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_R0_all_pairwise_comparisons.csv")),
  row.names = FALSE
)

# ==============================================================================
# 4. Endpoint contrast: last 5 stimuli
# ==============================================================================
last_n     <- 5
last_stims <- (max(df_resp$stimulus) - last_n + 1):max(df_resp$stimulus)

grid_end <- expand.grid(
  Genotype = levels(df_resp$Genotype),
  Block    = levels(df_resp$Block),
  stimulus = last_stims
) %>%
  mutate(
    stimulus0 = stimulus - 1,
    Genotype  = factor(Genotype, levels = levels(df_resp$Genotype)),
    Block     = factor(Block,    levels = levels(df_resp$Block))
  )

pred_end_draws <- fitted(
  fit_model,
  newdata    = grid_end,
  re_formula = NA,
  summary    = FALSE
)

n_draws <- nrow(pred_end_draws)

draws_end <- bind_rows(
  lapply(levels(df_resp$Genotype), function(g) {
    bind_rows(
      lapply(levels(df_resp$Block), function(b) {
        row_idx <- which(grid_end$Genotype == g & grid_end$Block == b)
        tibble(
          draw     = seq_len(n_draws),
          Genotype = factor(g, levels = levels(df_resp$Genotype)),
          Block    = factor(b, levels = levels(df_resp$Block)),
          value    = rowMeans(pred_end_draws[, row_idx])
        )
      })
    )
  })
)

end_summary <- draws_end %>%
  group_by(Genotype, Block) %>%
  summarise(
    median = median(value),
    low    = quantile(value, 0.025),
    high   = quantile(value, 0.975),
    .groups = "drop"
  )

print(end_summary)

write.csv(
  end_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_endpoint_summary.csv")),
  row.names = FALSE
)

end_comparison <- compare_nlpar(
  draws_end,
  value_name = "endpoint",
  rope       = 0.5,    # adjust to a meaningful distance unit for your data
  ratio      = TRUE    # ratios natural for a distance measure
)

end_comparison_summary <- end_comparison$summary
print(end_comparison_summary)

write.csv(
  end_comparison_summary,
  file.path(save_results_dir, paste0("nlme_", var_name, "_endpoint_all_pairwise_comparisons.csv")),
  row.names = FALSE
)

# ==============================================================================
# 5. Plots
# ==============================================================================

# k posterior density
p_k <- plot_posterior_density(
  k_draws_df,
  k,
  "Posterior distributions of habituation rate k",
  "Habituation rate k",
  paste0("nlme_", var_name, "_k_posteriors.png"),
  save_fig_dir = save_fig_dir
)

# k pairwise ratios
p_k_compare <- k_comparison_summary %>%
  ggplot(aes(
    x     = comparison,
    y     = median_ratio,
    ymin  = ratio_low,
    ymax  = ratio_high,
    color = Genotype_1
  )) +
  facet_wrap(~Block, scales = "free_x") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(linewidth = 0.8) +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x     = "Comparison",
    y     = "Habituation rate ratio",
    title = "All pairwise habituation rate comparisons",
    color = "Numerator genotype"
  )

print(p_k_compare)

ggsave(
  file.path(save_fig_dir, paste0("nlme_", var_name, "_k_all_pairwise_ratios.png")),
  p_k_compare,
  width = 10, height = 8, dpi = 300, bg = "white"
)

# Ainf posterior density
p_Ainf <- plot_posterior_density(
  Ainf_draws_df,
  Ainf,
  "Posterior distributions of asymptote Ainf (distance moved)",
  "Asymptote distance moved Ainf",
  paste0("nlme_", var_name, "_Ainf_posteriors.png"),
  save_fig_dir = save_fig_dir
)

# Ainf pairwise ratios
p_Ainf_compare <- Ainf_comparison_summary %>%
  ggplot(aes(
    x     = comparison,
    y     = median_ratio,
    ymin  = ratio_low,
    ymax  = ratio_high,
    color = Genotype_1
  )) +
  facet_wrap(~Block, scales = "free_x") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(linewidth = 0.8) +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x     = "Comparison",
    y     = "Ainf ratio",
    title = "All pairwise asymptote comparisons",
    color = "Numerator genotype"
  )

print(p_Ainf_compare)

ggsave(
  file.path(save_fig_dir, paste0("nlme_", var_name, "_Ainf_all_pairwise_ratios.png")),
  p_Ainf_compare,
  width = 10, height = 8, dpi = 300, bg = "white"
)

# R0 posterior density
p_R0 <- plot_posterior_density(
  R0_draws_df,
  R0,
  "Posterior distributions of initial response R0 (distance moved)",
  "Initial distance moved R0",
  paste0("nlme_", var_name, "_R0_posteriors.png"),
  save_fig_dir = save_fig_dir
)

# R0 pairwise ratios
p_R0_compare <- R0_comparison_summary %>%
  ggplot(aes(
    x     = comparison,
    y     = median_ratio,
    ymin  = ratio_low,
    ymax  = ratio_high,
    color = Genotype_1
  )) +
  facet_wrap(~Block, scales = "free_x") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(linewidth = 0.8) +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x     = "Comparison",
    y     = "R0 ratio",
    title = "All pairwise initial response comparisons",
    color = "Numerator genotype"
  )

print(p_R0_compare)

ggsave(
  file.path(save_fig_dir, paste0("nlme_", var_name, "_R0_all_pairwise_ratios.png")),
  p_R0_compare,
  width = 10, height = 8, dpi = 300, bg = "white"
)

# Endpoint posterior density
p_end_density <- ggplot(
  draws_end,
  aes(x = value, fill = Genotype)
) +
  facet_wrap(~Block, scales = "free_y") +
  geom_density(alpha = 0.35) +
  geom_vline(
    data = end_summary,
    aes(xintercept = median, color = Genotype),
    linewidth = 0.8, linetype = "dashed"
  ) +
  theme_pubr(base_size = 14) +
  labs(
    x     = paste0("Distance moved (mean of last ", last_n, " stimuli)"),
    y     = "Posterior density",
    title = paste0("Endpoint distance moved (last ", last_n, " stimuli)"),
    fill  = "Genotype",
    color = "Genotype"
  )

print(p_end_density)

ggsave(
  file.path(save_fig_dir, paste0("nlme_", var_name, "_endpoint_posteriors.png")),
  p_end_density,
  width = 10, height = 5, dpi = 300, bg = "white"
)

# Endpoint pointrange
p_end_summary <- end_summary %>%
  ggplot(aes(
    x     = Genotype,
    y     = median,
    ymin  = low,
    ymax  = high,
    color = Genotype
  )) +
  facet_wrap(~Block) +
  geom_pointrange(linewidth = 0.8, size = 0.8) +
  theme_pubr(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x     = "Genotype",
    y     = paste0("Distance moved (last ", last_n, " stimuli)"),
    title = paste0("Endpoint distance moved (last ", last_n, " stimuli)")
  )

print(p_end_summary)

ggsave(
  file.path(save_fig_dir, paste0("nlme_", var_name, "_endpoint_summary_pointrange.png")),
  p_end_summary,
  width = 8, height = 5, dpi = 300, bg = "white"
)

# Endpoint pairwise ratios
p_end_compare <- end_comparison_summary %>%
  ggplot(aes(
    x     = comparison,
    y     = median_ratio,
    ymin  = ratio_low,
    ymax  = ratio_high,
    color = Genotype_1
  )) +
  facet_wrap(~Block, scales = "free_x") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(linewidth = 0.8) +
  coord_flip() +
  theme_pubr(base_size = 14) +
  labs(
    x     = "Comparison",
    y     = paste0("Endpoint ratio (last ", last_n, " stimuli)"),
    title = paste0("All pairwise endpoint comparisons (last ", last_n, " stimuli)"),
    color = "Numerator genotype"
  )

print(p_end_compare)

ggsave(
  file.path(save_fig_dir, paste0("nlme_", var_name, "_endpoint_all_pairwise_ratios.png")),
  p_end_compare,
  width = 10, height = 8, dpi = 300, bg = "white"
)