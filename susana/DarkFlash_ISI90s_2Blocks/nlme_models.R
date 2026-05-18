# ==============================================================================
# 0. Packages
# ==============================================================================
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(nlme)
library(MASS)
library(emmeans)

source("C:/Users/NilsPC/Desktop/Susana/NLME/utils.R")

base_dir <- "C:/Users/NilsPC/Desktop/Susana/NLME/"
file_dir <- "C:/Users/NilsPC/Desktop/Susana/NLME/SPZ_ISI60_removed_non_responders_2stimuli.csv"


# ==============================================================================
# 1. Load data
# ==============================================================================
res <- load_data_darkflash_60s(file_dir, move_th = 0, take_peak = 0)

df_final     <- res$df_final
df_final_sub <- res$df_final_sub


# ==============================================================================
# 2. Prepare data for NLME max_peak model
# ==============================================================================
df_nlme <- df_final_sub %>%
  mutate(
    animal = interaction(Video, Well, drop = TRUE),
    stimulus = as.numeric(stimulus),
    Block = factor(Block),
    Genotype = factor(Genotype),
    y = log(max_peak + 1)
  ) %>%
  filter(
    is.finite(y),
    !is.na(stimulus),
    !is.na(Block),
    !is.na(Genotype),
    !is.na(animal)
  )

df_nlme_g <- groupedData(
  y ~ stimulus | animal,
  data = df_nlme
)

n_per_genotype <- df_nlme %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)

print(n_per_genotype)


# ==============================================================================
# 3. Plot raw habituation data
# ==============================================================================
p_raw <- ggplot(
  df_nlme,
  aes(
    x = stimulus,
    y = max_peak,
    color = Genotype,
    group = interaction(Genotype, animal)
  )
) +
  facet_wrap(~Block, scales = "free_x") +
  geom_line(alpha = 0.12) +
  stat_summary(
    aes(group = Genotype),
    fun = mean,
    geom = "line",
    linewidth = 1.4
  ) +
  theme_pubr(base_size = 14) +
  labs(
    x = "Stimulus number within block",
    y = "max_peak",
    color = "Genotype"
  )

print(p_raw)


# ==============================================================================
# 4. Get starting values from simple nonlinear model
# ==============================================================================
m_start <- nls(
  y ~ SSasymp(stimulus, Asym, R0, lrc),
  data = df_nlme
)

start_vals <- coef(m_start)
print(start_vals)


# ==============================================================================
# 5. Helper for fixed-effect starting values
# ==============================================================================
make_start <- function(data, formula_rhs = ~ Genotype * Block, start_vals) {
  X <- model.matrix(formula_rhs, data = data)
  n_coef <- ncol(X)
  
  c(
    Asym = c(start_vals["Asym"], rep(0, n_coef - 1)),
    R0   = c(start_vals["R0"],   rep(0, n_coef - 1)),
    lrc  = c(start_vals["lrc"],  rep(0, n_coef - 1))
  )
}

start_geno_block <- make_start(
  data = df_nlme,
  formula_rhs = ~ Genotype * Block,
  start_vals = start_vals
)


# ==============================================================================
# 6. Fit models
# ==============================================================================

# Null model: one global habituation curve
m_null <- nlme(
  y ~ SSasymp(stimulus, Asym, R0, lrc),
  
  data = df_nlme_g,
  
  fixed = Asym + R0 + lrc ~ 1,
  
  random = lrc ~ 1 | animal,
  
  start = start_vals,
  
  method = "ML",
  
  control = nlmeControl(
    maxIter = 500,
    pnlsMaxIter = 100,
    msMaxIter = 500,
    returnObject = TRUE
  )
)


# Block-only model
m_block <- nlme(
  y ~ SSasymp(stimulus, Asym, R0, lrc),
  
  data = df_nlme_g,
  
  fixed = list(
    Asym ~ Block,
    R0   ~ Block,
    lrc  ~ Block
  ),
  
  random = lrc ~ 1 | animal,
  
  start = c(
    Asym = c(start_vals["Asym"], 0),
    R0   = c(start_vals["R0"],   0),
    lrc  = c(start_vals["lrc"],  0)
  ),
  
  method = "ML",
  
  control = nlmeControl(
    maxIter = 500,
    pnlsMaxIter = 100,
    msMaxIter = 500,
    returnObject = TRUE
  )
)


# Genotype-only model
X_genotype <- model.matrix(~ Genotype, data = df_nlme)
start_genotype <- c(
  Asym = c(start_vals["Asym"], rep(0, ncol(X_genotype) - 1)),
  R0   = c(start_vals["R0"],   rep(0, ncol(X_genotype) - 1)),
  lrc  = c(start_vals["lrc"],  rep(0, ncol(X_genotype) - 1))
)

m_genotype <- nlme(
  y ~ SSasymp(stimulus, Asym, R0, lrc),
  
  data = df_nlme_g,
  
  fixed = list(
    Asym ~ Genotype,
    R0   ~ Genotype,
    lrc  ~ Genotype
  ),
  
  random = lrc ~ 1 | animal,
  
  start = start_genotype,
  
  method = "ML",
  
  control = nlmeControl(
    maxIter = 500,
    pnlsMaxIter = 100,
    msMaxIter = 500,
    returnObject = TRUE
  )
)


# Full experimental design model: Genotype * Block
m_genotype_block <- nlme(
  y ~ SSasymp(stimulus, Asym, R0, lrc),
  
  data = df_nlme_g,
  
  fixed = list(
    Asym ~ Genotype * Block,
    R0   ~ Genotype * Block,
    lrc  ~ Genotype * Block
  ),
  
  random = lrc ~ 1 | animal,
  
  start = start_geno_block,
  
  method = "ML",
  
  control = nlmeControl(
    maxIter = 500,
    pnlsMaxIter = 100,
    msMaxIter = 500,
    returnObject = TRUE
  )
)


# ==============================================================================
# 7. Model comparison
# ==============================================================================
anova(m_null, m_block, m_genotype, m_genotype_block)

summary(m_genotype_block)

warnings()


# ==============================================================================
# 8. Refit final model with REML if desired
# ==============================================================================
# Use ML for model comparison.
# Use REML for final estimates if model structure is fixed.

m_final <- update(m_genotype_block, method = "REML")

summary(m_final)


# ==============================================================================
# 9. Helper: fixed-effect prediction with approximate CI
# ==============================================================================
predict_nlme_fixed_ci <- function(model, new_data, formula_rhs = ~ Genotype * Block,
                                  n_sim = 1000, backtransform = TRUE) {
  
  beta <- fixed.effects(model)
  V <- vcov(model)
  
  X <- model.matrix(formula_rhs, data = new_data)
  
  n_coef <- ncol(X)
  
  beta_Asym <- beta[1:n_coef]
  beta_R0   <- beta[(n_coef + 1):(2 * n_coef)]
  beta_lrc  <- beta[(2 * n_coef + 1):(3 * n_coef)]
  
  Asym <- as.numeric(X %*% beta_Asym)
  R0   <- as.numeric(X %*% beta_R0)
  lrc  <- as.numeric(X %*% beta_lrc)
  
  fit_log <- Asym + (R0 - Asym) * exp(-exp(lrc) * new_data$stimulus)
  
  beta_sim <- MASS::mvrnorm(n_sim, mu = beta, Sigma = V)
  
  pred_sim <- apply(beta_sim, 1, function(b) {
    b_Asym <- b[1:n_coef]
    b_R0   <- b[(n_coef + 1):(2 * n_coef)]
    b_lrc  <- b[(2 * n_coef + 1):(3 * n_coef)]
    
    Asym_s <- as.numeric(X %*% b_Asym)
    R0_s   <- as.numeric(X %*% b_R0)
    lrc_s  <- as.numeric(X %*% b_lrc)
    
    Asym_s + (R0_s - Asym_s) * exp(-exp(lrc_s) * new_data$stimulus)
  })
  
  CI_low_log  <- apply(pred_sim, 1, quantile, probs = 0.025, na.rm = TRUE)
  CI_high_log <- apply(pred_sim, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  if (backtransform) {
    new_data$fit     <- exp(fit_log) - 1
    new_data$CI_low  <- exp(CI_low_log) - 1
    new_data$CI_high <- exp(CI_high_log) - 1
  } else {
    new_data$fit     <- fit_log
    new_data$CI_low  <- CI_low_log
    new_data$CI_high <- CI_high_log
  }
  
  new_data
}


# ==============================================================================
# 10. Plot function for NLME habituation model
# ==============================================================================
plot_habituation_nlme <- function(df_nlme, model,
                                  label = "max_peak",
                                  raw_var = "max_peak",
                                  formula_rhs = ~ Genotype * Block,
                                  backtransform = TRUE,
                                  n_sim = 1000) {
  
  raw_summary <- df_nlme %>%
    group_by(Block, Genotype, stimulus) %>%
    summarise(
      y_raw = mean(.data[[raw_var]], na.rm = TRUE),
      .groups = "drop"
    )
  
  blocks <- levels(df_nlme$Block)
  
  stim_ranges <- lapply(blocks, function(b) {
    range(df_nlme$stimulus[df_nlme$Block == b], na.rm = TRUE)
  })
  names(stim_ranges) <- blocks
  
  new_data <- bind_rows(lapply(blocks, function(b) {
    expand.grid(
      Block = b,
      stimulus = seq(
        stim_ranges[[b]][1],
        stim_ranges[[b]][2],
        length.out = 100
      )
    )
  })) %>%
    crossing(Genotype = levels(df_nlme$Genotype)) %>%
    mutate(
      Block = factor(Block, levels = levels(df_nlme$Block)),
      Genotype = factor(Genotype, levels = levels(df_nlme$Genotype))
    )
  
  pred_data <- predict_nlme_fixed_ci(
    model = model,
    new_data = new_data,
    formula_rhs = formula_rhs,
    n_sim = n_sim,
    backtransform = backtransform
  )
  
  labels <- setNames(
    paste0("Block ", blocks),
    blocks
  )
  
  ggplot(pred_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
    facet_wrap(
      ~Block,
      ncol = length(blocks),
      scales = "free_x",
      labeller = as_labeller(labels)
    ) +
    geom_point(
      data = raw_summary,
      aes(x = stimulus, y = y_raw, color = Genotype),
      size = 1,
      alpha = 0.26,
      inherit.aes = FALSE
    ) +
    geom_ribbon(
      aes(ymin = CI_low, ymax = CI_high),
      alpha = 0.10,
      color = NA
    ) +
    geom_line(aes(y = fit), linewidth = 1.5) +
    scale_x_continuous(n.breaks = 4) +
    labs(
      x = "Stimulus number within block",
      y = label,
      color = "Genotype",
      fill = "Genotype"
    ) +
    theme_pubr(base_size = 14) +
    theme(
      panel.spacing = unit(1.5, "lines")
    )
}


# ==============================================================================
# 11. Final plot
# ==============================================================================
p_hab <- plot_habituation_nlme(
  df_nlme = df_nlme,
  model = m_final,
  label = "max_peak",
  raw_var = "max_peak",
  formula_rhs = ~ Genotype * Block,
  backtransform = TRUE,
  n_sim = 1000
)

print(p_hab)

ggsave(
  filename = file.path(base_dir, "NLME_max_peak_habituation_by_Genotype_Block.png"),
  plot = p_hab,
  width = 12,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 12. Extract biologically useful parameters
# ==============================================================================

extract_curve_params <- function(model, df_nlme, formula_rhs = ~ Genotype * Block) {
  
  grid <- expand.grid(
    Genotype = levels(df_nlme$Genotype),
    Block = levels(df_nlme$Block)
  ) %>%
    mutate(
      Genotype = factor(Genotype, levels = levels(df_nlme$Genotype)),
      Block = factor(Block, levels = levels(df_nlme$Block))
    )
  
  beta <- fixed.effects(model)
  X <- model.matrix(formula_rhs, data = grid)
  n_coef <- ncol(X)
  
  beta_Asym <- beta[1:n_coef]
  beta_R0   <- beta[(n_coef + 1):(2 * n_coef)]
  beta_lrc  <- beta[(2 * n_coef + 1):(3 * n_coef)]
  
  grid %>%
    mutate(
      Asym_log = as.numeric(X %*% beta_Asym),
      R0_log   = as.numeric(X %*% beta_R0),
      lrc      = as.numeric(X %*% beta_lrc),
      k        = exp(lrc),
      half_life_stimuli = log(2) / k,
      Asym_max_peak = exp(Asym_log) - 1,
      R0_max_peak   = exp(R0_log) - 1
    )
}

curve_params <- extract_curve_params(
  model = m_final,
  df_nlme = df_nlme,
  formula_rhs = ~ Genotype * Block
)

print(curve_params)

write.csv(
  curve_params,
  file = file.path(base_dir, "NLME_max_peak_curve_parameters.csv"),
  row.names = FALSE
)


# ==============================================================================
# AR(1)
# ==============================================================================
df_nlme <- df_nlme %>%
  arrange(animal, Block, stimulus)

df_nlme_g <- groupedData(
  y ~ stimulus | animal,
  data = df_nlme
)

m_final_ar1 <- nlme(
  y ~ SSasymp(stimulus, Asym, R0, lrc),
  
  data = df_nlme_g,
  
  fixed = list(
    Asym ~ Genotype * Block,
    R0   ~ Genotype * Block,
    lrc  ~ Genotype * Block
  ),
  
  random = lrc ~ 1 | animal,
  
  correlation = corAR1(form = ~ stimulus | animal/Block),
  
  start = start_geno_block,
  
  method = "REML",
  
  control = nlmeControl(
    maxIter = 500,
    pnlsMaxIter = 100,
    msMaxIter = 500,
    returnObject = TRUE
  )
)

m_final_ml <- update(m_final, method = "ML")
m_final_ar1_ml <- update(m_final_ar1, method = "ML")

anova(m_final_ml, m_final_ar1_ml)

curve_params_ar1 <- extract_curve_params(
  model = m_final_ar1,
  df_nlme = df_nlme,
  formula_rhs = ~ Genotype * Block
)

print(curve_params_ar1)

p_hab_ar1 <- plot_habituation_nlme(
  df_nlme = df_nlme_g,
  model = m_final_ar1_ml,
  label = "max_peak",
  raw_var = "max_peak",
  formula_rhs = ~ Genotype * Block,
  backtransform = TRUE,
  n_sim = 1000
)

print(p_hab_ar1)
ggsave(
  filename = file.path(base_dir, "NLME_max_peak_habituation_by_Genotype_Block.png"),
  plot = p_hab_ar1,
  width = 12,
  height = 6,
  dpi = 300
)

########
# Habituation curves: separate panels for each Genotype × Block

plot_habituation_nlme_by_genotype_block <- function(df_nlme, model,
                                                    label = "max_peak",
                                                    raw_var = "max_peak",
                                                    formula_rhs = ~ Genotype * Block,
                                                    backtransform = TRUE,
                                                    n_sim = 1000) {
  
  raw_summary <- df_nlme %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      y_raw = mean(.data[[raw_var]], na.rm = TRUE),
      .groups = "drop"
    )
  
  new_data <- df_nlme %>%
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
      Genotype = factor(Genotype, levels = levels(df_nlme$Genotype)),
      Block = factor(Block, levels = levels(df_nlme$Block))
    )
  
  pred_data <- predict_nlme_fixed_ci(
    model = model,
    new_data = new_data,
    formula_rhs = formula_rhs,
    n_sim = n_sim,
    backtransform = backtransform
  )
  
  ggplot(pred_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
    facet_grid(Block ~ Genotype, scales = "free_y") +
    geom_point(
      data = raw_summary,
      aes(x = stimulus, y = y_raw, color = Genotype),
      alpha = 0.25,
      size = 1,
      inherit.aes = FALSE
    ) +
    geom_ribbon(
      aes(ymin = CI_low, ymax = CI_high),
      alpha = 0.12,
      color = NA
    ) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    theme_pubr(base_size = 14) +
    labs(
      x = "Stimulus number within block",
      y = label,
      color = "Genotype",
      fill = "Genotype"
    ) +
    theme(
      panel.spacing = unit(1.2, "lines")
    )
}

p_by_genotype_block <- plot_habituation_nlme_by_genotype_block(
  df_nlme = df_nlme,
  model = m_final_ar1,
  label = "Peak distance moved",
  raw_var = "max_peak",
  formula_rhs = ~ Genotype * Block,
  backtransform = TRUE,
  n_sim = 1000
)

print(p_by_genotype_block)

ggsave(
  filename = file.path(base_dir, "NLME_max_peak_habituation_each_genotype_each_block.png"),
  plot = p_by_genotype_block,
  width = 14,
  height = 7,
  dpi = 300
)