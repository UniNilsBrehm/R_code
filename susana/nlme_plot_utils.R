# ------------------------------------------------------------------------------
# Helper
# ------------------------------------------------------------------------------

.safe_logit <- function(p, eps = 1e-6) {
  qlogis(pmin(pmax(p, eps), 1 - eps))
}

# ==============================================================================
# RESPONSE PROBABILITY
# ==============================================================================
# ------------------------------------------------------------------------------
# Probability-scale plot using aggregated raw probabilities
# ------------------------------------------------------------------------------

plot_habituation_probability <- function(
    df_resp,
    fit_model,
    save_dir,
    var_name
) {

  new_resp <- df_resp %>%
    group_by(Genotype, Block) %>%
    summarise(
      stim_min = min(stimulus, na.rm = TRUE),
      stim_max = max(stimulus, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    reframe(
      stimulus = seq(stim_min, stim_max, length.out = 200),
      .by = c(Genotype, Block)
    ) %>%
    mutate(
      stimulus0 = stimulus - 1,
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block = factor(Block, levels = levels(df_resp$Block))
    )
  
  pred_resp <- fitted(
    fit_model,
    newdata = new_resp,
    re_formula = NA,
    summary = TRUE
  )
  
  pred_resp_data <- bind_cols(new_resp, as.data.frame(pred_resp)) %>%
    rename(
      fit = Estimate,
      CI_low = Q2.5,
      CI_high = Q97.5
    )
  
  raw_prob <- df_resp %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      response_prob = mean(move, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  p <- ggplot(pred_resp_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
    facet_grid(Block ~ Genotype, scales = "fixed") +
    geom_point(
      data = raw_prob,
      aes(x = stimulus, y = response_prob, size = n, color = Genotype),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_pubr(base_size = 14) +
    labs(
      x = "Stimulus number within block",
      y = "Response probability",
      title = "Bayesian nonlinear habituation curves on probability scale",
      color = "Genotype",
      fill = "Genotype",
      size = "N"
    ) +
    theme(
      legend.position = "top",
      panel.spacing = unit(1.2, "lines")
    )
  if (!is.null(save_dir)) {
    ggsave(
      file.path(save_dir, paste0("nlme_", var_name, "_habituation_curves_probability_scale.png")),
      p,
      width = 14,
      height = 7,
      dpi = 300,
      bg = "white"
    )
  }
  
  p
}

# ------------------------------------------------------------------------------
# Probability-scale plot with TRUE raw binary data
# ------------------------------------------------------------------------------

plot_habituation_probability_raw <- function(
    df_resp,
    fit_model,
    save_dir,
    var_name
) {

  new_resp <- df_resp %>%
    group_by(Genotype, Block) %>%
    summarise(
      stim_min = min(stimulus, na.rm = TRUE),
      stim_max = max(stimulus, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    reframe(
      stimulus = seq(stim_min, stim_max, length.out = 200),
      .by = c(Genotype, Block)
    ) %>%
    mutate(
      stimulus0 = stimulus - 1,
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block = factor(Block, levels = levels(df_resp$Block))
    )
  
  pred_resp <- fitted(
    fit_model,
    newdata = new_resp,
    re_formula = NA,
    summary = TRUE
  )
  
  pred_resp_data <- bind_cols(new_resp, as.data.frame(pred_resp)) %>%
    rename(
      fit = Estimate,
      CI_low = Q2.5,
      CI_high = Q97.5
    )
  
  p <- ggplot(pred_resp_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
    facet_grid(Block ~ Genotype, scales = "fixed") +
    geom_jitter(
      data = df_resp,
      aes(x = stimulus, y = move, color = Genotype),
      inherit.aes = FALSE,
      width = 0.15,
      height = 0.03,
      alpha = 0.12,
      size = 0.7
    ) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_pubr(base_size = 14) +
    labs(
      x = "Stimulus number within block",
      y = "Response probability",
      title = "Bayesian nonlinear habituation curves with raw binary data",
      color = "Genotype",
      fill = "Genotype"
    ) +
    theme(
      legend.position = "top",
      panel.spacing = unit(1.2, "lines")
    )
  if (!is.null(save_dir)) {
  ggsave(
    file.path(save_dir, paste0("nlme_", var_name, "_habituation_curves_probability_scale_true_raw_binary.png")),
    p,
    width = 14,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  }
  p
}

# ------------------------------------------------------------------------------
# Logit-scale plot with aggregated raw data
# ------------------------------------------------------------------------------

plot_habituation_logit_raw <- function(
    df_resp,
    fit_model,
    save_dir,
    var_name
) {
  
  new_resp <- df_resp %>%
    group_by(Genotype, Block) %>%
    summarise(
      stim_min = min(stimulus, na.rm = TRUE),
      stim_max = max(stimulus, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    reframe(
      stimulus = seq(stim_min, stim_max, length.out = 200),
      .by = c(Genotype, Block)
    ) %>%
    mutate(
      stimulus0 = stimulus - 1,
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block = factor(Block, levels = levels(df_resp$Block))
    )
  
  pred_prob <- fitted(
    fit_model,
    newdata = new_resp,
    re_formula = NA,
    summary = TRUE
  )
  
  pred_logit_data <- bind_cols(new_resp, as.data.frame(pred_prob)) %>%
    transmute(
      Genotype,
      Block,
      stimulus,
      stimulus0,
      fit = .safe_logit(Estimate),
      CI_low = .safe_logit(Q2.5),
      CI_high = .safe_logit(Q97.5)
    )
  
  raw_logit <- df_resp %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      n = n(),
      successes = sum(move, na.rm = TRUE),
      response_prob = successes / n,
      empirical_prob = (successes + 0.5) / (n + 1),
      response_logit = qlogis(empirical_prob),
      .groups = "drop"
    )
  
  p <- ggplot(pred_logit_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
    facet_grid(Block ~ Genotype, scales = "fixed") +
    geom_point(
      data = raw_logit,
      aes(x = stimulus, y = response_logit, size = n, color = Genotype),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    theme_pubr(base_size = 14) +
    labs(
      x = "Stimulus number within block",
      y = "logit(response probability)",
      title = "Bayesian nonlinear habituation curves on logit scale",
      subtitle = "Raw points use empirical logit: qlogis((successes + 0.5) / (n + 1))",
      color = "Genotype",
      fill = "Genotype",
      size = "N"
    ) +
    theme(
      legend.position = "top",
      panel.spacing = unit(1.2, "lines")
    )
  if (!is.null(save_dir)) {
  ggsave(
    file.path(save_dir, paste0("nlme_", var_name, "_habituation_curves_logit_scale_raw_aggregated.png")),
    p,
    width = 14,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  }
  p
}