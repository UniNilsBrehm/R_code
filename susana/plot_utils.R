# ==============================================================================
# plot_utils.R
# Utility plotting functions for habituation models
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(brms)

# ------------------------------------------------------------------------------
# Probability-scale plot using aggregated raw probabilities
# ------------------------------------------------------------------------------

plot_habituation_probability <- function(
    df_resp,
    fit_model,
    save_fig_dir,
    var_name
) {
  
  new_resp <- df_resp %>%
    group_by(Genotype, Block) %>%
    summarise(
      stim_min = min(stimulus, na.rm = TRUE),
      stim_max = max(stimulus, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Genotype, Block) %>%
    reframe(
      stimulus = seq(stim_min, stim_max, length.out = 200)
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
  
  pred_resp_data <- bind_cols(
    new_resp,
    as.data.frame(pred_resp)
  ) %>%
    rename(
      fit = Estimate,
      CI_low = Q2.5,
      CI_high = Q97.5
    )
  
  raw_prob <- df_resp %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      response_prob = mean(move, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(
    pred_resp_data,
    aes(x = stimulus, color = Genotype, fill = Genotype)
  ) +
    
    facet_grid(Block ~ Genotype, scales = "fixed") +
    
    geom_point(
      data = raw_prob,
      aes(
        x = stimulus,
        y = response_prob,
        color = Genotype
      ),
      inherit.aes = FALSE,
      alpha = 0.5,
      size = 1
    ) +
    
    geom_ribbon(
      aes(
        ymin = CI_low,
        ymax = CI_high
      ),
      alpha = 0.15,
      color = NA
    ) +
    
    geom_line(
      aes(y = fit),
      linewidth = 1.3
    ) +
    
    coord_cartesian(ylim = c(0, 1)) +
    
    theme_pubr(base_size = 14) +
    
    labs(
      x = "Stimulus number within block",
      y = "Response probability",
      title = "Bayesian nonlinear habituation curves on probability scale",
      color = "Genotype",
      fill = "Genotype"
    ) +
    
    theme(
      legend.position = "top",
      panel.spacing = unit(1.2, "lines")
    )
  
  if (!is.null(save_fig_dir)) {
    ggsave(
      filename = file.path(
        save_fig_dir,
        paste0(
          "nlme_",
          var_name,
          "_habituation_curves_probability_scale.png"
        )
      ),
      plot = p,
      width = 14,
      height = 7,
      dpi = 300
    )
  }
  return(p)
}


# ------------------------------------------------------------------------------
# Probability-scale plot with TRUE raw binary data
# ------------------------------------------------------------------------------

plot_habituation_probability_raw <- function(
    df_resp,
    fit_model,
    save_fig_dir,
    var_name
) {
  
  new_resp <- df_resp %>%
    group_by(Genotype, Block) %>%
    summarise(
      stim_min = min(stimulus, na.rm = TRUE),
      stim_max = max(stimulus, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Genotype, Block) %>%
    reframe(
      stimulus = seq(stim_min, stim_max, length.out = 200)
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
  
  pred_resp_data <- bind_cols(
    new_resp,
    as.data.frame(pred_resp)
  ) %>%
    rename(
      fit = Estimate,
      CI_low = Q2.5,
      CI_high = Q97.5
    )
  
  p <- ggplot(
    pred_resp_data,
    aes(x = stimulus, color = Genotype, fill = Genotype)
  ) +
    
    facet_grid(Block ~ Genotype, scales = "fixed") +
    
    geom_jitter(
      data = df_resp,
      aes(
        x = stimulus,
        y = move,
        color = Genotype
      ),
      inherit.aes = FALSE,
      width = 0.15,
      height = 0.03,
      alpha = 0.12,
      size = 0.7
    ) +
    
    geom_ribbon(
      aes(
        ymin = CI_low,
        ymax = CI_high
      ),
      alpha = 0.15,
      color = NA
    ) +
    
    geom_line(
      aes(y = fit),
      linewidth = 1.3
    ) +
    
    coord_cartesian(ylim = c(0, 1)) +
    
    theme_pubr(base_size = 14) +
    
    labs(
      x = "Stimulus number within block",
      y = "Response probability",
      title = "Bayesian nonlinear habituation curves on probability scale",
      color = "Genotype",
      fill = "Genotype"
    ) +
    
    theme(
      legend.position = "top",
      panel.spacing = unit(1.2, "lines")
    )
  if (!is.null(save_fig_dir)) {
    ggsave(
      filename = file.path(
        save_fig_dir,
        paste0(
          "nlme_",
          var_name,
          "_habituation_curves_probability_scale_true_raw_binary.png"
        )
      ),
      plot = p,
      width = 14,
      height = 7,
      dpi = 300
    )
  }
  
  return(p)
}


make_prediction_grid <- function(df_resp, n_points = 100) {
  df_resp %>%
    group_by(Genotype, Block) %>%
    summarise(
      stim_min = min(stimulus, na.rm = TRUE),
      stim_max = max(stimulus, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Genotype, Block) %>%
    reframe(
      stimulus = seq(stim_min, stim_max, length.out = n_points)
    ) %>%
    mutate(
      stimulus0 = stimulus - 1,
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block = factor(Block, levels = levels(df_resp$Block))
    )
}


get_model_predictions <- function(fit_model, newdata, scale = NULL) {
  if (is.null(scale)) {
    pred <- fitted(
      fit_model,
      newdata = newdata,
      re_formula = NA,
      summary = TRUE
    )
  } else {
    pred <- fitted(
      fit_model,
      newdata = newdata,
      re_formula = NA,
      scale = scale,
      summary = TRUE
    )
  }
  
  bind_cols(newdata, as.data.frame(pred)) %>%
    rename(
      fit = Estimate,
      CI_low = Q2.5,
      CI_high = Q97.5
    )
}


plot_habituation_response_scale <- function(
    df_resp,
    fit_model,
    var_name,
    save_fig_dir = NULL,
    true_raw_binary = FALSE,
    n_points = 100,
    width = 14,
    height = 7,
    dpi = 300
) {
  new_resp <- make_prediction_grid(df_resp, n_points = n_points)
  
  pred_resp_data <- get_model_predictions(
    fit_model = fit_model,
    newdata = new_resp
  )
  
  if (true_raw_binary) {
    raw_layer <- geom_jitter(
      data = df_resp,
      aes(x = stimulus, y = move, color = Genotype),
      inherit.aes = FALSE,
      width = 0.15,
      height = 0.03,
      alpha = 0.12,
      size = 0.7
    )
    
    file_suffix <- "_habituation_curves_true_raw_data.png"
    
  } else {
    raw_resp <- df_resp %>%
      group_by(Genotype, Block, stimulus) %>%
      summarise(
        response_prob = mean(move, na.rm = TRUE),
        .groups = "drop"
      )
    
    raw_layer <- geom_point(
      data = raw_resp,
      aes(x = stimulus, y = response_prob),
      alpha = 0.35,
      size = 1,
      inherit.aes = FALSE
    )
    
    file_suffix <- "_habituation_curves.png"
  }
  
  p <- ggplot(
    pred_resp_data,
    aes(x = stimulus, color = Genotype, fill = Genotype)
  ) +
    facet_grid(Block ~ Genotype, scales = "fixed") +
    raw_layer +
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
      y = "Response probability",
      color = "Genotype",
      fill = "Genotype"
    ) +
    theme(
      legend.position = "top",
      panel.spacing = unit(1.2, "lines")
    )
  
  if (!is.null(save_fig_dir)) {
    ggsave(
      filename = file.path(save_fig_dir, paste0("nlme_", var_name, file_suffix)),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  return(p)
}


plot_habituation_logit_scale <- function(
    df_resp,
    fit_model,
    var_name,
    save_fig_dir = NULL,
    true_raw_binary = FALSE,
    n_points = 200,
    binary_display_values = c(-4, 4),
    width = 14,
    height = 7,
    dpi = 300
) {
  new_resp_logit <- make_prediction_grid(df_resp, n_points = n_points)
  
  pred_logit_data <- get_model_predictions(
    fit_model = fit_model,
    newdata = new_resp_logit,
    scale = "linear"
  )
  
  if (true_raw_binary) {
    raw_binary_logit <- df_resp %>%
      mutate(
        logit_binary_display = ifelse(
          move == 1,
          binary_display_values[2],
          binary_display_values[1]
        )
      )
    
    raw_layer <- geom_jitter(
      data = raw_binary_logit,
      aes(x = stimulus, y = logit_binary_display, color = Genotype),
      inherit.aes = FALSE,
      width = 0.15,
      height = 0.15,
      alpha = 0.12,
      size = 0.7
    )
    
    file_suffix <- "_habituation_curves_logit_scale_true_raw_binary.png"
    
  } else {
    raw_prob <- df_resp %>%
      group_by(Genotype, Block, stimulus) %>%
      summarise(
        response_prob = mean(move, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(
        response_prob_clipped = pmin(pmax(response_prob, 0.001), 0.999),
        logit_prob = qlogis(response_prob_clipped)
      )
    
    raw_layer <- geom_point(
      data = raw_prob,
      aes(x = stimulus, y = logit_prob, color = Genotype),
      inherit.aes = FALSE,
      alpha = 0.5,
      size = 1
    )
    
    file_suffix <- "_habituation_curves_logit_scale_with_raw.png"
  }
  
  p <- ggplot(
    pred_logit_data,
    aes(x = stimulus, color = Genotype, fill = Genotype)
  ) +
    facet_grid(Block ~ Genotype, scales = "fixed") +
    raw_layer +
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
      y = "Latent response tendency (logit scale)",
      title = "Bayesian nonlinear habituation curves on latent logit scale",
      color = "Genotype",
      fill = "Genotype"
    ) +
    theme(
      legend.position = "top",
      panel.spacing = unit(1.2, "lines")
    )
  
  if (!is.null(save_fig_dir)) {
    ggsave(
      filename = file.path(save_fig_dir, paste0("nlme_", var_name, file_suffix)),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  return(p)
}
