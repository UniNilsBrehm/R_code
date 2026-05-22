# ------------------------------------------------------------------------------
# Probability-scale plot using aggregated raw probabilities
# One panel per Block x Genotype
# ------------------------------------------------------------------------------

plot_habituation_probability <- function(
    df_resp,
    fit_model,
    save_fig_dir = NULL,
    var_name = "response_prob",
    raw_data = c("aggregate", "binary"),
    binary_alpha = 0.06,
    binary_size = 0.35,
    jitter_width = 0.15,
    jitter_height = 0.025
) {
  
  raw_data <- match.arg(raw_data)
  
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
    # robust = TRUE,  # use median
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
      n = n(),
      .groups = "drop"
    )
  
  binary_data <- df_resp %>%
    mutate(
      move = as.numeric(move),
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block = factor(Block, levels = levels(df_resp$Block))
    )
  
  p <- ggplot(
    pred_resp_data,
    aes(
      x = stimulus,
      color = Genotype,
      fill = Genotype
    )
  ) +
    
    facet_grid(Block ~ Genotype, scales = "fixed")
  
  if (raw_data == "aggregate") {
    p <- p +
      geom_point(
        data = raw_prob,
        aes(
          x = stimulus,
          y = response_prob,
          # size = n,  # show number of animals as dot size
          color = Genotype
        ),
        inherit.aes = FALSE,
        alpha = 0.5,
        size = 1
      )
  }
  
  if (raw_data == "binary") {
    p <- p +
      geom_jitter(
        data = binary_data,
        aes(
          x = stimulus,
          y = move,
          color = Genotype
        ),
        inherit.aes = FALSE,
        width = jitter_width,
        height = jitter_height,
        alpha = binary_alpha,
        size = binary_size
      )
  }
  
  p <- p +
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
          "_habituation_curves_probability_scale_",
          raw_data,
          ".png"
        )
      ),
      plot = p,
      width = 14,
      height = 7,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(
        save_fig_dir,
        paste0(
          "nlme_",
          var_name,
          "_habituation_curves_probability_scale_",
          raw_data,
          ".pdf"
        )
      ),
      plot = p,
      width = 14,
      height = 7,
      useDingbats = FALSE
    )
  }
  
  return(p)
}


# ------------------------------------------------------------------------------
# Latent-scale plot using posterior_linpred()
# One panel per Block x Genotype
# ------------------------------------------------------------------------------

plot_habituation_latent <- function(
    df_resp,
    fit_model,
    save_fig_dir = NULL,
    var_name = "response_prob",
    eps = 0.5,
    integer_predictions = FALSE
) {
  
  # ---------------------------------------------------------------------------
  # Prediction grid
  # ---------------------------------------------------------------------------
  
  if (integer_predictions) {
    
    # Predict only at actually observed stimulus numbers
    # Useful if you do not want artificial-looking interpolation between flashes.
    new_resp <- df_resp %>%
      distinct(Genotype, Block, stimulus) %>%
      mutate(
        stimulus0 = stimulus - 1,
        Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
        Block = factor(Block, levels = levels(df_resp$Block))
      ) %>%
      arrange(Block, Genotype, stimulus)
    
  } else {
    
    # Smooth curve over stimulus range
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
  }
  
  # ---------------------------------------------------------------------------
  # Posterior predictions on latent logit scale
  # ---------------------------------------------------------------------------
  # posterior_linpred(..., transform = FALSE) returns the linear predictor.
  # For bernoulli(link = "logit"), this is logit(p).
  
  pred_latent_draws <- posterior_linpred(
    fit_model,
    newdata = new_resp,
    re_formula = NA,
    transform = FALSE
  )
  
  pred_latent_summary <- apply(
    pred_latent_draws,
    2,
    quantile,
    probs = c(0.025, 0.5, 0.975),
    na.rm = TRUE
  )
  
  pred_resp_data <- bind_cols(
    new_resp,
    as.data.frame(t(pred_latent_summary))
  ) %>%
    rename(
      CI_low = `2.5%`,
      fit = `50%`,
      CI_high = `97.5%`
    )
  
  # ---------------------------------------------------------------------------
  # Aggregated raw probabilities transformed to empirical logit
  # ---------------------------------------------------------------------------
  # eps avoids infinite logits when all animals respond or no animals respond.
  # response_logit = log((responders + eps) / (nonresponders + eps))
  
  raw_latent <- df_resp %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      y = sum(move, na.rm = TRUE),
      n = sum(!is.na(move)),
      response_prob = y / n,
      response_logit = log((y + eps) / (n - y + eps)),
      .groups = "drop"
    ) %>%
    mutate(
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block = factor(Block, levels = levels(df_resp$Block))
    )
  
  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  
  p <- ggplot(
    pred_resp_data,
    aes(
      x = stimulus,
      color = Genotype,
      fill = Genotype
    )
  ) +
    
    facet_grid(Block ~ Genotype, scales = "fixed") +
    
    geom_point(
      data = raw_latent,
      aes(
        x = stimulus,
        y = response_logit,
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
    
    theme_pubr(base_size = 14) +
    
    labs(
      x = "Stimulus number within block",
      y = "Latent response tendency, logit(p)",
      title = "Bayesian nonlinear habituation curves on latent logit scale",
      color = "Genotype",
      fill = "Genotype"
    ) +
    
    theme(
      legend.position = "top",
      panel.spacing = unit(1.2, "lines")
    )
  
  if (!is.null(save_fig_dir)) {
    suffix <- ifelse(
      integer_predictions,
      "_latent_logit_scale_integer_predictions",
      "_latent_logit_scale_smooth"
    )
    
    ggsave(
      filename = file.path(
        save_fig_dir,
        paste0("nlme_", var_name, "_habituation_curves", suffix, ".png")
      ),
      plot = p,
      width = 14,
      height = 7,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(
        save_fig_dir,
        paste0("nlme_", var_name, "_habituation_curves", suffix, ".pdf")
      ),
      plot = p,
      width = 14,
      height = 7,
      useDingbats = FALSE
    )
  }
  
  return(p)
}
plot_habituation_probability_animal_averaged <- function(
    df_resp,
    fit_model,
    save_fig_dir = NULL,
    var_name = "response_prob",
    ndraws = NULL
) {

  # --------------------------------------------------------------------------
  # 1. Animal information
  # --------------------------------------------------------------------------
  
  animal_info <- df_resp %>%
    distinct(Genotype, Block, animal) %>%
    mutate(
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block = factor(Block, levels = levels(df_resp$Block)),
      animal = factor(animal, levels = levels(df_resp$animal))
    )
  
  # --------------------------------------------------------------------------
  # 2. Stimulus grid per Genotype × Block
  # --------------------------------------------------------------------------
  
  stim_grid <- df_resp %>%
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
  
  # --------------------------------------------------------------------------
  # 3. Cross real animals with stimulus grid within their Genotype × Block
  # --------------------------------------------------------------------------
  
  new_resp_animals <- animal_info %>%
    left_join(
      stim_grid,
      by = c("Genotype", "Block"),
      relationship = "many-to-many"
    )
  
  # --------------------------------------------------------------------------
  # 4. Posterior expected response probabilities including animal effects
  # --------------------------------------------------------------------------
  
  if (is.null(ndraws)) {
    epred_animals <- posterior_epred(
      fit_model,
      newdata = new_resp_animals,
      re_formula = NULL,
      summary = FALSE
    )
  } else {
    epred_animals <- posterior_epred(
      fit_model,
      newdata = new_resp_animals,
      re_formula = NULL,
      summary = FALSE,
      ndraws = ndraws
    )
  }
  
  # epred_animals: draws × rows
  # Average over animals within each Genotype × Block × stimulus for each draw
  
  row_info <- new_resp_animals %>%
    mutate(.row = row_number())
  
  epred_long <- as.data.frame(epred_animals) %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(
      cols = -.draw,
      names_to = ".row",
      values_to = "response_prob_pred"
    ) %>%
    mutate(
      .row = as.integer(gsub("V", "", .row))
    ) %>%
    left_join(row_info, by = ".row")
  
  pred_animal_avg <- epred_long %>%
    group_by(.draw, Genotype, Block, stimulus, stimulus0) %>%
    summarise(
      response_prob_pred = mean(response_prob_pred, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Genotype, Block, stimulus, stimulus0) %>%
    summarise(
      fit = mean(response_prob_pred, na.rm = TRUE),
      CI_low = quantile(response_prob_pred, 0.025, na.rm = TRUE),
      CI_high = quantile(response_prob_pred, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --------------------------------------------------------------------------
  # 5. Raw aggregate response probabilities
  # --------------------------------------------------------------------------
  
  raw_prob <- df_resp %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      response_prob = mean(move, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # --------------------------------------------------------------------------
  # 6. Plot
  # --------------------------------------------------------------------------
  
  p <- ggplot(
    pred_animal_avg,
    aes(
      x = stimulus,
      color = Genotype,
      fill = Genotype
    )
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
    
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.25)
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    
    theme_pubr(base_size = 14) +
    
    labs(
      x = "Stimulus number within block",
      y = "Response probability",
      title = "Animal-averaged Bayesian nonlinear habituation curves",
      color = "Genotype",
      fill = "Genotype"
    ) +
    
    theme(
      legend.position = "top",
      panel.spacing = unit(1.2, "lines")
    )
  
  if (!is.null(save_fig_dir)) {
    
    dir.create(save_fig_dir, recursive = TRUE, showWarnings = FALSE)
    
    ggsave(
      filename = file.path(
        save_fig_dir,
        paste0(
          "nlme_",
          var_name,
          "_habituation_curves_animal_averaged.png"
        )
      ),
      plot = p,
      width = 14,
      height = 7,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(
        save_fig_dir,
        paste0(
          "nlme_",
          var_name,
          "_habituation_curves_animal_averaged.pdf"
        )
      ),
      plot = p,
      width = 14,
      height = 7,
      useDingbats = FALSE
    )
  }
  
  return(p)
}