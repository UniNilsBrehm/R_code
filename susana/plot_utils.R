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

# ------------------------------------------------------------------------------
# Plot on animal average (not marginal/conditional but with animal variation)
# One panel per Block x Genotype
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Response-scale plot for any habituation response variable
# (e.g., max_peak, max_cumsum, summed_distance)
# One panel per Block x Genotype
# ------------------------------------------------------------------------------
plot_habituation_response <- function(
    df_resp,
    fit_model,
    response_var = "max_peak",          # column name in df_resp matching the model's LHS
    y_label      = "Distance moved",    # axis label
    title        = NULL,                # plot title; auto-generated if NULL
    save_fig_dir = NULL,
    var_name     = NULL,                # filename component; defaults to response_var
    raw_data     = c("aggregate", "trials"),
    trials_alpha  = 0.05,
    trials_size   = 0.4,
    jitter_width  = 0.10,
    jitter_height = 0,
    y_limits      = NULL                # e.g. c(0, NA) or c(0, 50); NULL = auto
) {
  
  raw_data <- match.arg(raw_data)
  if (is.null(var_name)) var_name <- response_var
  if (is.null(title))    title    <- paste0("Bayesian nonlinear habituation curves (", response_var, ")")
  
  # Guard: column exists?
  if (!response_var %in% names(df_resp)) {
    stop("Column '", response_var, "' not found in df_resp.")
  }
  
  # ---- prediction grid --------------------------------------------------------
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
      Genotype  = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block     = factor(Block,    levels = levels(df_resp$Block))
    )
  
  # ---- posterior fitted values (population-level, response scale) -------------
  pred_resp <- fitted(
    fit_model,
    newdata    = new_resp,
    re_formula = NA,
    summary    = TRUE
  )
  
  pred_resp_data <- bind_cols(
    new_resp,
    as.data.frame(pred_resp)
  ) %>%
    rename(
      fit     = Estimate,
      CI_low  = Q2.5,
      CI_high = Q97.5
    )
  
  # ---- aggregated raw data: mean ± SE per stimulus ----------------------------
  raw_agg <- df_resp %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      y_mean = mean(.data[[response_var]], na.rm = TRUE),
      y_se   = sd(.data[[response_var]],  na.rm = TRUE) /
        sqrt(sum(!is.na(.data[[response_var]]))),
      n      = sum(!is.na(.data[[response_var]])),
      .groups = "drop"
    )
  
  # ---- per-trial raw data -----------------------------------------------------
  trial_data <- df_resp %>%
    mutate(
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block    = factor(Block,    levels = levels(df_resp$Block))
    )
  
  # ---- build plot -------------------------------------------------------------
  p <- ggplot(
    pred_resp_data,
    aes(x = stimulus, color = Genotype, fill = Genotype)
  ) +
    facet_grid(Block ~ Genotype, scales = "fixed")
  
  if (raw_data == "aggregate") {
    p <- p +
      geom_errorbar(
        data = raw_agg,
        aes(x = stimulus, ymin = y_mean - y_se, ymax = y_mean + y_se, color = Genotype),
        inherit.aes = FALSE,
        width = 0.2, alpha = 0.6
      ) +
      geom_point(
        data = raw_agg,
        aes(x = stimulus, y = y_mean, color = Genotype),
        inherit.aes = FALSE,
        alpha = 0.7, size = 1.2
      )
  }
  
  if (raw_data == "trials") {
    p <- p +
      geom_jitter(
        data = trial_data,
        aes(x = stimulus, y = .data[[response_var]], color = Genotype),
        inherit.aes = FALSE,
        width  = jitter_width,
        height = jitter_height,
        alpha  = trials_alpha,
        size   = trials_size
      )
  }
  
  p <- p +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    theme_pubr(base_size = 14) +
    labs(
      x     = "Stimulus number within block",
      y     = y_label,
      title = title,
      color = "Genotype",
      fill  = "Genotype"
    ) +
    theme(
      legend.position = "top",
      panel.spacing   = unit(1.2, "lines")
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  # ---- save -------------------------------------------------------------------
  if (!is.null(save_fig_dir)) {
    dir.create(save_fig_dir, recursive = TRUE, showWarnings = FALSE)
    
    fname <- paste0("nlme_", var_name,
                    "_habituation_curves_response_scale_", raw_data)
    
    ggsave(file.path(save_fig_dir, paste0(fname, ".png")),
           plot = p, width = 14, height = 7, dpi = 300)
    ggsave(file.path(save_fig_dir, paste0(fname, ".pdf")),
           plot = p, width = 14, height = 7, useDingbats = FALSE)
  }
  
  return(p)
}


# ------------------------------------------------------------------------------
# Animal-averaged version for any habituation response variable
# ------------------------------------------------------------------------------
plot_habituation_response_animal_averaged <- function(
    df_resp,
    fit_model,
    response_var = "max_peak",
    y_label      = "Distance moved",
    title        = NULL,
    save_fig_dir = NULL,
    var_name     = NULL,
    ndraws       = NULL,
    y_limits     = NULL,
    y_breaks     = waiver()
) {
  
  if (is.null(var_name)) var_name <- response_var
  if (is.null(title))    title    <- paste0("Animal-averaged Bayesian nonlinear habituation curves (", response_var, ")")
  
  if (!response_var %in% names(df_resp)) {
    stop("Column '", response_var, "' not found in df_resp.")
  }
  
  # --------------------------------------------------------------------------
  # 1. Animal information
  # --------------------------------------------------------------------------
  animal_info <- df_resp %>%
    distinct(Genotype, Block, animal) %>%
    mutate(
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block    = factor(Block,    levels = levels(df_resp$Block)),
      animal   = factor(animal,   levels = levels(df_resp$animal))
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
      Genotype  = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block     = factor(Block,    levels = levels(df_resp$Block))
    )
  
  # --------------------------------------------------------------------------
  # 3. Cross real animals with stimulus grid within their Genotype × Block
  # --------------------------------------------------------------------------
  new_resp_animals <- animal_info %>%
    left_join(stim_grid, by = c("Genotype", "Block"),
              relationship = "many-to-many")
  
  # --------------------------------------------------------------------------
  # 4. Posterior expected response including animal effects
  # --------------------------------------------------------------------------
  epred_args <- list(
    object     = fit_model,
    newdata    = new_resp_animals,
    re_formula = NULL,
    summary    = FALSE
  )
  if (!is.null(ndraws)) epred_args$ndraws <- ndraws
  
  epred_animals <- do.call(posterior_epred, epred_args)
  
  row_info <- new_resp_animals %>%
    mutate(.row = row_number())
  
  epred_long <- as.data.frame(epred_animals) %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(
      cols      = -.draw,
      names_to  = ".row",
      values_to = "y_pred"
    ) %>%
    mutate(.row = as.integer(gsub("V", "", .row))) %>%
    left_join(row_info, by = ".row")
  
  pred_animal_avg <- epred_long %>%
    group_by(.draw, Genotype, Block, stimulus, stimulus0) %>%
    summarise(y_pred = mean(y_pred, na.rm = TRUE), .groups = "drop") %>%
    group_by(Genotype, Block, stimulus, stimulus0) %>%
    summarise(
      fit     = mean(y_pred, na.rm = TRUE),
      CI_low  = quantile(y_pred, 0.025, na.rm = TRUE),
      CI_high = quantile(y_pred, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --------------------------------------------------------------------------
  # 5. Raw aggregate: mean ± SE per stimulus
  # --------------------------------------------------------------------------
  raw_agg <- df_resp %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      y_mean = mean(.data[[response_var]], na.rm = TRUE),
      y_se   = sd(.data[[response_var]],  na.rm = TRUE) /
        sqrt(sum(!is.na(.data[[response_var]]))),
      n      = sum(!is.na(.data[[response_var]])),
      .groups = "drop"
    )
  
  # --------------------------------------------------------------------------
  # 6. Plot
  # --------------------------------------------------------------------------
  p <- ggplot(
    pred_animal_avg,
    aes(x = stimulus, color = Genotype, fill = Genotype)
  ) +
    facet_grid(Block ~ Genotype, scales = "fixed") +
    
    geom_errorbar(
      data = raw_agg,
      aes(x = stimulus, ymin = y_mean - y_se, ymax = y_mean + y_se, color = Genotype),
      inherit.aes = FALSE,
      width = 0.2, alpha = 0.6
    ) +
    
    geom_point(
      data = raw_agg,
      aes(x = stimulus, y = y_mean, color = Genotype),
      inherit.aes = FALSE,
      alpha = 0.7, size = 1.2
    ) +
    
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    
    scale_y_continuous(breaks = y_breaks) +
    
    theme_pubr(base_size = 14) +
    
    labs(
      x     = "Stimulus number within block",
      y     = y_label,
      title = title,
      color = "Genotype",
      fill  = "Genotype"
    ) +
    
    theme(
      legend.position = "top",
      panel.spacing   = unit(1.2, "lines")
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  # --------------------------------------------------------------------------
  # 7. Save
  # --------------------------------------------------------------------------
  if (!is.null(save_fig_dir)) {
    dir.create(save_fig_dir, recursive = TRUE, showWarnings = FALSE)
    
    fname <- paste0("nlme_", var_name, "_habituation_curves_animal_averaged")
    
    ggsave(file.path(save_fig_dir, paste0(fname, ".png")),
           plot = p, width = 14, height = 7, dpi = 300)
    ggsave(file.path(save_fig_dir, paste0(fname, ".pdf")),
           plot = p, width = 14, height = 7, useDingbats = FALSE)
  }
  
  return(p)
}


# ==============================================================================
# Ordinal helpers
# ==============================================================================

# Compute expected category from category probability array (draws × N × K)
# Returns draws × N matrix of expected delay
.expected_category_from_epred <- function(epred_array, category_values) {
  # epred_array: array [draws, N, K]
  # category_values: numeric vector of length K, the values 0:4 (or any ordinal labels)
  n_dim <- dim(epred_array)
  stopifnot(length(n_dim) == 3, n_dim[3] == length(category_values))
  
  # Multiply each category slice by its value and sum across categories
  exp_cat <- array(0, dim = c(n_dim[1], n_dim[2]))
  for (k in seq_along(category_values)) {
    exp_cat <- exp_cat + epred_array[, , k] * category_values[k]
  }
  exp_cat  # draws × N
}

# Convert ordered factor or 0:4 column to numeric
.as_numeric_response <- function(x) {
  if (is.ordered(x) || is.factor(x)) {
    as.integer(as.character(x))
  } else {
    as.numeric(x)
  }
}


# ==============================================================================
# Population-level expected category plot (analog of plot_habituation_response)
# ==============================================================================
plot_habituation_ordinal <- function(
    df_resp,
    fit_model,
    response_var     = "delay",
    category_values  = 0:4,
    y_label          = "Expected delay category",
    title            = NULL,
    save_fig_dir     = NULL,
    var_name         = NULL,
    raw_data         = c("aggregate", "trials"),
    trials_alpha     = 0.05,
    trials_size      = 0.4,
    jitter_width     = 0.10,
    jitter_height    = 0.05,
    y_limits         = NULL
) {
  
  raw_data <- match.arg(raw_data)
  if (is.null(var_name)) var_name <- response_var
  if (is.null(title))    title    <- paste0("Bayesian nonlinear habituation (", response_var, ", expected category)")
  
  if (!response_var %in% names(df_resp)) {
    stop("Column '", response_var, "' not found in df_resp.")
  }
  
  # ---- prediction grid -------------------------------------------------------
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
      Genotype  = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block     = factor(Block,    levels = levels(df_resp$Block))
    )
  
  # ---- posterior_epred returns draws × N × K -----------------------------------
  epred_array <- posterior_epred(
    fit_model,
    newdata    = new_resp,
    re_formula = NA,
    summary    = FALSE
  )
  
  # Collapse over categories into expected delay (draws × N)
  exp_cat <- .expected_category_from_epred(epred_array, category_values)
  
  pred_resp_data <- new_resp %>%
    mutate(
      fit     = apply(exp_cat, 2, mean),
      CI_low  = apply(exp_cat, 2, quantile, 0.025),
      CI_high = apply(exp_cat, 2, quantile, 0.975)
    )
  
  # ---- aggregated raw data ---------------------------------------------------
  df_num <- df_resp %>%
    mutate(.y = .as_numeric_response(.data[[response_var]]))
  
  raw_agg <- df_num %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      y_mean = mean(.y, na.rm = TRUE),
      y_se   = sd(.y,  na.rm = TRUE) / sqrt(sum(!is.na(.y))),
      n      = sum(!is.na(.y)),
      .groups = "drop"
    )
  
  # ---- per-trial raw data ----------------------------------------------------
  trial_data <- df_num %>%
    mutate(
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block    = factor(Block,    levels = levels(df_resp$Block))
    )
  
  # ---- build plot ------------------------------------------------------------
  p <- ggplot(
    pred_resp_data,
    aes(x = stimulus, color = Genotype, fill = Genotype)
  ) +
    facet_grid(Block ~ Genotype, scales = "fixed")
  
  if (raw_data == "aggregate") {
    p <- p +
      geom_errorbar(
        data = raw_agg,
        aes(x = stimulus, ymin = y_mean - y_se, ymax = y_mean + y_se, color = Genotype),
        inherit.aes = FALSE,
        width = 0.2, alpha = 0.6
      ) +
      geom_point(
        data = raw_agg,
        aes(x = stimulus, y = y_mean, color = Genotype),
        inherit.aes = FALSE,
        alpha = 0.7, size = 1.2
      )
  }
  
  if (raw_data == "trials") {
    p <- p +
      geom_jitter(
        data = trial_data,
        aes(x = stimulus, y = .y, color = Genotype),
        inherit.aes = FALSE,
        width  = jitter_width,
        height = jitter_height,
        alpha  = trials_alpha,
        size   = trials_size
      )
  }
  
  p <- p +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    theme_pubr(base_size = 14) +
    labs(
      x     = "Stimulus number within block",
      y     = y_label,
      title = title,
      color = "Genotype",
      fill  = "Genotype"
    ) +
    theme(
      legend.position = "top",
      panel.spacing   = unit(1.2, "lines")
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  # ---- save ------------------------------------------------------------------
  if (!is.null(save_fig_dir)) {
    dir.create(save_fig_dir, recursive = TRUE, showWarnings = FALSE)
    fname <- paste0("nlme_", var_name,
                    "_habituation_curves_expected_category_", raw_data)
    ggsave(file.path(save_fig_dir, paste0(fname, ".png")),
           plot = p, width = 14, height = 7, dpi = 300)
    ggsave(file.path(save_fig_dir, paste0(fname, ".pdf")),
           plot = p, width = 14, height = 7, useDingbats = FALSE)
  }
  
  return(p)
}


# ==============================================================================
# Animal-averaged expected category plot
# ==============================================================================
plot_habituation_ordinal_animal_averaged <- function(
    df_resp,
    fit_model,
    response_var     = "delay",
    category_values  = 0:4,
    y_label          = "Expected delay category",
    title            = NULL,
    save_fig_dir     = NULL,
    var_name         = NULL,
    ndraws           = NULL,
    y_limits         = NULL,
    y_breaks         = waiver()
) {
  
  if (is.null(var_name)) var_name <- response_var
  if (is.null(title))    title    <- paste0("Animal-averaged Bayesian nonlinear habituation (", response_var, ", expected category)")
  
  if (!response_var %in% names(df_resp)) {
    stop("Column '", response_var, "' not found in df_resp.")
  }
  
  # ---- 1. Animals ------------------------------------------------------------
  animal_info <- df_resp %>%
    distinct(Genotype, Block, animal) %>%
    mutate(
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block    = factor(Block,    levels = levels(df_resp$Block)),
      animal   = factor(animal,   levels = levels(df_resp$animal))
    )
  
  # ---- 2. Stimulus grid ------------------------------------------------------
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
      Genotype  = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block     = factor(Block,    levels = levels(df_resp$Block))
    )
  
  # ---- 3. Cross animals with grid --------------------------------------------
  new_resp_animals <- animal_info %>%
    left_join(stim_grid, by = c("Genotype", "Block"),
              relationship = "many-to-many")
  
  # ---- 4. Posterior epred with animal effects --------------------------------
  epred_args <- list(
    object     = fit_model,
    newdata    = new_resp_animals,
    re_formula = NULL,
    summary    = FALSE
  )
  if (!is.null(ndraws)) epred_args$ndraws <- ndraws
  
  epred_array <- do.call(posterior_epred, epred_args)
  # epred_array: draws × N_rows × K
  
  # Collapse categories into expected delay: draws × N_rows
  exp_cat <- .expected_category_from_epred(epred_array, category_values)
  
  # Average over animals within each Genotype × Block × stimulus, per draw
  row_info <- new_resp_animals %>%
    mutate(.row = row_number())
  
  exp_cat_long <- as.data.frame(exp_cat) %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(
      cols      = -.draw,
      names_to  = ".row",
      values_to = "y_pred"
    ) %>%
    mutate(.row = as.integer(gsub("V", "", .row))) %>%
    left_join(row_info, by = ".row")
  
  pred_animal_avg <- exp_cat_long %>%
    group_by(.draw, Genotype, Block, stimulus, stimulus0) %>%
    summarise(y_pred = mean(y_pred, na.rm = TRUE), .groups = "drop") %>%
    group_by(Genotype, Block, stimulus, stimulus0) %>%
    summarise(
      fit     = mean(y_pred, na.rm = TRUE),
      CI_low  = quantile(y_pred, 0.025, na.rm = TRUE),
      CI_high = quantile(y_pred, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- 5. Raw aggregate ------------------------------------------------------
  df_num <- df_resp %>%
    mutate(.y = .as_numeric_response(.data[[response_var]]))
  
  raw_agg <- df_num %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      y_mean = mean(.y, na.rm = TRUE),
      y_se   = sd(.y,  na.rm = TRUE) / sqrt(sum(!is.na(.y))),
      n      = sum(!is.na(.y)),
      .groups = "drop"
    )
  
  # ---- 6. Plot ---------------------------------------------------------------
  p <- ggplot(
    pred_animal_avg,
    aes(x = stimulus, color = Genotype, fill = Genotype)
  ) +
    facet_grid(Block ~ Genotype, scales = "fixed") +
    
    geom_errorbar(
      data = raw_agg,
      aes(x = stimulus, ymin = y_mean - y_se, ymax = y_mean + y_se, color = Genotype),
      inherit.aes = FALSE,
      width = 0.2, alpha = 0.6
    ) +
    
    geom_point(
      data = raw_agg,
      aes(x = stimulus, y = y_mean, color = Genotype),
      inherit.aes = FALSE,
      alpha = 0.7, size = 1.2
    ) +
    
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    
    scale_y_continuous(breaks = y_breaks) +
    
    theme_pubr(base_size = 14) +
    
    labs(
      x     = "Stimulus number within block",
      y     = y_label,
      title = title,
      color = "Genotype",
      fill  = "Genotype"
    ) +
    
    theme(
      legend.position = "top",
      panel.spacing   = unit(1.2, "lines")
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  # ---- 7. Save ---------------------------------------------------------------
  if (!is.null(save_fig_dir)) {
    dir.create(save_fig_dir, recursive = TRUE, showWarnings = FALSE)
    fname <- paste0("nlme_", var_name, "_habituation_curves_animal_averaged_expected_category")
    ggsave(file.path(save_fig_dir, paste0(fname, ".png")),
           plot = p, width = 14, height = 7, dpi = 300)
    ggsave(file.path(save_fig_dir, paste0(fname, ".pdf")),
           plot = p, width = 14, height = 7, useDingbats = FALSE)
  }
  
  return(p)
}