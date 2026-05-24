plot_habituation_probability <- function(
    df_resp,
    fit_model,
    save_fig_dir = NULL,
    var_name = "response_prob",
    raw_data = c("aggregate", "binary"),
    binary_alpha = 0.06,
    binary_size = 0.35,
    jitter_width = 0.15,
    jitter_height = 0.025,
    stim0_offset = 0,                    # NEW: 0 for massed, 1 for dark-flash
    n_grid = 300,                        # NEW: more grid points
    facet_layout = c("genotype_rows",    # NEW: layout chooser
                     "block_rows",
                     "wrap")
) {
  
  raw_data     <- match.arg(raw_data)
  facet_layout <- match.arg(facet_layout)
  
  new_resp <- df_resp %>%
    group_by(Genotype, Block) %>%
    summarise(
      stim_min = min(stimulus, na.rm = TRUE),
      stim_max = max(stimulus, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    reframe(
      stimulus = seq(stim_min, stim_max, length.out = n_grid),
      .by = c(Genotype, Block)
    ) %>%
    mutate(
      stimulus0 = stimulus - stim0_offset,
      Genotype  = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block     = factor(Block,    levels = levels(df_resp$Block))
    )
  
  pred_resp <- fitted(
    fit_model,
    newdata    = new_resp,
    re_formula = NA,
    summary    = TRUE
  )
  
  pred_resp_data <- bind_cols(new_resp, as.data.frame(pred_resp)) %>%
    rename(fit = Estimate, CI_low = Q2.5, CI_high = Q97.5)
  
  raw_prob <- df_resp %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      response_prob = mean(move, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  binary_data <- df_resp %>%
    mutate(
      move     = as.numeric(move),
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block    = factor(Block,    levels = levels(df_resp$Block))
    )
  
  p <- ggplot(
    pred_resp_data,
    aes(x = stimulus, color = Genotype, fill = Genotype)
  )
  
  # ----- Facet layout ---------------------------------------------------------
  if (facet_layout == "genotype_rows") {
    p <- p + facet_grid(Genotype ~ Block, scales = "free_x")
  } else if (facet_layout == "block_rows") {
    p <- p + facet_grid(Block ~ Genotype, scales = "free_x")
  } else {
    p <- p + facet_wrap(Genotype ~ Block, scales = "free_x", ncol = 2)
  }
  
  # ----- Raw data layer -------------------------------------------------------
  if (raw_data == "aggregate") {
    p <- p + geom_point(
      data = raw_prob,
      aes(x = stimulus, y = response_prob, color = Genotype),
      inherit.aes = FALSE, alpha = 0.5, size = 1
    )
  }
  
  if (raw_data == "binary") {
    p <- p + geom_jitter(
      data = binary_data,
      aes(x = stimulus, y = move, color = Genotype),
      inherit.aes = FALSE,
      width  = jitter_width,
      height = jitter_height,
      alpha  = binary_alpha,
      size   = binary_size
    )
  }
  
  # ----- Model layer ----------------------------------------------------------
  p <- p +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.15, color = NA) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_pubr(base_size = 14) +
    labs(
      x     = "Stimulus number within block",
      y     = "Response probability",
      title = "Bayesian nonlinear habituation curves on probability scale",
      color = "Genotype",
      fill  = "Genotype"
    ) +
    theme(legend.position = "top", panel.spacing = unit(1.2, "lines"))
  
  if (!is.null(save_fig_dir)) {
    ggsave(
      filename = file.path(save_fig_dir,
                           paste0("nlme_", var_name, "_habituation_curves_probability_scale_", raw_data, ".png")),
      plot = p, width = 14, height = 7, dpi = 300
    )
    ggsave(
      filename = file.path(save_fig_dir,
                           paste0("nlme_", var_name, "_habituation_curves_probability_scale_", raw_data, ".pdf")),
      plot = p, width = 14, height = 7, useDingbats = FALSE
    )
  }
  
  return(p)
}