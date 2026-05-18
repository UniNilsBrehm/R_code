load_data_darkflash_60s <- function(file_dir, move_th = 1, keep = NULL, take_peak = 4) {
  # max_peak is distance moved
  # Load data from csv file
  df <- read_csv(file_dir)
  df_reduced <- df[, c("Block", "Well", "Video", "Peak", "Genotype", "Stimulus_New", "max_peak", "max_cumsum", "peak_maxdist")]
  
  # Print peaks
  message("Found peak numbers:")
  print(unique(df_reduced$Peak))
  
  # rename peak_maxdist to delay
  names(df_reduced)[names(df_reduced) == "peak_maxdist"] <- "delay"
  
  if (!is.null(keep)) {
    df_reduced <- df_reduced %>%
      dplyr::filter(Genotype %in% keep)
    message("Keeping only genotypes:")
    print(paste(keep, collapse = " -- "))
  }
  
  # Prepare for GLMM
  df_reduced <- df_reduced %>%
    mutate(
      Genotype = factor(Genotype),
      Genotype = relevel(Genotype, ref = "ABTL"),
    )
  
  df_reduced$Well <- as.factor(df_reduced$Well)
  df_reduced$Video <- as.factor(df_reduced$Video)
  df_reduced$Block <- as.factor(df_reduced$Block)
  df_reduced$Peak <- as.numeric(df_reduced$Peak)
  df_reduced$Stimulus_New <- as.numeric(df_reduced$Stimulus_New)
  
  df_reduced <- df_reduced %>%
    mutate(
      stimulus = Stimulus_New,
    )
  
  # For each Peak (-10 to 4) the distance move value (max_peak) is the same
  # so remove all duplicates
  df_final <- df_reduced %>%
    group_by(Block, Well, Video, Stimulus_New) %>%
    # slice(1) %>%          # or slice_head(n = 1), takes the first
    slice_tail(n = 1) %>%   # takes the last row in each group
    ungroup()
  
  df_final <- df_reduced %>%
    group_by(Block, Well, Video, Stimulus_New) %>%
    filter(Peak == take_peak) %>%
    ungroup()
  
  # Compute log of stimulus number (for exponential-like model)
  df_final <- df_final %>%
    mutate(
      stimulus_log = log(stimulus)
    )
  
  # Get responses and non-responses (response: max_peak > 0)
  df_final$move <- ifelse(df_final$max_peak > move_th, 1, 0)
  df_final_sub <- subset(df_final, move > 0)
  
  return(list(df_final = df_final, df_final_sub = df_final_sub))
}

model_residuals_check  <- function(model, df) {
  # Load required packages
  require(DHARMa)
  
  # Simulate residuals
  res <- simulateResiduals(model)
  
  # Diagnostic plots
  plot(res)
}

validate_model <- function(model, df) {
  # Load required packages
  require(DHARMa)
  require(performance)
  require(parameters)
  require(see)
  
  # Simulate residuals
  res <- simulateResiduals(model)
  
  # Diagnostic plots
  plot(res)
  plotResiduals(res, df$Genotype)
  plotResiduals(res, df$Block)
  # plotResiduals(res, df$stimulus_log)
  
  # Residual tests
  cat("\n--- DHARMa Tests ---\n")
  # print(testUniformity(res))
  # print(testOutliers(res, type = "bootstrap"))
  print(testDispersion(res))
  # print(testQuantiles(res))
  print(testZeroInflation(res))
  
  # Model performance metrics
  cat("\n--- Model Performance ---\n")
  print(model_performance(model))
  
  # Collinearity diagnostics
  cat("\n--- Collinearity ---\n")
  #print(check_collinearity(model))
  print(check_collinearity(update(model, . ~ Genotype + Block + stimulus_log)))
  
  # Return results invisibly
  invisible(list(
    residuals = res,
    performance = model_performance(model),
    collinearity = check_collinearity(model)
  ))
}

plot_habituation <- function(df_final, model, label,
                             transform = c("plogis", "exp", "none"),
                             raw_var = "max_cumsum") {
  # --- 0. Setup ----------------------------------------------
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggpubr)
  
  # Match the transformation argument
  transform <- match.arg(transform)
  
  # --- 1. Raw summary (depends on what you're plotting) -------
  # If you want probability: use mean(move)
  # Otherwise: use mean(raw_var)  (default: max_cumsum)
  raw_summary <- df_final %>%
    group_by(Block, Genotype, stimulus) %>%
    summarise(
      y_raw = if (transform == "plogis") {
        mean(move, na.rm = TRUE)
      } else {
        mean(.data[[raw_var]], na.rm = TRUE)
      },
      .groups = "drop"
    )
  
  # --- 2. Dynamic prediction grid ----------------------------
  blocks <- sort(unique(df_final$Block))
  
  stim_ranges <- lapply(blocks, function(b) {
    range(unique(df_final$stimulus[df_final$Block == b]))
  })
  names(stim_ranges) <- blocks
  
  new_data <- bind_rows(lapply(blocks, function(b) {
    expand.grid(
      Block = as.character(b),
      stimulus = stim_ranges[[as.character(b)]][1]:
        stim_ranges[[as.character(b)]][2]
    )
  })) %>%
    tidyr::crossing(Genotype = unique(df_final$Genotype)) %>%
    mutate(stimulus_log = log(stimulus))
  
  # Ensure no log(0)
  new_data <- new_data %>%
    mutate(
      stimulus = ifelse(stimulus <= 0, NA, stimulus),
      stimulus_log = log(stimulus)
    ) %>%
    filter(!is.na(stimulus_log))
  
  # Align factor levels
  new_data$Genotype <- factor(new_data$Genotype, levels = levels(df_final$Genotype))
  new_data$Block    <- factor(new_data$Block,    levels = levels(df_final$Block))
  
  # --- 3. Model predictions ----------------------------------
  pred <- predict(model, newdata = new_data, re.form = NA, se.fit = TRUE)
  
  transform_fun <- switch(transform,
                          plogis = plogis,
                          exp    = exp,
                          none   = identity)
  
  new_data <- new_data %>%
    mutate(
      fit     = transform_fun(pred$fit),
      CI_low  = transform_fun(pred$fit - 1.96 * pred$se.fit),
      CI_high = transform_fun(pred$fit + 1.96 * pred$se.fit)
    )
  
  # --- 4. Facet labels ---------------------------------------
  labels <- setNames(
    paste0("Block ", blocks, ": ", sapply(stim_ranges, `[`, 2), " flashes"),
    blocks
  )
  
  # --- 5. Plot -----------------------------------------------
  g <- ggplot(new_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
    facet_wrap(~Block, ncol = length(blocks), scales = "free_x",
               labeller = as_labeller(labels)) +
    geom_point(
      data = raw_summary,
      aes(x = stimulus, y = y_raw, color = Genotype),
      size = 1,
      alpha = 0.26
    ) +
    geom_line(aes(y = fit), size = 1.5) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
                alpha = 0.1, color = NA) +
    scale_x_continuous(n.breaks = 4) +
    labs(
      x = "Stimulus number (within block)",
      y = label,
      color = "Genotype",
      fill  = "Genotype"
    ) +
    theme_pubr(base_size = 14) +
    theme(
      panel.spacing = unit(1.5, "lines")
    )
  
  # Keep probability plots in [0,1]; otherwise don't constrain y
  if (transform == "plogis") {
    g <- g + expand_limits(y = c(0, 1))
  }
  
  return(g)
}


plot_habituation_glmm_by_genotype_block <- function(df_final, model,
                                                    label = "Peak distance moved (mm)",
                                                    raw_var = "max_peak",
                                                    transform = c("exp", "none", "plogis"),
                                                    n_points = 100) {
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggpubr)
  
  transform <- match.arg(transform)
  
  transform_fun <- switch(
    transform,
    exp = exp,
    none = identity,
    plogis = plogis
  )
  
  raw_summary <- df_final %>%
    group_by(Genotype, Block, stimulus) %>%
    summarise(
      y_raw = mean(.data[[raw_var]], na.rm = TRUE),
      .groups = "drop"
    )
  
  new_data <- df_final %>%
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
      stimulus_log = log(stimulus),
      Genotype = factor(Genotype, levels = levels(df_final$Genotype)),
      Block = factor(Block, levels = levels(df_final$Block))
    )
  
  pred <- predict(
    model,
    newdata = new_data,
    type = "link",
    # type = predict_type,
    se.fit = TRUE,
    re.form = NA
  )
  
  new_data <- new_data %>%
    mutate(
      fit = transform_fun(pred$fit),
      CI_low = transform_fun(pred$fit - 1.96 * pred$se.fit),
      CI_high = transform_fun(pred$fit + 1.96 * pred$se.fit)
    )
  
  ggplot(new_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
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
    scale_x_continuous(n.breaks = 4) +
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

save_plot <- function(g_plot, file_dir, width=10, height=7, dpi=600){
  ggsave(
    filename = paste0(file_dir, ".pdf"),
    plot = g_plot,
    width = width,
    height = height,
    units = "in"
  )
  ggsave(
    filename = paste0(file_dir, ".jpg"),
    dpi = dpi,
    plot = g_plot,
    width = width,
    height = height,
    units = "in"
  )
}