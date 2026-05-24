# PLOT FUNCTIONS
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

plot_habituation22 <- function(df_final, model, label,
                              transform = c("plogis", "exp", "none")) {
  # --- 0. Setup ----------------------------------------------
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggpubr)
  
  # Match the transformation argument
  transform <- match.arg(transform)
  
  # --- 1. Raw summary (depends on what you're plotting) -------
  # If you want probability: use mean(move)
  # If you want distance (max_cumsum): use mean(max_cumsum)
  raw_summary <- df_final %>%
    group_by(Block, Genotype, stimulus) %>%
    summarise(
      y_raw = if (transform == "plogis") mean(move, na.rm = TRUE) else mean(max_cumsum, na.rm = TRUE),
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
      # shape = 1,   # hollow circle
      size = 1,
      alpha = 0.26
    ) +
    geom_line(aes(y = fit), size = 1.5) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
                alpha = 0.1, color = NA) +
    scale_x_continuous(n.breaks = 4)+
    labs(
      x = "Stimulus number (within block)",
      y = label,
      color = "Genotype",
      fill  = "Genotype"
    ) +
    theme_pubr(base_size = 14) +
    theme(
      panel.spacing = unit(1.5, "lines")  # increase spacing
    )
  
  # Keep probability plots in [0,1]; otherwise don't constrain y
  if (transform == "plogis") {
    g <- g + expand_limits(y = c(0, 1))
  }
  
  return(g)
}

plot_habituation_old <- function(df_final, model, label, transform = c("plogis", "exp", "none")) {
  # --- 0. Setup ----------------------------------------------
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggpubr)
  
  raw_summary <- df_final %>%
    group_by(Block, Genotype, stimulus) %>%
    summarise(prob = mean(move), .groups = "drop")
  
  # Match the transformation argument
  transform <- match.arg(transform)
  
  # --- 1. Dynamic prediction grid ----------------------------
  blocks <- sort(unique(df_final$Block))
  
  # Determine stimulus range for each block dynamically
  stim_ranges <- lapply(blocks, function(b) {
    range(unique(df_final$stimulus[df_final$Block == b]))
  })
  names(stim_ranges) <- blocks
  
  # Build the prediction grid dynamically
  new_data <- bind_rows(lapply(blocks, function(b) {
    expand.grid(Block = as.character(b),
                stimulus = stim_ranges[[as.character(b)]][1]:
                  stim_ranges[[as.character(b)]][2])
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
  new_data$Block <- factor(new_data$Block, levels = levels(df_final$Block))
  
  # --- 2. Model predictions ----------------------------------
  pred <- predict(model, newdata = new_data, re.form = NA, se.fit = TRUE)
  
  # Choose transformation
  transform_fun <- switch(transform,
                          plogis = plogis,
                          exp    = exp,
                          none   = identity
  )
  
  new_data <- new_data %>%
    mutate(
      fit     = transform_fun(pred$fit),
      CI_low  = transform_fun(pred$fit - 1.96 * pred$se.fit),
      CI_high = transform_fun(pred$fit + 1.96 * pred$se.fit)
    )
  
  # --- 3. Facet labels ---------------------------------------
  labels <- setNames(
    paste0("Block ", blocks, ": ", sapply(stim_ranges, `[`, 2), " flashes"),
    blocks
  )
  
  # --- 4. Plot -----------------------------------------------
  g <- ggplot(new_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
    facet_wrap(~Block, ncol = length(blocks), scales = "free_x",
               labeller = as_labeller(labels)) +
    geom_point(data = raw_summary,
               aes(x = stimulus, y = prob, color = Genotype),
               alpha = 0.5, size = 1.5)+
    geom_line(aes(y = fit), linewidth = 1.2) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
                alpha = 0.1, color = NA) +
    expand_limits(y = c(0, 1)) +
    labs(
      x = "Stimulus number (within block)",
      y = label,
      color = "Genotype",
      fill  = "Genotype"
    ) +
    theme_pubr(base_size = 14)
  
  return(g)
}

# UTILITY FUNCTIONS
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

load_data_ISI60_DF <- function(file_dir, move_th = 1, keep = NULL) {
  
  # max_peak is distance moved
  
  # Load data from csv file
  df <- read_csv(file_dir)
  df_reduced <- df[, c("Block", "Well", "Video", "Peak", "Genotype", "Stimulus_New", "max_peak", "max_cumsum", "peak_maxdist")]
  
  names(df_reduced)[names(df_reduced) == "peak_maxdist"] <- "delay"
  
  if (!is.null(keep)) {
    # keep <- c("ABTL", "th, th2, tyr", "th, tyr")
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
      stimulus = Stimulus_New + 1,
    )
  
  df_final <- df_reduced %>%
    filter(Peak == 1)
  
  df_final <- df_final %>%
    mutate(
      stimulus_log = log(stimulus)
    )
  
  # Get responses and non-responses (response: max_peak > 0)
  df_final$move <- ifelse(df_final$max_peak > move_th, 1, 0)
  df_final_sub <- subset(df_final, move > 0)
  
  return(list(df_final = df_final, df_final_sub = df_final_sub))
}


load_data <- function(file_dir, move_th = 1, drop = NULL) {
  
  # max_peak is distance moved
  
  # Load data from csv file
  df <- read_csv(file_dir)
  
  df_reduced <- df[, c(
    "Block", "Well", "Video", "Peak", "Genotype",
    "Stimulus_New", "max_peak", "max_cumsum", "peak_maxdist"
  )]
  
  names(df_reduced)[names(df_reduced) == "peak_maxdist"] <- "delay"
  
  if (!is.null(drop)) {
    
    df_reduced <- df_reduced %>%
      dplyr::filter(!Genotype %in% drop)
    
    message("Dropping genotypes:")
    print(paste(drop, collapse = " -- "))
  }
  
  # Prepare for GLMM
  df_reduced <- df_reduced %>%
    mutate(
      Genotype = factor(Genotype),
      Genotype = relevel(Genotype, ref = "ABTL")
    )
  
  df_reduced$Well <- as.factor(df_reduced$Well)
  df_reduced$Video <- as.factor(df_reduced$Video)
  df_reduced$Block <- as.factor(df_reduced$Block)
  df_reduced$Peak <- as.numeric(df_reduced$Peak)
  df_reduced$Stimulus_New <- as.numeric(df_reduced$Stimulus_New)
  
  df_reduced <- df_reduced %>%
    mutate(
      stimulus = Stimulus_New + 1
    )
  
  df_final <- df_reduced %>%
    filter(Peak == 1)
  
  df_final <- df_final %>%
    mutate(
      stimulus_log = log(stimulus)
    )
  
  # Get responses and non-responses
  # response: max_peak > move_th
  df_final$move <- ifelse(df_final$max_peak > move_th, 1, 0)
  
  df_final_sub <- subset(df_final, move > 0)
  
  return(list(
    df_final = df_final,
    df_final_sub = df_final_sub
  ))
}


pretty_pairs <- function(emm) {
  df <- as.data.frame(pairs(emm))
  
  # Normalize column names
  if ("odds.ratio" %in% names(df)) {
    df <- dplyr::rename(df, estimate = odds.ratio)
  }
  
  df |>
    dplyr::mutate(
      estimate = round(estimate, 3),
      SE = round(SE, 3),
      z.ratio = round(z.ratio, 2),
      p.value = round(p.value, 4),
      sig = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE ~ ""
      )
    ) |>
    dplyr::select(
      dplyr::any_of(c("Block", "Genotype")),
      contrast, estimate, SE, z.ratio, p.value, sig
    )
}

pretty_blocks <- function(x) {
  # Convert to data.frame and add significance stars
  df <- as.data.frame(x)
  
  # Extract Genotype info if nested contrast object
  if ("Genotype" %in% names(attributes(x))) {
    df$Genotype <- rep(attr(x, "by.vars")[[1]], each = nrow(df) / length(attr(x, "by.vars")[[1]]))
  } else if (!"Genotype" %in% names(df)) {
    df$Genotype <- rep(unique(x@grid$Genotype), each = nrow(df) / length(unique(x@grid$Genotype)))
  }
  
  # Add significance stars
  df <- df %>%
    dplyr::mutate(
      sig = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(Genotype, contrast, odds.ratio, SE, z.ratio, p.value, sig)
  
  df
}

check_spread <- function(residuals_obj, group, data = NULL) {
  # Ensure required package
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package 'car' is required for Levene's test. Please install it with install.packages('car').")
  }
  
  # Extract residuals
  res <- residuals(residuals_obj)
  
  # Handle both vector or column name for group
  if (is.character(substitute(group))) {
    group_name <- deparse(substitute(group))
    group_var <- data[[group_name]]
  } else {
    group_var <- group
    group_name <- deparse(substitute(group))
  }
  
  # Run Levene's Test (center = median, more robust)
  lv <- car::leveneTest(res ~ group_var, center = "median")
  
  # Extract test info
  Fval <- lv$`F value`[1]
  df1  <- lv$Df[1]
  df2  <- lv$Df[2]
  pval <- lv$`Pr(>F)`[1]
  eta2_partial <- (Fval * df1) / (Fval * df1 + df2)
  
  # Print formatted summary
  cat("\n--- Levene's Test for Homogeneity of Variance ---\n")
  cat(sprintf("Grouping variable: %s\n", group_name))
  cat(sprintf("F(%d, %d) = %.3f, p = %.4f\n", df1, df2, Fval, pval))
  cat(sprintf("Partial eta² = %.5f\n", eta2_partial))
  
  # Interpret magnitude
  magnitude <- dplyr::case_when(
    eta2_partial < 0.01 ~ "very small",
    eta2_partial < 0.06 ~ "small",
    eta2_partial < 0.14 ~ "medium",
    TRUE ~ "large"
  )
  cat(sprintf("Interpretation: %s effect size\n", magnitude))
  
  # Return invisibly (for programmatic use)
  invisible(list(
    levene = lv,
    eta2_partial = eta2_partial,
    F = Fval,
    df = c(df1, df2),
    p.value = pval,
    magnitude = magnitude
  ))
}


get_emm_blocks <- function(model, df, n_middle = 10, n_last = 10, n_first = 10) {
  require(emmeans)
  
  # Detect blocks
  blocks <- sort(unique(df$Block))
  
  # Prepare results list
  emm_results <- list()
  
  for (b in blocks) {
    # Get all stimulus numbers for this block
    stim_block <- sort(unique(df$stimulus[df$Block == b]))
    n_stim <- length(stim_block)
    
    # --- Helper function to compute both emmeans and contrasts -------------
    make_emm_and_pairs <- function(idx, label) {
      emm <- emmeans(
        model,
        ~ Genotype,
        at = list(
          stimulus_log = log(idx),
          Block = as.character(b)
        ),
        cov.reduce = mean,
        type = "response"
      )
      
      pairwise <- pairs(emm)  # pairwise contrasts between genotypes
      list(emm = emm, pairs = pairwise)
    }
    
    # --- Full block --------------------------------------------------------
    emm_results[[paste0("Block", b, "_full")]] <- make_emm_and_pairs(stim_block, "full")
    
    # --- First stimuli -----------------------------------------------------
    first_idx <- stim_block[1:min(n_first, n_stim)]
    emm_results[[paste0("Block", b, "_first")]] <- make_emm_and_pairs(first_idx, "first")
    
    # --- Middle stimuli ----------------------------------------------------
    mid_start <- floor(n_stim / 2 - n_middle / 2)
    mid_idx <- stim_block[
      pmax(1, mid_start):(pmin(n_stim, mid_start + n_middle))
    ]
    emm_results[[paste0("Block", b, "_middle")]] <- make_emm_and_pairs(mid_idx, "middle")
    
    # --- Last stimuli ------------------------------------------------------
    last_idx <- tail(stim_block, n_last)
    emm_results[[paste0("Block", b, "_last")]] <- make_emm_and_pairs(last_idx, "last")
  }
  
  return(emm_results)
}


get_emm_consecutive_blocks <- function(model, df, n_stim = 9, compare_single = TRUE) {
  # Ensure required package
  require(emmeans)
  
  # Detect all block levels and sort them
  blocks <- sort(unique(df$Block))
  
  # Prepare result container
  emm_results <- list()
  
  # Iterate over consecutive block pairs: (1 vs 2), (2 vs 3), etc.
  for (i in seq_len(length(blocks) - 1)) {
    b1 <- blocks[i]
    b2 <- blocks[i + 1]
    
    # Get all stimuli per block
    stim_b1 <- sort(unique(df$stimulus[df$Block == b1]))
    stim_b2 <- sort(unique(df$stimulus[df$Block == b2]))
    
    # Determine how many stimuli to use for comparison
    n1 <- min(length(stim_b1), n_stim)
    n2 <- min(length(stim_b2), n_stim)
    
    # Last N stimuli of block 1 and first N of block 2
    last_idx_b1  <- tail(stim_b1, n1)
    first_idx_b2 <- head(stim_b2, n2)
    
    # --- Compute EMMs for each subset ----------------------------------------
    emm_b1 <- emmeans(
      model, ~ Block | Genotype,
      at = list(stimulus_log = log(last_idx_b1), Block = as.character(b1)),
      cov.reduce = mean, type = "response"
    )
    
    emm_b2 <- emmeans(
      model, ~ Block | Genotype,
      at = list(stimulus_log = log(first_idx_b2), Block = as.character(b2)),
      cov.reduce = mean, type = "response"
    )
    
    # Combine and contrast (recovery-style comparison)
    combined_emms <- rbind(emm_b1, emm_b2)
    block_contrast <- contrast(
      combined_emms,
      method = "revpairwise",
      by = "Genotype",
      adjust = "Tukey"
    )
    
    # Store results
    name_base <- paste0("Block", b1, "_vs_Block", b2)
    emm_results[[paste0(name_base, "_combined")]]  <- combined_emms
    emm_results[[paste0(name_base, "_contrast")]] <- block_contrast
    
    # --- Optional: single-stimulus comparison ---------------------------------
    if (compare_single) {
      emm_b1_single <- emmeans(
        model, ~ Block | Genotype,
        at = list(stimulus_log = log(max(stim_b1)), Block = as.character(b1)),
        type = "response"
      )
      emm_b2_single <- emmeans(
        model, ~ Block | Genotype,
        at = list(stimulus_log = log(min(stim_b2)), Block = as.character(b2)),
        type = "response"
      )
      
      combined_single <- rbind(emm_b1_single, emm_b2_single)
      block_single_contrast <- contrast(
        combined_single,
        method = "revpairwise",
        by = "Genotype",
        adjust = "Tukey"
      )
      
      emm_results[[paste0(name_base, "_single_combined")]]  <- combined_single
      emm_results[[paste0(name_base, "_single_contrast")]] <- block_single_contrast
    }
  }
  
  return(emm_results)
}

get_emm_consecutive_blocks_first <- function(model, df, n_stim = 9, compare_single = TRUE) {
  # Ensure required package
  require(emmeans)
  
  # Detect all block levels and sort them
  blocks <- sort(unique(df$Block))
  
  # Prepare result container
  emm_results <- list()
  
  # Iterate over consecutive block pairs: (1 vs 2), (2 vs 3), etc.
  for (i in seq_len(length(blocks) - 1)) {
    b1 <- blocks[i]
    b2 <- blocks[i + 1]
    
    # Get all stimuli per block
    stim_b1 <- sort(unique(df$stimulus[df$Block == b1]))
    stim_b2 <- sort(unique(df$stimulus[df$Block == b2]))
    
    # Determine how many stimuli to use for comparison
    n1 <- min(length(stim_b1), n_stim)
    n2 <- min(length(stim_b2), n_stim)
    
    # First N stimuli of block 1 and first N of block 2
    first_idx_b1  <- head(stim_b1, n1)
    first_idx_b2 <- head(stim_b2, n2)
    
    # --- Compute EMMs for each subset ----------------------------------------
    emm_b1 <- emmeans(
      model, ~ Block | Genotype,
      at = list(stimulus_log = log(first_idx_b1), Block = as.character(b1)),
      cov.reduce = mean, type = "response"
    )
    
    emm_b2 <- emmeans(
      model, ~ Block | Genotype,
      at = list(stimulus_log = log(first_idx_b2), Block = as.character(b2)),
      cov.reduce = mean, type = "response"
    )
    
    # Combine and contrast (recovery-style comparison)
    combined_emms <- rbind(emm_b1, emm_b2)
    block_contrast <- contrast(
      combined_emms,
      method = "revpairwise",
      by = "Genotype",
      adjust = "Tukey"
    )
    
    # Store results
    name_base <- paste0("Block", b1, "_vs_Block", b2)
    emm_results[[paste0(name_base, "_combined")]]  <- combined_emms
    emm_results[[paste0(name_base, "_contrast")]] <- block_contrast
    
    # --- Optional: single-stimulus comparison ---------------------------------
    if (compare_single) {
      emm_b1_single <- emmeans(
        model, ~ Block | Genotype,
        at = list(stimulus_log = log(min(stim_b1)), Block = as.character(b1)),
        type = "response"
      )
      emm_b2_single <- emmeans(
        model, ~ Block | Genotype,
        at = list(stimulus_log = log(min(stim_b2)), Block = as.character(b2)),
        type = "response"
      )
      
      combined_single <- rbind(emm_b1_single, emm_b2_single)
      block_single_contrast <- contrast(
        combined_single,
        method = "revpairwise",
        by = "Genotype",
        adjust = "Tukey"
      )
      
      emm_results[[paste0(name_base, "_single_combined")]]  <- combined_single
      emm_results[[paste0(name_base, "_single_contrast")]] <- block_single_contrast
    }
  }
  
  return(emm_results)
}


write_section <- function(title, desc = NULL) {
  cat("\n", strrep("-", 80), "\n", sep = "")
  cat("###", title, "###\n")
  if (!is.null(desc)) cat(desc, "\n")
  cat(strrep("-", 80), "\n\n", sep = "")
}


write_emm_report <- function(
    model,
    emm_blocks,
    emm_slopes,
    emm_slopes_pairs,
    emm_between,
    emm_between_first,
    emm_between_blocks_slopes,
    emm_between_blocks_slopes_pairs,
    outfile,
    model_name = "GLMM",
    note_about_interactions = TRUE
) {
  # Ensure output directory exists
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  
  # Helper function for nice section headers
  write_section <- function(title, desc = NULL) {
    cat("\n", strrep("-", 80), "\n", sep = "")
    cat("###", title, "###\n")
    if (!is.null(desc)) cat(desc, "\n")
    cat(strrep("-", 80), "\n\n", sep = "")
  }
  
  # Comput Model Performance
  m_perf <- model_performance(model)
  
  # Start writing output
  sink(outfile)
  
  # --- Header -----------------------------------------------------------------
  cat("###", model_name, "Estimated Marginal Means & Contrasts ###\n\n")
  cat("Model summary:\n")
  print(summary(model))
  cat("\n\n")
  
  if (note_about_interactions) {
    cat("NOTE: Results may be misleading due to involvement in interactions.\n")
    cat("      This message is expected when factors participate in interaction terms.\n")
    cat("      EMMs shown here are conditional on those variables (cov.reduce = mean).\n\n")
  }
  
  # --- Model Performance ------------------------------------------------------
  write_section("Model Performance")
  print(m_perf)
  cat("\n\n")
  
  # --- 1. Block-level EMMs ----------------------------------------------------
  write_section("Block-wise EMMs and Contrasts")
  
  info <- "
  Interpretation for estimate (odds ratio):\n
  A ratio of 1.00 → no difference between genotypes.\n
  A ratio < 1.00 → the first genotype (left of “/”) has a lower predicted response than the second.\n
  A ratio > 1.00 → the first genotype has a higher predicted response.\n
  Habituation Rate (Slope): \n
  Interpretation: Negative slope -> faster habituation.\n\n
  "
  cat(info)
  cat("\n")
  
  for (nm in names(emm_blocks)) {
    cat("##", nm, "\n\n")
    cat("Estimated Marginal Means:\n")
    print(summary(emm_blocks[[nm]]$emm), row.names = FALSE)
    
    cat("\nPairwise Contrasts:\n")
    print(summary(emm_blocks[[nm]]$pairs), row.names = FALSE)
    cat("\n\n")
  }
  
  # --- 2. Habituation Slopes (within blocks) ---------------------------------
  write_section("Habituation Slopes (within blocks)")
  cat("EMM Trends (slope per Block):\n\n")
  print(summary(emm_slopes), row.names = FALSE)
  cat("\nPairwise Contrasts between Genotypes:\n\n")
  print(summary(emm_slopes_pairs), row.names = FALSE)
  
  # --- 3. Between-block comparisons ------------------------------------------
  write_section("Between-block Comparisons")
  cat("Last n vs first n\n")
  for (nm in names(emm_between)) {
    cat("##", nm, "\n\n")
    if (grepl("_contrast$", nm)) {
      cat("Contrast results (last vs first):\n")
      print(summary(emm_between[[nm]]), row.names = FALSE)
    } else if (grepl("_combined$", nm)) {
      cat("Combined EMMs (last vs first):\n")
      print(summary(emm_between[[nm]]), row.names = FALSE)
    }
    cat("\n\n")
  }
  
  cat("First n vs first n\n")
  for (nm in names(emm_between_first)) {
    cat("##", nm, "\n\n")
    if (grepl("_contrast$", nm)) {
      cat("Contrast results (first vs first):\n")
      print(summary(emm_between_first[[nm]]), row.names = FALSE)
    } else if (grepl("_combined$", nm)) {
      cat("Combined EMMs (first vs first):\n")
      print(summary(emm_between_first[[nm]]), row.names = FALSE)
    }
    cat("\n\n")
  }
  
  # --- 4. Between-block slopes -----------------------------------------------
  write_section("Between-block Habituation Slopes")
  cat("EMM Trends (slope between blocks):\n\n")
  print(summary(emm_between_blocks_slopes), row.names = FALSE)
  cat("\nPairwise Contrasts:\n\n")
  print(summary(emm_between_blocks_slopes_pairs), row.names = FALSE)
  
  # --- Wrap up ---------------------------------------------------------------
  sink()
  cat("All EMM results written to:", outfile, "\n")
}
