# UTILITY FUNCTIONS

compute_emm_response_prob <- function(m_prob, df, base_dir){
  # --- Block-level EMMs ---------------------------------------------------------
  emm_prob <- get_emm_blocks(m_prob, df)
  
  # --- Habituation slopes (within blocks) --------------------------------------
  emm_prob_slopes <- emtrends(m_prob, ~ drug_condition | Block, var = "stimulus_log")
  emm_prob_slopes_pairs <- pairs(emm_prob_slopes)
  
  # --- Between-block comparisons (1 vs 2, 2 vs 3, etc.) -------------------------
  emm_prob_between <- get_emm_consecutive_blocks(m_prob, df, n_stim = 20)
  
  # --- Between-block slopes -----------------------------------------------------
  emm_prob_between_blocks_slopes <- emtrends(m_prob, ~ Block | drug_condition, var = "stimulus_log")
  emm_prob_between_blocks_slopes_pairs <- pairs(emm_prob_between_blocks_slopes)
  
  write_emm_report(
    model = m_prob,
    emm_blocks = emm_prob,
    emm_slopes = emm_prob_slopes,
    emm_slopes_pairs = emm_prob_slopes_pairs,
    emm_between = emm_prob_between,
    emm_between_blocks_slopes = emm_prob_between_blocks_slopes,
    emm_between_blocks_slopes_pairs = emm_prob_between_blocks_slopes_pairs,
    outfile = file.path(base_dir, "response_prob", "glmm_response_prob_comparisons.txt"),
    model_name = "GLMM Response Probability Model"
  )
}

compute_emm_distance <- function(m, df, base_dir){
  # --- Block-level EMMs ---------------------------------------------------------
  emm <- get_emm_blocks(m, df)
  
  # --- Habituation slopes (within blocks) --------------------------------------
  emm_slopes <- emtrends(m, ~ drug_condition | Block, var = "stimulus_log")
  emm_slopes_pairs <- pairs(emm_slopes)
  
  # --- Between-block comparisons (1 vs 2, 2 vs 3, etc.) -------------------------
  emm_between <- get_emm_consecutive_blocks(m, df, n_stim = 20)
  
  # --- Between-block slopes -----------------------------------------------------
  emm_between_blocks_slopes <- emtrends(m, ~ Block | drug_condition, var = "stimulus_log")
  emm_between_blocks_slopes_pairs <- pairs(emm_between_blocks_slopes)
  
  write_emm_report(
    model = m,
    emm_blocks = emm,
    emm_slopes = emm_slopes,
    emm_slopes_pairs = emm_slopes_pairs,
    emm_between = emm_between,
    emm_between_blocks_slopes = emm_between_blocks_slopes,
    emm_between_blocks_slopes_pairs = emm_between_blocks_slopes_pairs,
    outfile = file.path(base_dir, "distance", "glmm_distance_comparisons.txt"),
    model_name = "GLMM Distance Model"
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

get_effects_logit <- function(emm, category = "treatment", p_th = 0.05) {
  emm_tbl <- as.data.frame(emm)
  contr <- pairs(emm, by = category)
  contr_tbl <- as.data.frame(contr)
  sig_contr <- subset(contr_tbl, p.value < p_th)
  
  # Split contrasts (e.g., "control - dk")
  sig_contr2 <- sig_contr %>%
    tidyr::separate(contrast, into = c("drug1", "drug2"), sep = " - ")
  
  # Join with emmeans table
  left <- sig_contr2 %>%
    dplyr::left_join(emm_tbl,
                     by = c("drug1" = "drug_condition", category))
  
  right <- left %>%
    dplyr::left_join(emm_tbl,
                     by = c("drug2" = "drug_condition", category),
                     suffix = c("_1", "_2"))
  
  # Compute probability effect difference
  results <- right %>%
    dplyr::mutate(
      prob_1 = plogis(emmean_1),
      prob_2 = plogis(emmean_2),
      prob_diff = prob_1 - prob_2
    ) %>%
    dplyr::select(
      treatment,
      drug1,
      drug2,
      prob_1,
      prob_2,
      prob_diff,
      p.value   # ← include p-value
    )
  
  return(results)
}

get_effects_log <- function(emm, category = "treatment", p_th = 0.05) {
  emm_tbl <- as.data.frame(emm)
  contr <- pairs(emm, by = category)
  contr_tbl <- as.data.frame(contr)
  sig_contr <- subset(contr_tbl, p.value < p_th)
  
  if (nrow(sig_contr) == 0) return(NULL)
  
  # Split contrasts (e.g., "control - dk")
  sig_contr2 <- sig_contr %>%
    tidyr::separate(contrast, into = c("drug1", "drug2"), sep = " - ")
  
  # Join with emmeans table
  left <- sig_contr2 %>%
    dplyr::left_join(
      emm_tbl,
      by = c("drug1" = "drug_condition", category)
    )
  
  right <- left %>%
    dplyr::left_join(
      emm_tbl,
      by = c("drug2" = "drug_condition", category),
      suffix = c("_1", "_2")
    )
  
  # Compute effects on response scale
  results <- right %>%
    dplyr::mutate(
      val_1 = exp(emmean_1),
      val_2 = exp(emmean_2),
      diff = val_2 - val_1,
      ratio = val_2 / val_1,
      pct_diff = (ratio - 1) * 100
    ) %>%
    dplyr::select(
      !!category, drug1, drug2,
      val_1, val_2,
      diff,
      ratio,
      pct_diff,
      p.value
    )
  
  return(results)
}

get_slope_effects_log <- function(slopes, category = "treatment", p_th = 0.05) {
  
  slope_tbl <- as.data.frame(slopes)
  contr <- pairs(slopes, by = category)
  contr_tbl <- as.data.frame(contr)
  sig_contr <- subset(contr_tbl, p.value < p_th)
  
  if (nrow(sig_contr) == 0) return(NULL)
  
  # Split contrasts (e.g., "control - dk")
  sig_contr2 <- sig_contr %>%
    tidyr::separate(contrast, into = c("drug1", "drug2"), sep = " - ")
  
  # Join with slope table
  left <- sig_contr2 %>%
    dplyr::left_join(
      slope_tbl,
      by = c("drug1" = "drug_condition", category)
    )
  
  right <- left %>%
    dplyr::left_join(
      slope_tbl,
      by = c("drug2" = "drug_condition", category),
      suffix = c("_1", "_2")
    )
  
  # Compute slope effects on the response scale
  results <- right %>%
    dplyr::mutate(
      slope_1 = exp(stimulus_log.trend_1),
      slope_2 = exp(stimulus_log.trend_2),
      diff = slope_2 - slope_1,
      ratio   = slope_2 / slope_1,      # relative slope difference
      pct_diff = (ratio - 1) * 100
    ) %>%
    dplyr::select(
      !!category,
      drug1, drug2,
      slope_1, slope_2,
      diff,
      ratio,
      pct_diff,
      p.value
    )
  
  return(results)
}

pl <- function(val1, val2){
  require(glue)
  v1_plogis <- plogis(val1)
  v2_plogis <- plogis(val2)
  print("Back-transfomred values from logit scale and the difference:")
  print(glue("val1: {v1_plogis}"))
  print(glue("val2: {v2_plogis}"))
  print(glue("difference: {v1_plogis - v2_plogis}"))
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
  plotResiduals(res, df$drug_condition)
  plotResiduals(res, df$Block)
  
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
  print(check_collinearity(update(model, . ~ drug_condition + Block + stimulus_log)))
}

load_48_well_plates <- function(file_dir){
  # Load data from csv file
  df <- read_csv(file_dir)
  df$response <- as.numeric(df$response)
  df <- df %>%
    group_by(ISI_block, fish_id, plate_id) %>%
    mutate(tap_in_block = tap_num - min(tap_num) + 1) %>%
    ungroup()
  
  df <- df %>%
    rename(
      stimulus = tap_in_block,
      Block = ISI_block
    )
  
  df$fish_id <- as.factor(df$fish_id)
  df$plate_id <- as.factor(df$plate_id)
  df$drug_condition <- as.factor(df$drug_condition)
  df$Block <- as.factor(df$Block)
  
  # Log transform tap number (for exponential decay)
  df$stimulus_log<- log(df$stimulus)
  
  # Set References
  df$drug_condition <- factor(df$drug_condition, levels = c("control", "acute", "dk"))
  # df$Block <- factor(df$Block, levels = c("ISI_90s", "ISI_5s_block"))
    
  # Make plate ID unqiue
  df <- df %>%
    mutate(plate_uid = paste(treatment, plate_id, sep = "_"))
  
  # plate_ids renumbered independently within each treatment:
  # df <- df %>%
  #   group_by(treatment, plate_id) %>%
  #   mutate(plate_uid = cur_group_id()) %>%
  #   ungroup()
  
  df_sub <- subset(df, response > 0)
  return(list(df = df, df_sub = df_sub))
}

load_data <- function(file_dir) {
  
  df <- read_csv(file_dir)
  df$prob <- as.numeric(df$prob)       # 0/1
  df$fish_id <- as.factor(df$fish_id)
  df$plate_id <- as.factor(df$plate_id)
  df$drug_condition <- as.factor(df$drug_condition)
  df$ISI_block <- as.factor(df$ISI_block)
  
  # Define tap number within block
  df <- df %>%
    group_by(ISI_block, fish_id, plate_id) %>%
    mutate(tap_in_block = tap_num - min(tap_num) + 1) %>%
    ungroup()
  
  df <- df %>%
    rename(
      stimulus = tap_in_block,
      Block = ISI_block
    )
  
  # Log transform tap number (for exponential decay)
  df$stimulus_log<- log(df$stimulus)
  
  # Set References
  df$drug_condition <- factor(df$drug_condition, levels = c("control", "acute", "dk"))
  df$Block <- factor(df$Block, levels = c("ISI_90s", "ISI_5s_block1", "ISI_5s_block2"))
  
  df_sub <- subset(df, prob > 0)
  
  return(list(df = df, df_sub = df_sub))
}

plot_habituation <- function(df_final, model, label, Ymin, Ymax, response_var, 
                             transform = c("plogis", "exp", "none")) {
  # --- 0. Setup ----------------------------------------------
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggpubr)
  
  transform <- match.arg(transform)
  
  # --- 1. Dynamic prediction grid ----------------------------
  blocks <- sort(unique(df_final$Block))
  
  stim_ranges <- lapply(blocks, function(b) {
    range(unique(df_final$stimulus[df_final$Block == b]))
  })
  names(stim_ranges) <- blocks
  
  new_data <- bind_rows(lapply(blocks, function(b) {
    expand.grid(Block = as.character(b),
                stimulus = stim_ranges[[as.character(b)]][1]:
                  stim_ranges[[as.character(b)]][2])
  })) %>%
    tidyr::crossing(drug_condition = unique(df_final$drug_condition)) %>%
    mutate(stimulus_log = log(stimulus)) %>%
    filter(stimulus > 0)
  
  # Align factor levels
  new_data$drug_condition <- factor(new_data$drug_condition, levels = levels(df_final$drug_condition))
  new_data$Block <- factor(new_data$Block, levels = levels(df_final$Block))
  
  # --- 2. Model predictions ----------------------------------
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
  
  # --- 3. Compute raw data means per stimulus -----------------
  raw_means <- df_final %>%
    filter(stimulus > 0) %>%
    group_by(Block, drug_condition, stimulus) %>%
    summarise(
      raw_mean = mean(.data[[response_var]], na.rm = TRUE),
      .groups = "drop"
    )
  
  raw_means$drug_condition <- factor(raw_means$drug_condition, 
                                     levels = levels(df_final$drug_condition))
  raw_means$Block <- factor(raw_means$Block, 
                            levels = levels(df_final$Block))
  
  # --- 4. Facet labels ----------------------------------------
  labels <- setNames(
    paste0("Block ", blocks, ": ", sapply(stim_ranges, `[`, 2), " taps"),
    blocks
  )
  
  # --- 5. Plot -------------------------------------------------
  g <- ggplot(new_data, aes(x = stimulus, color = drug_condition, fill = drug_condition)) +
    facet_wrap(~Block, ncol = length(blocks), scales = "free_x",
               labeller = as_labeller(labels)) +
    
    # Raw mean data points
    geom_point(
      data = raw_means,
      aes(x = stimulus, y = raw_mean, color = drug_condition),
      size = 2, alpha = 0.6,
      position = position_dodge(width = 0.4)
    ) +
    
    # Model predictions
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
                alpha = 0.1, color = NA) +
    geom_line(aes(y = fit), linewidth = 1.2) +
    coord_cartesian(ylim = c(Ymin, Ymax)) +
    labs(
      x = "Stimulus number (within block)",
      y = label,
      color = "drug_condition",
      fill  = "drug_condition"
    ) +
    theme_pubr(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2)
    )
  
  return(g)
}


plot_habituation_old <- function(df_final, model, label, Ymin, Ymax, transform = c("plogis", "exp", "none")) {
  # --- 0. Setup ----------------------------------------------
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggpubr)
  
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
    tidyr::crossing(drug_condition = unique(df_final$drug_condition)) %>%
    mutate(stimulus_log = log(stimulus))
  
  # Ensure no log(0)
  new_data <- new_data %>%
    mutate(
      stimulus = ifelse(stimulus <= 0, NA, stimulus),
      stimulus_log = log(stimulus)
    ) %>%
    filter(!is.na(stimulus_log))
  
  # Align factor levels
  new_data$drug_condition <- factor(new_data$drug_condition, levels = levels(df_final$drug_condition))
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
    paste0("Block ", blocks, ": ", sapply(stim_ranges, `[`, 2), " taps"),
    blocks
  )
  
  # --- 4. Plot -----------------------------------------------
  g <- ggplot(new_data, aes(x = stimulus, color = drug_condition, fill = drug_condition)) +
    facet_wrap(~Block, ncol = length(blocks), scales = "free_x",
               labeller = as_labeller(labels)) +
    geom_line(aes(y = fit), linewidth = 1.2) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
                alpha = 0.1, color = NA) +
    expand_limits(y = c(0, 1)) +
    coord_cartesian(ylim = c(Ymin, Ymax)) +
    labs(
      x = "Stimulus number (within block)",
      y = label,
      color = "drug_condition",
      fill  = "drug_condition"
    ) +
    theme_pubr(base_size = 14)+
    theme(
      panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2)
    )
  
  return(g)
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
        ~ drug_condition,
        at = list(
          stimulus_log = log(idx),
          Block = as.character(b)
        ),
        cov.reduce = mean,
        type = "response"
      )
      
      pairwise <- pairs(emm)  # pairwise contrasts between drug_conditions
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
    # either use n_stim or less
    n1 <- min(length(stim_b1), n_stim)
    n2 <- min(length(stim_b2), n_stim)
    
    # Last N stimuli of block 1 and first N of block 2
    last_idx_b1  <- tail(stim_b1, n1)
    first_idx_b2 <- head(stim_b2, n2)
    
    # --- Compute EMMs for each subset ----------------------------------------
    emm_b1 <- emmeans(
      model, ~ Block | drug_condition,
      at = list(stimulus_log = log(last_idx_b1), Block = as.character(b1)),
      cov.reduce = mean, type = "response"
    )
    
    emm_b2 <- emmeans(
      model, ~ Block | drug_condition,
      at = list(stimulus_log = log(first_idx_b2), Block = as.character(b2)),
      cov.reduce = mean, type = "response"
    )
    
    # Combine and contrast (recovery-style comparison)
    combined_emms <- rbind(emm_b1, emm_b2)
    block_contrast <- contrast(
      combined_emms,
      method = "revpairwise",
      by = "drug_condition",
      adjust = "Tukey"
    )
    
    # Store results
    name_base <- paste0("Block", b1, "_vs_Block", b2)
    emm_results[[paste0(name_base, "_combined")]]  <- combined_emms
    emm_results[[paste0(name_base, "_contrast")]] <- block_contrast
    
    # --- Optional: single-stimulus comparison ---------------------------------
    # This compares the first stimulus of Block x and the last of Block x+1
    if (compare_single) {
      emm_b1_single <- emmeans(
        model, ~ Block | drug_condition,
        at = list(stimulus_log = log(max(stim_b1)), Block = as.character(b1)),
        type = "response"
      )
      emm_b2_single <- emmeans(
        model, ~ Block | drug_condition,
        at = list(stimulus_log = log(min(stim_b2)), Block = as.character(b2)),
        type = "response"
      )
      
      combined_single <- rbind(emm_b1_single, emm_b2_single)
      block_single_contrast <- contrast(
        combined_single,
        method = "revpairwise",
        by = "drug_condition",
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
    # either use n_stim or less
    n1 <- min(length(stim_b1), n_stim)
    n2 <- min(length(stim_b2), n_stim)
    
    # First N stimuli of block 1 and first N of block 2
    first_idx_b1  <- head(stim_b1, n1)
    first_idx_b2 <- head(stim_b2, n2)
    
    # --- Compute EMMs for each subset ----------------------------------------
    emm_b1 <- emmeans(
      model, ~ Block | drug_condition,
      at = list(stimulus_log = log(first_idx_b1), Block = as.character(b1)),
      cov.reduce = mean, type = "response"
    )
    
    emm_b2 <- emmeans(
      model, ~ Block | drug_condition,
      at = list(stimulus_log = log(first_idx_b2), Block = as.character(b2)),
      cov.reduce = mean, type = "response"
    )
    
    # Combine and contrast (recovery-style comparison)
    combined_emms <- rbind(emm_b1, emm_b2)
    block_contrast <- contrast(
      combined_emms,
      method = "revpairwise",
      by = "drug_condition",
      adjust = "Tukey"
    )
    
    # Store results
    name_base <- paste0("Block", b1, "_vs_Block", b2)
    emm_results[[paste0(name_base, "_combined")]]  <- combined_emms
    emm_results[[paste0(name_base, "_contrast")]] <- block_contrast
    
    # --- Optional: single-stimulus comparison ---------------------------------
    # This compares the first stimulus of Block x and the first of Block x+1
    if (compare_single) {
      emm_b1_single <- emmeans(
        model, ~ Block | drug_condition,
        at = list(stimulus_log = log(min(stim_b1)), Block = as.character(b1)),
        type = "response"
      )
      emm_b2_single <- emmeans(
        model, ~ Block | drug_condition,
        at = list(stimulus_log = log(min(stim_b2)), Block = as.character(b2)),
        type = "response"
      )
      
      combined_single <- rbind(emm_b1_single, emm_b2_single)
      block_single_contrast <- contrast(
        combined_single,
        method = "revpairwise",
        by = "drug_condition",
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
  
  for (nm in names(emm_between)) {
    cat("##", nm, "\n\n")
    if (grepl("_contrast$", nm)) {
      cat("Contrast results:\n")
      print(summary(emm_between[[nm]]), row.names = FALSE)
    } else if (grepl("_combined$", nm)) {
      cat("Combined EMMs:\n")
      print(summary(emm_between[[nm]]), row.names = FALSE)
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
