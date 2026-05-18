# UTILITY FUNCTIONS

extract_random_effects_and_icc <- function(glmm) {
  # 1. Check if the model is a glmmTMB object
  if (!inherits(glmm, "glmmTMB")) {
    stop("Input must be a fitted glmmTMB model object.")
  }
  
  # 2. Extract Random Effects Variances (Conditional model)
  vc <- glmmTMB::VarCorr(glmm)$cond
  
  # Extract specific variance values
  var_fish <- vc$fish[1, 1]
  var_roi  <- vc$roi[1, 1]
  
  # 3. Extract the Residual Variance Proxy (Dispersion Estimate for Gamma(log))
  # This is needed for the ICC approximation
  var_dispersion <- glmmTMB::sigma(glmm) 
  
  # 4. Calculate Total Variances
  var_total_random <- var_fish + var_roi
  var_total_approx <- var_total_random + var_dispersion # Total variance on the log scale
  
  # 5. Calculate Proportions and ICCs
  # Proportion of Conditional Random Variance
  prop_fish_random <- var_fish / var_total_random
  prop_roi_random  <- var_roi / var_total_random
  
  # ICC (Intra-Class Correlation) - Proportion of Total Variance
  icc_fish <- var_fish / var_total_approx
  icc_roi  <- var_roi / var_total_approx
  
  # ICC for Residuals (The remaining unexplained variance)
  icc_residual <- var_dispersion / var_total_approx
  
  # 6. Display the results in a formatted way
  cat("\n--- Random Effects & ICC Analysis (Gamma GLMM, Log Link) ---\n")
  cat(sprintf("Total Variance (Approx, Log-Scale): %.5f\n", var_total_approx))
  cat(sprintf("Dispersion Estimate (Residual Proxy): %.5f\n\n", var_dispersion))
  
  # Display Proportions of Conditional Random Variance (Your original output)
  cat("A. Proportions of Conditional Random Variance:\n")
  cat(sprintf("   > Variance (fish): %.5f (%.2f%% of total RANDOM variance)\n", var_fish, prop_fish_random * 100))
  cat(sprintf("   > Variance (roi):  %.5f (%.2f%% of total RANDOM variance)\n", var_roi, prop_roi_random * 100))
  cat(sprintf("Total Conditional Random Variance: %.5f\n\n", var_total_random))
  
  # Display ICCs (New output)
  cat("\nB. Intra-Class Correlation (ICC):\n")
  cat(sprintf("   > ICC (fish): %.3f (%.1f%% of total log-scale variance)\n", icc_fish, icc_fish * 100))
  cat(sprintf("   > ICC (roi):  %.3f (%.1f%% of total log-scale variance)\n", icc_roi, icc_roi * 100))
  cat(sprintf("   > Unexplained (Residual): %.3f (%.1f%% of total log-scale variance)\n", icc_residual, icc_residual * 100))
  
  # Optional: Return a data frame with all the values
  return(
    invisible(
      data.frame(
        Group = c("fish", "roi", "residual"),
        Variance = c(var_fish, var_roi, var_dispersion),
        ICC = c(icc_fish, icc_roi, icc_residual)
      )
    )
  )
}


simple_validate <- function(model){
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
  print(check_collinearity(model))
  print(check_collinearity(update(model, . ~ condition)))
  
  # Diagnostic plots
  plot(res)
  
}


load_data <- function(file_dir) {
  
  # Load data from csv file
  df <- read_csv(file_dir)
  
  # Remove baseline condition (pseudo)
  # df <- subset(df_org, condition != "base_line")
  
  df$fish <- as.factor(df$fish)
  df$roi <- as.factor(df$roi)
  df$stim_name <- as.factor(df$stim_name)
  
  df <- df %>%
    mutate(fish = recode(fish,
                         "1" = "1",
                         "3" = "2",
                         "4" = "3",
                         "5" = "4",
                         "6" = "5",
                         "9" = "6"))
  df <- df %>%
    mutate(fish = factor(fish, levels = c("1", "2", "3", "4", "5", "6")))
  
  # Categorical GLMM
  df$condition <- factor(df$condition, 
                         levels = c("stim_only", "motor_only", "stim_motor"))
  
  df$score <- df$score + 1e-6  # avoid zeros
  df$score_norm99 <- ifelse(df$score_norm99 < 0, 0, df$score_norm99)
  df$score_norm99 <- df$score_norm99 + 1e-6
  df$log1p_score <- log1p(df$score)
  df$condition <- factor(df$condition, 
                         levels = c("stim_only", "motor_only", "stim_motor"))
  
  
  return(df)
}

plot_synergy_dual_scale <- function(model, data, label, cond_var, order) {
  require(emmeans)
  require(ggplot2)
  require(ggpubr)
  require(dplyr)
  
  palette <- list(
    "motor_only" = c("#D95F02"),
    "stim_only" = c("#1B9E77"),
    "stim+motor" = c("#7570B3")
  )
  
  # --- 1. Estimated marginal means on log scale ---
  emm_log  <- emmeans(model, as.formula(paste("~", cond_var))) 
  emm_resp <- emmeans(model, as.formula(paste("~", cond_var)), type = "response")
  
  emm_log_df  <- as.data.frame(emm_log)
  emm_resp_df <- as.data.frame(emm_resp)
  
  # --- 2. Condition order ---
  levels_order <- order
  
  # if (!is.null(order)) {
  #   levels_order <- order
  # } else {
  #   cond_levels <- levels(as.factor(data[[cond_var]]))
  #   if (all(c("motor_only", "stim_only", "stim+motor") %in% cond_levels)) {
  #     levels_order <- c("motor_only", "stim_only", "stim+motor")
  #   } else {
  #     levels_order <- cond_levels
  #   }
  # }
  # 
  emm_log_df[[cond_var]]  <- factor(emm_log_df[[cond_var]],  levels = levels_order)
  emm_resp_df[[cond_var]] <- factor(emm_resp_df[[cond_var]], levels = levels_order)
  
  # --- 3. Compute additive line parameters (first two → extrapolate to last) ---
  pts <- emm_log_df %>%
    filter(.data[[cond_var]] %in% levels_order[1:2]) %>%
    arrange(match(.data[[cond_var]], levels_order[1:2]))
  
  x_vals <- as.numeric(pts[[cond_var]])
  slope <- diff(pts$emmean) / diff(x_vals)
  intercept <- pts$emmean[1] - slope * x_vals[1]
  x_target <- which(levels_order == tail(levels_order, 1))
  
  # --- 4. Compute additive predicted mean and CI at extrapolated x ---
  V <- vcov(emm_log)
  L <- rep(0, length(levels_order))
  L[1:2] <- c(-1, 2)
  
  add_SE_log <- sqrt(as.numeric(t(L) %*% V %*% L))
  add_pred_log <- intercept + slope * x_target
  add_CI_log <- add_pred_log + c(-1.96, 1.96) * add_SE_log
  
  add_pred_resp <- exp(add_pred_log)
  add_CI_resp <- exp(add_CI_log)
  
  # --- 5. Data frames for plotting ---
  x_ext <- seq(0.5, length(levels_order) + 0.5, length.out = 200)
  line_df_log  <- data.frame(x = x_ext, y = intercept + slope * x_ext)
  line_df_resp <- data.frame(x = x_ext, y = exp(intercept + slope * x_ext))
  
  dot_log  <- data.frame(x = x_target, y = add_pred_log,
                         ymin = add_CI_log[1], ymax = add_CI_log[2])
  dot_resp <- data.frame(x = x_target, y = add_pred_resp,
                         ymin = add_CI_resp[1], ymax = add_CI_resp[2])
  
  # --- 6. Plot log scale ---
  p_log <- ggplot(emm_log_df, aes(x = as.numeric(.data[[cond_var]]),
                                  y = emmean, color = .data[[cond_var]])) +
    geom_line(data = line_df_log, aes(x = x, y = y),
              inherit.aes = FALSE, color = "gray40", linetype = "dotted", linewidth = 1) +
    geom_errorbar(data = dot_log, aes(x = x, ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE, width = 0.1, color = "gray20", linewidth = 0.8) +
    geom_point(data = dot_log, aes(x = x, y = y),
               shape = 21, size = 4.5, fill = "gray20", color = "white", stroke = 1.3,
               inherit.aes = FALSE) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                  width = 0.1, linewidth = 1) +
    scale_x_continuous(breaks = 1:length(levels_order), labels = levels_order) +
    scale_color_manual(values = palette) +
    labs(title = "A. Log scale (model space)",
         subtitle = "Gray dot = additive expectation (with 95% CI)",
         x = cond_var, y = paste0("Predicted log(", label, ")")) +
    theme_ssi(base_size = 10) +  # <- apply same theme
    theme(legend.position = "none")
  
  # --- 7. Plot response scale ---
  p_resp <- ggplot(emm_resp_df, aes(x = as.numeric(.data[[cond_var]]),
                                    y = response, color = .data[[cond_var]])) +
    geom_line(data = line_df_resp, aes(x = x, y = y),
              inherit.aes = FALSE, color = "gray40", linetype = "dotted", linewidth = 1.2) +
    geom_errorbar(data = dot_resp, aes(x = x, ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE, width = 0.1, color = "gray20", linewidth = 0.8) +
    geom_point(data = dot_resp, aes(x = x, y = y),
               shape = 21, size = 4.5, fill = "gray20", color = "white", stroke = 1.3,
               inherit.aes = FALSE) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                  width = 0.1, linewidth = 1) +
    scale_x_continuous(breaks = 1:length(levels_order), labels = levels_order) +
    scale_color_manual(values = palette) +
    labs(title = "B. Response scale (back-transformed)",
         subtitle = "Gray dot = additive expectation (with 95% CI)",
         x = cond_var, y = paste0("Predicted ", label)) +
    theme_ssi(base_size = 10) +  # <- apply same theme
    theme(legend.position = "none")
  
  ggarrange(p_log, p_resp, ncol = 2, align = "hv")
}


plot_synergy_summary <- function(model, cond_var = "condition",
                                 palette = c("#1B9E77", "#D95F02", "#7570B3")) {
  require(emmeans)
  require(ggplot2)
  require(dplyr)
  
  # --- 1. Get pairwise ratios and synergy ratio ---
  emm_resp <- emmeans(model, as.formula(paste("~", cond_var)), type = "response")
  
  pairs_df <- as.data.frame(pairs(emm_resp)) %>%
    mutate(type = "pairwise")
  
  synergy_df <- as.data.frame(
    contrast(emm_resp, list("stim+motor_vs_additive" = c(-0.5, -0.5, 1)))
  ) %>%
    mutate(contrast = "synergy_vs_additive",
           type = "synergy")
  
  # Combine
  df_all <- bind_rows(pairs_df, synergy_df) %>%
    mutate(
      contrast = factor(contrast, levels = rev(unique(contrast))),
      CI_low  = ratio / exp(1.96 * SE / ratio * 0 + 0),  # just placeholder to avoid warnings
      CI_low  = ratio / exp(1.96 * SE),
      CI_high = ratio * exp(1.96 * SE)
    )
  
  # --- 2. Plot ---
  ggplot(df_all, aes(x = contrast, y = ratio, color = type, fill = type)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    geom_point(size = 3.5) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15, linewidth = 1) +
    scale_y_log10() +
    scale_color_manual(values = c(pairwise = palette[1], synergy = "#E41A1C")) +
    scale_fill_manual(values = c(pairwise = palette[1], synergy = "#E41A1C")) +
    coord_flip() +
    labs(
      title = "Synergy summary (response scale)",
      subtitle = "Ratios > 1 indicate higher response in numerator condition",
      x = "Contrast",
      y = "Response ratio (log scale)"
    ) +
    theme_pubr(base_size = 14)
}


rescale_to_01 <- function(x, eps = 1e-6) {
  y <- (x + 1) / 2
  y <- pmin(pmax(y, eps), 1 - eps)
  return(y)
}


extract_model_summary <- function(model) {
  # For lmer models → emmean
  if (inherits(model, "lmerMod")) {
    emm <- emmeans(model, ~1)
    out <- as.data.frame(summary(emm))
    tibble(
      Estimate = out$emmean[1],
      CI_low   = out$lower.CL[1],
      CI_high  = out$upper.CL[1]
    )
    
    # For glmmTMB models (Gamma/log, Beta/logit, etc.) → response scale
  } else if (inherits(model, "glmmTMB")) {
    emm <- emmeans(model, ~1, type = "response")
    out <- as.data.frame(summary(emm))
    tibble(
      Estimate = out$response[1],
      CI_low   = out$asymp.LCL[1],
      CI_high  = out$asymp.UCL[1]
    )
    
  } else {
    stop("Unsupported model type: ", class(model)[1])
  }
}


plot_fish_clouds_dualpanel <- function(
    df,
    models,
    palette = "Dark2",
    spacing = 6,        # larger = more white gap between fish clouds
    jitter_width = 0.03,
    point_size = 1.8
) {
  require(ggplot2)
  require(dplyr)
  require(ggpubr)
  require(tidyr)
  require(emmeans)
  
  # --- Prepare data (long format) ---
  df_long <- df %>%
    select(fish,
           vis_spont_index,
           mixed_spont_index,
           mixed_vis_index,
           synergy_index,
           synergy) %>%
    pivot_longer(
      cols = -fish,
      names_to = "index",
      values_to = "value"
    ) %>%
    mutate(
      index = factor(
        index,
        levels = c(
          "vis_spont_index", "mixed_spont_index", "mixed_vis_index",
          "synergy_index", "synergy"
        ),
        labels = c(
          "Vis vs Spont", "Mixed vs Spont", "Mixed vs Vis",
          "Synergy Index", "Synergy Ratio"
        )
      )
    )
  
  # --- Compute per-fish offsets for cloud separation ---
  df_long <- df_long %>%
    group_by(index) %>%
    mutate(
      fish = factor(fish),
      fish_n = as.numeric(fish),
      x = as.numeric(index) + (fish_n - mean(fish_n)) / spacing
    ) %>%
    ungroup()
  
  # --- Extract model summaries (means + CIs) ---
  results <- bind_rows(
    cbind(index = "Vis vs Spont",   extract_model_summary(models$m_vis_spont)),
    cbind(index = "Mixed vs Spont", extract_model_summary(models$m_mixed_spont)),
    cbind(index = "Mixed vs Vis",   extract_model_summary(models$m_mixed_vis)),
    cbind(index = "Synergy Index",  extract_model_summary(models$m_synergy_index)),
    cbind(index = "Synergy Ratio",  extract_model_summary(models$m_synergy_gamma))
  )
  
  results$index <- factor(results$index, levels = levels(df_long$index))
  
  # --- Split data for dual-panel plotting ---
  df_idx   <- df_long  %>% filter(index != "Synergy Ratio")
  df_syn   <- df_long  %>% filter(index == "Synergy Ratio")
  res_idx  <- results  %>% filter(index != "Synergy Ratio")
  res_syn  <- results  %>% filter(index == "Synergy Ratio")
  
  # --- Panel A: Indices (−1 → 1) ---
  p_idx <- ggplot() +
    geom_jitter(
      data = df_idx,
      aes(x = x, y = value, color = fish),
      size = point_size, alpha = 0.6,
      width = jitter_width, height = 0
    ) +
    geom_errorbar(
      data = res_idx,
      aes(x = as.numeric(index), ymin = CI_low, ymax = CI_high),
      color = "black", linewidth = 1, width = 0
    ) +
    geom_point(
      data = res_idx,
      aes(x = as.numeric(index), y = Estimate),
      size = 3, color = "black", fill = "gray60",
      shape = 21, stroke = 1.1
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    scale_color_brewer(palette = palette) +
    scale_x_continuous(
      breaks = 1:length(levels(df_idx$index)),
      labels = levels(df_idx$index)
    ) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      title = "A. Selectivity Indices",
      # subtitle = "Colored dots = neurons grouped by fish; black = model mean ± 95% CI",
      x = "",
      y = "Index (−1 to 1)",
      color = "Fish"
    )+
    theme_ssi(base_size = 10) +  # <- apply same theme
    theme(legend.position = "none")  # optionally override
  
  # --- Panel B: Synergy Ratio (>0) ---
  p_syn <- ggplot() +
    geom_jitter(
      data = df_syn,
      aes(x = x, y = value, color = fish),
      size = point_size, alpha = 0.6,
      width = jitter_width, height = 0
    ) +
    geom_errorbar(
      data = res_syn,
      aes(x = as.numeric(index), ymin = CI_low, ymax = CI_high),
      color = "black", linewidth = 1, width = 0
    ) +
    geom_point(
      data = res_syn,
      aes(x = as.numeric(index), y = Estimate),
      size = 3, color = "black", fill = "gray60",
      shape = 21, stroke = 1.1
    ) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "gray60") +
    scale_color_brewer(palette = palette) +
    scale_x_continuous(
      breaks = 1,
      labels = "Synergy Ratio"
    ) +
    # coord_cartesian(ylim = c(0, max(df_syn$value, na.rm = TRUE) * 1.2)) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(
      title = "B. Synergy Ratio",
      # subtitle = "Dotted line = additivity (ratio = 1)",
      x = "",
      y = "Synergy Ratio",
      color = "Fish"
    )+
    theme_ssi(base_size = 10) +  # <- apply same theme
    theme(legend.position = "none")  # optionally override
  
  
  # --- Combine both panels side-by-side ---
  # ggarrange(p_idx, p_syn, ncol = 2, widths = c(3, 1), align = "hv")
  ggarrange(
    p_idx + theme(plot.margin = margin(5, 15, 5, 5)),  # add right margin
    p_syn + theme(plot.margin = margin(5, 5, 5, 10)),  # add left margin
    ncol = 2, widths = c(3, 1), align = "hv"
  )
}


# =====================================================================
# test_ssi_significance()
# Runs emmeans tests for all SSI models (lmer + Gamma)
# =====================================================================

test_ssi_significance <- function(models) {
  library(emmeans)
  library(dplyr)
  library(purrr)
  
  # Define which model compares to 0 vs 1
  tests <- tibble(
    name = c("Vis vs Spont", "Mixed vs Spont", "Mixed vs Vis", 
             "Synergy Index", "Synergy Ratio"),
    model = c("m_vis_spont", "m_mixed_spont", "m_mixed_vis", 
              "m_synergy_index", "m_synergy_gamma"),
    null  = c(0, 0, 0, 0, 1)
  )
  
  # Helper for a single model
  run_test <- function(model, null_value) {
    if (inherits(model, "lmerMod")) {
      emm <- emmeans(model, ~1)
      test_out <- test(emm, null = null_value)
      est <- summary(emm)$emmean
      lower <- summary(emm)$lower.CL
      upper <- summary(emm)$upper.CL
      pval <- test_out$p.value
      tibble(Estimate = est, CI_low = lower, CI_high = upper, p = pval)
      
    } else if (inherits(model, "glmmTMB")) {
      emm <- emmeans(model, ~1, type = "response")
      test_out <- test(emm, null = null_value)
      est <- summary(emm)$response
      lower <- summary(emm)$asymp.LCL
      upper <- summary(emm)$asymp.UCL
      pval <- test_out$p.value
      tibble(Estimate = est, CI_low = lower, CI_high = upper, p = pval)
      
    } else {
      tibble(Estimate = NA, CI_low = NA, CI_high = NA, p = NA)
    }
  }
  
  # Loop over all models
  results <- tests %>%
    mutate(
      res = pmap(list(model, null), function(model, null) {
        run_test(models[[model]], null)
      })
    ) %>%
    unnest(res) %>%
    mutate(
      Signif = case_when(
        is.na(p) ~ NA_character_,
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE      ~ "n.s."
      )
    )
  
  results %>%
    select(name, Estimate, CI_low, CI_high, null, p, Signif)
}

theme_ssi <- function(base_size = 10) {
  theme_pubr(base_size = base_size) +
    theme(
      # --- Axis text and titles ---
      axis.title   = element_text(size = base_size + 1, face = "bold"),
      axis.text    = element_text(size = base_size),
      axis.line    = element_line(linewidth = 1.0, color = "black"),
      axis.ticks   = element_line(linewidth = 0.5, color = "black"),
      axis.ticks.length = unit(0.12, "cm"),
      
      # --- Panel and grid ---
      # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      # panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      # panel.grid.minor = element_blank(),
      
      # --- Titles ---
      plot.title   = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle= element_text(size = base_size, hjust = 0),
      
      # --- Legend ---
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text  = element_text(size = base_size - 1),
      legend.key.size = unit(0.5, "lines"),
      legend.spacing.x = unit(0.4, "lines"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box.spacing = unit(0.2, "lines"),
      
      # --- Margins ---
      plot.margin = margin(5, 5, 5, 5)
    )
}


# Parallel Hierarchical Bootstrap 
parallel_hier_boot <- function(data, fish_var, value_var, 
                               B = 10000, workers = 4,
                               null_value = 0) {
  
  require(dplyr)
  require(future.apply)
  
  # 1. Pre-split values by fish (fastest representation)
  groups <- split(data[[value_var]], data[[fish_var]])
  fish_ids <- names(groups)
  n_fish   <- length(groups)
  
  # 2. Start parallel workers
  plan(multisession, workers = workers)
  
  # 3. Parallel bootstrap
  boot_stats <- future_sapply(seq_len(B), function(i) {
    # sample fish with replacement
    fish_sample <- sample(fish_ids, n_fish, replace = TRUE)
    
    # for each fish, resample ROIs
    sampled_values <- unlist(
      lapply(fish_sample, function(f) {
        g <- groups[[f]]
        g[sample.int(length(g), size = length(g), replace = TRUE)]
      }),
      use.names = FALSE
    )
    
    mean(sampled_values, na.rm = TRUE)
  }, future.seed = TRUE)
  
  # 4. Summaries
  tibble(
    mean_boot = mean(boot_stats),
    ci_lower  = quantile(boot_stats, 0.025),
    ci_upper  = quantile(boot_stats, 0.975),
    
    # p-value for null: mean = null_value
    p_twoside = 2 * min(
      mean(boot_stats <= null_value),
      mean(boot_stats >= null_value)
    ),
    
    null_value = null_value,
    distribution = list(boot_stats)
  )
}

make_ssi_boot_table <- function(ssi_boot_results) {
  tibble::tibble(
    name = c(
      "Vis vs Spont",
      "Mixed vs Spont",
      "Mixed vs Vis",
      "Synergy Index",
      "Synergy Ratio"
    ),
    Estimate = c(
      ssi_boot_results$vis_spont$mean_boot,
      ssi_boot_results$mixed_spont$mean_boot,
      ssi_boot_results$mixed_vis$mean_boot,
      ssi_boot_results$synergy_index$mean_boot,
      ssi_boot_results$synergy_ratio$mean_boot
    ),
    CI_low = c(
      ssi_boot_results$vis_spont$ci_lower,
      ssi_boot_results$mixed_spont$ci_lower,
      ssi_boot_results$mixed_vis$ci_lower,
      ssi_boot_results$synergy_index$ci_lower,
      ssi_boot_results$synergy_ratio$ci_lower
    ),
    CI_high = c(
      ssi_boot_results$vis_spont$ci_upper,
      ssi_boot_results$mixed_spont$ci_upper,
      ssi_boot_results$mixed_vis$ci_upper,
      ssi_boot_results$synergy_index$ci_upper,
      ssi_boot_results$synergy_ratio$ci_upper
    ),
    null = c(
      ssi_boot_results$vis_spont$null_value,
      ssi_boot_results$mixed_spont$null_value,
      ssi_boot_results$mixed_vis$null_value,
      ssi_boot_results$synergy_index$null_value,
      ssi_boot_results$synergy_ratio$null_value
    ),
    p = c(
      ssi_boot_results$vis_spont$p_twoside,
      ssi_boot_results$mixed_spont$p_twoside,
      ssi_boot_results$mixed_vis$p_twoside,
      ssi_boot_results$synergy_index$p_twoside,
      ssi_boot_results$synergy_ratio$p_twoside
    )
  ) %>%
    dplyr::mutate(
      Signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE      ~ "n.s."
      )
    )
}

plot_boot <- function(boot_res, title = "Bootstrap Distribution") {
  library(ggplot2)
  
  null <- boot_res$null_value   # automatically get the correct null
  
  tibble(x = unlist(boot_res$distribution)) %>%
    ggplot(aes(x)) +
    geom_histogram(bins = 50, color="white", fill="steelblue") +
    geom_vline(xintercept = null, color="red", linewidth = 1) +
    labs(
      title = title,
      subtitle = paste("Null =", null),
      x = "Bootstrap mean",
      y = "Count"
    ) +
    theme_bw()
}


plot_all_boot <- function(ssi_boot_results) {
  library(ggplot2)
  library(patchwork)
  
  plot_list <- lapply(names(ssi_boot_results), function(name) {
    plot_boot(ssi_boot_results[[name]], title = name)
  })
  
  wrap_plots(plot_list, ncol = 2)
}

plot_fish_clouds_dualpanel_boot <- function(
    df,
    ssi_boot_results,
    palette = "Dark2",
    spacing = 10,        # SAME default as LMM version
    jitter_width = 0.02, # SAME default as LMM version
    point_size = 1       # SAME default as LMM version
) {
  require(ggplot2)
  require(dplyr)
  require(ggpubr)
  require(tidyr)
  
  # ---- PREPARE DATA ----
  df_long <- df %>%
    select(
      fish,
      vis_spont_index,
      mixed_spont_index,
      mixed_vis_index,
      synergy_index,
      synergy
    ) %>%
    pivot_longer(
      cols = -fish,
      names_to = "index",
      values_to = "value"
    ) %>%
    mutate(
      index = factor(
        index,
        levels = c(
          "vis_spont_index",
          "mixed_spont_index",
          "mixed_vis_index",
          "synergy_index",
          "synergy"
        ),
        labels = c(
          "Vis vs Spont",
          "Mixed vs Spont",
          "Mixed vs Vis",
          "Synergy Index",
          "Synergy Ratio"
        )
      )
    )
  
  # ---- jitter position using same spacing logic ----
  df_long <- df_long %>%
    group_by(index) %>%
    mutate(
      fish = factor(fish),
      fish_n = as.numeric(fish),
      x = as.numeric(index) + (fish_n - mean(fish_n)) / spacing
    ) %>%
    ungroup()
  
  # ---- BOOTSTRAP summaries ----
  boot_summary <- make_ssi_boot_table(ssi_boot_results)
  boot_summary$index <- factor(boot_summary$name,
                               levels = levels(df_long$index))
  
  # ---- split ----
  df_idx <- df_long %>%
    filter(index != "Synergy Ratio") %>%
    mutate(index = droplevels(index))
  
  df_syn <- df_long %>%
    filter(index == "Synergy Ratio") %>%
    mutate(index = droplevels(index))
  
  res_idx <- boot_summary %>%
    filter(index != "Synergy Ratio") %>%
    mutate(index = droplevels(index))
  
  res_syn <- boot_summary %>%
    filter(index == "Synergy Ratio") %>%
    mutate(index = droplevels(index))
  
  # ---- Panel A ----
  p_idx <- ggplot() +
    geom_jitter(
      data = df_idx,
      aes(x = x, y = value, color = fish),
      size = point_size,
      alpha = 0.6,
      width = jitter_width
    ) +
    geom_errorbar(
      data = res_idx,
      aes(x = as.numeric(index), ymin = CI_low, ymax = CI_high),
      color = "black",
      linewidth = 1,
      width = 0
    ) +
    geom_point(
      data = res_idx,
      aes(x = as.numeric(index), y = Estimate),
      size = 3,
      shape = 21,
      fill = "gray60",
      color = "black",
      stroke = 1.1
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    scale_color_brewer(palette = palette) +
    scale_x_continuous(
      breaks = seq_along(levels(df_idx$index)),
      labels = levels(df_idx$index)
    ) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      title = "A. Selectivity Indices (Bootstrap)",
      x = "",
      y = "Index (−1 to 1)",
      color = "Fish"
    ) +
    theme_ssi(base_size = 10) +
    theme(legend.position = "none")
  
  # ---- Panel B: Synergy Ratio ----
  
  # force both layers to use the SAME x-position
  x_syn <- 1  
  
  p_syn <- ggplot() +
    geom_jitter(
      data = df_syn,
      aes(x = x_syn + (fish_n - mean(fish_n)) / spacing, y = value, color = fish),
      size = point_size,
      alpha = 0.6,
      width = jitter_width
    ) +
    geom_errorbar(
      data = res_syn,
      aes(x = x_syn, ymin = CI_low, ymax = CI_high),
      color = "black",
      linewidth = 1,
      width = 0
    ) +
    geom_point(
      data = res_syn,
      aes(x = x_syn, y = Estimate),
      size = 3,
      color = "black",
      fill = "gray60",
      shape = 21,
      stroke = 1.1
    ) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "gray60") +
    scale_color_brewer(palette = palette) +
    scale_x_continuous(
      breaks = 1,
      labels = "Synergy Ratio"
    ) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(
      title = "B. Synergy Ratio (Bootstrap)",
      x = "",
      y = "Synergy Ratio",
      color = "Fish"
    ) +
    theme_ssi(base_size = 10) +
    theme(legend.position = "none")
  
  # ---- combine panels ----
  ggarrange(
    p_idx + theme(plot.margin = margin(5, 15, 5, 5)),
    p_syn + theme(plot.margin = margin(5, 5, 5, 10)),
    ncol = 2,
    widths = c(3, 1)
  )
}

save_boot_distributions_csv <- function(boot_results, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (name in names(boot_results)) {
    dist <- boot_results[[name]]$distribution[[1]]
    out_file <- file.path(out_dir, paste0("boot_dist_", name, ".csv"))
    
    write.csv(
      data.frame(value = dist),
      out_file,
      row.names = FALSE
    )
  }
}

save_boot_distributions_combined <- function(boot_results, filename) {
  # Extract all distributions, names become column names
  dist_list <- lapply(boot_results, function(x) x$distribution[[1]])
  
  # Combine into data frame
  df <- as.data.frame(dist_list)
  
  # Write CSV
  write.csv(df, filename, row.names = FALSE)
}

write_table_txt <- function(tbl, filename) {
  txt <- capture.output(print(tbl))
  writeLines(txt, filename)
}

save_glmm_full_report <- function(
    glmm_summary,
    emm_log,
    emm_resp,
    pairs_log,
    pairs_resp,
    synergy,
    filename
) {
  # Clean delete if exists — no warnings
  if (file.exists(filename)) {
    try(file.remove(filename), silent = TRUE)
  }
  
  con <- file(filename, open = "wt")
  on.exit(close(con))
  
  write_section <- function(title, obj) {
    writeLines(paste0("\n=== ", title, " ===\n"), con)
    writeLines(capture.output(print(obj)), con)
  }
  
  # ----- helper to add significance stars -----
  add_stars <- function(df) {
    df <- as.data.frame(df)
    if (!("p.value" %in% names(df))) return(df)
    
    df$stars <- cut(df$p.value,
                    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", ".", " ")
    )
    df
  }
  
  # ---- write each section ----
  write_section("Summary of glmm_score", glmm_summary)
  write_section("EMMEANS (log scale)", emm_log)
  write_section("EMMEANS (response scale)", emm_resp)
  write_section("Pairwise comparisons (log scale)", add_stars(pairs_log))
  write_section("Pairwise comparisons (response scale)", add_stars(pairs_resp))
  write_section("Additive Synergy", add_stars(synergy))
  
  message("Saved report to: ", filename)
}


save_plot <- function(plot, base_path, width, height, units = "cm", dpi = 300) {
  ggsave(paste0(base_path, ".pdf"), plot, width = width, height = height, units = units)
  # ggsave(paste0(base_path, ".svg"), plot, width = width, height = height, units = units, dpi = dpi)
  ggsave(paste0(base_path, ".png"), plot, width = width, height = height, units = units, dpi = dpi)
  # ggsave(paste0(base_path, ".jpg"), plot, width = width, height = height, units = units, dpi = dpi)
}
