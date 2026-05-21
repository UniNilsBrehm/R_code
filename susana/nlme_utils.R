# I/O

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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RESPONSE PROBABILITY
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DIAGNOSTIC
# ------------------------------------------------------------------------------
diagnose_brms_nlme_model <- function(
    fit_model,
    df = NULL,
    save_results_dir,
    var_name = "response_prob",
    response_col = "move",
    ndraws_ppc = 100
) {
  library(brms)
  library(posterior)
  library(bayesplot)
  library(ggplot2)
  library(dplyr)
  library(loo)
  
  dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
  
  diag_txt <- file.path(save_results_dir, paste0("nlme_", var_name, "_diagnostics.txt"))
  
  sink(diag_txt)
  
  cat("============================================================\n")
  cat("MODEL SUMMARY\n")
  cat("============================================================\n\n")
  print(summary(fit_model))
  
  cat("\n\n============================================================\n")
  cat("SAMPLER DIAGNOSTICS\n")
  cat("============================================================\n\n")
  
  np <- nuts_params(fit_model)
  
  n_div <- sum(np$Parameter == "divergent__" & np$Value == 1)
  n_tree <- sum(np$Parameter == "treedepth__" & np$Value >= 12)
  
  cat("Number of divergent transitions:", n_div, "\n")
  cat("Number of max treedepth hits:", n_tree, "\n\n")
  
  cat("Rule of thumb:\n")
  cat("- Divergences should be 0.\n")
  cat("- Max treedepth hits should ideally be 0 or very rare.\n")
  cat("- Rhat should be close to 1.00, usually < 1.01.\n")
  cat("- Bulk ESS and tail ESS should not be very small.\n\n")
  
  cat("\n\n============================================================\n")
  cat("RHAT / ESS CHECK\n")
  cat("============================================================\n\n")
  
  draws <- as_draws_df(fit_model)
  summ <- summarise_draws(draws)
  
  print(
    summ %>%
      arrange(desc(rhat)) %>%
      select(variable, mean, sd, rhat, ess_bulk, ess_tail) %>%
      head(50)
  )
  
  cat("\nWorst Rhat:\n")
  print(max(summ$rhat, na.rm = TRUE))
  
  cat("\nSmallest bulk ESS:\n")
  print(min(summ$ess_bulk, na.rm = TRUE))
  
  cat("\nSmallest tail ESS:\n")
  print(min(summ$ess_tail, na.rm = TRUE))
  
  cat("\n\n============================================================\n")
  cat("LOO\n")
  cat("============================================================\n\n")
  
  loo_res <- tryCatch(
    loo(fit_model),
    error = function(e) e
  )
  
  print(loo_res)
  
  cat("\n\n============================================================\n")
  cat("BAYES R2\n")
  cat("============================================================\n\n")
  
  r2_res <- tryCatch(
    bayes_R2(fit_model),
    error = function(e) e
  )
  
  print(r2_res)
  
  sink()
  
  # ---------------------------------------------------------------------------
  # Trace plots for key population parameters
  # ---------------------------------------------------------------------------
  key_pars <- variables(draws)
  key_pars <- key_pars[grepl("^b_", key_pars)]
  
  if (length(key_pars) > 0) {
    p_trace <- mcmc_trace(
      as_draws_array(fit_model),
      pars = key_pars[seq_len(min(12, length(key_pars)))]
    )
    
    ggsave(
      file.path(save_results_dir, paste0("nlme_", var_name, "_traceplots.png")),
      p_trace,
      width = 12,
      height = 8,
      dpi = 300,
      bg = "white"
    )
  }
  
  # ---------------------------------------------------------------------------
  # Posterior predictive checks
  # ---------------------------------------------------------------------------
  p_ppc_bars <- pp_check(fit_model, type = "bars", ndraws = ndraws_ppc)
  
  ggsave(
    file.path(save_results_dir, paste0("nlme_", var_name, "_ppcheck_bars.png")),
    p_ppc_bars,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  p_ppc_stat <- pp_check(fit_model, type = "stat", stat = "mean", ndraws = ndraws_ppc)
  
  ggsave(
    file.path(save_results_dir, paste0("nlme_", var_name, "_ppcheck_mean.png")),
    p_ppc_stat,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  # ---------------------------------------------------------------------------
  # Calibration plot: observed response rate vs fitted probability
  # ---------------------------------------------------------------------------
  if (!is.null(df)) {
    fitted_df <- df %>%
      mutate(
        fitted_prob = fitted(fit_model, summary = TRUE)[, "Estimate"],
        fitted_bin = cut(
          fitted_prob,
          breaks = quantile(fitted_prob, probs = seq(0, 1, 0.1), na.rm = TRUE),
          include.lowest = TRUE
        )
      ) %>%
      group_by(fitted_bin) %>%
      summarise(
        mean_predicted = mean(fitted_prob, na.rm = TRUE),
        mean_observed = mean(.data[[response_col]], na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
    
    p_cal <- ggplot(fitted_df, aes(x = mean_predicted, y = mean_observed)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_point(aes(size = n)) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      theme_bw() +
      labs(
        x = "Mean fitted probability",
        y = "Observed response rate",
        size = "N",
        title = "Calibration plot"
      )
    
    ggsave(
      file.path(save_results_dir, paste0("nlme_", var_name, "_calibration.png")),
      p_cal,
      width = 6,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
  
  message("Diagnostics written to: ", diag_txt)
  invisible(TRUE)
}

# ==============================================================================
# Helper functions for response prob model comparison
# ==============================================================================
make_nlpar_draws <- function(fit_model, df_resp, nlpar, transform_fun = identity) {
  
  grid <- expand.grid(
    Genotype = levels(df_resp$Genotype),
    Block = levels(df_resp$Block)
  ) %>%
    mutate(
      stimulus0 = 0,
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block = factor(Block, levels = levels(df_resp$Block))
    )
  
  draws <- posterior_linpred(
    fit_model,
    newdata = grid,
    nlpar = nlpar,
    re_formula = NA,
    transform = FALSE
  )
  
  bind_rows(
    lapply(seq_len(nrow(grid)), function(i) {
      tibble(
        draw = seq_len(nrow(draws)),
        Genotype = grid$Genotype[i],
        Block = grid$Block[i],
        value_raw = draws[, i],
        value = transform_fun(draws[, i])
      )
    })
  )
}


summarise_nlpar <- function(draws_df, value_name) {
  
  draws_df %>%
    group_by(Genotype, Block) %>%
    summarise(
      median = median(value),
      low = quantile(value, 0.025),
      high = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    rename(
      !!paste0(value_name, "_median") := median,
      !!paste0(value_name, "_low") := low,
      !!paste0(value_name, "_high") := high
    )
}


compare_nlpar <- function(draws_df, value_name, rope = 0.02, ratio = FALSE) {
  
  comp <- draws_df %>%
    rename(
      Genotype_1 = Genotype,
      value_1 = value
    ) %>%
    inner_join(
      draws_df %>%
        rename(
          Genotype_2 = Genotype,
          value_2 = value
        ),
      by = c("draw", "Block"),
      relationship = "many-to-many"
    ) %>%
    filter(as.integer(Genotype_1) < as.integer(Genotype_2)) %>%
    mutate(
      comparison = paste(Genotype_1, "vs", Genotype_2),
      difference = value_1 - value_2
    )
  
  if (ratio) {
    
    comp <- comp %>%
      mutate(ratio_value = value_1 / value_2)
    
    summary <- comp %>%
      group_by(Block, Genotype_1, Genotype_2, comparison) %>%
      summarise(
        median_difference = median(difference),
        diff_low = quantile(difference, 0.025),
        diff_high = quantile(difference, 0.975),
        
        median_ratio = median(ratio_value),
        ratio_low = quantile(ratio_value, 0.025),
        ratio_high = quantile(ratio_value, 0.975),
        
        prob_Genotype_1_higher = mean(value_1 > value_2),
        prob_Genotype_2_higher = mean(value_1 < value_2),
        
        ROPE_prob_10percent = mean(abs(log(ratio_value)) < log(1.10)),
        
        evidence_strength = case_when(
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.995 ~ "extreme",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.97  ~ "very strong",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.90  ~ "strong",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.75  ~ "moderate",
          TRUE ~ "weak"
        ),
        .groups = "drop"
      )
    
  } else {
    
    summary <- comp %>%
      group_by(Block, Genotype_1, Genotype_2, comparison) %>%
      summarise(
        median_difference = median(difference),
        diff_low = quantile(difference, 0.025),
        diff_high = quantile(difference, 0.975),
        
        prob_Genotype_1_higher = mean(value_1 > value_2),
        prob_Genotype_2_higher = mean(value_1 < value_2),
        
        ROPE_prob_small_effect = mean(abs(difference) < rope),
        
        evidence_strength = case_when(
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.995 ~ "extreme",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.97  ~ "very strong",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.90  ~ "strong",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.75  ~ "moderate",
          TRUE ~ "weak"
        ),
        .groups = "drop"
      )
  }
  
  list(draws = comp, summary = summary)
}

# ==============================================================================
# Helper functions for response probability model comparison
# ==============================================================================

make_nlpar_draws <- function(fit_model, df_resp, nlpar, transform_fun = identity) {
  
  grid <- expand.grid(
    Genotype = levels(df_resp$Genotype),
    Block = levels(df_resp$Block)
  ) %>%
    mutate(
      stimulus0 = 0,
      Genotype = factor(Genotype, levels = levels(df_resp$Genotype)),
      Block = factor(Block, levels = levels(df_resp$Block))
    )
  
  draws <- posterior_linpred(
    fit_model,
    newdata = grid,
    nlpar = nlpar,
    re_formula = NA,
    transform = FALSE
  )
  
  bind_rows(
    lapply(seq_len(nrow(grid)), function(i) {
      tibble(
        draw = seq_len(nrow(draws)),
        Genotype = grid$Genotype[i],
        Block = grid$Block[i],
        value_raw = draws[, i],
        value = transform_fun(draws[, i])
      )
    })
  )
}


summarise_nlpar <- function(draws_df, value_name) {
  
  draws_df %>%
    group_by(Genotype, Block) %>%
    summarise(
      median = median(value),
      low = quantile(value, 0.025),
      high = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    rename(
      !!paste0(value_name, "_median") := median,
      !!paste0(value_name, "_low") := low,
      !!paste0(value_name, "_high") := high
    )
}


compare_nlpar <- function(draws_df, value_name, rope = 0.02, ratio = FALSE) {
  
  comp <- draws_df %>%
    rename(
      Genotype_1 = Genotype,
      value_1 = value
    ) %>%
    inner_join(
      draws_df %>%
        rename(
          Genotype_2 = Genotype,
          value_2 = value
        ),
      by = c("draw", "Block"),
      relationship = "many-to-many"
    ) %>%
    filter(as.integer(Genotype_1) < as.integer(Genotype_2)) %>%
    mutate(
      comparison = paste(Genotype_1, "vs", Genotype_2),
      difference = value_1 - value_2
    )
  
  if (ratio) {
    
    comp <- comp %>%
      mutate(ratio_value = value_1 / value_2)
    
    summary <- comp %>%
      group_by(Block, Genotype_1, Genotype_2, comparison) %>%
      summarise(
        median_difference = median(difference),
        diff_low = quantile(difference, 0.025),
        diff_high = quantile(difference, 0.975),
        
        median_ratio = median(ratio_value),
        ratio_low = quantile(ratio_value, 0.025),
        ratio_high = quantile(ratio_value, 0.975),
        
        prob_Genotype_1_higher = mean(value_1 > value_2),
        prob_Genotype_2_higher = mean(value_1 < value_2),
        
        ROPE_prob_10percent = mean(abs(log(ratio_value)) < log(1.10)),
        
        evidence_strength = case_when(
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.995 ~ "extreme",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.97  ~ "very strong",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.90  ~ "strong",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.75  ~ "moderate",
          TRUE ~ "weak"
        ),
        .groups = "drop"
      )
    
  } else {
    
    summary <- comp %>%
      group_by(Block, Genotype_1, Genotype_2, comparison) %>%
      summarise(
        median_difference = median(difference),
        diff_low = quantile(difference, 0.025),
        diff_high = quantile(difference, 0.975),
        
        prob_Genotype_1_higher = mean(value_1 > value_2),
        prob_Genotype_2_higher = mean(value_1 < value_2),
        
        ROPE_prob_small_effect = mean(abs(difference) < rope),
        
        evidence_strength = case_when(
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.995 ~ "extreme",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.97  ~ "very strong",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.90  ~ "strong",
          pmax(prob_Genotype_1_higher, prob_Genotype_2_higher) > 0.75  ~ "moderate",
          TRUE ~ "weak"
        ),
        .groups = "drop"
      )
  }
  
  list(draws = comp, summary = summary)
}


plot_posterior_density <- function(
    draws_df,
    xvar,
    title,
    xlab,
    filename,
    save_fig_dir,
    block_limits = NULL
) {
  
  p <- ggplot(
    draws_df,
    aes(x = {{ xvar }}, fill = Genotype)
  ) +
    facet_wrap(~Block, scales = "free") +
    geom_density(alpha = 0.35) +
    theme_pubr(base_size = 14) +
    labs(
      x = xlab,
      y = "Posterior density",
      title = title,
      fill = "Genotype"
    )
  
  if (!is.null(block_limits)) {
    
    scale_list <- lapply(
      names(block_limits),
      function(b) {
        as.formula(
          paste0(
            'Block == "', b,
            '" ~ scale_x_continuous(limits = c(',
            block_limits[[b]][1], ", ",
            block_limits[[b]][2], "))"
          )
        )
      }
    )
    
    p <- p +
      ggh4x::facetted_pos_scales(x = scale_list)
  }
  
  print(p)
  
  ggsave(
    file.path(save_fig_dir, filename),
    p,
    width = 10,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  p
}


plot_pairwise_differences <- function(
    summary_df,
    ylab,
    title,
    filename,
    save_fig_dir
) {
  
  p <- summary_df %>%
    ggplot(
      aes(
        x = comparison,
        y = median_difference,
        ymin = diff_low,
        ymax = diff_high,
        color = Genotype_1
      )
    ) +
    facet_wrap(~Block, scales = "free_x") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(linewidth = 0.8) +
    coord_flip() +
    theme_pubr(base_size = 14) +
    labs(
      x = "Comparison",
      y = ylab,
      title = title,
      color = "Genotype 1"
    )
  
  print(p)
  
  ggsave(
    file.path(save_fig_dir, filename),
    p,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  p
}