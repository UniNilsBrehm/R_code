  # ==============================================================================
  # Aggregated response probability
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
  library(readr)
  library(tidyr)
  
  # ==============================================================================
  # Paths
  # ==============================================================================
  # source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_Massed_vs_Spaced//utils.R")
  source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/NLME/DarkFlash_Massed_vs_Spaced//utils.R")
  
  # base_dir <- "D:/WorkingData/Susana/NLME/DarkFlash_Joint_SpacedVsMassed/"
  base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/NLME/DarkFlash_Joint_SpacedVsMassed/"
  
  file_massed <- file.path(
    base_dir,
    "data_files",
    "SPZ_Massed_Training_7Nov2025.csv"
  )
  
  file_spaced <- file.path(
    base_dir,
    "data_files",
    "SPZ_Spaced_Training_Nov2025.csv"    
  )
  
  var_name <- "response_prob"
  col_name <- "move"
  
  save_fig_dir     <- file.path(base_dir, "figs",    "nlme_joint_response_prob", var_name)
  save_results_dir <- file.path(base_dir, "results", "nlme_joint_response_prob", var_name)
  save_model_dir   <- file.path(base_dir, "models")
  
  dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)
  dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(save_model_dir,   recursive = TRUE, showWarnings = FALSE)
  
  
  # ==============================================================================
  # Load and prepare both datasets
  # ==============================================================================
  res_massed <- load_data(file_massed, move_th = 1, drop = c('th2, tyr', 'tyr'))
  res_spaced <- load_data(file_spaced, move_th = 1, drop = c('th2, tyr', 'tyr'))
  
  # Use df_final (NOT df_final_sub) -- df_final_sub is responders-only.
  df_massed <- res_massed$df_final
  df_spaced <- res_spaced$df_final
  
  # ------------------------------------------------------------------------------
  # Tag each row with its Training condition and define BlockRole
  # ------------------------------------------------------------------------------
  # In each experiment, the LAST block is the memory test. All preceding blocks
  # are training blocks.
  massed_blocks_train <- "1"            # massed: block 1 = training
  massed_block_test   <- "2"            # massed: block 2 = test
  
  spaced_blocks_train <- c("1","2","3","4")   # spaced: blocks 1-4 = training
  spaced_block_test   <- "5"                  # spaced: block 5 = test
  
  
  df_massed_tagged <- df_massed %>%
    mutate(
      Training  = "massed",
      BlockRole = ifelse(as.character(Block) == massed_block_test, "test", "training")
    )
  
  df_spaced_tagged <- df_spaced %>%
    mutate(
      Training  = "spaced",
      BlockRole = ifelse(as.character(Block) == spaced_block_test, "test", "training")
    )
  
  
  # ------------------------------------------------------------------------------
  # CRITICAL: animal IDs must be unique across experiments.
  # Same Video x Well combo from different experiments would otherwise collide.
  # We prefix the animal label with the Training condition.
  # ------------------------------------------------------------------------------
  df_all <- bind_rows(df_massed_tagged, df_spaced_tagged) %>%
    mutate(
      stimulus    = as.numeric(stimulus),
      stimulus0   = stimulus - 1,
      stimulus_log   = log(stimulus),
      Training    = factor(Training,  levels = c("massed", "spaced")),
      Block       = factor(Block),
      BlockRole   = factor(BlockRole, levels = c("training", "test")),
      Genotype    = factor(Genotype),
      Video       = factor(Video),
      Well        = factor(Well),
      animal      = factor(paste0(Training, "_", Video, ".", Well))
    )
  
  # # Remove all non-responders to stimulus 1 in Block 1
  # df_filtered <- df_all %>%
  #   group_by(animal) %>%
  #   filter(!any(Block == 1 & stimulus == 1 & move == 0)) %>%
  #   ungroup()
  # 
  # summary_compare <- bind_rows(
  #   df_all %>%
  #     distinct(animal, Genotype, Training) %>%
  #     mutate(dataset = "before"),
  #   
  #   df_filtered %>%
  #     distinct(animal, Genotype, Training) %>%
  #     mutate(dataset = "after")
  # ) %>%
  #   dplyr::count(dataset, Genotype, Training)
  # 
  # summary_compare
  # df_all <- df_filtered

  
  # ============================================================================
  # Thresholds
  ggplot(df_spaced, aes(x = max_cumsum)) +
    geom_histogram(bins = 100)
  
  # Log-scale view, useful because movement amplitudes are usually right-skewed
  ggplot(df_all %>% filter(max_peak > 0), aes(x = log10(max_peak))) +
    geom_histogram(bins = 100) +
    theme_pubr() +
    labs(
      x = "log10(max_peak)",
      y = "Count",
      title = "Distribution of non-zero max_peak values"
    )
  
  # Test effect of different thresholds on overall response prob.
  thresholds <- seq(0, 2, by = 0.1)
  
  threshold_summary <- purrr::map_dfr(thresholds, function(th) {
    
    df_all %>%
      mutate(move_tmp = as.numeric(max_peak > th)) %>%
      group_by(Training, BlockRole, Genotype) %>%
      summarise(
        p_move = mean(move_tmp, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(threshold = th)
  })
  
  ggplot(
    threshold_summary,
    aes(x = threshold, y = p_move, color = Genotype)
  ) +
    geom_line(linewidth = 1) +
    facet_grid(BlockRole ~ Training) +
    theme_pubr() +
    labs(
      x = "Threshold on max_peak",
      y = "Response probability",
      title = "Sensitivity of response probability to movement threshold"
    )
  
  # Use mixture model to find threshold
  library(mclust)
  
  peak_dat <- df_all %>%
    filter(max_peak > 0) %>%
    mutate(log_peak = log10(max_peak)) %>%
    pull(log_peak)
  
  mix_fit <- Mclust(peak_dat, G = 2)
  
  summary(mix_fit)
  plot(mix_fit, what = 'density')
  
  df_mix <- tibble(
    log_peak = peak_dat,
    class = factor(mix_fit$classification)
  )
  
  ggplot(df_mix, aes(x = log_peak, fill = class)) +
    geom_histogram(bins = 100, alpha = 0.6, position = "identity") +
    theme_pubr() +
    labs(
      x = "log10(max_peak)",
      y = "Count",
      title = "Mixture-model classification of response amplitudes"
    )
  
  
  # ==============================================================================
  # Aggregated response probability + no-asymptote fit using stimulus_log
  # ==============================================================================
  
  raw_summary_joint <- df_all %>%
    group_by(
      Training,
      Block,
      BlockRole,
      Genotype,
      stimulus,
      stimulus0,
      stimulus_log
    ) %>%
    summarise(
      n_trials = sum(!is.na(move)),
      n_move   = sum(move, na.rm = TRUE),
      p_move   = mean(move, na.rm = TRUE),
      .groups  = "drop"
    )
  
  
  # ==============================================================================
  # Fit no-asymptote model using stimulus_log
  # ==============================================================================
  
  # Model:
  #
  #   p_move = amplitude * exp(-rate * stimulus_log)
  #
  # Since stimulus_log = log(stimulus), this is equivalent to:
  #
  #   p_move = amplitude * stimulus^(-rate)
  #
  # Parameters:
  #   amplitude = predicted response probability at stimulus 1
  #   rate      = decay strength
  #
  # No asymptote is included.
  
  fit_exp_no_asym_logstim_one_group <- function(dat) {
    
    dat <- dat %>%
      filter(
        is.finite(stimulus_log),
        is.finite(p_move)
      ) %>%
      arrange(stimulus_log)
    
    n_stim <- length(unique(dat$stimulus_log))
    
    if (nrow(dat) < 3 || n_stim < 3) {
      return(tibble(
        amplitude = NA_real_,
        rate = NA_real_,
        fit_type = "not_fitted",
        converged = FALSE
      ))
    }
    
    start_p <- dat$p_move[which.min(dat$stimulus_log)]
    
    fit <- tryCatch(
      nls(
        p_move ~ amplitude * exp(-rate * stimulus_log),
        data = dat,
        start = list(
          amplitude = min(max(start_p, 0.001), 0.999),
          rate = 0.5
        ),
        algorithm = "port",
        lower = c(
          amplitude = 0,
          rate = 1e-6
        ),
        upper = c(
          amplitude = 1,
          rate = 20
        ),
        control = nls.control(
          maxiter = 500,
          warnOnly = TRUE
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      return(tibble(
        amplitude = NA_real_,
        rate = NA_real_,
        fit_type = "failed",
        converged = FALSE
      ))
    }
    
    coefs <- coef(fit)
    
    tibble(
      amplitude = unname(coefs["amplitude"]),
      rate = unname(coefs["rate"]),
      fit_type = "no_asym_logstim",
      converged = TRUE
    )
  }
  
  
  exp_params <- raw_summary_joint %>%
    group_by(Training, Block, BlockRole, Genotype) %>%
    group_modify(~ fit_exp_no_asym_logstim_one_group(.x)) %>%
    ungroup()
  
  print(exp_params)
  
  
  # ==============================================================================
  # Build prediction grid
  # ==============================================================================
  
  new_data_exp <- raw_summary_joint %>%
    group_by(Training, Block, BlockRole, Genotype) %>%
    summarise(
      stim_min = min(stimulus, na.rm = TRUE),
      stim_max = max(stimulus, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      exp_params,
      by = c("Training", "Block", "BlockRole", "Genotype")
    ) %>%
    rowwise() %>%
    mutate(
      grid = list(tibble(
        stimulus = seq(stim_min, stim_max, length.out = 200)
      ))
    ) %>%
    unnest(grid) %>%
    ungroup() %>%
    mutate(
      stimulus0 = stimulus - 1,
      stimulus_log = log(stimulus),
      fit = amplitude * exp(-rate * stimulus_log),
      fit = pmin(pmax(fit, 0), 1),
      Training = factor(Training, levels = levels(df_all$Training)),
      Block    = factor(Block, levels = levels(df_all$Block)),
      Genotype = factor(Genotype, levels = levels(df_all$Genotype))
    )
  
  
  # ==============================================================================
  # Massed plot
  # ==============================================================================
  
  p_massed_curves_exp <- ggplot(
    new_data_exp %>% filter(Training == "massed"),
    aes(x = stimulus, color = Genotype)
  ) +
    facet_grid(Genotype ~ Block, scales = "free_x") +
    
    geom_point(
      data = raw_summary_joint %>% filter(Training == "massed"),
      aes(x = stimulus, y = p_move, color = Genotype),
      inherit.aes = FALSE,
      alpha = 0.45,
      size = 1.0
    ) +
    
    geom_line(
      aes(y = fit),
      linewidth = 1.1,
      na.rm = TRUE
    ) +
    
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_pubr(base_size = 12) +
    labs(
      x = "Stimulus number within block",
      y = "Response probability",
      title = "Aggregated response probability: Massed training",
      subtitle = "Points = aggregated means; lines = no-asymptote fit using log(stimulus)"
    ) +
    theme(
      legend.position = "none"
    )
  
  
  # ==============================================================================
  # Spaced plot
  # ==============================================================================
  
  p_spaced_curves_exp <- ggplot(
    new_data_exp %>% filter(Training == "spaced"),
    aes(x = stimulus, color = Genotype)
  ) +
    facet_grid(Genotype ~ Block, scales = "free_x") +
    
    geom_point(
      data = raw_summary_joint %>% filter(Training == "spaced"),
      aes(x = stimulus, y = p_move, color = Genotype),
      inherit.aes = FALSE,
      alpha = 0.45,
      size = 1.0
    ) +
    
    geom_line(
      aes(y = fit),
      linewidth = 1.1,
      na.rm = TRUE
    ) +
    
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_pubr(base_size = 12) +
    labs(
      x = "Stimulus number within block",
      y = "Response probability",
      title = "Aggregated response probability: Spaced training",
      subtitle = "Points = aggregated means; lines = no-asymptote fit using log(stimulus)"
    ) +
    theme(
      legend.position = "top"
    )
  
  
  print(p_spaced_curves_exp)
  print(p_massed_curves_exp)
  
