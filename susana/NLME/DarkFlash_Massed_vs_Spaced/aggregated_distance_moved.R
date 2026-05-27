  # ==============================================================================
  # Aggregated max_peak and max_cumsum — exponential decay fit using log(stimulus)
  # ==============================================================================
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
  library(readr)
  library(minpack.lm)
  
  # ==============================================================================
  # Paths
  # ==============================================================================
  source("C:/UniFreiburg/Code/R_code/susana/NLME/DarkFlash_Massed_vs_Spaced//utils.R")
  # source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/NLME/DarkFlash_Massed_vs_Spaced//utils.R")
  
  base_dir <- "D:/WorkingData/Susana/NLME/DarkFlash_Joint_SpacedVsMassed/"
  # base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/NLME/DarkFlash_Joint_SpacedVsMassed/"
  
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
  
  # ----------------------------------------------------------------------------
  # Settings
  select_first_responders <- FALSE
  conditional_on_response <- TRUE
  th_for_move <- 1
  
  
  # ==============================================================================
  # Load and prepare both data sets
  # ==============================================================================
  res_massed <- load_data(file_massed, move_th = th_for_move, drop = c('th2, tyr', 'tyr'))
  res_spaced <- load_data(file_spaced, move_th = th_for_move, drop = c('th2, tyr', 'tyr'))
  
  res_massed <- load_data(file_massed, move_th = th_for_move)
  res_spaced <- load_data(file_spaced, move_th = th_for_move)
  
  
  if (conditional_on_response) {
    df_massed <- res_massed$df_final_sub
    df_spaced <- res_spaced$df_final_sub
  }else{
    df_massed <- res_massed$df_final
    df_spaced <- res_spaced$df_final
  }
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
  
  # Remove all non-responders to stimulus 1 in Block 1
  if (select_first_responders) {
    message("Removing all non-responders to first stimulus")
    df_filtered <- df_all %>%
      group_by(animal) %>%
      filter(!any(Block == 1 & stimulus == 1 & move == 0)) %>%
      ungroup()
    
    summary_compare <- bind_rows(
      df_all %>%
        distinct(animal, Genotype, Training) %>%
        mutate(dataset = "before"),
      
      df_filtered %>%
        distinct(animal, Genotype, Training) %>%
        mutate(dataset = "after")
    ) %>%
      dplyr::count(dataset, Genotype, Training)
    
    summary_compare
    df_all <- df_filtered
  }
  
  # ------------------------------------------------------------------------------
  # Variables to loop over
  # col_name   : column in df_all to average
  # var_name   : used for directory names and plot labels
  # y_label    : y-axis label in plots
  # ------------------------------------------------------------------------------
  targets <- list(
    list(col_name = "max_peak",   var_name = "max_peak",   y_label = "Mean peak distance moved"),
    list(col_name = "max_cumsum", var_name = "max_cumsum", y_label = "Mean summed distance moved")
  )
  
  # ==============================================================================
  # Switches — change these only
  # ==============================================================================
  use_log_stimulus  <- TRUE   # TRUE = log(stimulus),  FALSE = raw stimulus
  use_asymptote     <- FALSE  # TRUE = add asymptote,  FALSE = no asymptote
  use_double_exp    <- FALSE   # TRUE = double exponential, FALSE = single exponential
  fit_max_stimulus  <- Inf   # Inf = use all;  e.g. 100 = fit only stimulus 1:100
  extrapolate_curve <- TRUE   # TRUE = draw curve over full range, FALSE = clip to fit range
  
  # ==============================================================================
  # Fitting function
  # ==============================================================================
  fit_exp_logstim_continuous <- function(dat, response_col, use_log, use_asym, use_double) {
    
    dat <- dat %>%
      rename(y_resp = all_of(response_col)) %>%
      mutate(x = if (use_log) stimulus_log else stimulus) %>%
      filter(is.finite(x), is.finite(y_resp), !is.na(y_resp)) %>%
      arrange(x)
    
    n_stim <- length(unique(dat$x))
    if (nrow(dat) < 3 || n_stim < 3) {
      return(tibble(
        amplitude_fast = NA_real_, rate_fast = NA_real_,
        amplitude_slow = NA_real_, rate_slow = NA_real_,
        asymptote      = NA_real_,
        fit_type       = "not_fitted",
        converged      = FALSE,
        AIC            = NA_real_,
        BIC            = NA_real_,
        RSS            = NA_real_,
        n_obs          = nrow(dat)
      ))
    }
    
    # --- starting values via log-linear approximation ---
    lm_init <- tryCatch(
      lm(log(y_resp + 1e-6) ~ x, data = dat),
      error = function(e) NULL
    )
    if (!is.null(lm_init)) {
      start_amp  <- max(exp(coef(lm_init)[1]), 1e-6)
      start_rate <- max(-coef(lm_init)[2],     1e-6)
    } else {
      start_amp  <- max(dat$y_resp[which.min(dat$x)], 1e-6)
      start_rate <- 0.5
    }
    start_asym <- max(min(dat$y_resp, na.rm = TRUE), 0)
    
    rate_candidates <- unique(c(start_rate, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0))
    
    # --- inner fitter: tries all rate candidates for a given formula/bounds ---
    try_fit <- function(formula_nls, starts, lower_bounds, upper_bounds, is_double) {
      fit <- NULL
      for (rate_start in rate_candidates) {
        starts$rate_fast <- rate_start * ifelse(is_double, 5, 1)
        if (is_double && "rate_slow" %in% names(starts)) {
          starts$rate_slow <- rate_start * 0.2
        }
        
        fit <- tryCatch(
          withCallingHandlers(
            nlsLM(
              formula_nls,
              data    = dat,
              start   = starts,
              lower   = lower_bounds,
              upper   = upper_bounds,
              control = nls.lm.control(maxiter = 500, ftol = 1e-6, ptol = 1e-6)
            ),
            warning = function(w) invokeRestart("muffleWarning")
          ),
          error = function(e) NULL
        )
        if (!is.null(fit) && any(!is.finite(coef(fit)))) fit <- NULL
        if (!is.null(fit)) break
      }
      fit
    }
    
    # --- attempt double exponential ---
    fit      <- NULL
    fit_used <- NULL
    
    if (use_double) {
      if (use_asym) {
        formula_nls  <- y_resp ~ asymptote +
          amplitude_fast * exp(-rate_fast * x) +
          amplitude_slow * exp(-rate_slow * x)
        starts       <- list(amplitude_fast = start_amp * 0.5, rate_fast = start_rate * 5,
                             amplitude_slow = start_amp * 0.5, rate_slow = start_rate * 0.2,
                             asymptote      = start_asym)
        lower_bounds <- c(amplitude_fast = 0, rate_fast = 1e-6,
                          amplitude_slow = 0, rate_slow = 1e-6, asymptote = 0)
        upper_bounds <- c(amplitude_fast = Inf, rate_fast = 20,
                          amplitude_slow = Inf, rate_slow = 20,  asymptote = Inf)
      } else {
        formula_nls  <- y_resp ~ amplitude_fast * exp(-rate_fast * x) +
          amplitude_slow * exp(-rate_slow * x)
        starts       <- list(amplitude_fast = start_amp * 0.5, rate_fast = start_rate * 5,
                             amplitude_slow = start_amp * 0.5, rate_slow = start_rate * 0.2)
        lower_bounds <- c(amplitude_fast = 0, rate_fast = 1e-6,
                          amplitude_slow = 0, rate_slow = 1e-6)
        upper_bounds <- c(amplitude_fast = Inf, rate_fast = 20,
                          amplitude_slow = Inf, rate_slow = 20)
      }
      fit      <- try_fit(formula_nls, starts, lower_bounds, upper_bounds, is_double = TRUE)
      fit_used <- if (!is.null(fit)) "double" else NULL
    }
    
    # --- fallback to single exponential if double failed or not requested ---
    if (is.null(fit)) {
      if (use_asym) {
        formula_nls  <- y_resp ~ asymptote + amplitude_fast * exp(-rate_fast * x)
        starts       <- list(amplitude_fast = start_amp, rate_fast = start_rate,
                             asymptote = start_asym)
        lower_bounds <- c(amplitude_fast = 0, rate_fast = 1e-6, asymptote = 0)
        upper_bounds <- c(amplitude_fast = Inf, rate_fast = 20,  asymptote = Inf)
      } else {
        formula_nls  <- y_resp ~ amplitude_fast * exp(-rate_fast * x)
        starts       <- list(amplitude_fast = start_amp, rate_fast = start_rate)
        lower_bounds <- c(amplitude_fast = 0, rate_fast = 1e-6)
        upper_bounds <- c(amplitude_fast = Inf, rate_fast = 20)
      }
      fit      <- try_fit(formula_nls, starts, lower_bounds, upper_bounds, is_double = FALSE)
      fit_used <- if (!is.null(fit)) "single_fallback" else NULL
    }
    
    # --- total failure ---
    if (is.null(fit)) {
      return(tibble(
        amplitude_fast = NA_real_, rate_fast = NA_real_,
        amplitude_slow = NA_real_, rate_slow = NA_real_,
        asymptote      = NA_real_,
        fit_type       = "failed",
        converged      = FALSE,
        AIC            = NA_real_,
        BIC            = NA_real_,
        RSS            = NA_real_,
        n_obs          = nrow(dat)
      ))
    }
    
    coefs    <- coef(fit)
    get_coef <- function(name) if (name %in% names(coefs)) unname(coefs[name]) else NA_real_
    
    tibble(
      amplitude_fast = get_coef("amplitude_fast"),
      rate_fast      = get_coef("rate_fast"),
      amplitude_slow = get_coef("amplitude_slow"),
      rate_slow      = get_coef("rate_slow"),
      asymptote      = get_coef("asymptote"),
      fit_type       = paste0(
        fit_used, "_",
        if (use_asym) "asym" else "no_asym", "_",
        if (use_log)  "logstim" else "rawstim"
      ),
      converged = TRUE,
      AIC       = AIC(fit),
      BIC       = BIC(fit),
      RSS       = sum(residuals(fit)^2),
      n_obs     = nrow(dat)
    )
  }
  # ==============================================================================
  # Main loop
  # ==============================================================================
  for (target in targets) {
    
    col_name <- target$col_name
    var_name <- target$var_name
    y_label  <- target$y_label
    
    message("\n==============================")
    message("Processing: ", var_name)
    message("==============================")
    
    save_fig_dir     <- file.path(base_dir, "figs",    "aggregated", var_name)
    save_results_dir <- file.path(base_dir, "results", "aggregated", var_name)
    dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)
    dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # --------------------------------------------------------------------------
    # Aggregate
    # --------------------------------------------------------------------------
    raw_summary_joint <- df_all %>%
      group_by(Training, Block, BlockRole, Genotype, stimulus, stimulus0, stimulus_log) %>%
      summarise(
        n_trials  = sum(!is.na(.data[[col_name]])),
        mean_resp = mean(.data[[col_name]], na.rm = TRUE),
        .groups   = "drop"
      )
    
    # --------------------------------------------------------------------------
    # Fit
    # --------------------------------------------------------------------------
    exp_params <- raw_summary_joint %>%
      filter(stimulus <= fit_max_stimulus) %>%
      group_by(Training, Block, BlockRole, Genotype) %>%
      group_modify(~ fit_exp_logstim_continuous(
        dat          = .x,
        response_col = "mean_resp",
        use_log      = use_log_stimulus,
        use_asym     = use_asymptote,
        use_double   = use_double_exp       # <-- new argument
      )) %>%
      ungroup()
    
    print(exp_params)
    readr::write_csv(
      exp_params,
      file.path(save_results_dir, paste0("exp_params_", var_name, ".csv"))
    )
    
    # --------------------------------------------------------------------------
    # Prediction grid
    # --------------------------------------------------------------------------
    new_data_exp <- raw_summary_joint %>%
      group_by(Training, Block, BlockRole, Genotype) %>%
      summarise(
        stim_min = min(stimulus),
        stim_max = if (extrapolate_curve) max(stimulus)
        else min(max(stimulus), fit_max_stimulus),
        .groups  = "drop"
      ) %>%
      left_join(exp_params, by = c("Training", "Block", "BlockRole", "Genotype")) %>%
      rowwise() %>%
      mutate(grid = list(tibble(stimulus = seq(stim_min, stim_max, length.out = 200)))) %>%
      unnest(grid) %>%
      ungroup() %>%
      mutate(
        stimulus0    = stimulus - 1,
        stimulus_log = log(stimulus),
        x_pred       = if (use_log_stimulus) stimulus_log else stimulus,
        fit = case_when(
          use_double_exp & use_asymptote  ~ asymptote +
            amplitude_fast * exp(-rate_fast * x_pred) +
            amplitude_slow * exp(-rate_slow * x_pred),
          use_double_exp & !use_asymptote ~ amplitude_fast * exp(-rate_fast * x_pred) +
            amplitude_slow * exp(-rate_slow * x_pred),
          !use_double_exp & use_asymptote  ~ asymptote +
            amplitude_fast * exp(-rate_fast * x_pred),
          !use_double_exp & !use_asymptote ~ amplitude_fast * exp(-rate_fast * x_pred)
        ),
        Training = factor(Training, levels = levels(df_all$Training)),
        Block    = factor(Block,    levels = levels(df_all$Block)),
        Genotype = factor(Genotype, levels = levels(df_all$Genotype))
      ) %>%
      select(-x_pred)
    
    # --------------------------------------------------------------------------
    # Plot helper
    # --------------------------------------------------------------------------
    make_plot <- function(training_label, legend_pos = "none") {
      ggplot(
        new_data_exp %>% filter(Training == training_label),
        aes(x = stimulus, color = Genotype)
      ) +
        facet_grid(Genotype ~ Block, scales = "free_x") +
        geom_point(
          data = raw_summary_joint %>% filter(Training == training_label),
          aes(x = stimulus, y = mean_resp, color = Genotype),
          inherit.aes = FALSE, alpha = 0.45, size = 1.0
        ) +
        geom_line(aes(y = fit), linewidth = 1.1, na.rm = TRUE) +
        theme_pubr(base_size = 12) +
        labs(
          x        = "Stimulus number within block",
          y        = y_label,
          title    = paste0("Aggregated ", var_name, ": ",
                            tools::toTitleCase(training_label), " training"),
          subtitle = paste0(
            "Points = group means; lines = ",
            if (use_asymptote) "asymptote + " else "",
            "exp decay on ",
            if (use_log_stimulus) "log(stimulus)" else "stimulus"
          )
        ) +
        theme(legend.position = legend_pos)
    }
    
    p_massed <- make_plot("massed", legend_pos = "none")
    p_spaced <- make_plot("spaced", legend_pos = "top")
    
    print(p_massed)
    print(p_spaced)
    ggsave(
      file.path(save_fig_dir, paste0(var_name, "_aggregated_massed.png")),
      p_massed , width = 10, height = 8, dpi = 300, bg = "white"
    )
    ggsave(
      file.path(save_fig_dir, paste0(var_name, "_aggregated_spaced.png")),
      p_spaced , width = 10, height = 8, dpi = 300, bg = "white"
    )
    # --------------------------------------------------------------------------
    # Per-animal raw traces
    # --------------------------------------------------------------------------
    raw_per_animal <- df_all %>%
      group_by(Training, Block, BlockRole, Genotype, animal, stimulus) %>%
      summarise(
        mean_resp = mean(.data[[col_name]], na.rm = TRUE),
        .groups   = "drop"
      )
    
    make_individual_plot <- function(training_label) {
      dat <- raw_per_animal %>% filter(Training == training_label)
      
      ggplot(dat, aes(x = stimulus, y = mean_resp, group = animal, color = Genotype)) +
        facet_grid(Genotype ~ Block, scales = "free_x") +
        coord_cartesian(ylim = c(0, 30)) +
        geom_line(alpha = 0.25, linewidth = 0.5) +
        theme_pubr(base_size = 11) +
        labs(
          x        = "Stimulus number within block",
          y        = y_label,
          title    = paste0("Individual traces — ", var_name, ": ",
                            tools::toTitleCase(training_label), " training"),
          subtitle = "One line per fish"
        ) +
        geom_line(
          data = new_data_exp %>% filter(Training == training_label),
          aes(x = stimulus, y = fit, group = Genotype),
          color = "black", linewidth = 1.0, inherit.aes = FALSE
        )+
        theme(legend.position = "none")
    }
    
    p_ind_massed <- make_individual_plot("massed")
    p_ind_spaced <- make_individual_plot("spaced")
    
    print(p_ind_massed)
    print(p_ind_spaced)    
    ggsave(
      file.path(save_fig_dir, paste0(var_name, "_ind_massed.png")),
      p_ind_massed , width = 10, height = 8, dpi = 300, bg = "white"
    )
    ggsave(
      file.path(save_fig_dir, paste0(var_name, "_ind_spaced.png")),
      p_ind_spaced , width = 10, height = 8, dpi = 300, bg = "white"
    )

    
    message("Done: ", var_name)
  }
  
  # ============================================================================
  # COMPARE POWER LAW TO EXPONENTIL FIT  
  # # --- run power law ---
  # use_log_stimulus <- TRUE
  # exp_params_log <- raw_summary_joint %>%
  #   filter(stimulus <= fit_max_stimulus) %>%
  #   group_by(Training, Block, BlockRole, Genotype) %>%
  #   group_modify(~ fit_exp_logstim_continuous(
  #     dat          = .x,
  #     response_col = "mean_resp",
  #     use_log      = use_log_stimulus,
  #     use_asym     = use_asymptote,
  #     use_double   = use_double_exp
  #   )) %>%
  #   ungroup()
  # 
  # # --- run raw exponential ---
  # use_log_stimulus <- FALSE
  # exp_params_raw <- raw_summary_joint %>%
  #   filter(stimulus <= fit_max_stimulus) %>%
  #   group_by(Training, Block, BlockRole, Genotype) %>%
  #   group_modify(~ fit_exp_logstim_continuous(
  #     dat          = .x,
  #     response_col = "mean_resp",
  #     use_log      = use_log_stimulus,
  #     use_asym     = use_asymptote,
  #     use_double   = use_double_exp
  #   )) %>%
  #   ungroup()
  # 
  # # --- compare ---
  # model_comparison <- bind_rows(exp_params_log, exp_params_raw) %>%
  #   filter(converged) %>%
  #   select(Training, Block, BlockRole, Genotype, fit_type, AIC, BIC, RSS, n_obs) %>%
  #   arrange(Training, Block, BlockRole, Genotype, AIC)
  # 
  # print(model_comparison)
  # 
  # readr::write_csv(
  #   model_comparison,
  #   file.path(save_results_dir, paste0("model_comparison_", var_name, ".csv"))
  # )
  

  
