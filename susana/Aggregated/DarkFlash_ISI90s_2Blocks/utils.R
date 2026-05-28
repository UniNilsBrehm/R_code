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

fit_aggregate_exp <- function(
    df,
    outcome = c("response_prob", "distance", "cumsum_distance", "delay"),
    genotype_col = "Genotype",
    block_col    = "Block",
    stimulus_col = "stimulus",
    move_col     = "move",
    peak_col     = "max_peak",
    cumsum_col   = "max_cumsum",
    delay_col    = "delay",
    genotype_order = NULL,
    colors         = NULL,
    y_limits       = NULL,
    y_break        = NULL
) {
  outcome <- match.arg(outcome)
  
  # ── 1. Prepare data ──────────────────────────────────────────────────────────
  geno_levels <- if (!is.null(genotype_order)) genotype_order else levels(factor(df[[genotype_col]]))
  
  df <- df %>%
    mutate(
      stimulus  = as.numeric(.data[[stimulus_col]]),
      stimulus0 = stimulus - 1,
      Genotype  = factor(.data[[genotype_col]], levels = geno_levels),
      Block     = factor(.data[[block_col]])
    )
  
  outcome_var <- switch(outcome,
                        response_prob   = move_col,
                        distance        = peak_col,
                        cumsum_distance = cumsum_col,
                        delay           = delay_col
  )
  
  df <- df %>%
    mutate(.outcome_raw = as.numeric(.data[[outcome_var]]))
  
  df_agg <- df %>%
    group_by(Genotype, Block, stimulus, stimulus0) %>%
    summarise(
      n_animals = n(),
      y_mean    = mean(.outcome_raw, na.rm = TRUE),
      .groups   = "drop"
    )
  
  y_label <- switch(outcome,
                    response_prob   = "Response probability",
                    distance        = "Distance moved (peak)",
                    cumsum_distance = "Cumulative distance moved",
                    delay           = "Delay (ordinal)"
  )
  
  # ── 2. Define model ──────────────────────────────────────────────────────────
  get_bounds <- function(dat, outcome) {
    y      <- dat$y_mean
    y_min  <- min(y, na.rm = TRUE)
    y_max  <- max(y, na.rm = TRUE)
    y_rng  <- max(y_max - y_min, 1e-6)
    
    switch(outcome,
           response_prob = list(
             start = list(Asym = max(0.01, y_min), R0 = min(0.99, y_max), k = 0.1),
             lower = c(Asym = 0.001, R0 = 0.001, k = 0.0001),
             upper = c(Asym = 0.999, R0 = 0.999, k = 10)
           ),
           delay = list(
             start = list(Asym = max(0, y_min), R0 = min(4, y_max), k = 0.1),
             lower = c(Asym = 0, R0 = 0, k = 0.0001),
             upper = c(Asym = 4, R0 = 4, k = 10)
           ),
           {
             list(
               start = list(Asym = max(0, y_min), R0 = max(y_min + 1e-3, y_max), k = 0.1),
               lower = c(Asym = 0,         R0 = 0,         k = 0.0001),
               upper = c(Asym = y_max * 5, R0 = y_max * 5, k = 10)
             )
           }
    )
  }
  
  fit_exp <- function(dat, outcome) {
    bounds <- get_bounds(dat, outcome)
    tryCatch(
      nls(
        y_mean ~ Asym + (R0 - Asym) * exp(-k * stimulus0),
        data      = dat,
        weights   = n_animals,
        start     = bounds$start,
        algorithm = "port",
        lower     = bounds$lower,
        upper     = bounds$upper,
        control   = nls.control(maxiter = 500, warnOnly = TRUE)
      ),
      error = function(e) {
        message("NLS failed for ", unique(dat$Genotype), " ", unique(dat$Block), ": ", e$message)
        NULL
      }
    )
  }
  
  # ── 3. Fit per Genotype × Block ──────────────────────────────────────────────
  fit_keys  <- df_agg %>% distinct(Genotype, Block) %>% arrange(Genotype, Block)
  split_dat <- df_agg %>% group_by(Genotype, Block) %>% group_split()
  fits      <- map(split_dat, fit_exp, outcome = outcome)
  
  # ── 4. Extract parameters ────────────────────────────────────────────────────
  params <- map2_dfr(fits, seq_along(fits), function(mod, i) {
    base <- tibble(Genotype = fit_keys$Genotype[i], Block = fit_keys$Block[i])
    if (is.null(mod)) {
      bind_cols(base, tibble(Asym = NA_real_, R0 = NA_real_, k = NA_real_))
    } else {
      cc <- coef(mod)
      bind_cols(base, tibble(
        Asym = unname(cc["Asym"]),
        R0   = unname(cc["R0"]),
        k    = unname(cc["k"])
      ))
    }
  }) %>%
    mutate(half_life_stimuli = log(2) / k)
  
  # ── 5. Prediction curves ─────────────────────────────────────────────────────
  pred <- df_agg %>%
    group_by(Genotype, Block) %>%
    summarise(stim_min = min(stimulus), stim_max = max(stimulus), .groups = "drop") %>%
    group_by(Genotype, Block) %>%
    reframe(stimulus = seq(stim_min, stim_max, length.out = 100)) %>%
    mutate(stimulus0 = stimulus - 1) %>%
    left_join(params, by = c("Genotype", "Block")) %>%
    mutate(fit = Asym + (R0 - Asym) * exp(-k * stimulus0))
  
  # ── 6. Color scales ──────────────────────────────────────────────────────────
  color_scales <- if (!is.null(colors)) {
    list(
      scale_color_manual(values = colors),
      scale_fill_manual(values = colors)
    )
  } else {
    list()
  }
  
  # ── 7. Plot ──────────────────────────────────────────────────────────────────
  p <- ggplot(pred, aes(x = stimulus, color = Genotype, fill = Genotype)) +
    facet_grid(Block ~ Genotype, scales = "fixed") +
    geom_point(
      data        = df_agg,
      aes(x = stimulus, y = y_mean, color = Genotype),
      alpha       = 0.45,
      size        = 1,
      inherit.aes = FALSE
    ) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    color_scales +
    scale_x_continuous(n.breaks = 4) +
    coord_cartesian(ylim = y_limits) +
    scale_y_continuous(
      breaks = if (!is.null(y_limits) && !is.null(y_break))
        seq(y_limits[1], y_limits[2], by = y_break)
      else
        waiver()
    ) +
    theme_pubr(base_size = 14) +
    labs(
      x     = "Stimulus number within block",
      y     = y_label,
      color = "Genotype",
      fill  = "Genotype",
      title = paste("Aggregate exponential fit —", y_label)
    ) +
    theme(
      legend.position = "top",
      panel.spacing   = unit(1.2, "lines"),
      panel.grid.major = element_line(color = "grey95", linewidth = 0.5)
    )
  
  list(plot = p, params = params)
}


# ==============================================================================
# BOOTSTRAP CIs
bootstrap_exp_ci <- function(
    df,
    outcome      = c("response_prob", "distance", "cumsum_distance", "delay"),
    n_boot       = 1000,
    ci_level     = 0.95,
    genotype_col  = "Genotype",
    block_col     = "Block",
    stimulus_col  = "stimulus",
    move_col      = "move",
    peak_col      = "max_peak",
    cumsum_col    = "max_cumsum",
    delay_col     = "delay",
    animal_col    = NULL,   # if NULL, built from Video × Well
    seed         = 42
) {
  outcome <- match.arg(outcome)
  set.seed(seed)
  
  alpha <- 1 - ci_level
  
  outcome_var <- switch(outcome,
                        response_prob   = move_col,
                        distance        = peak_col,
                        cumsum_distance = cumsum_col,
                        delay           = delay_col
  )
  
  y_label <- switch(outcome,
                    response_prob   = "Response probability",
                    distance        = "Distance moved (peak)",
                    cumsum_distance = "Cumulative distance moved",
                    delay           = "Delay (ordinal)"
  )
  
  # ── 1. Prep: create stimulus0 and animal ID ──────────────────────────────────
  df <- df %>%
    mutate(
      stimulus  = as.numeric(.data[[stimulus_col]]),
      stimulus0 = stimulus - 1,
      Genotype  = factor(.data[[genotype_col]]),
      Block     = factor(.data[[block_col]]),
      .y        = as.numeric(.data[[outcome_var]]),
      .animal   = if (!is.null(animal_col)) {
        as.character(.data[[animal_col]])
      } else {
        as.character(interaction(Video, Well, drop = TRUE))
      }
    )
  
  # ── 2. Bounds helper (same as fit_aggregate_exp) ─────────────────────────────
  get_bounds <- function(dat, outcome) {
    y     <- dat$y_mean
    y_max <- max(y, na.rm = TRUE)
    switch(outcome,
           response_prob = list(
             start = list(Asym = max(0.01, min(y)), R0 = min(0.99, max(y)), k = 0.1),
             lower = c(Asym = 0.001, R0 = 0.001, k = 0.0001),
             upper = c(Asym = 0.999, R0 = 0.999, k = 10)
           ),
           delay = list(
             start = list(Asym = max(0, min(y)), R0 = min(4, max(y)), k = 0.1),
             lower = c(Asym = 0, R0 = 0,      k = 0.0001),
             upper = c(Asym = 4, R0 = 4,      k = 10)
           ),
           list(
             start = list(Asym = max(0, min(y)), R0 = max(min(y) + 1e-3, max(y)), k = 0.1),
             lower = c(Asym = 0,         R0 = 0,         k = 0.0001),
             upper = c(Asym = y_max * 5, R0 = y_max * 5, k = 10)
           )
    )
  }
  
  fit_nls <- function(dat, outcome) {
    bounds <- get_bounds(dat, outcome)
    tryCatch(
      nls(
        y_mean ~ Asym + (R0 - Asym) * exp(-k * stimulus0),
        data      = dat,
        weights   = n_animals,
        start     = bounds$start,
        algorithm = "port",
        lower     = bounds$lower,
        upper     = bounds$upper,
        control   = nls.control(maxiter = 500, warnOnly = TRUE)
      ),
      error = function(e) NULL
    )
  }
  
  # ── 3. Bootstrap per Genotype × Block ────────────────────────────────────────
  groups <- df %>% distinct(Genotype, Block) %>% arrange(Genotype, Block)
  
  boot_ci <- pmap_dfr(groups, function(Genotype, Block) {
    dat_grp <- df %>% filter(Genotype == !!Genotype, Block == !!Block)
    animals  <- unique(dat_grp$.animal)
    
    stim_vals  <- sort(unique(dat_grp$stimulus))
    stim_grid  <- tibble(
      stimulus  = seq(min(stim_vals), max(stim_vals), length.out = 100),
      stimulus0 = seq(min(stim_vals) - 1, max(stim_vals) - 1, length.out = 100)
    )
    
    boot_curves <- map(seq_len(n_boot), function(b) {
      sampled  <- sample(animals, length(animals), replace = TRUE)
      
      dat_boot <- map_dfr(sampled, function(a) {
        dat_grp %>% filter(.animal == a)
      }) %>%
        group_by(stimulus, stimulus0) %>%
        summarise(
          n_animals = n(),
          y_mean    = mean(.y, na.rm = TRUE),
          .groups   = "drop"
        )
      
      mod <- fit_nls(dat_boot, outcome)
      if (is.null(mod)) return(NULL)
      
      cc <- coef(mod)
      tibble(
        stimulus = stim_grid$stimulus,
        fit      = cc["Asym"] + (cc["R0"] - cc["Asym"]) * exp(-cc["k"] * stim_grid$stimulus0)
      )
    }) %>%
      bind_rows()
    
    boot_curves %>%
      group_by(stimulus) %>%
      summarise(
        ci_lo = quantile(fit, alpha / 2,     na.rm = TRUE),
        ci_hi = quantile(fit, 1 - alpha / 2, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Genotype = Genotype, Block = Block)
  })
  
  # ── 4. Point estimates via fit_aggregate_exp ──────────────────────────────────
  res <- fit_aggregate_exp(
    df, outcome = outcome,
    genotype_col = genotype_col, block_col = block_col,
    stimulus_col = stimulus_col, move_col = move_col,
    peak_col = peak_col, cumsum_col = cumsum_col, delay_col = delay_col
  )
  
  # Extract the smooth pred data from the plot layers
  pred <- res$plot$layers[[2]]$data  # geom_line layer
  
  pred_ci <- pred %>%
    left_join(boot_ci, by = c("Genotype", "Block", "stimulus"))
  
  # ── 5. Observed aggregated means ─────────────────────────────────────────────
  df_agg <- df %>%
    group_by(Genotype, Block, stimulus, stimulus0) %>%
    summarise(n_animals = n(), y_mean = mean(.y, na.rm = TRUE), .groups = "drop")
  
  # ── 6. Plot ───────────────────────────────────────────────────────────────────
  p <- ggplot(pred_ci, aes(x = stimulus, color = Genotype, fill = Genotype)) +
    facet_grid(Block ~ Genotype, scales = "fixed") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, color = NA) +
    geom_point(
      data        = df_agg,
      aes(x = stimulus, y = y_mean, color = Genotype),
      alpha       = 0.45, size = 1,
      inherit.aes = FALSE
    ) +
    geom_line(aes(y = fit), linewidth = 1.3) +
    theme_pubr(base_size = 14) +
    labs(
      x     = "Stimulus number within block",
      y     = y_label,
      color = "Genotype",
      fill  = "Genotype",
      title = paste0("Aggregate exponential fit — ", y_label,
                     " (", ci_level * 100, "% bootstrap CI, n = ", n_boot, ")")
    ) +
    theme(legend.position = "top", panel.spacing = unit(1.2, "lines"))
  
  list(plot = p, params = res$params, boot_ci = boot_ci)
}