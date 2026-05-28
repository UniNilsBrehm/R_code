# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
library(readr)        # Reading CSV files
library(ggplot2)      # Visualization
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(ggpubr)       # Publication-ready plots
library(purrr)
library(coin)         # permutation tests
library(rstatix)      # Kruskal-Wallis + Dunn's post-hoc
library(patchwork)    # plot composition

source("C:/UniFreiburg/Code/R_code/susana/Aggregated//DarkFlash_ISI90s_2Blocks/utils.R")
base_dir <- "D:/WorkingData/Susana/Aggregated/DarkFlash_ISI90s_2Blocks"

file_dir <- file.path(
  base_dir,
  "data_files",
  "SPZ_ISI60_removed_non_responders_2stimuli.csv"
)

# ==============================================================================
# Load Data
# ==============================================================================
message("Loading data...")
res <- load_data_darkflash_60s(file_dir, move_th = 0, take_peak = 0)

df_final     <- res$df_final
df_final_sub <- res$df_final_sub
# Replace exact-zero delays with epsilon to avoid safe_ratio NA (zero denominator)
df_final_sub$delay[df_final_sub$delay == 0] <- 0.0001
df_final_sub$delay_ord <- factor(df_final_sub$delay, ordered = TRUE)

# Number of animals per genotype
n_per_genotype <- df_final %>%
  distinct(Video, Well, Genotype) %>%
  count(Genotype)
print(n_per_genotype)

save_fig_dir     <- file.path(base_dir, "figs")
save_results_dir <- file.path(base_dir, "results")
models_dir       <- file.path(base_dir, "models")

dir.create(save_fig_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(save_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir,       recursive = TRUE, showWarnings = FALSE)

geno_colors <- c(
  "ABTL"         = "#440154",
  "tyr"          = "#3B528B",
  "th2, tyr"     = "#21908C",
  "th, tyr"      = "#5DC863",
  "th, th2, tyr" = "#FDE725"
)
geno_order <- names(geno_colors)   # order for x-axis: derived from color map

# ==============================================================================
# Window settings — adjust here
# ==============================================================================
FIRST_WIN_END  <- 3    # stimuli 1:FIRST_WIN_END   = "first" window
LAST_WIN_START <- 57   # stimuli LAST_WIN_START:59  = "last"  window

# Total stimuli per block (used to compute center window)
N_STIM <- 59
CENTER_STIM <- round(N_STIM / 2)   # = 30
CENTER_WIN_HALF <- 1                # center window = (CENTER_STIM - 1) : (CENTER_STIM + 1)
CENTER_WIN_START <- CENTER_STIM - CENTER_WIN_HALF   # 29
CENTER_WIN_END   <- CENTER_STIM + CENTER_WIN_HALF   # 31

# ==============================================================================
# 1. Per-animal block averages
# ==============================================================================

# --- Response probability (from df_final — all trials) ---
df_block_prob <- df_final %>%
  mutate(stimulus = as.numeric(stimulus)) %>%
  group_by(Video, Well, Genotype, Block) %>%
  summarise(
    mean_prob    = mean(as.numeric(move), na.rm = TRUE),
    first_prob   = mean(as.numeric(move)[stimulus <= FIRST_WIN_END],  na.rm = TRUE),
    center_prob  = mean(as.numeric(move)[stimulus >= CENTER_WIN_START &
                                           stimulus <= CENTER_WIN_END], na.rm = TRUE),
    last_prob    = mean(as.numeric(move)[stimulus >= LAST_WIN_START],  na.rm = TRUE),
    .groups      = "drop"
  )

# --- Peak distance + cumulative distance + delay (from df_final_sub — responders only) ---
df_block_sub <- df_final_sub %>%
  mutate(stimulus = as.numeric(stimulus)) %>%
  group_by(Video, Well, Genotype, Block) %>%
  summarise(
    # Peak distance
    mean_peak    = mean(max_peak, na.rm = TRUE),
    first_peak   = mean(max_peak[stimulus <= FIRST_WIN_END],  na.rm = TRUE),
    center_peak  = mean(max_peak[stimulus >= CENTER_WIN_START &
                                   stimulus <= CENTER_WIN_END], na.rm = TRUE),
    last_peak    = mean(max_peak[stimulus >= LAST_WIN_START],  na.rm = TRUE),
    # Cumulative distance
    mean_cumsum  = mean(max_cumsum, na.rm = TRUE),
    first_cumsum = mean(max_cumsum[stimulus <= FIRST_WIN_END],  na.rm = TRUE),
    center_cumsum= mean(max_cumsum[stimulus >= CENTER_WIN_START &
                                     stimulus <= CENTER_WIN_END], na.rm = TRUE),
    last_cumsum  = mean(max_cumsum[stimulus >= LAST_WIN_START],  na.rm = TRUE),
    # Delay
    mean_delay   = mean(delay, na.rm = TRUE),
    first_delay  = mean(delay[stimulus <= FIRST_WIN_END],  na.rm = TRUE),
    center_delay = mean(delay[stimulus >= CENTER_WIN_START &
                                stimulus <= CENTER_WIN_END], na.rm = TRUE),
    last_delay   = mean(delay[stimulus >= LAST_WIN_START],  na.rm = TRUE),
    .groups      = "drop"
  )

# ==============================================================================
# 2. Helper functions
# ==============================================================================

safe_ratio <- function(a, b) {
  case_when(
    is.na(a) | is.na(b)   ~ NA_real_,
    is.nan(a) | is.nan(b) ~ NA_real_,
    (a + b) == 0          ~ NA_real_,
    TRUE                  ~ (a - b) / (a + b)
  )
}

# HI_between = (B1 - B2) / (B1 + B2) for all outcomes
# Positive = B1 > B2 (higher in Block1 than Block2)
hi_between <- function(b1, b2) {
  safe_ratio(b1, b2)
}

report_nas <- function(df, label) {
  df %>%
    summarise(
      outcome = label,
      across(where(is.numeric), ~ sum(is.na(.x))),
      n_total = n()
    )
}

# ==============================================================================
# 3. Compute indices — Response probability
# ==============================================================================

df_indices_prob <- df_block_prob %>%
  pivot_wider(
    names_from  = Block,
    values_from = c(mean_prob, first_prob, center_prob, last_prob)
  ) %>%
  mutate(
    # --- Between-block HI (4 window variants) ---
    HI_between_mean   = hi_between(mean_prob_Block1,   mean_prob_Block2),
    HI_between_first  = hi_between(first_prob_Block1,  first_prob_Block2),
    HI_between_center = hi_between(center_prob_Block1, center_prob_Block2),
    HI_between_last   = hi_between(last_prob_Block1,   last_prob_Block2),
    # --- Within-block SI ---
    SI_block1         = safe_ratio(first_prob_Block1,  last_prob_Block1),
    SI_block2         = safe_ratio(first_prob_Block2,  last_prob_Block2),
    # --- Absolute drop ---
    abs_drop_b1       = first_prob_Block1 - last_prob_Block1,
    abs_drop_b2       = first_prob_Block2 - last_prob_Block2,
    # --- Proportional drop ---
    prop_drop_b1      = case_when(
      is.na(first_prob_Block1) | first_prob_Block1 == 0 ~ NA_real_,
      TRUE ~ (first_prob_Block1 - last_prob_Block1) / first_prob_Block1
    ),
    prop_drop_b2      = case_when(
      is.na(first_prob_Block2) | first_prob_Block2 == 0 ~ NA_real_,
      TRUE ~ (first_prob_Block2 - last_prob_Block2) / first_prob_Block2
    ),
    # --- Asymptote & initial reactivity ---
    asymptote_b1      = last_prob_Block1,
    asymptote_b2      = last_prob_Block2,
    initial_b1        = first_prob_Block1,
    initial_b2        = first_prob_Block2,
    outcome           = "Response probability"
  ) %>%
  select(Video, Well, Genotype, outcome,
         HI_between_mean, HI_between_first, HI_between_center, HI_between_last,
         SI_block1, SI_block2,
         abs_drop_b1, abs_drop_b2,
         prop_drop_b1, prop_drop_b2,
         asymptote_b1, asymptote_b2,
         initial_b1, initial_b2)

# ==============================================================================
# 4. Compute indices — Peak distance
# ==============================================================================

df_indices_peak <- df_block_sub %>%
  pivot_wider(
    names_from  = Block,
    values_from = c(mean_peak,   first_peak,   center_peak,   last_peak,
                    mean_cumsum, first_cumsum, center_cumsum, last_cumsum,
                    mean_delay,  first_delay,  center_delay,  last_delay)
  ) %>%
  mutate(
    HI_between_mean   = hi_between(mean_peak_Block1,   mean_peak_Block2),
    HI_between_first  = hi_between(first_peak_Block1,  first_peak_Block2),
    HI_between_center = hi_between(center_peak_Block1, center_peak_Block2),
    HI_between_last   = hi_between(last_peak_Block1,   last_peak_Block2),
    SI_block1         = safe_ratio(first_peak_Block1,  last_peak_Block1),
    SI_block2         = safe_ratio(first_peak_Block2,  last_peak_Block2),
    abs_drop_b1       = first_peak_Block1 - last_peak_Block1,
    abs_drop_b2       = first_peak_Block2 - last_peak_Block2,
    prop_drop_b1      = case_when(
      is.na(first_peak_Block1) | first_peak_Block1 == 0 ~ NA_real_,
      TRUE ~ (first_peak_Block1 - last_peak_Block1) / first_peak_Block1
    ),
    prop_drop_b2      = case_when(
      is.na(first_peak_Block2) | first_peak_Block2 == 0 ~ NA_real_,
      TRUE ~ (first_peak_Block2 - last_peak_Block2) / first_peak_Block2
    ),
    asymptote_b1      = last_peak_Block1,
    asymptote_b2      = last_peak_Block2,
    initial_b1        = first_peak_Block1,
    initial_b2        = first_peak_Block2,
    outcome           = "Peak distance"
  ) %>%
  select(Video, Well, Genotype, outcome,
         HI_between_mean, HI_between_first, HI_between_center, HI_between_last,
         SI_block1, SI_block2,
         abs_drop_b1, abs_drop_b2,
         prop_drop_b1, prop_drop_b2,
         asymptote_b1, asymptote_b2,
         initial_b1, initial_b2)

# ==============================================================================
# 5. Compute indices — Cumulative distance
# ==============================================================================

df_indices_cumsum <- df_block_sub %>%
  pivot_wider(
    names_from  = Block,
    values_from = c(mean_peak,   first_peak,   center_peak,   last_peak,
                    mean_cumsum, first_cumsum, center_cumsum, last_cumsum,
                    mean_delay,  first_delay,  center_delay,  last_delay)
  ) %>%
  mutate(
    HI_between_mean   = hi_between(mean_cumsum_Block1,   mean_cumsum_Block2),
    HI_between_first  = hi_between(first_cumsum_Block1,  first_cumsum_Block2),
    HI_between_center = hi_between(center_cumsum_Block1, center_cumsum_Block2),
    HI_between_last   = hi_between(last_cumsum_Block1,   last_cumsum_Block2),
    SI_block1         = safe_ratio(first_cumsum_Block1,  last_cumsum_Block1),
    SI_block2         = safe_ratio(first_cumsum_Block2,  last_cumsum_Block2),
    abs_drop_b1       = first_cumsum_Block1 - last_cumsum_Block1,
    abs_drop_b2       = first_cumsum_Block2 - last_cumsum_Block2,
    prop_drop_b1      = case_when(
      is.na(first_cumsum_Block1) | first_cumsum_Block1 == 0 ~ NA_real_,
      TRUE ~ (first_cumsum_Block1 - last_cumsum_Block1) / first_cumsum_Block1
    ),
    prop_drop_b2      = case_when(
      is.na(first_cumsum_Block2) | first_cumsum_Block2 == 0 ~ NA_real_,
      TRUE ~ (first_cumsum_Block2 - last_cumsum_Block2) / first_cumsum_Block2
    ),
    asymptote_b1      = last_cumsum_Block1,
    asymptote_b2      = last_cumsum_Block2,
    initial_b1        = first_cumsum_Block1,
    initial_b2        = first_cumsum_Block2,
    outcome           = "Cumulative distance"
  ) %>%
  select(Video, Well, Genotype, outcome,
         HI_between_mean, HI_between_first, HI_between_center, HI_between_last,
         SI_block1, SI_block2,
         abs_drop_b1, abs_drop_b2,
         prop_drop_b1, prop_drop_b2,
         asymptote_b1, asymptote_b2,
         initial_b1, initial_b2)

# ==============================================================================
# 6. Compute indices — Delay
# ==============================================================================

df_indices_delay <- df_block_sub %>%
  pivot_wider(
    names_from  = Block,
    values_from = c(mean_peak,   first_peak,   center_peak,   last_peak,
                    mean_cumsum, first_cumsum, center_cumsum, last_cumsum,
                    mean_delay,  first_delay,  center_delay,  last_delay)
  ) %>%
  mutate(
    # Delay: same as all outcomes — HI = (B1 - B2) / (B1 + B2)
    HI_between_mean   = hi_between(mean_delay_Block1,   mean_delay_Block2),
    HI_between_first  = hi_between(first_delay_Block1,  first_delay_Block2),
    HI_between_center = hi_between(center_delay_Block1, center_delay_Block2),
    HI_between_last   = hi_between(last_delay_Block1,   last_delay_Block2),
    # Delay SI: same as all outcomes — (first - last) / (first + last)
    SI_block1         = safe_ratio(first_delay_Block1,  last_delay_Block1),
    SI_block2         = safe_ratio(first_delay_Block2,  last_delay_Block2),
    # Absolute drop: same as all outcomes — first - last
    abs_drop_b1       = first_delay_Block1 - last_delay_Block1,
    abs_drop_b2       = first_delay_Block2 - last_delay_Block2,
    prop_drop_b1      = case_when(
      is.na(first_delay_Block1) | first_delay_Block1 == 0 ~ NA_real_,
      TRUE ~ (first_delay_Block1 - last_delay_Block1) / first_delay_Block1
    ),
    prop_drop_b2      = case_when(
      is.na(first_delay_Block2) | first_delay_Block2 == 0 ~ NA_real_,
      TRUE ~ (first_delay_Block2 - last_delay_Block2) / first_delay_Block2
    ),
    asymptote_b1      = last_delay_Block1,
    asymptote_b2      = last_delay_Block2,
    initial_b1        = first_delay_Block1,
    initial_b2        = first_delay_Block2,
    outcome           = "Delay"
  ) %>%
  select(Video, Well, Genotype, outcome,
         HI_between_mean, HI_between_first, HI_between_center, HI_between_last,
         SI_block1, SI_block2,
         abs_drop_b1, abs_drop_b2,
         prop_drop_b1, prop_drop_b2,
         asymptote_b1, asymptote_b2,
         initial_b1, initial_b2)

# --- NA diagnostics ---
na_report <- bind_rows(
  df_indices_prob   %>% summarise(outcome = "Response probability",
                                  across(where(is.numeric), ~ sum(is.na(.x))), n = n()),
  df_indices_peak   %>% summarise(outcome = "Peak distance",
                                  across(where(is.numeric), ~ sum(is.na(.x))), n = n()),
  df_indices_cumsum %>% summarise(outcome = "Cumulative distance",
                                  across(where(is.numeric), ~ sum(is.na(.x))), n = n()),
  df_indices_delay  %>% summarise(outcome = "Delay",
                                  across(where(is.numeric), ~ sum(is.na(.x))), n = n())
)
print(na_report)

# ==============================================================================
# 7. Statistical tests
# ==============================================================================

run_nonparam_tests <- function(df, index_var, label = "") {
  df_clean <- df %>%
    filter(!is.na(.data[[index_var]]), !is.nan(.data[[index_var]]))
  
  n_removed <- nrow(df) - nrow(df_clean)
  if (n_removed > 0) message(label, ": removed ", n_removed, " NA animals from test")
  
  # Drop unused factor levels — empty levels crash kruskal_test internally
  df_clean <- df_clean %>% mutate(Genotype = droplevels(factor(Genotype)))
  
  # Guard: need at least 2 genotypes with >0 observations to run any test
  group_counts <- table(df_clean$Genotype)
  n_groups     <- sum(group_counts > 0)
  if (n_groups < 2) {
    message(label, ": fewer than 2 groups after NA removal — skipping tests")
    empty_dunn <- tibble(group1 = character(), group2 = character(),
                         p = numeric(), p.adj = numeric(), p.adj.signif = character())
    return(list(kruskal = NULL, dunn = empty_dunn, data = df_clean))
  }
  
  # Use base R kruskal.test to avoid rstatix add_column bug on edge cases
  kw_base <- kruskal.test(
    as.formula(paste(index_var, "~ Genotype")), data = df_clean
  )
  kw <- tibble(
    .y.        = index_var,
    n          = nrow(df_clean),
    statistic  = kw_base$statistic,
    df         = kw_base$parameter,
    p          = kw_base$p.value,
    method     = "Kruskal-Wallis"
  )
  
  # Filter to ABTL-only contrasts before p-adjustment:
  # 4 comparisons instead of 10 => more power after BH correction
  # Use tryCatch to handle rstatix edge-case crashes
  dunn_raw <- tryCatch({
    df_clean %>%
      dunn_test(as.formula(paste(index_var, "~ Genotype")), p.adjust.method = "none") %>%
      filter(group1 == "ABTL" | group2 == "ABTL")
  }, error = function(e) {
    message(label, ": dunn_test failed (", conditionMessage(e), ") — skipping post-hoc")
    tibble(group1 = character(), group2 = character(), p = numeric())
  })
  
  if (nrow(dunn_raw) == 0) {
    message(label, ": no ABTL contrasts possible after filtering")
    dunn <- tibble(group1 = character(), group2 = character(),
                   p = numeric(), p.adj = numeric(), p.adj.signif = character())
  } else {
    dunn <- dunn_raw %>%
      mutate(p.adj = p.adjust(p, method = "BH")) %>%
      add_significance("p.adj")
  }
  
  list(kruskal = kw, dunn = dunn, data = df_clean)
}

# Helper: run a given index across all four outcome dataframes
run_all_outcomes <- function(index_var, label_suffix) {
  list(
    prob   = run_nonparam_tests(df_indices_prob,   index_var, paste("Prob",   label_suffix)),
    peak   = run_nonparam_tests(df_indices_peak,   index_var, paste("Peak",   label_suffix)),
    cumsum = run_nonparam_tests(df_indices_cumsum, index_var, paste("Cumsum", label_suffix)),
    delay  = run_nonparam_tests(df_indices_delay,  index_var, paste("Delay",  label_suffix))
  )
}

# --- HI_between (four window variants) ---
stats_hi_mean   <- run_all_outcomes("HI_between_mean",   "HI_between_mean")
stats_hi_first  <- run_all_outcomes("HI_between_first",  "HI_between_first")
stats_hi_center <- run_all_outcomes("HI_between_center", "HI_between_center")
stats_hi_last   <- run_all_outcomes("HI_between_last",   "HI_between_last")

# --- Within-block SI ---
stats_si_b1  <- run_all_outcomes("SI_block1",    "SI_block1")
stats_si_b2  <- run_all_outcomes("SI_block2",    "SI_block2")

# --- Additional indices ---
stats_abs_b1  <- run_all_outcomes("abs_drop_b1",  "abs_drop_b1")
stats_abs_b2  <- run_all_outcomes("abs_drop_b2",  "abs_drop_b2")
stats_prop_b1 <- run_all_outcomes("prop_drop_b1", "prop_drop_b1")
stats_prop_b2 <- run_all_outcomes("prop_drop_b2", "prop_drop_b2")
stats_asym_b1 <- run_all_outcomes("asymptote_b1", "asymptote_b1")
stats_asym_b2 <- run_all_outcomes("asymptote_b2", "asymptote_b2")
stats_init_b1 <- run_all_outcomes("initial_b1",   "initial_b1")
stats_init_b2 <- run_all_outcomes("initial_b2",   "initial_b2")

# ==============================================================================
# 8. Save results
# ==============================================================================

save_stats <- function(stats_list, filename) {
  sink(file.path(save_results_dir, filename))
  if (is.null(stats_list$kruskal)) {
    cat("Kruskal-Wallis: skipped (fewer than 2 groups)\n")
  } else {
    cat("Kruskal-Wallis:\n");  print(stats_list$kruskal)
  }
  cat("\nDunn post-hoc:\n"); print(stats_list$dunn)
  sink()
}

save_all <- function(stats_group, prefix) {
  save_stats(stats_group$prob,   paste0(prefix, "_prob.txt"))
  save_stats(stats_group$peak,   paste0(prefix, "_peak.txt"))
  save_stats(stats_group$cumsum, paste0(prefix, "_cumsum.txt"))
  save_stats(stats_group$delay,  paste0(prefix, "_delay.txt"))
}

save_all(stats_hi_mean,   "nonparam_HI_between_mean")
save_all(stats_hi_first,  "nonparam_HI_between_first")
save_all(stats_hi_center, "nonparam_HI_between_center")
save_all(stats_hi_last,   "nonparam_HI_between_last")
save_all(stats_si_b1,     "nonparam_SI_block1")
save_all(stats_si_b2,     "nonparam_SI_block2")
save_all(stats_abs_b1,    "nonparam_abs_drop_b1")
save_all(stats_abs_b2,    "nonparam_abs_drop_b2")
save_all(stats_prop_b1,   "nonparam_prop_drop_b1")
save_all(stats_prop_b2,   "nonparam_prop_drop_b2")
save_all(stats_asym_b1,   "nonparam_asymptote_b1")
save_all(stats_asym_b2,   "nonparam_asymptote_b2")
save_all(stats_init_b1,   "nonparam_initial_b1")
save_all(stats_init_b2,   "nonparam_initial_b2")

# ==============================================================================
# 9. Plot helpers
# ==============================================================================

# Label helpers — avoids "59-59" when last window is a single stimulus
last_win_label  <- if (LAST_WIN_START == N_STIM) {
  paste0("stim ", N_STIM, " only")
} else {
  paste0("stim ", LAST_WIN_START, "\u2013", N_STIM)
}
first_win_label  <- paste0("stim 1\u2013", FIRST_WIN_END)
center_win_label <- paste0("stim ", CENTER_WIN_START, "\u2013", CENTER_WIN_END)

get_dunn_annot <- function(stats_obj, df, index_var) {
  annot <- stats_obj$dunn %>%
    filter(p.adj < 0.05) %>%
    select(group1, group2, p.adj, p.adj.signif)
  
  if (nrow(annot) == 0) return(annot)
  
  annot <- annot %>%
    add_xy_position(
      data          = df %>%
        filter(!is.na(.data[[index_var]])) %>%
        mutate(Genotype = factor(Genotype, levels = geno_order)),
      formula       = as.formula(paste(index_var, "~ Genotype")),
      fun           = "max",
      step.increase = 0.08
    )
  
  # Drop rows where y.position could not be computed (all-NA group, etc.)
  if ("y.position" %in% names(annot)) {
    annot <- annot %>% filter(!is.na(y.position))
  }
  annot
}

plot_index <- function(stats_obj, index_var, title,
                       y_label, expand_top = 0.35) {
  df_clean <- stats_obj$data %>%
    mutate(Genotype = factor(Genotype, levels = geno_order))
  annot    <- get_dunn_annot(stats_obj, df_clean, index_var)
  
  # Determine y limits: extend ceiling to accommodate significance brackets
  y_vals  <- df_clean[[index_var]]
  y_min   <- min(y_vals, na.rm = TRUE)
  y_max   <- max(y_vals, na.rm = TRUE)
  y_range <- if (is.finite(y_max - y_min)) y_max - y_min else 1
  
  has_brackets <- nrow(annot) > 0 &&
    "y.position" %in% names(annot) &&
    any(is.finite(annot$y.position))
  
  if (has_brackets) {
    bracket_top <- max(annot$y.position[is.finite(annot$y.position)])
    y_upper     <- max(y_max, bracket_top) + y_range * 0.08
  } else {
    y_upper <- y_max + y_range * expand_top
  }
  y_lower <- y_min - y_range * 0.05
  
  p <- ggplot(df_clean, aes(x = Genotype, y = .data[[index_var]], color = Genotype)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
    scale_color_manual(values = geno_colors) +
    coord_cartesian(ylim = c(y_lower, y_upper)) +
    labs(title = title, y = y_label, x = NULL) +
    theme_pubr(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 30, hjust = 1),
      plot.title      = element_text(size = 11, face = "bold")
    )
  
  if (nrow(annot) > 0) {
    p <- p +
      stat_pvalue_manual(
        annot,
        label         = "p.adj.signif",
        tip.length    = 0.01,
        step.increase = 0.08,
        hide.ns       = TRUE
      )
  }
  p
}

# Helper: build one row of four panels (prob | peak | cumsum | delay)
plot_row <- function(stats_list, index_var, y_label) {
  plot_index(stats_list$prob,   index_var, "Response probability", y_label) |
    plot_index(stats_list$peak,   index_var, "Peak distance",        "") |
    plot_index(stats_list$cumsum, index_var, "Cumulative distance",  "") |
    plot_index(stats_list$delay,  index_var, "Delay",                "")
}

# ==============================================================================
# 10. HI_between — four window variants stacked in one figure
# ==============================================================================

p_hi_combined <-
  (plot_row(stats_hi_mean,   "HI_between_mean",   "HI (mean)") /
     plot_row(stats_hi_first,  "HI_between_first",  "HI (first)") /
     plot_row(stats_hi_center, "HI_between_center", "HI (center)") /
     plot_row(stats_hi_last,   "HI_between_last",   "HI (last)")) +
  plot_annotation(
    title    = "Between-block habituation index — by stimulus window",
    subtitle = paste0(
      "HI = (B1 - B2) / (B1 + B2)  |  positive = B1 > B2  |  ",
      "first = ", first_win_label,
      "  |  center = ", center_win_label,
      "  |  last = ", last_win_label
    ),
    tag_levels = list(c(
      "Mean",   "", "", "",
      "First",  "", "", "",
      "Center", "", "", "",
      "Last",   "", "", ""
    )),
    theme = theme(
      plot.title        = element_text(size = 15, face = "bold"),
      plot.subtitle     = element_text(size = 10),
      plot.tag          = element_text(size = 10, face = "bold", angle = 90),
      plot.tag.position = "left"
    )
  )

print(p_hi_combined)
ggsave(
  file.path(save_fig_dir, "nonparam_HI_between_all_windows.png"),
  p_hi_combined,
  width = 18, height = 22, dpi = 300, bg = "white"
)

# ==============================================================================
# 11. SI_within — Block1 and Block2 rows
# ==============================================================================

p_si_combined <-
  (plot_row(stats_si_b1, "SI_block1", "SI_within") /
     plot_row(stats_si_b2, "SI_block2", "")) +
  plot_annotation(
    title    = "Within-block habituation index",
    subtitle = paste0(
      "SI = (first - last) / (first + last)  |  positive = habituation  |  ",
      "first = ", first_win_label, "  |  last = ", last_win_label
    ),
    tag_levels = list(c("Block 1", "", "", "", "Block 2", "", "", "")),
    theme = theme(
      plot.title        = element_text(size = 15, face = "bold"),
      plot.subtitle     = element_text(size = 11),
      plot.tag          = element_text(size = 11, face = "bold", angle = 90),
      plot.tag.position = "left"
    )
  )

print(p_si_combined)
ggsave(
  file.path(save_fig_dir, "nonparam_SI_within_all_outcomes.png"),
  p_si_combined,
  width = 18, height = 10, dpi = 300, bg = "white"
)

# ==============================================================================
# 12. Additional indices — absolute drop, proportional drop, asymptote, initial
# ==============================================================================

p_abs <- (plot_row(stats_abs_b1, "abs_drop_b1", "Absolute drop") /
            plot_row(stats_abs_b2, "abs_drop_b2", "")) +
  plot_annotation(
    title    = "Absolute drop (first \u2212 last)",
    subtitle = paste0("first = ", first_win_label, "  |  last = ", last_win_label),
    tag_levels = list(c("Block 1", "", "", "", "Block 2", "", "", "")),
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10),
                  plot.tag = element_text(size = 10, face = "bold", angle = 90),
                  plot.tag.position = "left")
  )

p_prop <- (plot_row(stats_prop_b1, "prop_drop_b1", "Proportional drop") /
             plot_row(stats_prop_b2, "prop_drop_b2", "")) +
  plot_annotation(
    title    = "Proportional drop  (first \u2212 last) / first",
    subtitle = paste0("first = ", first_win_label, "  |  last = ", last_win_label),
    tag_levels = list(c("Block 1", "", "", "", "Block 2", "", "", "")),
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10),
                  plot.tag = element_text(size = 10, face = "bold", angle = 90),
                  plot.tag.position = "left")
  )

p_asym <- (plot_row(stats_asym_b1, "asymptote_b1", "Asymptote") /
             plot_row(stats_asym_b2, "asymptote_b2", "")) +
  plot_annotation(
    title    = "Habituated asymptote (mean of last window)",
    subtitle = paste0("last = ", last_win_label),
    tag_levels = list(c("Block 1", "", "", "", "Block 2", "", "", "")),
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10),
                  plot.tag = element_text(size = 10, face = "bold", angle = 90),
                  plot.tag.position = "left")
  )

p_init <- (plot_row(stats_init_b1, "initial_b1", "Initial reactivity") /
             plot_row(stats_init_b2, "initial_b2", "")) +
  plot_annotation(
    title    = "Initial reactivity (mean of first window)",
    subtitle = paste0("first = ", first_win_label),
    tag_levels = list(c("Block 1", "", "", "", "Block 2", "", "", "")),
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10),
                  plot.tag = element_text(size = 10, face = "bold", angle = 90),
                  plot.tag.position = "left")
  )

print(p_abs);  print(p_prop); print(p_asym); print(p_init)

ggsave(file.path(save_fig_dir, "nonparam_abs_drop.png"),
       p_abs,  width = 18, height = 10, dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "nonparam_prop_drop.png"),
       p_prop, width = 18, height = 10, dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "nonparam_asymptote.png"),
       p_asym, width = 18, height = 10, dpi = 300, bg = "white")
ggsave(file.path(save_fig_dir, "nonparam_initial_reactivity.png"),
       p_init, width = 18, height = 10, dpi = 300, bg = "white")

message("========== FINISHED ==========")