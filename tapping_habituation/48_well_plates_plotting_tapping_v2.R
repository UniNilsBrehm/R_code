# ==============================================================================
# Full script: One figure per drug/treatment with 4–5 panels (2x2 + optional LD)
# Layout per treatment:
#   col 1: ISI 90s (row1 prob, row2 distance)
#   col 2: ISI 5s  (row1 prob, row2 distance)
#   col 3: Light/Dark (spans rows 1–2) if available for that drug
#
# Notes:
# - It writes ONE image per treatment to: base_dir/figures_by_treatment/
# - Light/Dark is matched to the same "treatment" labels (Control/Acute/DK) and
#   uses the LD "drug" field to decide which LD curve belongs to each treatment.
#   If LD "drug" names differ from df$treatment names, edit `treatment_to_ld_drug`.
# ==============================================================================

# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
library(readr)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(tidyr)
library(performance)
library(ggpubr)
library(scales)
library(data.table)

library(patchwork)
library(purrr)
library(stringr)

# Load helper functions (validation, plotting, EMM utilities, reporting)
source("C:/UniFreiburg/Code/R_code/tapping_habituation/utils_tapping.R")

# Base directory for saving results
base_dir <- "D:/WorkingData/NoldusBehavior/Tapping_Habituation/R_data/48wp"

# Output directory
out_dir <- file.path(base_dir, "figures_by_treatment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load data
# ==============================================================================
res <- load_48_well_plates(file.path(base_dir, "all_treatments_48wellplates.csv"))
df <- res$df
df$auc <- df$auc + 0.0001
df_sub <- res$df_sub

# Standardize treatment names
df <- df %>%
  mutate(
    treatment = ifelse(treatment == "SCH", "SCH-23390", treatment),
    treatment = ifelse(treatment == "SKF_50", "SKF-38393-50uM", treatment),
    treatment = ifelse(treatment == "SKF_10", "SKF-38393-10uM", treatment),
    treatment = ifelse(treatment == "DA", "Dopamine", treatment)
  )

df_sub <- df_sub %>%
  mutate(
    treatment = ifelse(treatment == "SCH", "SCH-23390", treatment),
    treatment = ifelse(treatment == "SKF_50", "SKF-38393-50uM", treatment),
    treatment = ifelse(treatment == "SKF_10", "SKF-38393-10uM", treatment),
    treatment = ifelse(treatment == "DA", "Dopamine", treatment)
  )

df_5s  <- subset(df, Block %in% c("ISI_5s_block"))
df_90s <- subset(df, Block %in% c("ISI_90s"))

# Fish counts (optional)
fish_counts <- df %>%
  distinct(fish_id, treatment, drug_condition) %>%
  count(treatment, drug_condition, name = "n_fish")
print(fish_counts)

# Treatments present
treatments_all <- sort(unique(df$treatment))

# ==============================================================================
# 2. Fit models
# ==============================================================================

# ISI 5s: Response probability
m_prob_5s <- glmmTMB(
  response ~ drug_condition * stimulus_log * treatment + (1 | plate_uid/fish_id),
  data = df_5s,
  family = binomial
)

# ISI 5s: Distance moved (AUC)
m_dist_5s <- glmmTMB(
  auc ~ drug_condition * stimulus_log * treatment + (1 | plate_uid/fish_id),
  data = df_5s,
  family = Gamma(link = "log")
)

# ISI 90s: Response probability
m_prob_90s <- glmmTMB(
  response ~ drug_condition * treatment + stimulus_log + (1 | plate_id/fish_id),
  data = df_90s,
  family = binomial
)

# ISI 90s: Distance moved (AUC)
df_90s$auc <- df_90s$auc + 0.0001
m_dist_90s <- glmmTMB(
  auc ~ drug_condition * stimulus_log * treatment + (1 | plate_id/fish_id),
  data = df_90s,
  family = Gamma(link = "log")
)

# ==============================================================================
# 3. Build prediction grids + observed summaries (store separately)
# ==============================================================================

# ---------------------------
# ISI 5s: PROB
# ---------------------------
valid_combos_5s <- df_5s %>% distinct(drug_condition, treatment)

pred_grid_prob_5s <- expand.grid(stimulus = 1:20) %>%
  mutate(stimulus_log = log(stimulus)) %>%
  tidyr::crossing(valid_combos_5s)

pred_prob_5s <- predict(m_prob_5s, newdata = pred_grid_prob_5s, re.form = NA, se.fit = TRUE)

pred_grid_prob_5s <- pred_grid_prob_5s %>%
  mutate(
    fit     = plogis(pred_prob_5s$fit),
    CI_low  = plogis(pred_prob_5s$fit - 1.96 * pred_prob_5s$se.fit),
    CI_high = plogis(pred_prob_5s$fit + 1.96 * pred_prob_5s$se.fit)
  )

obs_prob_5s <- as.data.table(df_5s)[, .(
  response = mean(response == 1)
), by = .(treatment, stimulus, drug_condition)]


# ---------------------------
# ISI 5s: DIST
# ---------------------------
pred_grid_dist_5s <- expand.grid(stimulus = 1:20) %>%
  mutate(stimulus_log = log(stimulus)) %>%
  tidyr::crossing(valid_combos_5s)

pred_dist_5s <- predict(m_dist_5s, newdata = pred_grid_dist_5s, re.form = NA, se.fit = TRUE)

pred_grid_dist_5s <- pred_grid_dist_5s %>%
  mutate(
    fit     = exp(pred_dist_5s$fit),
    CI_low  = exp(pred_dist_5s$fit - 1.96 * pred_dist_5s$se.fit),
    CI_high = exp(pred_dist_5s$fit + 1.96 * pred_dist_5s$se.fit)
  )

obs_dist_5s <- as.data.table(df_5s)[, .(
  response = mean(auc)
), by = .(treatment, stimulus, drug_condition)]


# ---------------------------
# ISI 90s: PROB
# ---------------------------
valid_combos_90s <- df_90s %>% distinct(drug_condition, treatment)

pred_grid_prob_90s <- expand.grid(stimulus = 1:3) %>%
  mutate(stimulus_log = log(stimulus)) %>%
  tidyr::crossing(valid_combos_90s)

pred_prob_90s <- predict(m_prob_90s, newdata = pred_grid_prob_90s, re.form = NA, se.fit = TRUE)

pred_grid_prob_90s <- pred_grid_prob_90s %>%
  mutate(
    fit     = plogis(pred_prob_90s$fit),
    CI_low  = plogis(pred_prob_90s$fit - 1.96 * pred_prob_90s$se.fit),
    CI_high = plogis(pred_prob_90s$fit + 1.96 * pred_prob_90s$se.fit)
  )

obs_prob_90s <- as.data.table(df_90s)[, .(
  response = mean(response == 1)
), by = .(treatment, stimulus, drug_condition)]


# ---------------------------
# ISI 90s: DIST
# ---------------------------
pred_grid_dist_90s <- expand.grid(stimulus = 1:3) %>%
  mutate(stimulus_log = log(stimulus)) %>%
  tidyr::crossing(valid_combos_90s)

pred_dist_90s <- predict(m_dist_90s, newdata = pred_grid_dist_90s, re.form = NA, se.fit = TRUE)

pred_grid_dist_90s <- pred_grid_dist_90s %>%
  mutate(
    fit     = exp(pred_dist_90s$fit),
    CI_low  = exp(pred_dist_90s$fit - 1.96 * pred_dist_90s$se.fit),
    CI_high = exp(pred_dist_90s$fit + 1.96 * pred_dist_90s$se.fit)
  )

obs_dist_90s <- as.data.table(df_90s)[, .(
  response = mean(auc)
), by = .(treatment, stimulus, drug_condition)]

# ==============================================================================
# 4. Light/Dark: load + summarize + helper plot function (optional panel)
# ==============================================================================
ld_file <- file.path(base_dir, "light_dark_all_treatments_48wellplates.csv")
df_ld <- read_csv(ld_file)

# Fish counts (optional)
fish_counts_ld <- df_ld %>%
  distinct(larva_id, drug, treatment) %>%
  count(drug, treatment, name = "n_fish")
print(fish_counts_ld)


ld_df <- df_ld %>%
  select(-`...1`) %>%
  mutate(
    treatment     = factor(treatment, levels = c("Control","Acute","DK")),
    phase         = factor(phase, levels = c("Light","Dark")),
    block         = factor(block, levels = c("L1","D1","L2","D2")),
    plate         = factor(plate),
    larva_id      = factor(larva_id),
    drug          = factor(drug),
    concentration = factor(concentration),
    time_c        = as.numeric(scale(time_sec, center = TRUE, scale = TRUE))
  ) %>%
  filter(!is.na(block))

sum_df <- ld_df %>%
  group_by(drug, concentration, treatment, time_bin, time_sec, phase, block) %>%
  summarise(
    n    = n_distinct(larva_id),
    mean = mean(distance, na.rm = TRUE),
    sd   = sd(distance, na.rm = TRUE),
    se   = sd / sqrt(n),
    ci95 = 1.96 * se,
    .groups = "drop"
  )

dark_rects <- data.frame(
  xmin = c(20*60, 60*60),
  xmax = c(40*60, 80*60)
)

# If LD "drug" strings don't match df$treatment strings, edit this mapping.
# By default, we assume LD drug == treatment.
treatment_to_ld_drug <- setNames(treatments_all, treatments_all)

# ==============================================================================
# 5. Helper functions for the 4 main panels
# ==============================================================================
drug_cols <- c(
  control = "#F8766D",  # red
  dk      = "#619CFF",  # blue
  acute   = "#00BA38"   # green
)

# A true blank/white panel that still occupies space
blank_panel <- function() {
  ggplot() +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position  = "none"
    )
}

# Optional: blank panel WITH a title "Light/Dark" if you want consistent labeling
blank_panel_titled <- function(title = "Light/Dark") {
  ggplot() +
    theme_void() +
    ggtitle(title) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position  = "none"
    )
}

base_theme_small <- function() {
  theme_pubr(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
}

plot_ld_for_drug <- function(treatment_name) {
  ld_drug <- unname(treatment_to_ld_drug[[treatment_name]])
  if (is.null(ld_drug) || is.na(ld_drug)) return(NULL)
  
  d <- sum_df %>% 
    filter(as.character(drug) == as.character(ld_drug)) %>%
    mutate(
      treatment_key = tolower(as.character(treatment)),
      treatment_key = dplyr::recode(
        treatment_key,
        "control" = "control",
        "dk"      = "dk",
        "acute"   = "acute",
        .default  = NA_character_
      )
    )
  
  if (nrow(d) == 0) return(NULL)
  
  ggplot(d, aes(x = time_sec, y = mean, color = treatment_key)) +
    geom_rect(
      data = dark_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      alpha = 0.10
    ) +
    geom_ribbon(
      aes(ymin = mean - ci95, ymax = mean + ci95, fill = treatment_key),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = drug_cols, na.value = "grey50") +
    scale_fill_manual(values = drug_cols, na.value = "grey50") +
    guides(
      color = guide_legend(title = NULL),
      fill  = "none"
    ) +
    scale_x_continuous(
      name   = "Time (min)",
      breaks = seq(0, 80, by = 20) * 60,
      labels = seq(0, 80, by = 20)
    ) +
    ylab("Distance (mm/30s)") +
    ggtitle("Light/Dark") +
    theme_bw(base_size = 12) +
    theme(
      # Remove all grid lines
      panel.grid = element_blank(),
      
      # Remove top and right borders
      panel.border = element_blank(),
      
      # Add only left and bottom axis lines (thicker)
      axis.line.x = element_line(size = 0.5, color = "black"),
      axis.line.y = element_line(size = 0.5, color = "black"),
      
      legend.position = NULL,
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )
}

plot_prob_panel <- function(pred_df, obs_df, treatment_name, title_txt,
                            x_breaks, x_limits, y_limits = c(0, 1), y_breaks = seq(0, 1, 0.2)) {
  
  p_pred <- pred_df %>% filter(treatment == treatment_name)
  p_obs  <- obs_df  %>% filter(treatment == treatment_name)
  
  ggplot(p_pred, aes(x = stimulus, color = drug_condition)) +
    geom_line(aes(y = fit), linewidth = 1.0) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = drug_condition),
                alpha = 0.18, color = NA) +
    geom_point(data = p_obs,
               aes(x = stimulus, y = response, color = drug_condition),
               alpha = 0.6, size = 1.3) +
    
    # FIXED COLORS HERE
    scale_color_manual(values = drug_cols) +
    scale_fill_manual(values = drug_cols) +
    guides(
      color = guide_legend(title = NULL),
      fill  = "none"
    ) +
    scale_x_continuous(breaks = x_breaks, limits = x_limits,
                       expand = expansion(mult = c(0.02, 0.05))) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks,
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = "Tap #", y = "Response prob", color = NULL, fill = NULL, title = title_txt) +
    base_theme_small()+
    theme(plot.title = element_text(hjust = 0.5))
}

plot_dist_panel <- function(pred_df, obs_df, treatment_name, title_txt,
                            x_breaks, x_limits, y_breaks_fun = scales::breaks_extended(n = 5)) {
  
  p_pred <- pred_df %>% filter(treatment == treatment_name)
  p_obs  <- obs_df  %>% filter(treatment == treatment_name)
  
  ggplot(p_pred, aes(x = stimulus, color = drug_condition)) +
    geom_line(aes(y = fit), linewidth = 1.0) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = drug_condition),
                alpha = 0.18, color = NA) +
    geom_point(data = p_obs,
               aes(x = stimulus, y = response, color = drug_condition),
               alpha = 0.6, size = 1.3) +
    
    # FIXED COLORS
    scale_color_manual(values = drug_cols) +
    scale_fill_manual(values = drug_cols) +
    guides(
      color = guide_legend(title = NULL),
      fill  = "none"
    ) +
    scale_x_continuous(breaks = x_breaks, limits = x_limits,
                       expand = expansion(mult = c(0.02, 0.05))) +
    scale_y_continuous(breaks = y_breaks_fun,
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = "Tap #", y = "Distance moved (mm)", color = NULL, fill = NULL, title = title_txt) +
    base_theme_small()+
    theme(plot.title = element_text(hjust = 0.5))
}

# ==============================================================================
# 6. Build one figure per treatment (2x2 + optional LD column)
# ==============================================================================

# Control how wide/tall output is
fig_width_in  <- 14
fig_height_in <- 7
fig_dpi <- 600

for (tr in treatments_all) {
  
  # Panels:
  p90_prob <- plot_prob_panel(
    pred_df = pred_grid_prob_90s,
    obs_df  = obs_prob_90s,
    treatment_name = tr,
    title_txt = "ISI 90s",
    x_breaks = c(1, 2, 3),
    x_limits = c(1, 3),
    y_limits = c(0, 1),
    y_breaks = seq(0, 1, 0.2)
  )
  
  p90_dist <- plot_dist_panel(
    pred_df = pred_grid_dist_90s,
    obs_df  = obs_dist_90s,
    treatment_name = tr,
    title_txt = NULL,
    x_breaks = c(1, 2, 3),
    x_limits = c(1, 3)
  )
  
  p5_prob <- plot_prob_panel(
    pred_df = pred_grid_prob_5s,
    obs_df  = obs_prob_5s,
    treatment_name = tr,
    title_txt = "ISI 5s",
    x_breaks = c(1, 5, 10, 15, 20),
    x_limits = c(1, 20),
    y_limits = c(0, 1),
    y_breaks = seq(0, 1, 0.2)
  )
  
  p5_dist <- plot_dist_panel(
    pred_df = pred_grid_dist_5s,
    obs_df  = obs_dist_5s,
    treatment_name = tr,
    title_txt = NULL,
    x_breaks = c(1, 5, 10, 15, 20),
    x_limits = c(1, 20)
  )
  
  # Optional LD (spans both rows)
  p_ld <- NULL
  if (!is.null(plot_ld_for_drug)) {
    p_ld <- plot_ld_for_drug(tr)
  }
  
  # Assemble layout:
  left_block <- (p90_prob / p90_dist) | (p5_prob / p5_dist)
  
  # Always define the 3rd panel
  p_ld <- blank_panel()
  if (!is.null(plot_ld_for_drug)) {
    tmp <- plot_ld_for_drug(tr)
    if (!is.null(tmp)) p_ld <- tmp
  }
  
  # Always use the same widths (so nothing stretches)
  full_plot <- (left_block | p_ld) +
    plot_layout(widths = c(1, 1, 1.8), guides = "collect") +   # <- key line
    plot_annotation(title = tr) &
    theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      legend.position = "top"                                  # or "bottom"
    )
  # Save
  out_file <- file.path(out_dir, paste0("Panels_", str_replace_all(tr, "[^A-Za-z0-9_-]", "_"), ".jpg"))
  ggsave(out_file, plot = full_plot, width = fig_width_in, height = fig_height_in, units = "in", dpi = fig_dpi)
  
  message("Saved: ", out_file)
}

# ==============================================================================
# End
# ==============================================================================