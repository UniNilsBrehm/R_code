# ==============================================================================
# Combined summary figures across response variables
# ==============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(forcats)
library(patchwork)
library(scales)

source("C:/Users/NilsPC/Desktop/Susana/R_code/susana/GLMM/DarkFlash_Massed_vs_Spaced/utils.R")

base_dir <- "C:/Users/NilsPC/Desktop/Susana/Susana/GLMM/DarkFlash_Massed_vs_Spaced"

fig_dir <- file.path(base_dir, "figs", "combined_readout_summary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dirs <- list(
  response_prob = file.path(base_dir, "results", "glmm_joint_response_prob"),
  delay         = file.path(base_dir, "results", "ordinal_joint_delay"),
  max_peak      = file.path(base_dir, "results", "glmm_joint_max_peak"),
  max_cumsum    = file.path(base_dir, "results", "glmm_joint_summed_distance")
)

plot_genotypes <- c("ABTL", "th, th2, tyr", "th, tyr")
plot_colors <- genotype_colors[plot_genotypes]

read_contrast_file <- function(path, readout, direction, unit_label) {
  
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      readout = readout,
      direction = direction,
      unit_label = unit_label,
      
      benefit_est_raw  = direction * diff_est,
      benefit_low_raw  = direction * diff_low,
      benefit_high_raw = direction * diff_high,
      
      benefit_est  = benefit_est_raw,
      benefit_low  = pmin(benefit_low_raw, benefit_high_raw),
      benefit_high = pmax(benefit_low_raw, benefit_high_raw),
      
      benefit_z = benefit_est / diff_se,
      benefit_p = p,
      sig = case_when(
        benefit_p < 0.001 ~ "***",
        benefit_p < 0.01  ~ "**",
        benefit_p < 0.05  ~ "*",
        benefit_p < 0.1   ~ ".",
        TRUE              ~ ""
      )
    ) %>%
    select(-benefit_est_raw, -benefit_low_raw, -benefit_high_raw)
}

# direction:
#   +1 means larger spaced-massed difference = better spaced memory
#   -1 means smaller spaced-massed difference = better spaced memory
#

all_contrasts <- bind_rows(
  read_contrast_file(
    file.path(dirs$response_prob, "ALL_contrasts_spaced_vs_massed.csv"),
    readout = "Response probability",
    direction = -1,
    unit_label = "Lower P(response) = better retention"
  ),
  read_contrast_file(
    file.path(dirs$delay, "ALL_contrasts_spaced_vs_massed_expected_delay.csv"),
    readout = "Response delay",
    direction = 1,
    unit_label = "Longer delay = better retention"
  ),
  read_contrast_file(
    file.path(dirs$max_peak, "ALL_contrasts_spaced_vs_massed_max_peak.csv"),
    readout = "Peak distance",
    direction = -1,
    unit_label = "Lower peak movement = better retention"
  ),
  read_contrast_file(
    file.path(dirs$max_cumsum, "ALL_contrasts_spaced_vs_massed_max_cumsum.csv"),
    readout = "Summed distance",
    direction = -1,
    unit_label = "Lower summed movement = better retention"
  )
) %>%
  mutate(
    Genotype = factor(Genotype, levels = plot_genotypes),
    readout = factor(
      readout,
      levels = c(
        "Response probability",
        "Response delay",
        "Peak distance",
        "Summed distance"
      )
    ),
    contrast = factor(
      contrast,
      levels = c(
        "A_stim1",
        "B_mean_stim1_to_3",
        "C_mean_all_test_stim",
        "D_raw_recovery",
        "D_normalized_recovery",
        "D_delay_change"
      )
    ),
    contrast_label = recode(
      as.character(contrast),
      "A_stim1" = "A: test stim 1",
      "B_mean_stim1_to_3" = "B: mean stim 1–3",
      "C_mean_all_test_stim" = "C: mean all test stim",
      "D_raw_recovery" = "D: raw recovery",
      "D_normalized_recovery" = "D: normalized recovery",
      "D_delay_change" = "D: delay change"
    )
  )

# ==============================================================================
# FIGURE 1: Main memory contrast C across all readouts
# ==============================================================================

main_C <- all_contrasts %>%
  filter(contrast == "C_mean_all_test_stim")

p_heat_C <- main_C %>%
  ggplot(aes(x = readout, y = Genotype, fill = benefit_z)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig), size = 6, color = "black") +
  scale_fill_gradient2(
    low = "#b2182b",
    mid = "white",
    high = "#2166ac",
    midpoint = 0,
    name = "Spaced benefit\nstandardized z"
  ) +
  theme_pub(base_size = 13) +
  labs(
    x = NULL,
    y = NULL,
    title = "Main memory contrast C: spaced benefit across behavioral readouts",
    subtitle = "Blue = spaced better than massed; red = spaced worse than massed"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(fig_dir, "FIG1_main_contrastC_heatmap_all_readouts.png"),
  p_heat_C, width = 8.5, height = 4.5, dpi = 300, bg = "white"
)

# ==============================================================================
# FIGURE 2: Natural-unit forest plots for contrast C
# ==============================================================================

p_forest_C <- main_C %>%
  ggplot(aes(x = benefit_est, y = Genotype, color = Genotype)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(
    aes(xmin = benefit_low, xmax = benefit_high),
    linewidth = 0.8,
    size = 0.7
  ) +
  geom_text(
    aes(label = sig),
    nudge_y = 0.22,
    color = "black",
    size = 5
  ) +
  facet_wrap(~ readout, scales = "free_x", ncol = 2) +
  scale_color_manual(values = plot_colors, drop = FALSE) +
  theme_pub(base_size = 13) +
  labs(
    x = "Spaced benefit in natural units\npositive = spaced better than massed",
    y = NULL,
    title = "Main test-block effect across readouts",
    subtitle = "Contrast C: mean across all shared test stimuli"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave(
  file.path(fig_dir, "FIG2_main_contrastC_forest_natural_units.png"),
  p_forest_C, width = 10, height = 7, dpi = 300, bg = "white"
)

# ==============================================================================
# FIGURE 3: Contrast C as standardized effects
# ==============================================================================

p_bar_C <- main_C %>%
  ggplot(aes(x = readout, y = benefit_z, fill = Genotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = sig),
    position = position_dodge(width = 0.75),
    vjust = ifelse(main_C$benefit_z >= 0, -0.4, 1.3),
    size = 5
  ) +
  scale_fill_manual(values = plot_colors, drop = FALSE) +
  theme_pub(base_size = 13) +
  labs(
    x = NULL,
    y = "Standardized spaced benefit, z",
    title = "Dissociation across readouts",
    subtitle = "Positive values indicate stronger spaced retention"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave(
  file.path(fig_dir, "FIG3_contrastC_standardized_effects.png"),
  p_bar_C, width = 9, height = 5.5, dpi = 300, bg = "white"
)

# ==============================================================================
# FIGURE 4: Gating/timing vs movement
# ==============================================================================

gating_vs_motor <- main_C %>%
  mutate(
    domain = case_when(
      readout %in% c("Response probability", "Response delay") ~ "Response gating / timing",
      readout %in% c("Peak distance", "Summed distance") ~ "Movement magnitude"
    )
  )

p_domain <- gating_vs_motor %>%
  ggplot(aes(x = readout, y = benefit_z, color = Genotype, group = Genotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ domain, scales = "free_x") +
  scale_color_manual(values = plot_colors, drop = FALSE) +
  theme_pub(base_size = 13) +
  labs(
    x = NULL,
    y = "Standardized spaced benefit, z",
    title = "Spaced-memory phenotype separates response gating from movement magnitude",
    subtitle = "th, tyr shows strongest dissociation: impaired gating/timing but retained movement suppression"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave(
  file.path(fig_dir, "FIG4_gating_vs_movement_dissociation.png"),
  p_domain, width = 10, height = 5.5, dpi = 300, bg = "white"
)

# ==============================================================================
# FIGURE 5: All contrasts A/B/C/D across readouts
# ==============================================================================

p_all_heat <- all_contrasts %>%
  filter(!is.na(benefit_z)) %>%
  ggplot(aes(x = contrast_label, y = readout, fill = benefit_z)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sig), size = 4.5, color = "black") +
  facet_wrap(~ Genotype, ncol = 1) +
  scale_fill_gradient2(
    low = "#b2182b",
    mid = "white",
    high = "#2166ac",
    midpoint = 0,
    name = "Spaced benefit\nstandardized z"
  ) +
  theme_pub(base_size = 12) +
  labs(
    x = NULL,
    y = NULL,
    title = "Full contrast summary across readouts",
    subtitle = "Blue = spaced better than massed; red = spaced worse than massed"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(fig_dir, "FIG5_all_contrasts_heatmap.png"),
  p_all_heat, width = 10, height = 9, dpi = 300, bg = "white"
)

# ==============================================================================
# FIGURE 6: Multi-panel summary
# ==============================================================================

combined_main <- (p_heat_C / p_domain) +
  plot_annotation(
    title = "Spaced vs massed habituation memory across behavioral dimensions",
    subtitle = "Positive/blue values indicate a spaced-training benefit after direction alignment",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12)
    )
  )

ggsave(
  file.path(fig_dir, "FIG6_combined_main_summary.png"),
  combined_main, width = 11, height = 10, dpi = 300, bg = "white"
)



# ==============================================================================
# Export combined plotting table
# ==============================================================================

write.csv(
  all_contrasts,
  file.path(fig_dir, "combined_all_contrasts_direction_aligned.csv"),
  row.names = FALSE
)

message("Combined figures saved to: ", fig_dir)