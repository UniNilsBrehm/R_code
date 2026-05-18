reshape_data <- function(file_dir){
  df <- read.csv(file_dir)
  
  df <- df %>%
    group_by(Trial, Well) %>%
    mutate(Time_s = (row_number() - 1) * 10) %>%
    ungroup()
  
  # Clean + process
  df <- df %>%
    mutate(
      Distance = ifelse(Distance == "-", 0, Distance),
      Distance = as.numeric(Distance)
    )
  
  # Wide format with time axis
  df_wide <- df %>%
    select(Time_s, Well, Distance) %>%
    pivot_wider(
      names_from = Well,
      values_from = Distance,
      values_fill = 0
    ) %>%
    arrange(Time_s)
  
  return(df_wide)
}

# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
library(readr)        # Reading CSV files
library(DHARMa)       # Residual diagnostics for (G)LMMs
library(emmeans)      # Estimated marginal means (EMMs) and contrasts
library(glmmTMB)      # Generalized linear mixed models
library(ggplot2)      # Visualization
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(performance)  # Model diagnostics (AIC, R², etc.)
library(ggpubr)       # Publication-ready plots
library(scales)
library(stringr)
library(pracma)  # for trapz()

# Load helper functions (validation, plotting, EMM utilities, reporting)
source("C:/UniFreiburg/Code/R_code/tapping_habituation/utils_tapping.R")

# Base directory for saving results
base_dir<- "D:/WorkingData/Susana/"
df01 <- reshape_data(file.path(base_dir, "Statistics-Habituationv3_23_Photokinesis.csv"))
df02 <- reshape_data(file.path(base_dir, "Statistics-Habituationv3_24_Photokinesis.csv"))


df01_long <- df01 %>%
  pivot_longer(
    cols = -Time_s,
    names_to = "Well",
    values_to = "Distance"
  )

df01_long <- df01_long %>%
  mutate(Genotype = case_when(
    Well %in% c("A5", "B2", "B6", "D1") ~ "th, th2, tyr",
    Well %in% c("D2", "B4") ~ "tyr",
    TRUE ~ NA_character_
  ))

df02_long <- df02 %>%
  pivot_longer(
    cols = -Time_s,
    names_to = "Well",
    values_to = "Distance"
  )

df02_long <- df02_long %>%
  mutate(Genotype = case_when(
    Well %in% c("D3", "D5", "D6") ~ "th, th2, tyr",
    Well %in% c("A4", "C1") ~ "tyr",
    TRUE ~ NA_character_
  ))

df_long <- bind_rows(
  df01_long %>% mutate(Plate = "01"),
  df02_long %>% mutate(Plate = "02")
)

df_long <- df_long %>%
  mutate(Time_min = Time_s / 60)

phases <- data.frame(
  xmin = c(0, 20, 40, 60),
  xmax = c(20, 40, 60, 80),
  phase = c("light", "dark", "light", "dark")
)

df_long <- df_long %>%
  mutate(
    phase_block = floor(Time_min / 20),
    phase = ifelse(phase_block %% 2 == 0, "light", "dark")
  )

ggplot(df_long, aes(x = Time_min, y = Distance, color = Genotype)) +
  geom_rect(data = phases,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = phase),
            inherit.aes = FALSE,
            alpha = 0.15) +
  scale_fill_manual(values = c(light = "yellow", dark = "grey30")) +
  geom_line() +
  facet_wrap(~ Well) +
  labs(x = "Time (min)", y = "Distance", title = "Time vs Traces with Stimuli") +
  theme_minimal()


# Integrate (distance moved per 30 s):
df_30s <- df_long %>%
  mutate(time_bin30 = floor(Time_s / 30)) %>%
  group_by(Plate, Well, Genotype, time_bin30) %>%
  summarise(
    Time_s = min(Time_s),
    Time_min = min(Time_min),
    Distance = sum(Distance, na.rm = TRUE),
    phase = first(phase),
    .groups = "drop"
  )

sum_df <- df_30s %>%
  group_by(Genotype, Time_min, phase) %>%   # or Time_s if you prefer seconds
  summarise(
    n    = n_distinct(Well),
    mean = mean(Distance, na.rm = TRUE),
    sd   = sd(Distance, na.rm = TRUE),
    se   = sd / sqrt(n),
    ci95 = 1.96 * se,
    .groups = "drop"
  )


# Plot Mean Traces
ggplot(sum_df, aes(x = Time_min, y = mean, color = Genotype, fill = Genotype)) +
  geom_rect(
    data = phases,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = c("yellow", "grey30", "yellow", "grey30"),
    alpha = 0.12
  ) +
  geom_ribbon(
    aes(ymin = mean - ci95, ymax = mean + ci95),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "th, th2, tyr" = "red",
    "tyr" = "blue"
  )) +
  scale_fill_manual(values = c(
    "th, th2, tyr" = "#f4a3a3",  # light red
    "tyr" = "#9ec9ff"            # light blue
  )) +
  labs(x = "Time (min)", y = "Mean distance", title = "Mean traces by genotype") +
  theme_minimal()



