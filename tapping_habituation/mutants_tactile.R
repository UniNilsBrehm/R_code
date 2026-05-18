reshape_data <- function(file_dir){
  df <- read.csv(file_dir)
  
  # helper: convert "H:MM:SS.sss" to seconds
  to_seconds <- function(x) {
    parts <- str_split_fixed(x, ":", 3)
    h <- as.numeric(parts[,1])
    m <- as.numeric(parts[,2])
    s <- as.numeric(parts[,3])
    h * 3600 + m * 60 + s
  }
  
  # Clean + process
  df <- df %>%
    mutate(
      Distance = ifelse(Distance == "-", 0, Distance),
      Distance = as.numeric(Distance)
    ) %>%
    separate(Time, into = c("t_start", "t_end"), sep = "-", fill = "right") %>%
    mutate(
      t_start = ifelse(t_start == "Start", "0:00:00.000", t_start),
      t_start_s = to_seconds(t_start),
      t_end_s   = to_seconds(t_end),
      Time_s    = (t_start_s + t_end_s) / 2   # midpoint time
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
df_01_mutants <- reshape_data(file.path(base_dir, "Statistics-SPZ_Behavior_48wells_Nils1_17ms.csv"))
df_01_control <- reshape_data(file.path(base_dir, "Statistics-SPZ_Behavior_48wells_Nils1_17ms_tyr.csv"))
df_02_mutants <- reshape_data(file.path(base_dir, "Statistics-SPZ_Behavior_48wells_Nils2_17ms.csv"))
df_02_control <- reshape_data(file.path(base_dir, "Statistics-SPZ_Behavior_48wells_Nils2_17ms_tyr.csv"))



df_long <- df_02_control %>%
  pivot_longer(
    cols = -Time_s,
    names_to = "Well",
    values_to = "Distance"
  )

stim_times <- c(
  30,120,210,300,305,310,315,320,325,330,335,340,345,350,
  355,360,365,370,375,380,385,390,395,
  695,700,705,710,715,720,725,730,735,740,745,750,755,
  760,765,770,775,780,785,790
)

df_long <- df_long %>%
  mutate(time_bin = floor(Time_s / 1) * 1)

df_auc <- df_long %>%
  arrange(Well, Time_s) %>%
  group_by(Well, time_bin) %>%
  summarise(
    AUC = trapz(Time_s, Distance),
    .groups = "drop"
  )

ggplot(df_long, aes(x = Time_s, y = Distance)) +
  geom_line() +
  geom_vline(xintercept = stim_times, color = "red", alpha = 0.2) +
  facet_wrap(~ Well, scales = "free_y") +
  labs(x = "Time (s)", y = "Distance", title = "Time vs Traces with Stimuli") +
  theme_minimal()

ggplot(df_auc, aes(x = time_bin, y = AUC)) +
  geom_line() +
  geom_vline(xintercept = stim_times, color = "red", alpha = 0.2) +
  facet_wrap(~ Well, scales = "free_y") +
  labs(x = "Time (s)", y = "Distance", title = "Time vs Traces with Stimuli") +
  theme_minimal()

