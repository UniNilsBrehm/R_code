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
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# Helper Functions
# ==============================================================================
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

plot_model_response <- function(model, data, y_var, y_label = y_var) {
  library(ggeffects)
  library(ggplot2)
  library(dplyr)
  
  pred <- ggpredict(
    model,
    terms = c("block_tap [all]", "genotype", "phase")
  )
  
  ggplot(pred, aes(x = x, y = predicted, color = group, fill = group)) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.2,
      color = NA
    ) +
    geom_line(linewidth = 1.2) +
    geom_point(
      data = data %>% filter(phase != "isolated"),
      aes(x = block_tap, y = .data[[y_var]], color = genotype),
      inherit.aes = FALSE,
      alpha = 0.5,
      position = position_jitter(width = 0.15, height = 0)
    ) +
    facet_wrap(~facet) +
    labs(
      x = "Tap number",
      y = y_label,
      color = "Genotype",
      fill = "Genotype"
    ) +
    theme_classic(base_size = 14)
}

# ==============================================================================
# Load Data
# ==============================================================================
# Base directory for saving results
base_dir<- "C:/Users/NilsPC/Desktop/Susana/Susana/Tapping"
df_01_mutants <- reshape_data(file.path(base_dir, "data_files", "Statistics-SPZ_Behavior_48wells_Nils1_17ms.csv"))
df_01_control <- reshape_data(file.path(base_dir, "data_files", "Statistics-SPZ_Behavior_48wells_Nils1_17ms_tyr.csv"))
df_02_mutants <- reshape_data(file.path(base_dir, "data_files", "Statistics-SPZ_Behavior_48wells_Nils2_17ms.csv"))
df_02_control <- reshape_data(file.path(base_dir, "data_files", "Statistics-SPZ_Behavior_48wells_Nils2_17ms_tyr.csv"))


tap_times <- c(
  30,120,210,300,305,310,315,320,325,330,335,340,345,350,
  355,360,365,370,375,380,385,390,395,
  695,700,705,710,715,720,725,730,735,740,745,750,755,
  760,765,770,775,780,785,790
)


stim_info <- tibble(
  tap_id = seq_along(tap_times),
  tap_time = tap_times
) %>%
  mutate(
    phase = case_when(
      tap_id <= 3 ~ "isolated",
      tap_time < 500 ~ "block1",
      TRUE ~ "block2"
    ),
    
    isi = c(NA, diff(tap_times)),
    
    block_tap = case_when(
      phase == "isolated" ~ tap_id,
      phase == "block1" ~ row_number() - 3,
      phase == "block2" ~ row_number() - 23
    )
  )

df_mutants <- full_join(
  df_01_mutants,
  df_02_mutants,
  by = "Time_s"
)

df_controls <- full_join(
  df_01_control,
  df_02_control,
  by = "Time_s"
)

df_long_controls <- df_controls %>%
  pivot_longer(
    cols = -Time_s,
    names_to = "Well",
    values_to = "Distance"
  ) %>%
  mutate(
    time_bin = floor(Time_s),
    genotype = "control"
  )

df_long_mutants <- df_mutants %>%
  pivot_longer(
    cols = -Time_s,
    names_to = "Well",
    values_to = "Distance"
  ) %>%
  mutate(
    time_bin = floor(Time_s),
    genotype = "mutant"
  )

df_long_all <- bind_rows(df_long_controls, df_long_mutants)


ggplot(df_long_all, aes(x = Time_s, y = Distance)) +
  geom_line() +
  geom_vline(xintercept = tap_times, color = "red", alpha = 0.2) +
  facet_wrap(~ Well, scales = "free_y") +
  labs(x = "Time (s)", y = "Distance", title = "Time vs Traces with Stimuli") +
  theme_minimal()


# Extract response variables
response_window <- 2
threshold_peak <- 1

responses <- map_dfr(seq_along(tap_times), function(i) {
  tap_time <- tap_times[i]
  
  df_long_all %>%
    filter(
      Time_s >= tap_time,
      Time_s <= tap_time + response_window
    ) %>%
    group_by(genotype, Well) %>%
    summarise(
      tap_id = i,
      tap_time = tap_time,
      peak_distance = max(Distance, na.rm = TRUE),
      summed_distance = sum(Distance, na.rm = TRUE),
      response_delay = {
        above <- which(Distance >= threshold_peak)
        if (length(above) == 0) NA_real_ else Time_s[above[1]] - tap_time
      },
      response_yes_no = peak_distance >= threshold_peak,
      .groups = "drop"
    )
})

responses <- responses %>%
  left_join(stim_info, by = c("tap_id", "tap_time"))


responses_summary <- responses %>%
  group_by(genotype, tap_id) %>%
  summarise(
    mean_peak = mean(peak_distance, na.rm = TRUE),
    sem_peak = sd(peak_distance, na.rm = TRUE) / sqrt(n()),
    mean_sum = mean(summed_distance, na.rm = TRUE),
    sem_sum = sd(summed_distance, na.rm = TRUE) / sqrt(n()),
    mean_delay = mean(response_delay, na.rm = TRUE),
    sem_delay = sd(response_delay, na.rm = TRUE) / sqrt(n()),
    response_rate = mean(response_yes_no, na.rm = TRUE),
    .groups = "drop"
  )


plot_df_sem <- responses %>%
  group_by(genotype, phase, block_tap) %>%
  summarise(
    mean_peak = mean(peak_distance, na.rm = TRUE),
    sem_peak = sd(peak_distance, na.rm = TRUE) / sqrt(sum(!is.na(peak_distance))),
    
    mean_sum = mean(summed_distance, na.rm = TRUE),
    sem_sum = sd(summed_distance, na.rm = TRUE) / sqrt(sum(!is.na(summed_distance))),
    
    mean_delay = mean(response_delay, na.rm = TRUE),
    sem_delay = sd(response_delay, na.rm = TRUE) / sqrt(sum(!is.na(response_delay))),
    
    mean_prob = mean(response_yes_no, na.rm = TRUE),
    sem_prob = sd(as.numeric(response_yes_no), na.rm = TRUE) / sqrt(sum(!is.na(response_yes_no))),
    
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(
      mean_peak, sem_peak,
      mean_sum, sem_sum,
      mean_delay, sem_delay,
      mean_prob, sem_prob
    ),
    names_to = c(".value", "response_variable"),
    names_pattern = "(mean|sem)_(.*)"
  )

ggplot(
  plot_df_sem,
  aes(
    x = block_tap,
    y = mean,
    color = genotype,
    group = genotype
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = mean - sem,
      ymax = mean + sem
    ),
    width = 0.2
  ) +
  facet_grid(
    response_variable ~ phase,
    scales = "free"
  ) +
  labs(
    x = "Stimulus number within phase",
    y = "Mean ± SEM",
    color = "Genotype",
    title = "Response variables across stimulus phases"
  ) +
  theme_classic()

# Model
# Response Probability
responses <- responses %>%
  mutate(response_yes_no_num = as.numeric(response_yes_no))
responses_sub <- responses %>%
  filter(response_yes_no == TRUE)

m_response_prob <- glmmTMB(
  response_yes_no_num ~ genotype * phase * log(block_tap) + (1 | Well),
  data = responses %>% filter(phase != "isolated"),
  family = binomial()
)
summary(m_response_prob)

# Peak distance moved
m_peak <- glmmTMB(
  peak_distance  ~ genotype * phase * log(block_tap) + (1 | Well),
  data = responses_sub %>% filter(phase != "isolated"),
  family = Gamma(link="log")
)
summary(m_peak)

# Summed distance moved
m_sum <- glmmTMB(
  summed_distance  ~ genotype * phase * log(block_tap) + (1 | Well),
  data = responses_sub %>% filter(phase != "isolated"),
  family = Gamma(link="log")
)
summary(m_sum)

# delay
m_delay <- glmmTMB(
  response_delay  ~ genotype * phase * log(block_tap) + (1 | Well),
  data = responses_sub %>% filter(phase != "isolated"),
  family = gaussian(link='identity')
)
summary(m_delay)


# Check Residuals
library(DHARMa)
sim <- simulateResiduals(m_response_prob)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

# Model Plots
plot_model_response(m_peak, responses_sub, "peak_distance", "Peak distance moved")

plot_model_response(m_sum, responses_sub, "summed_distance", "Summed distance moved")

plot_model_response(m_delay, responses_sub, "response_delay", "Response delay")

plot_model_response(m_response_prob, responses, "response_yes_no_num", "Response probability")
