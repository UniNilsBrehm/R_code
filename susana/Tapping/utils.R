load_data_tapping <- function(file_dir, move_th = 1, keep = NULL, take_peak = 0, ref = "ABTL") {
  # max_peak is distance moved
  # Load data from csv file
  df <- read_csv(file_dir)
  df_reduced <- df[, c("Well", "Video", "Peak", "Genotype", "Stimulus_New", "max_peak", "max_cumsum", "peak_maxdist")]
  
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
      Genotype = relevel(Genotype, ref = ref),
    )
  
  df_reduced$Well <- as.factor(df_reduced$Well)
  df_reduced$Video <- as.factor(df_reduced$Video)
  df_reduced$Peak <- as.numeric(df_reduced$Peak)
  df_reduced$Stimulus_New <- as.numeric(df_reduced$Stimulus_New)
  
  df_reduced <- df_reduced %>%
    mutate(
      stimulus = Stimulus_New,
    )
  
  # For each Peak (-10 to 4) the distance move value (max_peak) is the same
  df_final <- df_reduced %>%
    group_by(Well, Video, Stimulus_New) %>%
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
