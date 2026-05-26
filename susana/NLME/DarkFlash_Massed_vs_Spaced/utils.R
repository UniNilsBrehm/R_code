load_data <- function(file_dir, move_th = 1, drop = NULL) {
  
  # max_peak is distance moved
  
  # Load data from csv file
  df <- read_csv(file_dir)
  
  df_reduced <- df[, c(
    "Block", "Well", "Video", "Peak", "Genotype",
    "Stimulus_New", "max_peak", "max_cumsum", "peak_maxdist"
  )]
  
  names(df_reduced)[names(df_reduced) == "peak_maxdist"] <- "delay"
  
  if (!is.null(drop)) {
    
    df_reduced <- df_reduced %>%
      dplyr::filter(!Genotype %in% drop)
    
    message("Dropping genotypes:")
    print(paste(drop, collapse = " -- "))
  }
  
  # Prepare for GLMM
  df_reduced <- df_reduced %>%
    mutate(
      Genotype = factor(Genotype),
      Genotype = relevel(Genotype, ref = "ABTL")
    )
  
  df_reduced$Well <- as.factor(df_reduced$Well)
  df_reduced$Video <- as.factor(df_reduced$Video)
  df_reduced$Block <- as.factor(df_reduced$Block)
  df_reduced$Peak <- as.numeric(df_reduced$Peak)
  df_reduced$Stimulus_New <- as.numeric(df_reduced$Stimulus_New)
  
  df_reduced <- df_reduced %>%
    mutate(
      stimulus = Stimulus_New + 1
    )
  
  df_final <- df_reduced %>%
    filter(Peak == 1)
  
  df_final <- df_final %>%
    mutate(
      stimulus_log = log(stimulus)
    )
  
  # Get responses and non-responses
  # response: max_peak > move_th
  df_final$move <- ifelse(df_final$max_peak > move_th, 1, 0)
  
  df_final_sub <- subset(df_final, move > 0)
  
  return(list(
    df_final = df_final,
    df_final_sub = df_final_sub
  ))
}
