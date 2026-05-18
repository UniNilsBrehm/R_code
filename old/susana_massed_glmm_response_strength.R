library(readr)
library(lme4)
library(lmerTest)  # for p-values
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(scales)
library(data.table)
library(tibble)
library(tidyr)
library(purrr)
library(furrr)
library(performance)
library(sjPlot)
library(ggpubr)

# ==============================================================================
# HELPER FUNCTIONS
check_spread <- function(residuals_obj, group, data = NULL) {
  # Ensure required package
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package 'car' is required for Levene's test. Please install it with install.packages('car').")
  }
  
  # Extract residuals
  res <- residuals(residuals_obj)
  
  # Handle both vector or column name for group
  if (is.character(substitute(group))) {
    group_name <- deparse(substitute(group))
    group_var <- data[[group_name]]
  } else {
    group_var <- group
    group_name <- deparse(substitute(group))
  }
  
  # Run Levene's Test (center = median, more robust)
  lv <- car::leveneTest(res ~ group_var, center = "median")
  
  # Extract test info
  Fval <- lv$`F value`[1]
  df1  <- lv$Df[1]
  df2  <- lv$Df[2]
  pval <- lv$`Pr(>F)`[1]
  eta2_partial <- (Fval * df1) / (Fval * df1 + df2)
  
  # Print formatted summary
  cat("\n--- Levene's Test for Homogeneity of Variance ---\n")
  cat(sprintf("Grouping variable: %s\n", group_name))
  cat(sprintf("F(%d, %d) = %.3f, p = %.4f\n", df1, df2, Fval, pval))
  cat(sprintf("Partial eta² = %.5f\n", eta2_partial))
  
  # Interpret magnitude
  magnitude <- dplyr::case_when(
    eta2_partial < 0.01 ~ "very small",
    eta2_partial < 0.06 ~ "small",
    eta2_partial < 0.14 ~ "medium",
    TRUE ~ "large"
  )
  cat(sprintf("Interpretation: %s effect size\n", magnitude))
  
  # Return invisibly (for programmatic use)
  invisible(list(
    levene = lv,
    eta2_partial = eta2_partial,
    F = Fval,
    df = c(df1, df2),
    p.value = pval,
    magnitude = magnitude
  ))
}

# ==============================================================================
# Load data from csv file
df <- read_csv("D:/WorkingData/Susana/SPZ_Massed_Training_7Nov2025.csv")
df_reduced <- df[, c("Block", "Well", "Video", "Peak", "Genotype", "Stimulus_New", "max_peak", "max_cumsum", "peak_maxdist")]

df_reduced <- df_reduced %>%
  rename(
    "delay" = "peak_maxdist"
  )

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
    stimulus = Stimulus_New + 1,
  )

df_final <- df_reduced %>%
  filter(Peak == 1)

df_final <- df_final %>%
  mutate(
    stimulus_log = log(stimulus)
  )

# Get responses and non-responses (response: max_peak > 0)
move_th <- 1
df_final$move <- ifelse(df_final$max_peak > move_th, 1, 0)
df_final_sub <- subset(df_final, move > 0)
df_final_sub$delay_ord <- ordered(df_final_sub$delay) 

# Get number of animals
n_per_genotype <- df_final %>%
  dplyr::distinct(Video, Well, Genotype) %>%
  dplyr::count(Genotype)

# # Check data set
# df_final %>%
#   group_by(Block) %>%
#   summarise(
#     min_stimulus = min(stimulus, na.rm = TRUE),
#     max_stimulus = max(stimulus, na.rm = TRUE),
#     n_unique = n_distinct(stimulus)
#   )
# df_final %>%
#   group_by(Block) %>%
#   summarise(
#     missing = list(setdiff(
#       seq(min(stimulus, na.rm = TRUE), max(stimulus, na.rm = TRUE)),
#       unique(stimulus)
#     ))
#   )
# 
# df_final %>%
#   filter(Block == 1) %>%
#   count(stimulus) %>%
#   ggplot(aes(x = stimulus, y = n)) +
#   geom_col(fill = "steelblue") +
#   labs(
#     title = "Counts per stimulus (Block 1)",
#     x = "Stimulus number",
#     y = "Number of observations"
#   ) +
#   theme_minimal()

# Histograms
h1 <- ggplot(df_final, aes(x = max_peak)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_grid(col=vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "All Data",
    x = "Peak Distance Moved",
    y = "Count"
  )

ggsave(
  filename = "D:/WorkingData/Susana/results/massed/figures/distance_moved/peak_dist_all_dat.pdf",
  plot = h1,
  width = 6, height = 4
)

h2 <- ggplot(df_final_sub, aes(x = max_peak)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_grid(col=vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Only responsive",
    x = "Peak Distance Moved",
    y = "Count"
  )

ggsave(
  filename = "D:/WorkingData/Susana/results/massed/figures/distance_moved/peak_dist_selected.pdf",
  plot = h2,
  width = 6, height = 4
)

# ==============================================================================
# 2. Magnitude of movement (max_peak)
# Distance moved at Peak of Response
m_peak <- glmmTMB(
  max_peak  ~ Genotype * stimulus_log  * Block + (1|Video/Well),
  # dispformula = ~ Genotype,
  family = Gamma(link = "log"),
  data = df_final_sub
)

# Distance moved summed after stimulus onset
m_sum <- glmmTMB(
  max_cumsum  ~ Genotype * stimulus_log  * Block + (1|Video/Well),
  # dispformula = ~ Genotype,
  family = Gamma(link = "log"),
  data = df_final_sub
)

# Delay to stimulus onset
# m_delay_poisson <- glmmTMB(
#   delay ~ Genotype * stimulus_log * Block + (1 | Video/Well),
#   # family = nbinom2(link = "log"),
#   family = poisson(link = "log"),
#   data = df_final_sub
# )

m_delay_gaussian <- glmmTMB(
  delay ~ Genotype * stimulus * Block + (1 | Video/Well),
  family = gaussian(link = "identity"),
  data = df_final_sub
)


# compare_performance(m_delay_poisson, m_delay_gaussian)
summary(m_peak)
summary(m_sum)
summary(m_delay_gaussian)

# Validation
# Peak
res_m_peak <- simulateResiduals(m_peak)
plot(res_m_peak)
plotResiduals(res_m_peak, df_final_sub$Genotype)
plotResiduals(res_m_peak, df_final_sub$Block)

testDispersion(res_m_peak)
testZeroInflation(res_m_peak)
testUniformity(res_m_peak) 
testOutliers(res_m_peak, type = "bootstrap")

check_spread(res_m_peak, df_final_sub$Genotype)
check_spread(res_m_peak, df_final_sub$Block)


m_peak_performance <- model_performance(m_peak)  # AIC, R2, RMSE, ICC, etc.
check_collinearity(m_peak) # Multicollinearity (VIFs)
plot_model(m_peak, type = "re", sort.est = TRUE)

# Summed
res_m_sum <- simulateResiduals(m_sum)
plot(res_m_sum)

testDispersion(res_m_sum)
testZeroInflation(res_m_sum)
testUniformity(res_m_sum)
testOutliers(res_m_sum, type = "bootstrap")

plotResiduals(res_m_sum, df_final_sub$Genotype)
plotResiduals(res_m_sum, df_final_sub$Block)

m_summed_performance <-model_performance(m_sum)  # AIC, R2, RMSE, ICC, etc.
check_collinearity(m_sum) # Multicollinearity (VIFs)

# Delays
# Summed
res_m_delay <- simulateResiduals(m_delay_gaussian)
plot(res_m_delay)

testUniformity(res_m_delay)
testOutliers(res_m_delay, type = "bootstrap")
testDispersion(res_m_delay)
testZeroInflation(res_m_delay)

plotResiduals(res_m_delay, df_final_sub$Genotype)
plotResiduals(res_m_delay, df_final_sub$Block)

m_delayed_performance <-model_performance(m_delay_gaussian)  # AIC, R2, RMSE, ICC, etc.
check_collinearity(m_delay_gaussian) # Multicollinearity (VIFs)


# check_heteroscedasticity(m2)
# check_distribution(m2)
# check_collinearity(update(m2, . ~ Genotype + Block + stimulus_log))


# Plot Peak
# --- 1. Prediction grid ------------------------------------
new_data_2 <- bind_rows(
  expand.grid(Block = "1", stimulus = 1:477),
  expand.grid(Block = "2", stimulus = 1:9)
) %>%
  tidyr::crossing(Genotype = unique(df_final$Genotype)) %>%
  mutate(stimulus_log = log(stimulus))

# --- 2. Model predictions ----------------------------------
pred <- predict(m_peak, newdata = new_data_2, re.form = NA, se.fit = TRUE)

# log link → back-transform with exp()
new_data_2 <- new_data_2 %>%
  mutate(
    fit     = exp(pred$fit),
    CI_low  = exp(pred$fit - 1.96 * pred$se.fit),
    CI_high = exp(pred$fit + 1.96 * pred$se.fit)
  )

# --- 3. Plot -----------------------------------------------
g_peak <- ggplot(new_data_2, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_wrap(
    ~Block, ncol = 2, scales = "free_x",
    labeller = as_labeller(c(
      `1` = "Block 1: 477 flashes",
      `2` = "Block 2: 8 flashes"
    ))
  ) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.1, color = NA) +
  labs(
    x = "Stimulus number (within block)",
    y = "Peak distance moved",
    color = "Genotype",
    fill  = "Genotype"
  ) +
  theme_pubr(base_size = 14)

g_peak
ggsave(
  filename = "D:/WorkingData/Susana/results/massed/figures/distance_moved/peak_distance_habituation_curves.pdf",
  plot = g_peak,
  width = 6, height = 4
)


# Plot Sum
# --- 1. Prediction grid ------------------------------------
new_data_2 <- bind_rows(
  expand.grid(Block = "1", stimulus = 1:477),
  expand.grid(Block = "2", stimulus = 1:9)
) %>%
  tidyr::crossing(Genotype = unique(df_final$Genotype)) %>%
  mutate(stimulus_log = log(stimulus))

# --- 2. Model predictions ----------------------------------
pred <- predict(m_sum, newdata = new_data_2, re.form = NA, se.fit = TRUE)

# log link → back-transform with exp()
new_data_2 <- new_data_2 %>%
  mutate(
    fit     = exp(pred$fit),
    CI_low  = exp(pred$fit - 1.96 * pred$se.fit),
    CI_high = exp(pred$fit + 1.96 * pred$se.fit)
  )

# --- 3. Plot -----------------------------------------------
g_summed <- ggplot(new_data_2, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_wrap(
    ~Block, ncol = 2, scales = "free_x",
    labeller = as_labeller(c(
      `1` = "Block 1: 477 flashes",
      `2` = "Block 2: 8 flashes"
    ))
  ) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.1, color = NA) +
  labs(
    x = "Stimulus number (within block)",
    y = "Summed distance moved",
    color = "Genotype",
    fill  = "Genotype"
  ) +
  theme_pubr(base_size = 14)

g_summed
ggsave(
  filename = "D:/WorkingData/Susana/results/massed/figures/distance_moved/summed_distance_habituation_curves.pdf",
  plot = g_summed,
  width = 6, height = 4
)


# Plot Delay
# --- 1. Prediction grid ------------------------------------
new_data_delay <- bind_rows(
  expand.grid(Block = "1", stimulus = 1:477),
  expand.grid(Block = "2", stimulus = 1:9)
) %>%
  tidyr::crossing(Genotype = unique(df_final$Genotype)) %>%
  mutate(stimulus_log = log(stimulus))

# --- 2. Model predictions ----------------------------------
pred_delay <- predict(m_delay_gaussian, newdata = new_data_delay,
                      re.form = NA, se.fit = TRUE)

# Gaussian link = identity → no back-transformation
new_data_delay <- new_data_delay %>%
  mutate(
    fit     = pred_delay$fit,
    CI_low  = pred_delay$fit - 1.96 * pred_delay$se.fit,
    CI_high = pred_delay$fit + 1.96 * pred_delay$se.fit
  )

# --- 3. Plot -----------------------------------------------
g_delay <- ggplot(new_data_delay,
       aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_wrap(
    ~Block, ncol = 2, scales = "free_x",
    labeller = as_labeller(c(
      `1` = "Block 1: 477 flashes",
      `2` = "Block 2: 8 flashes"
    ))
  ) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.1, color = NA) +
  labs(
    x = "Stimulus number (within block)",
    y = "Response delay (s)",
    color = "Genotype",
    fill  = "Genotype"
  ) +
  theme_pubr(base_size = 14)

g_delay
ggsave(
  filename = "D:/WorkingData/Susana/results/massed/figures/distance_moved/delay_distance_habituation_curves.pdf",
  plot = g_delay,
  width = 6, height = 4
)


# ------------------------------------------------------------------------------
# Estimated Marginal Means
# Helper Function
pretty_pairs <- function(emm) {
  df <- as.data.frame(pairs(emm))
  
  # Normalize estimate column name (handle both "ratio" or "estimate")
  if ("ratio" %in% names(df)) {
    df <- dplyr::rename(df, estimate = ratio)
  } else if (!"estimate" %in% names(df)) {
    stop("Expected a column named 'estimate' or 'ratio' in pairs(emm).")
  }
  
  df |>
    dplyr::mutate(
      estimate = round(estimate, 3),
      SE = round(SE, 3),
      z.ratio = round(z.ratio, 2),
      p.value = round(p.value, 4),
      sig = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE ~ ""
      )
    ) |>
    dplyr::select(
      dplyr::any_of(c("Block", "Genotype")),
      contrast, estimate, SE, z.ratio, p.value, sig
    )
}

pretty_blocks <- function(x) {
  df <- as.data.frame(x)
  
  # Try to recover Genotype info if emmeans object was by = "Genotype"
  if (!"Genotype" %in% names(df)) {
    by_vars <- attr(x, "by.vars")
    if (!is.null(by_vars) && "Genotype" %in% by_vars) {
      df$Genotype <- rep(x@by[[1]], each = nrow(df) / length(x@by[[1]]))
    } else if (!is.null(x@grid$Genotype)) {
      df$Genotype <- rep(unique(x@grid$Genotype), each = nrow(df) / length(unique(x@grid$Genotype)))
    } else {
      df$Genotype <- NA
    }
  }
  
  # --- normalize estimate column name ---
  if ("odds.ratio" %in% names(df)) {
    df <- dplyr::rename(df, estimate = odds.ratio)
  } else if ("ratio" %in% names(df)) {
    df <- dplyr::rename(df, estimate = ratio)
  } else if (!"estimate" %in% names(df)) {
    stop("No column named 'estimate', 'ratio', or 'odds.ratio' found in the input.")
  }
  
  # --- rounding and significance stars ---
  df |>
    dplyr::mutate(
      estimate = round(estimate, 3),
      SE = round(SE, 3),
      z.ratio = round(z.ratio, 2),
      p.value = round(as.numeric(p.value), 4),
      sig = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE ~ ""
      )
    ) |>
    dplyr::select(Genotype, contrast, estimate, SE, z.ratio, p.value, sig)
}


# ------------------------------------------------------------------------------
# OVERALL RESPONSIVNESS
# ------------------------------------------------------------------------------
# PEAK Distance
# ------------------------------------------------------------------------------
# Block 1 (stimuli 1–477)
emm_peak_block1 <- emmeans(
  m_peak,
  ~ Genotype,
  at = list(stimulus_log = log(1:477), Block = "1"),
  cov.reduce = mean,
  type = "response"
)

# Block 2 (stimuli 1–9)
emm_peak_block2 <- emmeans(
  m_peak,
  ~ Genotype,
  at = list(stimulus_log = log(1:9), Block = "2"),
  cov.reduce = mean,
  type = "response"
)

# First Dark Flashes
emm_peak_block1_first_stim <- emmeans(
  m_peak,
  ~ Genotype,
  at = list(stimulus_log = log(1:10), Block = "1"),   
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)

# Middle Dark Flashes
emm_peak_block1_middle_stim <- emmeans(
  m_peak,
  ~ Genotype,
  at = list(stimulus_log = log(230:240), Block = "1"),   
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)

# Last Dark Flashes
emm_peak_block1_last_stim <- emmeans(
  m_peak,
  ~ Genotype,
  at = list(stimulus_log = log(466:477), Block = "1"),  
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)


# ------------------------------------------------------------------------------
# Habituation slopes
emm_peak_slopes <- emtrends(m_peak,
                       ~ Genotype | Block,
                       var = "stimulus_log")

# ------------------------------------------------------------------------------
# BLOCK 1 vs 2

# First Stimuli in Block 1 vs Block 2
emm_peak_between_blocks_first_stimulus <- emmeans(
  m_peak, ~ Block | Genotype, 
  at = list(stimulus_log = log(1)), 
  type = "response" 
)

# --- Compare Blocks per Genotype (recovery test) --- 
# The Last 9 Stimuli in Block 1 vs. The First 9 Stimuli in Block 2
emm_peak_block1_recovery <- emmeans( 
  m_peak, ~ Block | Genotype, 
  at = list(stimulus_log = log(466:477), Block="1"), 
  type = "response" ) 

emm_peak_block2_recovery <- emmeans( 
  m_peak, ~ Block | Genotype, 
  at = list(stimulus_log = log(1:9), Block="2"), 
  type = "response" ) 

combined_peak_emms <- rbind(emm_peak_block1_recovery, emm_peak_block2_recovery) 
# Use the contrast function to compare Block levels within each Genotype level. +
# The 'by = "Genotype"' argument ensures the comparisons are done separately for each genotype. 
block_peak_comparisons <- contrast(
  combined_peak_emms, 
  method = "revpairwise",
  by = "Genotype", 
  adjust = "Tukey") # You can choose your desired p-value adjustment 

# ------------------------------------------------------------------------------
# The Last Stim in Block 1 vs. The First Stim in Block 2
emm_peak_block1_last_recovery <- emmeans( 
  m_peak, ~ Block | Genotype, 
  at = list(stimulus_log = log(477), Block="1"), 
  type = "response" ) 

emm_peak_block2_first_recovery <- emmeans( 
  m_peak, ~ Block | Genotype, 
  at = list(stimulus_log = log(1), Block="2"), 
  type = "response" ) 

combined_peak_last_first_emms <- rbind(emm_peak_block1_last_recovery, emm_peak_block2_first_recovery) 
# Use the contrast function to compare Block levels within each Genotype level. +
# The 'by = "Genotype"' argument ensures the comparisons are done separately for each genotype. 
block_peak_last_first_comparisons <- contrast(
  combined_peak_last_first_emms, 
  method = "revpairwise",
  by = "Genotype", 
  adjust = "Tukey") # You can choose your desired p-value adjustment 

# ------------------------------------------------------------------------------
# Habituation slope between blocks
emm_peak_between_blocks_slopes <- emtrends(m_peak,
                                      ~ Block | Genotype,
                                      var = "stimulus_log")

# ------------------------------------------------------------------------------
# Create nice data frames
# Create a list of all results
# Define a named list of descriptions for each results section
descriptions <- list(
  responsiveness_block1 = "Pairwise contrasts of responsiveness within Block 1.",
  responsiveness_block2 = "Pairwise contrasts of responsiveness within Block 2.",
  responsiveness_block1_first = "Responsiveness contrasts within Block 1 (first stimulus subset).",
  responsiveness_block1_middle = "Responsiveness contrasts within Block 1 (middle stimuli subset).",
  responsiveness_block1_last = "Responsiveness contrasts within Block 1 (last stimulus subset).",
  habituation_slope_block1_and_block2 = "Comparisons of habituation slopes between Block 1 and Block 2.",
  between_blocks_responsiveness_first = "Between-block comparison for the first stimulus.",
  between_blocks_habituation_slope = "Between-block comparison for habituation slopes.",
  between_block_10_stimuli = "Block-wise comparison: 9 last stimuli of Block 1 vs. 9 first stimuli of Block 2.",
  between_block_last_first = "Comparison between last stimulus of Block 1 and first stimulus of Block 2."
)

info <- "
Interpretation for estimate (odds ratio):\n
A ratio of 1.00 → no difference between genotypes.\n
A ratio < 1.00 → the first genotype (left of “/”) has a lower predicted response than the second.\n
A ratio > 1.00 → the first genotype has a higher predicted response.\n
Habituation Rate (Slope): \n
Interpretation: Negative slope -> faster habituation.\n\n
"

results <- list(
  responsiveness_block1 = pretty_pairs(emm_peak_block1),
  responsiveness_block2 = pretty_pairs(emm_peak_block2 ),
  responsiveness_block1_first = pretty_pairs(emm_peak_block1_first_stim),
  responsiveness_block1_middle = pretty_pairs(emm_peak_block1_middle_stim ),
  responsiveness_block1_last = pretty_pairs(emm_peak_block1_last_stim),
  habituation_slope_block1_and_block2 = pretty_pairs(emm_peak_slopes),
  between_blocks_responsiveness_first = pretty_pairs(emm_peak_between_blocks_first_stimulus),
  between_blocks_habituation_slope = pretty_pairs(emm_peak_between_blocks_slopes),
  between_block_10_stimuli = pretty_blocks(block_peak_comparisons),
  between_block_last_first = pretty_blocks(block_peak_last_first_comparisons)
)


# Define output path
outfile <- "D:/WorkingData/Susana/results/massed/glmm_peak_distance_comparisons.txt"

# Write everything to the file
sink(outfile)

cat("### GLMM Peak Distance Moved (Massed) ###\n\n")
print(summary(m_peak), row.names = FALSE)
cat("\n")
print(n_per_genotype)
cat("\n\n")

cat("### Estimated Marginal Means Comparisons ###\n\n")
cat(info)

for (nm in names(results)) {
  cat("----", nm, "----\n")
  if (!is.null(descriptions[[nm]])) {
    cat(descriptions[[nm]], "\n\n")
  }
  print(results[[nm]], row.names = FALSE)
  cat("\n\n")
}

sink()
cat("All results written to:", outfile, "\n")


# Summed Distance
# ------------------------------------------------------------------------------
# Block 1 (stimuli 1–477)
emm_sum_block1 <- emmeans(
  m_sum,
  ~ Genotype,
  at = list(stimulus_log = log(1:477), Block = "1"),
  cov.reduce = mean,
  type = "response"
)

# Block 2 (stimuli 1–9)
emm_sum_block2 <- emmeans(
  m_sum,
  ~ Genotype,
  at = list(stimulus_log = log(1:9), Block = "2"),
  cov.reduce = mean,
  type = "response"
)

# First Dark Flashes
emm_sum_block1_first_stim <- emmeans(
  m_sum,
  ~ Genotype,
  at = list(stimulus_log = log(1:10), Block = "1"),   
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)

# Middle Dark Flashes
emm_sum_block1_middle_stim <- emmeans(
  m_sum,
  ~ Genotype,
  at = list(stimulus_log = log(230:240), Block = "1"),   
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)

# Last Dark Flashes
emm_sum_block1_last_stim <- emmeans(
  m_sum,
  ~ Genotype,
  at = list(stimulus_log = log(466:477), Block = "1"),  
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)


# ------------------------------------------------------------------------------
# Habituation slopes
emm_sum_slopes <- emtrends(m_sum,
                            ~ Genotype | Block,
                            var = "stimulus_log")

# ------------------------------------------------------------------------------
# BLOCK 1 vs 2

# First Stimuli in Block 1 vs Block 2
emm_sum_between_blocks_first_stimulus <- emmeans(
  m_sum, ~ Block | Genotype, 
  at = list(stimulus_log = log(1)), 
  type = "response" 
)

# --- Compare Blocks per Genotype (recovery test) --- 
# The Last 9 Stimuli in Block 1 vs. The First 9 Stimuli in Block 2
emm_sum_block1_recovery <- emmeans( 
  m_sum, ~ Block | Genotype, 
  at = list(stimulus_log = log(466:477), Block="1"), 
  type = "response" ) 

emm_sum_block2_recovery <- emmeans( 
  m_sum, ~ Block | Genotype, 
  at = list(stimulus_log = log(1:9), Block="2"), 
  type = "response" ) 

combined_sum_emms <- rbind(emm_sum_block1_recovery, emm_sum_block2_recovery) 
# Use the contrast function to compare Block levels within each Genotype level. +
# The 'by = "Genotype"' argument ensures the comparisons are done separately for each genotype. 
block_sum_comparisons <- contrast(
  combined_sum_emms, 
  method = "revpairwise",
  by = "Genotype", 
  adjust = "Tukey") # You can choose your desired p-value adjustment 

# ------------------------------------------------------------------------------
# The Last Stim in Block 1 vs. The First Stim in Block 2
emm_sum_block1_last_recovery <- emmeans( 
  m_sum, ~ Block | Genotype, 
  at = list(stimulus_log = log(477), Block="1"), 
  type = "response" ) 

emm_sum_block2_first_recovery <- emmeans( 
  m_sum, ~ Block | Genotype, 
  at = list(stimulus_log = log(1), Block="2"), 
  type = "response" ) 

combined_sum_last_first_emms <- rbind(emm_sum_block1_last_recovery, emm_sum_block2_first_recovery) 
# Use the contrast function to compare Block levels within each Genotype level. +
# The 'by = "Genotype"' argument ensures the comparisons are done separately for each genotype. 
block_sum_last_first_comparisons <- contrast(
  combined_sum_last_first_emms, 
  method = "revpairwise",
  by = "Genotype", 
  adjust = "Tukey") # You can choose your desired p-value adjustment 

# ------------------------------------------------------------------------------
# Habituation slope between blocks
emm_sum_between_blocks_slopes <- emtrends(m_sum,
                                           ~ Block | Genotype,
                                           var = "stimulus_log")

# ------------------------------------------------------------------------------
# Create nice data frames

results <- list(
  responsiveness_block1 = pretty_pairs(emm_sum_block1),
  responsiveness_block2 = pretty_pairs(emm_sum_block2 ),
  responsiveness_block1_first = pretty_pairs(emm_sum_block1_first_stim),
  responsiveness_block1_middle = pretty_pairs(emm_sum_block1_middle_stim ),
  responsiveness_block1_last = pretty_pairs(emm_sum_block1_last_stim),
  habituation_slope_block1_and_block2 = pretty_pairs(emm_sum_slopes),
  between_blocks_responsiveness_first = pretty_pairs(emm_sum_between_blocks_first_stimulus),
  between_blocks_habituation_slope = pretty_pairs(emm_sum_between_blocks_slopes),
  between_block_10_stimuli = pretty_blocks(block_sum_comparisons),
  between_block_last_first = pretty_blocks(block_sum_last_first_comparisons)
)


# Define output path
outfile <- "D:/WorkingData/Susana/results/massed/glmm_summed_distance_comparisons.txt"

# Write everything to the file
sink(outfile)

cat("### GLMM Summed Distance Moved (Massed) ###\n\n")
print(summary(m_sum), row.names = FALSE)
cat("\n")
print(n_per_genotype)
cat("\n\n")

cat("### Estimated Marginal Means Comparisons ###\n\n")
cat(info)
for (nm in names(results)) {
  cat("----", nm, "----\n")
  if (!is.null(descriptions[[nm]])) {
    cat(descriptions[[nm]], "\n\n")
  }
  print(results[[nm]], row.names = FALSE)
  cat("\n\n")
}

sink()
cat("All results written to:", outfile, "\n")


# ------------------------------------------------------------------------------
# Response Delay
# ------------------------------------------------------------------------------
# Block 1 (stimuli 1–477)
emm_delay_block1 <- emmeans(
  m_delay_gaussian,
  ~ Genotype,
  at = list(stimulus_log = log(1:477), Block = "1"),
  cov.reduce = mean,
  type = "response"
)

# Block 2 (stimuli 1–9)
emm_delay_block2 <- emmeans(
  m_delay_gaussian,
  ~ Genotype,
  at = list(stimulus_log = log(1:9), Block = "2"),
  cov.reduce = mean,
  type = "response"
)

# First Dark Flashes
emm_delay_block1_first_stim <- emmeans(
  m_delay_gaussian,
  ~ Genotype,
  at = list(stimulus_log = log(1:10), Block = "1"),   
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)

# Middle Dark Flashes
emm_delay_block1_middle_stim <- emmeans(
  m_delay_gaussian,
  ~ Genotype,
  at = list(stimulus_log = log(230:240), Block = "1"),   
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)

# Last Dark Flashes
emm_delay_block1_last_stim <- emmeans(
  m_delay_gaussian,
  ~ Genotype,
  at = list(stimulus_log = log(466:477), Block = "1"),  
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)


# ------------------------------------------------------------------------------
# Habituation slopes
emm_delay_slopes <- emtrends(m_delay_gaussian,
                            ~ Genotype | Block,
                            var = "stimulus")

# ------------------------------------------------------------------------------
# BLOCK 1 vs 2

# First Stimuli in Block 1 vs Block 2
emm_delay_between_blocks_first_stimulus <- emmeans(
  m_delay_gaussian, ~ Block | Genotype, 
  at = list(stimulus_log = log(1)), 
  type = "response" 
)

# --- Compare Blocks per Genotype (recovery test) --- 
# The Last 9 Stimuli in Block 1 vs. The First 9 Stimuli in Block 2
emm_delay_block1_recovery <- emmeans( 
  m_delay_gaussian, ~ Block | Genotype, 
  at = list(stimulus_log = log(466:477), Block="1"), 
  type = "response" ) 

emm_delay_block2_recovery <- emmeans( 
  m_delay_gaussian, ~ Block | Genotype, 
  at = list(stimulus_log = log(1:9), Block="2"), 
  type = "response" ) 

combined_delay_emms <- rbind(emm_delay_block1_recovery, emm_delay_block2_recovery) 
# Use the contrast function to compare Block levels within each Genotype level. +
# The 'by = "Genotype"' argument ensures the comparisons are done separately for each genotype. 
block_delay_comparisons <- contrast(
  combined_delay_emms, 
  method = "revpairwise",
  by = "Genotype", 
  adjust = "Tukey") # You can choose your desired p-value adjustment 

# ------------------------------------------------------------------------------
# The Last Stim in Block 1 vs. The First Stim in Block 2
emm_delay_block1_last_recovery <- emmeans( 
  m_delay_gaussian, ~ Block | Genotype, 
  at = list(stimulus_log = log(477), Block="1"), 
  type = "response" ) 

emm_delay_block2_first_recovery <- emmeans( 
  m_delay_gaussian, ~ Block | Genotype, 
  at = list(stimulus_log = log(1), Block="2"), 
  type = "response" ) 

combined_delay_last_first_emms <- rbind(emm_delay_block1_last_recovery, emm_delay_block2_first_recovery) 
# Use the contrast function to compare Block levels within each Genotype level. +
# The 'by = "Genotype"' argument ensures the comparisons are done separately for each genotype. 
block_delay_last_first_comparisons <- contrast(
  combined_delay_last_first_emms, 
  method = "revpairwise",
  by = "Genotype", 
  adjust = "Tukey") # You can choose your desired p-value adjustment 

# ------------------------------------------------------------------------------
# Habituation slope between blocks
emm_delay_between_blocks_slopes <- emtrends(m_delay_gaussian,
                                           ~ Block | Genotype,
                                           var = "stimulus")

# ------------------------------------------------------------------------------
# Create nice data frames
results <- list(
  responsiveness_block1 = pretty_pairs(emm_delay_block1),
  responsiveness_block2 = pretty_pairs(emm_delay_block2 ),
  responsiveness_block1_first = pretty_pairs(emm_delay_block1_first_stim),
  responsiveness_block1_middle = pretty_pairs(emm_delay_block1_middle_stim ),
  responsiveness_block1_last = pretty_pairs(emm_delay_block1_last_stim),
  habituation_slope_block1_and_block2 = pretty_pairs(emm_delay_slopes),
  between_blocks_responsiveness_first = pretty_pairs(emm_delay_between_blocks_first_stimulus),
  between_blocks_habituation_slope = pretty_pairs(emm_delay_between_blocks_slopes),
  between_block_10_stimuli = pretty_blocks(block_delay_comparisons),
  between_block_last_first = pretty_blocks(block_delay_last_first_comparisons)
)


# Define output path
outfile <- "D:/WorkingData/Susana/results/massed/glmm_delay_distance_comparisons.txt"

# Write everything to the file
sink(outfile)

cat("### GLMM Response Delay (Massed) ###\n\n")
print(summary(m_delay_gaussian), row.names = FALSE)
cat("\n")
print(n_per_genotype)
cat("\n\n")

cat("### Estimated Marginal Means Comparisons ###\n\n")
cat(info)
for (nm in names(results)) {
  cat("----", nm, "----\n")
  if (!is.null(descriptions[[nm]])) {
    cat(descriptions[[nm]], "\n\n")
  }
  print(results[[nm]], row.names = FALSE)
  cat("\n\n")
}

sink()
cat("All results written to:", outfile, "\n")

# ==============================================================================
# ==============================================================================
# TEST PLOTS
# --- 1. Filter prediction data to ABTL only ---
pred_abtl <- new_data_delay %>%
  filter(Genotype == "ABTL")

# --- 2. Compute mean observed delay per stimulus (and Block) ---
raw_means_abtl <- df_final_sub %>%
  filter(Genotype == "ABTL") %>%
  group_by(Block, stimulus) %>%
  summarise(
    mean_delay = mean(delay, na.rm = TRUE),
    se_delay   = sd(delay, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# --- 3. Plot ---
ggplot(pred_abtl, aes(x = stimulus)) +
  facet_wrap(
    ~Block, ncol = 2, scales = "free_x",
    labeller = as_labeller(c(
      `1` = "Block 1: 477 flashes",
      `2` = "Block 2: 8 flashes"
    ))
  ) +
  # Observed means (points only)
  geom_point(
    data = raw_means_abtl,
    aes(y = mean_delay),
    color = "grey30",
    size = 2,
    alpha = 0.8
  ) +
  # Predicted fit line and CI ribbon
  geom_line(aes(y = fit), color = "#0072B2", size = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              fill = "#0072B2", alpha = 0.1, color = NA) +
  labs(
    x = "Stimulus number (within block)",
    y = "Delay (s)",
    title = "Predicted and observed mean delay for ABTL"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )


# Delays Violin Plots
ggplot(df_final_sub, aes(x = as.factor(delay), y = stimulus, fill = as.factor(delay))) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "grey30") +
  facet_wrap(
    ~Block, ncol = 2, scales = "free_y",
    labeller = as_labeller(c(
      `1` = "Block 1: 477 flashes",
      `2` = "Block 2: 8 flashes"
    ))
  ) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1, name = "Delay (s)") +
  labs(
    x = "Delay category (s)",
    y = "Stimulus number",
    title = "Stimulus distribution per delay category and block"
  ) +
  theme_pubr(base_size = 14)
