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
library(ggpubr)


# Load data from csv file
df <- read_csv("D:/WorkingData/Susana/SPZ_Massed_Training_7Nov2025.csv")
df_reduced <- df[, c("Block", "Well", "Video", "Peak", "Genotype", "Stimulus_New", "max_peak", "max_cumsum")]

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

# Get number of animals
n_per_genotype <- df_final %>%
  dplyr::distinct(Video, Well, Genotype) %>%
  dplyr::count(Genotype)

# Histograms
h1 <- ggplot(df_final, aes(x = move)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_grid(col=vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Distribution",
    x = "Response (No/Yes)",
    y = "Count"
  )

ggsave(
  filename = "D:/WorkingData/Susana/results/massed/figures/response_prob/response_binary_dist.pdf",
  plot = h1,
  width = 6, height = 4
)

# ==============================================================================
# 1. Probability of movement
m1 <- glmmTMB(
  #move ~ Genotype * Block * stimulus_log + (1 | Video/Well) +     
  #(1 + stimulus_log | Video/Well),  # random intercept + slope,
  move ~ Genotype * Block * stimulus_log + (1 | Video/Well),
  family = binomial(link = "logit"), 
  # family = binomial(link = "probit"), 
  data = df_final,
  # control = glmmTMBControl(optimizer = "nloptwrap")  # Default
  # control = glmmTMBControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))  # Bound optimization by quadratic approximation
  # control = glmmTMBControl(optimizer = "Nelder_Mead")
)

# Validation
res_m1 <- simulateResiduals(m1)
plot(res_m1)
plotResiduals(res_m1, df_final$Genotype)
plotResiduals(res_m1, df_final$Block)
plotResiduals(res_m1, df_final$stimulus_log)

testUniformity(res_m1) 
testOutliers(res_m1, type = "bootstrap")
testDispersion(res_m1)
testQuantiles(res_m1)
testCategorical(res_m1, catPred = df_final$Genotype)
testCategorical(res_m1, catPred = df_final$Block)
testZeroInflation(res_m1)

model_performance(m1)  # AIC, R2, RMSE, ICC, etc.
check_collinearity(m1) # Multicollinearity (VIFs)
# check_collinearity(update(m2, . ~ Genotype + Block + stimulus_log))

# ------------------------------------------------------------------------------
# Plot Habituation Curves
# --- 1. Prediction grid ------------------------------------
new_data <- bind_rows(
  expand.grid(Block = "1", stimulus = 1:477),
  expand.grid(Block = "2", stimulus = 1:9)
) %>%
  tidyr::crossing(Genotype = unique(df_final$Genotype)) %>%
  mutate(stimulus_log = log(stimulus))

# --- 2. Model predictions ----------------------------------
pred <- predict(m1, newdata = new_data, re.form = NA, se.fit = TRUE)

new_data <- new_data %>%
  mutate(
    fit     = plogis(pred$fit),
    CI_low  = plogis(pred$fit - 1.96 * pred$se.fit),
    CI_high = plogis(pred$fit + 1.96 * pred$se.fit)
  )

# --- 3. Plot -----------------------------------------------
g <- ggplot(new_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_wrap(~Block, ncol = 2, scales = "free_x",
             labeller = as_labeller(c(`1` = "Block 1: 477 flashes",
                                      `2` = "Block 2: 8 flashes"))) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.1, color = NA) +
  expand_limits(y = c(0, 1)) +
  labs(
    x = "Stimulus number (within block)",
    y = "Response probability",
    color = "Genotype",
    fill  = "Genotype"
  ) +
  theme_pubr(base_size = 14)

g
ggsave(
  filename = "D:/WorkingData/Susana/results/massed/figures/response_prob/response_prob_habituation_curves.pdf",
  plot = g,
  width = 6, height = 4
)

# Plot only taps 1 to 9
new_data_short <- new_data %>% dplyr::filter(stimulus <= 9)

ggplot(new_data_short, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_wrap(
    ~Block, ncol = 2, scales = "free_x",
    labeller = as_labeller(c(`1` = "Block 1: first 9 flashes",
                             `2` = "Block 2: 9 flashes"))
  ) +
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.1, color = NA) +
  expand_limits(y = c(0, 1)) +
  labs(
    x = "Stimulus number (within block)",
    y = "Predicted probability of movement",
    color = "Genotype",
    fill  = "Genotype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

# ------------------------------------------------------------------------------
# Estimated Marginal Means
# Helper Function
pretty_pairs <- function(emm) {
  df <- as.data.frame(pairs(emm))
  
  # Normalize column names
  if ("odds.ratio" %in% names(df)) {
    df <- dplyr::rename(df, estimate = odds.ratio)
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
  # Convert to data.frame and add significance stars
  df <- as.data.frame(x)
  
  # Extract Genotype info if nested contrast object
  if ("Genotype" %in% names(attributes(x))) {
    df$Genotype <- rep(attr(x, "by.vars")[[1]], each = nrow(df) / length(attr(x, "by.vars")[[1]]))
  } else if (!"Genotype" %in% names(df)) {
    df$Genotype <- rep(unique(x@grid$Genotype), each = nrow(df) / length(unique(x@grid$Genotype)))
  }
  
  # Add significance stars
  df <- df %>%
    dplyr::mutate(
      sig = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(Genotype, contrast, odds.ratio, SE, z.ratio, p.value, sig)
  
  df
}

# ------------------------------------------------------------------------------
# OVERALL RESPONSIVNESS
# All taps (mean y shift of curve): mean response probability across the entire block
# Block 1 (stimuli 1–477)
emm_block1 <- emmeans(
  m1,
  ~ Genotype,
  at = list(stimulus_log = log(1:477), Block = "1"),
  cov.reduce = mean,
  type = "response"
)

# Block 2 (stimuli 1–9)
emm_block2 <- emmeans(
  m1,
  ~ Genotype,
  at = list(stimulus_log = log(1:9), Block = "2"),
  cov.reduce = mean,
  type = "response"
)

# First Dark Flashes
emm_block1_first_stim <- emmeans(
  m1,
  ~ Genotype,
  at = list(stimulus_log = log(1:10), Block = "1"),   
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)

# Middle Dark Flashes
emm_block1_middle_stim <- emmeans(
  m1,
  ~ Genotype,
  at = list(stimulus_log = log(230:240), Block = "1"),   
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)

# Last Dark Flashes
emm_block1_last_stim <- emmeans(
  m1,
  ~ Genotype,
  at = list(stimulus_log = log(466:477), Block = "1"),  
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)


# ------------------------------------------------------------------------------
# Habituation slopes
emm_slopes <- emtrends(m1,
                       ~ Genotype | Block,
                       var = "stimulus_log")

# ------------------------------------------------------------------------------
# BLOCK 1 vs 2

# First Stimuli in Block 1 vs Block 2
emm_between_blocks_first_stimulus <- emmeans(
  m1, ~ Block | Genotype, 
  at = list(stimulus_log = log(1)), 
  type = "response" 
)

# --- Compare Blocks per Genotype (recovery test) --- 
# The Last 9 Stimuli in Block 1 vs. The First 9 Stimuli in Block 2
emm_block1_recovery <- emmeans( 
  m1, ~ Block | Genotype, 
  at = list(stimulus_log = log(466:477), Block="1"), 
  type = "response" ) 

emm_block2_recovery <- emmeans( 
  m1, ~ Block | Genotype, 
  at = list(stimulus_log = log(1:9), Block="2"), 
  type = "response" ) 

combined_emms <- rbind(emm_block1_recovery, emm_block2_recovery) 
# Use the contrast function to compare Block levels within each Genotype level. +
# The 'by = "Genotype"' argument ensures the comparisons are done separately for each genotype. 
block_comparisons <- contrast(
  combined_emms, 
  method = "revpairwise",
  by = "Genotype", 
  adjust = "Tukey") # You can choose your desired p-value adjustment 

# ------------------------------------------------------------------------------
# The Last Stim in Block 1 vs. The First Stim in Block 2
emm_block1_last_recovery <- emmeans( 
  m1, ~ Block | Genotype, 
  at = list(stimulus_log = log(477), Block="1"), 
  type = "response" ) 

emm_block2_first_recovery <- emmeans( 
  m1, ~ Block | Genotype, 
  at = list(stimulus_log = log(1), Block="2"), 
  type = "response" ) 

combined_last_first_emms <- rbind(emm_block1_last_recovery, emm_block2_first_recovery) 
# Use the contrast function to compare Block levels within each Genotype level. +
# The 'by = "Genotype"' argument ensures the comparisons are done separately for each genotype. 
block_last_first_comparisons <- contrast(
  combined_last_first_emms, 
  method = "revpairwise",
  by = "Genotype", 
  adjust = "Tukey") # You can choose your desired p-value adjustment 

# ------------------------------------------------------------------------------
# Habituation slope between blocks
emm_between_blocks_slopes <- emtrends(m1,
                                      ~ Block | Genotype,
                                      var = "stimulus_log")

# ------------------------------------------------------------------------------
# Create nice data frames
# Create a list of all results
results <- list(
  responsiveness_block1 = pretty_pairs(emm_block1),
  responsiveness_block2 = pretty_pairs(emm_block2),
  responsiveness_block1_first = pretty_pairs(emm_block1_first_stim),
  responsiveness_block1_middle = pretty_pairs(emm_block1_middle_stim),
  responsiveness_block1_last = pretty_pairs(emm_block1_last_stim),
  habituation_slope_block1_and_block2 = pretty_pairs(emm_slopes),
  between_blocks_responsiveness_first = pretty_pairs(emm_between_blocks_first_stimulus),
  between_blocks_habituation_slope = pretty_pairs(emm_between_blocks_slopes),
  between_block_10_stimuli = pretty_blocks(block_comparisons),
  between_block_last_first = pretty_blocks(block_last_first_comparisons)
)

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

# Define output path
outfile <- "D:/WorkingData/Susana/results/massed/glmm_response_prob_comparisons.txt"

# Write everything to the file
sink(outfile)

cat("### GLMM Response Probability (Massed) ###\n\n")
print(summary(m1), row.names = FALSE)
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
# Predict stuff (extrapolate)
# Create prediction grid: same range for both blocks
new_data_extrap <- bind_rows(
  expand.grid(Block = "1", stimulus = 1:477),
  expand.grid(Block = "2", stimulus = 1:477)   # extend Block 2 artificially
) %>%
  tidyr::crossing(Genotype = unique(df_final$Genotype)) %>%
  mutate(stimulus_log = log(stimulus))

# Predict from model
pred <- predict(m1, newdata = new_data_extrap, re.form = NA, se.fit = TRUE)

new_data_extrap <- new_data_extrap %>%
  mutate(
    fit     = plogis(pred$fit),
    CI_low  = plogis(pred$fit - 1.96 * pred$se.fit),
    CI_high = plogis(pred$fit + 1.96 * pred$se.fit),
    extrapolated = ifelse(Block == "2" & stimulus > 9, TRUE, FALSE)
  )

ggplot(new_data_extrap, aes(x = stimulus, y = fit, color = Genotype)) +
  facet_wrap(
    ~Block,
    ncol = 2,
    labeller = as_labeller(c(`1` = "Block 1 (observed)", `2` = "Block 2 (partial extrapolation)"))
  ) +
  geom_line(aes(linetype = extrapolated), linewidth = 1.1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = Genotype), alpha = 0.1, color = NA) +
  scale_linetype_manual(
    values = c("FALSE" = "solid", "TRUE" = "dotted"),  # ← change "dotted" here
    guide = "none"
  ) +
  expand_limits(y = c(0, 1)) +
  labs(
    x = "Stimulus number (within block)",
    y = "Predicted probability of movement",
    color = "Genotype",
    fill  = "Genotype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )
