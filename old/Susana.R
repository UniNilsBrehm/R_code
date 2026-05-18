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


# Load data from csv file
df <- read_csv("D:/WorkingData/Susana/SPZ_Massed_Training_Nov2025.csv")
# df <- read_csv("D:/WorkingData/Susana/SPZ_Spaced_Training_Nov2025.csv")
# response = max_peak
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
  filter(Peak == 0)

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
ggplot(df_final, aes(x = max_peak)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_grid(col=vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Distribution",
    x = "X",
    y = "Count"
  )

# df_final_sub$max_peak_log <- log(df_final_sub$max_peak)

ggplot(df_final_sub, aes(x = max_peak)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_grid(col=vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Distribution",
    x = "X",
    y = "Count"
  )

ggplot(df_final, aes(x = move)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_grid(col=vars(Genotype)) +
  theme_minimal() +
  labs(
    title = "Distribution",
    x = "X",
    y = "Count"
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


# Habituation slopes
emm_slopes <- emtrends(m1,
                       ~ Genotype | Block,
                       var = "stimulus_log")

# ------------------------------------------------------------------------------
# BLOCK 1 vs 2
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
pretty_pairs(emm_block1)
pretty_pairs(emm_block2)
pretty_pairs(emm_block1_first_stim)
pretty_pairs(emm_block1_middle_stim)
pretty_pairs(emm_block1_last_stim)

pretty_pairs(emm_slopes)

pretty_pairs(emm_between_blocks_tap1)
pretty_pairs(emm_between_blocks_slopes)
pretty_blocks(block_comparisons)
pretty_blocks(block_last_first_comparisons)

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
ggplot(new_data, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_wrap(~Block, ncol = 2, scales = "free_x",
             labeller = as_labeller(c(`1` = "Block 1: 477 flashes",
                                      `2` = "Block 2: 8 flashes"))) +
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

# ==============================================================================
# ==============================================================================
# 2. Magnitude of movement (max_peak)
# df_final$max_peak_gamma <- df_final$max_peak + 0.0001
# df_final_sub$max_peak_gamma <- df_final_sub$max_peak + 0.0001
# df_final_sub$stimulus_log_c <- scale(df_final_sub$stimulus_log, center = TRUE, scale = TRUE)

m_peak <- glmmTMB(
  max_peak  ~ Genotype * stimulus_log  * Block + (1|Video/Well),
  # dispformula = ~ Genotype,
  family = Gamma(link = "log"),
  data = df_final_sub
)

m_sum <- glmmTMB(
  max_cumsum  ~ Genotype * stimulus_log  * Block + (1|Video/Well),
  # dispformula = ~ Genotype,
  family = Gamma(link = "log"),
  data = df_final_sub
)

summary(m_peak)
summary(m_sum)

# Validation
res_m_peak <- simulateResiduals(m_peak)
plot(res_m_peak)
testDispersion(res_m_peak)
testZeroInflation(res_m_peak)
plotResiduals(res_m_peak, df_final_sub$Genotype)
plotResiduals(res_m_peak, df_final_sub$Block)

res_m_sum <- simulateResiduals(m_sum)
plot(res_m_sum)
testDispersion(res_m_sum)
testZeroInflation(res_m_sum)
plotResiduals(res_m_sum, df_final_sub$Genotype)
plotResiduals(res_m_sum, df_final_sub$Block)

model_performance(m_peak)  # AIC, R2, RMSE, ICC, etc.
check_collinearity(m_peak) # Multicollinearity (VIFs)
model_performance(m_sum)  # AIC, R2, RMSE, ICC, etc.
check_collinearity(m_sum) # Multicollinearity (VIFs)
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
ggplot(new_data_2, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_wrap(
    ~Block, ncol = 2, scales = "free",
    labeller = as_labeller(c(
      `1` = "Block 1: 477 flashes",
      `2` = "Block 2: 8 flashes"
    ))
  ) +
  geom_line(aes(y = fit), size = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.1, color = NA) +
  labs(
    x = "Stimulus number (within block)",
    y = "Predicted movement peak",
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
ggplot(new_data_2, aes(x = stimulus, color = Genotype, fill = Genotype)) +
  facet_wrap(
    ~Block, ncol = 2, scales = "free",
    labeller = as_labeller(c(
      `1` = "Block 1: 477 flashes",
      `2` = "Block 2: 8 flashes"
    ))
  ) +
  geom_line(aes(y = fit), size = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.1, color = NA) +
  labs(
    x = "Stimulus number (within block)",
    y = "Predicted Movement Summed",
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
# All taps (mean y shift of curve): mean response probability across the entire block
emm_overall_2 <- emmeans(
  m2,
  ~ Genotype | Block,
  at = list(stimulus_log = log(1:477)),   # all taps within block
  cov.reduce = mean,                         # average over those predictions
  type = "response"
)

# Habituation slopes
emm_slopes_2 <- emtrends(m2,
                       ~ Genotype | Block,
                       var = "stimulus_log",
                       type = "response")

# BLOCK 1 vs 2
# Overall responsiveness between blocks: All taps (mean y shift of curve): mean response probability across the entire block
emm_block_overall_2 <- emmeans(m2,
                             ~ Block | Genotype,
                             at = list(stimulus_log= log(1:477)),   # all taps within block
                             cov.reduce = mean,                         # average over those predictions
                             type = "response")


# Habituation slope between blocks
emm_block_slopes_2 <- emtrends(m2,
                             ~ Block | Genotype,
                             var = "stimulus_log",
                             type = "response")

# --- Block 1: last 9 stimuli ---
emm_last_block1 <- emmeans(
  m_peak,
  ~ Genotype,
  at = list(stimulus_log = log(469:477), Block = "1"),
  cov.reduce = mean,
  type = "response"
)

# --- Block 2: all 9 stimuli ---
emm_block2 <- emmeans(
  m_peak,
  ~ Genotype,
  at = list(stimulus_log = log(1:9), Block = "2"),
  cov.reduce = mean,
  type = "response"
)

# --- Combine into a single comparison ---
emm_recovery <- rbind(
  as.data.frame(emm_last_block1) %>% mutate(Block = "1_last9"),
  as.data.frame(emm_block2) %>% mutate(Block = "2_first9")
)

# --- Compare Blocks per Genotype (recovery test) ---
contrast(
  emmeans(
    m_peak,
    ~ Block | Genotype,
    at = list(stimulus_log = c(log(469:477), log(1:9))),
    cov.reduce = mean,
    type = "response"
  ),
  method = "revpairwise",
  adjust = "tukey"
)




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

