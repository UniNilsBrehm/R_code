# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
library(readr)        # Reading CSV files
library(DHARMa)       # Residual diagnostics for (G)LMMs
library(emmeans)      # Estimated marginal means (EMMs) and contrasts
library(glmmTMB)      # Generalized linear mixed models
library(lme4)        # Linear Mixed Model
library(ggplot2)      # Visualization
library(ggpubr)
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(performance)  # Model diagnostics (AIC, R², etc.)

source("prda_synergy/utils_prda.R")
base_dir = "D:/WorkingData/PrTecDA_Data/PrDA_somas_Ca_imaging/R_data"

# ==============================================================================
# Load Data
# ==============================================================================
# Load data from csv file
df <- load_data("D:/WorkingData/PrTecDA_Data/PrDA_somas_Ca_imaging/R_data/stim_motor_data_synergy.csv")
n_neurons <- n_distinct(df$roi)
n_neurons

# ==============================================================================
# Exploratory Distributions
# ==============================================================================
# Score Distribution
h_scores <- ggplot(df, aes(x = score)) +
  geom_histogram(binwidth = 0.2, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(condition)) +
  labs(
    title = "Scores",
    x = "Score",
    y = "Count"
  )
h_scores


# Score Norm. 99th Distribution
h_scores_norm99 <- ggplot(df, aes(x = score_norm99)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(condition)) +
  labs(
    title = "Scores",
    x = "Score",
    y = "Count"
  )
h_scores_norm99

# Ca AUC Distribution
h_auc <- ggplot(df, aes(x = ca_auc)) +
  geom_histogram(binwidth = 3, fill = "skyblue", color = "black") +
  facet_grid(cols = vars(condition)) +
  labs(
    title = "Ca area under the curve",
    x = "AUC",
    y = "Count"
  )
h_auc

# ==============================================================================
# GLMM
# ==============================================================================
# Scores
# df_clean <- df %>%
#   filter(!stim_name %in% c("grating_0", "dark_loom"))
glmm_score <- glmmTMB(
  score ~ condition + (1 | fish) + (1 | roi),
  data = df,
  family = Gamma(link = "log")
)

# glmm_score <- glmmTMB(
#   score ~ condition + (1 | fish) + (1 | roi),
#   data = df,
#   family = gaussian(link = "identity")
# )

# Ca AUC
glmm_auc <- glmmTMB(
  ca_auc ~ condition + (1 | fish) + (1 | roi),
  data = df,
  family = Gamma(link = "log")
)

# glmm_auc <- glmmTMB(
#   ca_auc ~ condition + (1 | fish) + (1 | roi),
#   data = df,
#   family = gaussian(link = "identity")
# )

# Norm. Scores
glmm_score_norm99 <- glmmTMB(
  score_norm99 ~ condition + (1 | fish) + (1 | roi),
  data = df,
  family = Gamma(link = "log")
)

# ==============================================================================
# Validation
# ==============================================================================
validate_model(glmm_score, df)
validate_model(glmm_score_norm99, df)
validate_model(glmm_auc, df)

# Residual Hist
hist(resid(glmm_score, type="pearson"), breaks="FD", main="Score Residuals")
hist(resid(glmm_score_norm99, type="pearson"), breaks="FD", main=" Norm. Score Residuals")
hist(resid(glmm_auc, type="pearson"), breaks="FD", main="Ca AUC Residuals")

# ==============================================================================
# Contrasts
# ==============================================================================
# SCORES
# ------------------------------------------------------------------------------
# Estimated marginal means
emm_log  <- emmeans(glmm_score, ~ condition) 
emm_resp <- emmeans(glmm_score, ~ condition, type = "response")

# Synergy in log difference of scores
# synergy_log_scale <- contrast(emm_log, list("synergy" = c(-0.5, -0.5, 1)))

# Synergy in response scale of scores as ratio (x times larger than additive)
# synergy_response_scale <- contrast(emm_resp, list("synergy" = c(-0.5, -0.5, 1)))

# additive synergy contrast: 1*both - 1*stim - 1*motor
emm_resp_rg <- regrid(emm_resp)
additive_synergy <- contrast(
  emm_resp_rg,
  list(additive = c(-1, -1, 1))
)

pairs_log <- pairs(emm_log)
pairs_resp <- pairs(emm_resp)

# ADDITIVE SYN
emm_resp <- emmeans(glmm_score, ~ condition, type = "response")
# Convert to a data.frame-like summary
emm_tab <- summary(emm_resp)
# Extract the vector of means and covariance matrix
b   <- emm_tab$response   # or emm_tab$emmean, see names(emm_tab)
lev <- emm_tab$condition
V   <- vcov(emm_resp)

L <- numeric(length(lev))
L[lev == "stim_only"]   <- -1
L[lev == "motor_only"]  <- -1
L[lev == "stim_motor"]  <-  1

est  <- sum(L * b)              # estimated difference on original scale
se   <- sqrt( t(L) %*% V %*% L )# delta-method SE
z    <- as.numeric(est / se)
pval <- 2 * pnorm(-abs(z))

c(estimate = est, SE = se, z = z, p = pval)

# Bootstrap
set.seed(123)          # for reproducibility
B <- 100               # number of bootstrap replicates (increase in real analysis)
# original estimate
orig_est <- get_synergy(glmm_score)

# storage for bootstrap contrasts
boot_contr <- numeric(B)

for (b in seq_len(B)) {
  # simulate new response from fitted model
  y_sim <- simulate(glmm_score, nsim = 1)[[1]]
  
  # refit model with simulated response
  df_boot <- df
  df_boot$score <- y_sim
  
  fit_boot <- try(
    glmmTMB(
      score ~ condition + (1 | fish) + (1 | roi),
      data = df_boot,
      family = Gamma(link = "log")
    ),
    silent = TRUE
  )
  
  if (inherits(fit_boot, "try-error")) {
    boot_contr[b] <- NA
  } else {
    boot_contr[b] <- get_synergy(fit_boot)
  }
}

# remove failed fits
boot_contr <- boot_contr[!is.na(boot_contr)]

# bootstrap SE
boot_se <- sd(boot_contr)

# percentile 95% CI
ci_low  <- quantile(boot_contr, 0.025)
ci_high <- quantile(boot_contr, 0.975)

# one-sided p-value for H1: synergy > 0
p_one_sided <- mean(boot_contr <= 0)
p_two_sided <- mean(abs(boot_contr) >= abs(orig_est))

list(
  original_estimate = orig_est,
  bootstrap_SE      = boot_se,
  CI_95             = c(lower = ci_low, upper = ci_high),
  p_one_sided       = p_one_sided,
  p_two_sided       = p_two_sided
)
####

# plot_synergy_summary(glmm_score)
g_scores <- plot_synergy_dual_scale(
  glmm_score, df, "Score", "condition",
  order = c("motor_only", "stim_only", "stim_motor")
  )
g_scores
# save_plot(g_scores, file.path(base_dir, "GLMM_scores_synergy"), width=20, height=10, units = "cm", dpi = 300)

save_glmm_full_report(
  glmm_summary          = summary(glmm_score),
  emm_log               = emm_log,
  emm_resp              = emm_resp,
  pairs_log             = pairs_log,
  pairs_resp            = pairs_resp,
  synergy               = additive_synergy, 
  filename              = file.path(base_dir, "results", "glmm_score_results.txt")
)

# Export EMMEANS (response scale)
emm_resp_df <- as.data.frame(emm_resp)
write.csv(emm_resp_df, file.path(base_dir, "results", "glmm_score_emm_response.csv"), row.names = FALSE)

# Export synergy contrast
synergy_df <- as.data.frame(additive_synergy)
write.csv(synergy_df, file.path(base_dir, "results", "glmm_score_synergy.csv"), row.names = FALSE)

# SCORES Norm99
# ------------------------------------------------------------------------------
# Estimated marginal means
emm_norm99_log  <- emmeans(glmm_score_norm99, ~ condition) 
emm_norm99_resp <- emmeans(glmm_score_norm99, ~ condition, type = "response")

pairs_norm99_log <- pairs(emm_norm99_log)
pairs_norm99_resp <- pairs(emm_norm99_resp)

# additive synergy contrast: 1*both - 1*stim - 1*motor
emm_norm99_resp_rg <- regrid(emm_norm99_resp)
additive_synergy_norm99 <- contrast(
  emm_norm99_resp_rg,
  list(additive = c(-1, -1, 1))
)

# plot_synergy_summary(glmm_score_norm99)
g_score_norm99 <- plot_synergy_dual_scale(
  glmm_score_norm99, df, "Norm. Score", "condition",
  order = c("motor_only", "stim_only", "stim+motor")
)

g_score_norm99
# save_plot(g_score_norm99, file.path(base_dir, "GLMM_score_norm99_synergy"), width=20, height=10, units = "cm", dpi = 300)


save_glmm_full_report(
  glmm_summary          = summary(glmm_score_norm99),
  emm_log               = emm_norm99_log,
  emm_resp              = emm_norm99_resp,
  pairs_log             = pairs_norm99_log,
  pairs_resp            = pairs_norm99_resp,
  synergy               = additive_synergy_norm99, 
  filename              = file.path(base_dir, "results", "glmm_score_norm99_results.txt")
)

# Export EMMEANS (response scale)
emm_norm99_resp_df <- as.data.frame(emm_norm99_resp)
write.csv(emm_norm99_resp_df, file.path(base_dir, "results", "glmm_norm99_score_emm_response.csv"), row.names = FALSE)

# Export synergy contrast
synergy_norm99_df <- as.data.frame(additive_synergy_norm99)
write.csv(synergy_df, file.path(base_dir, "results", "glmm_norm99_score_synergy.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# CA AUC
# ------------------------------------------------------------------------------
# Estimated marginal means
emm_auc_log  <- emmeans(glmm_auc, ~ condition) 
emm_auc_resp <- emmeans(glmm_auc, ~ condition, type = "response")

# additive synergy contrast: 1*both - 1*stim - 1*motor
emm_auc_resp_rg <- regrid(emm_auc_resp)
additive_synergy_auc <- contrast(
  emm_auc_resp_rg,
  list(additive = c(-1, -1, 1)),
  type="response"
)


pairs_auc_log <- pairs(emm_auc_log)
pairs_auc_resp <- pairs(emm_auc_resp)

# plot_synergy_summary(glmm_auc)
g_auc <- plot_synergy_dual_scale(
  glmm_auc, df, "Ca AUC", "condition",
  order = c("stim_only", "motor_only", "stim+motor")
)
g_auc
# save_plot(g_auc, file.path(base_dir, "GLMM_ca_auc_synergy"), width=20, height=10, units = "cm", dpi = 300)


save_glmm_full_report(
  glmm_summary          = summary(glmm_auc),
  emm_log               = emm_auc_log,
  emm_resp              = emm_auc_resp,
  pairs_log             = pairs_auc_log,
  pairs_resp            = pairs_auc_resp,
  synergy               = additive_synergy_auc, 
  filename              = file.path(base_dir, "results", "glmm_ca_auc_results.txt")
)

# Export EMMEANS (response scale)
emm_auc_resp_df <- as.data.frame(emm_auc_resp)
write.csv(emm_auc_resp_df, file.path(base_dir, "results", "glmm_auc_emm_response.csv"), row.names = FALSE)

# Export synergy contrast
synergy_auc_df <- as.data.frame(additive_synergy_auc)
write.csv(synergy_auc_df, file.path(base_dir, "results", "glmm_auc_synergy.csv"), row.names = FALSE)


# ==============================================================================
# Selectivity Indices
# ==============================================================================
df_ssi_org <- read_csv("D:/WorkingData/PrTecDA_Data/PrDA_somas_Ca_imaging/R_data/score_metrics_per_neuron.csv")
df_ssi_org <- df_ssi_org %>%
  mutate(fish = recode(fish,
                       "1" = "1",
                       "3" = "2",
                       "4" = "3",
                       "5" = "4",
                       "6" = "5",
                       "9" = "6"))
df_ssi_org <- df_ssi_org %>%
  mutate(fish = factor(fish, levels = c("1", "2", "3", "4", "5", "6")))

# Remove bad data
df_ssi <- df_ssi_org %>%
  filter(
    between(vis_spont_index, -1, 1),
    between(mixed_spont_index, -1, 1),
    between(mixed_vis_index, -1, 1),
    between(synergy_index, -1, 1)
  )

tibble(
  before = nrow(df_ssi_org),
  after  = nrow(df_ssi),
  removed = nrow(df_ssi_org) - nrow(df_ssi)
)
n_neurons <- n_distinct(df_ssi$roi)
n_neurons

# Histograms
hist(df_ssi$vis_spont_index, breaks=50, main="vis_spont_index")
hist(df_ssi$mixed_spont_index , breaks=50, main="mixed_spont_index")
hist(df_ssi$mixed_vis_index , breaks=50, main="mixed_vis_index")
hist(df_ssi$synergy_index , breaks=50, main="synergy_index")
hist(df_ssi$synergy , breaks=50, main="synergy")


# 1. vis_spont_index
m_vis_spont <- lmer(
  vis_spont_index ~ 1 + (1 | fish),
  data = df_ssi
)
simple_validate(m_vis_spont)

# 2. mixed_spont_index
m_mixed_spont <- lmer(
  mixed_spont_index ~ 1 + (1 | fish),
  data = df_ssi
)
simple_validate(m_mixed_spont)

# 3. mixed_vis_index
m_mixed_vis <- lmer(
  mixed_vis_index ~ 1 + (1 | fish),
  data = df_ssi
)
simple_validate(m_mixed_vis)

# 4. synergy_index
m_synergy_index <- lmer(
  synergy_index ~ 1 + (1 | fish),
  data = df_ssi
)
simple_validate(m_synergy_index)

# 5. synergy ratio (positive continuous)
m_synergy_gamma <- glmmTMB(
  synergy ~ 1 + (1 | fish),
  data = df_ssi,
  family = Gamma(link = "log")
)
simple_validate(m_synergy_gamma)


models <- list(
  m_vis_spont      = m_vis_spont,
  m_mixed_spont    = m_mixed_spont,
  m_mixed_vis      = m_mixed_vis,
  m_synergy_index  = m_synergy_index,
  m_synergy_gamma  = m_synergy_gamma
)

# Test if SSI are different from 0 (1)
ssi_tests <- test_ssi_significance(models)
ssi_tests
write_table_txt(ssi_tests,
                file.path(base_dir, "results", "ssi_tests.txt"))

# Plot
g_ssi <- plot_fish_clouds_dualpanel(df_ssi, models, spacing = 10, jitter_width = 0.02, point_size = 1)
g_ssi
# save_plot(g_ssi, file.path(base_dir, "LMM_Selectivity_Indices"), width=20, height=14, units = "cm", dpi = 300)

# ==============================================================================
# Hierarchical Bootstrap
# ==============================================================================

ssi_boot_results <- list(
  vis_spont     = parallel_hier_boot(df_ssi, "fish", "vis_spont_index", null_value = 0, B=100000, workers=12),
  mixed_spont   = parallel_hier_boot(df_ssi, "fish", "mixed_spont_index", null_value = 0, B=100000, workers=12),
  mixed_vis     = parallel_hier_boot(df_ssi, "fish", "mixed_vis_index", null_value = 0, B=100000, workers=12),
  synergy_index = parallel_hier_boot(df_ssi, "fish", "synergy_index", null_value = 0, B=100000, workers=12),
  synergy_ratio = parallel_hier_boot(df_ssi, "fish", "synergy", null_value = 1, B=100000, workers=12)
)

ssi_boot_table <- make_ssi_boot_table(ssi_boot_results)
ssi_boot_table
write_table_txt(ssi_boot_table,
                file.path(base_dir, "results", "ssi_boot_table.txt"))

save_boot_distributions_combined(
  ssi_boot_results,
  file.path(base_dir, "results", "bootstrap_dist.csv")
  )

write.csv(as.data.frame(ssi_boot_table), file.path(base_dir, "results", "selectivity_index_bootstrap.csv"), row.names = FALSE)

# Plot the boostrapped distributions
g_boot <- plot_all_boot(ssi_boot_results)
g_boot
# save_plot(g_boot, file.path(base_dir, "Bootstrap_Distributions_Selectivity_Indices"), width=20, height=14, units = "cm", dpi = 300)

g_ssi_boot <- plot_fish_clouds_dualpanel_boot(
  df_ssi,
  ssi_boot_results,
  spacing = 10,
  jitter_width = 0.02,
  point_size = 1
)
g_ssi_boot
# save_plot(g_ssi_boot, file.path(base_dir, "Bootstrap_Selectivity_Indices"), width=20, height=14, units = "cm", dpi = 300)