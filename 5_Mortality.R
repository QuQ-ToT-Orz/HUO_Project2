#### 7 Mortality Analysis with Dispersion Index ####
gc()
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(purrr)
library(survey)
library(compareC)

#### 1 Load Data ####
dir_path <- "./data/"

load(file = paste0(dir_path, "dispersion/dispersion_df_old.RData"))
load(file = paste0(dir_path, "dispersion/dispersion_summary_old.RData"))
load(file = paste0(dir_path, "count/data_analysis_old.RData"))
load(file = paste0(dir_path, "runlength/event_analysis_old.RData"))
load(file = paste0(dir_path, "count/Act_Analysis_old.RData"))

load(file = paste0(dir_path, "dispersion/dispersion_df_new.RData"))
load(file = paste0(dir_path, "dispersion/dispersion_summary_new.RData"))
load(file = paste0(dir_path, "mims/data_analysis_new.RData"))
load(file = paste0(dir_path, "runlength/event_analysis_new.RData"))
load(file = paste0(dir_path, "mims/Act_Analysis_new.RData"))

# Aggregate activity measures from daily to subject level
activity_summary <- Act_Analysis %>%
  group_by(SEQN) %>%
  summarise(
    TAC = mean(TAC, na.rm = TRUE),
    Peak30 = mean(Peak30, na.rm = TRUE),
    ASTP = mean(ASTP, na.rm = TRUE),
    SATP = mean(SATP, na.rm = TRUE),
    SBout = mean(SBout, na.rm = TRUE),
    ABout = mean(ABout, na.rm = TRUE),
    ST = mean(ST, na.rm = TRUE),
    WT = mean(WT, na.rm = TRUE),
    MVPA = mean(MVPA, na.rm = TRUE),
    .groups = "drop"
  )

# Extract night regularity metrics
max_cv <- max(event_analysis$mean_sleep_cv, na.rm = TRUE)
night_metrics <- event_analysis %>%
  select(SEQN, days_early_activer, days_late_activer, consolidated_sleep_days, mean_sleep_cv) %>%
  distinct() %>%
  mutate(mean_sleep_cv = ifelse(is.nan(mean_sleep_cv), max_cv, mean_sleep_cv))

sleep_event_summary <- event_analysis %>%                     
    mutate(SEQN = as.numeric(as.character(SEQN))) %>%  
    # convert factor numeric                                            
    group_by(SEQN, WEEKDAY) %>%                                 
    summarise(sleep_events = sum(start >= 60 & start < 300), .groups = "drop")%>%                                                       
    group_by(SEQN) %>%                                          
    summarise(mean_sleep_events = mean(sleep_events, na.rm = TRUE), .groups = "drop")     

#### 2 Validation Against Hawkes MLE ####
source("./function_process.R")

# Load marked MLE
mle_marked_file <- paste0("./data/hawkes/", "fits_marked_all_old.RData")
mle_marked_file <- paste0("./data/hawkes/", "fits_marked_all_new.RData")
load(mle_marked_file)

mle_marked_df <- process_fits_data(fits_list_all) %>%
  mutate(
    SEQN = as.numeric(as.character(id)),
    n_mle_marked = alpha/beta,
    mu_mle_marked = mu * 60
  ) %>%
  select(SEQN, mu_mle_marked, n_mle_marked)

# Load unmarked MLE
mle_unmarked_file <- paste0("./data/hawkes/", "fits_unmarked_all_old.RData")
mle_unmarked_file <- paste0("./data/hawkes/", "fits_unmarked_all_new.RData")
load(mle_unmarked_file)

mle_unmarked_df <- process_fits_data(fits_list_all) %>%
  mutate(
    SEQN = as.numeric(as.character(id)),
    n_mle_unmarked = alpha/beta,
    mu_mle_unmarked = mu * 60
  ) %>%
  select(SEQN, mu_mle_unmarked, n_mle_unmarked)

# Merge MLE estimates
mle_df <- mle_marked_df %>%
  inner_join(mle_unmarked_df, by = "SEQN")

# Merge with dispersion estimates
validation_df <- dispersion_summary %>%
  mutate(SEQN = as.numeric(as.character(SEQN))) %>%
  inner_join(mle_df, by = "SEQN")

analysis_df <- data_analysis %>%
  mutate(SEQN = as.numeric(as.character(SEQN))) %>%
  inner_join(validation_df, by = "SEQN") %>%
  left_join(activity_summary %>% mutate(SEQN = as.numeric(as.character(SEQN))), by = "SEQN") %>%
  left_join(night_metrics %>% mutate(SEQN = as.numeric(as.character(SEQN))), by = "SEQN") %>%
  left_join(sleep_event_summary, by = "SEQN") %>%
  mutate(
    # MLE: derive lambda from mu/(1-n)
    lambda_mle_marked = mu_mle_marked / (1 - n_mle_marked),
    lambda_mle_unmarked = mu_mle_unmarked / (1 - n_mle_unmarked)
  )

 analysis_df <- analysis_df %>%
    mutate(
      TAC = sqrt(TAC),
      Peak30 = sqrt(Peak30)
    )
  cor_vars <- c("n_raw", "n_adj", "n_marks", "n_adj_marks",
  "n_mle_unmarked", "n_mle_marked",

  "mu_star_raw", "mu_star_adj", "mu_star_marks", "mu_star_adj_marks",
  "mu_mle_unmarked", "mu_mle_marked",

  "lambda_count", "lambda_marks", "lambda_mle_unmarked", "lambda_mle_marked",

  "TAC", "Peak30", "ASTP", "SATP", "WT", "MVPA", "ST")
  library(corrplot)
  cor_matrix <- cor(analysis_df[, cor_vars], use = "pairwise.complete.obs")
  corrplot(cor_matrix, method = "color", type = "upper", addCoef.col = "black",
           number.cex = 0.5, tl.cex = 0.6, tl.srt = 45)

  pdf("Output/dispersion/variables/hip_cor.pdf", width = 12, height = 10)
  pdf("Output/dispersion/variables/wrist_cor.pdf", width = 12, height = 10)
  corrplot(cor_matrix, method = "color", type = "upper", addCoef.col = "black",
           number.cex = 0.5, tl.cex = 0.6, tl.srt = 45)
  dev.off()

#### 3 Merge with Covariates for Mortality Analysis ####
analysis_df <- analysis_df %>%
  mutate(
    # Survival
    surv_time = permth_exm / 12,
    event = mortstat
  )

#### Survey-Weighted Cox Models ####
cat("\n=== Survey-Weighted Cox Models (NHANES Design) ===\n")

# Create survey design object
survey_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  # weights = ~WTMEC2YR,
  weights = ~I(as.numeric(WTMEC2YR) / 2),
  nest = TRUE,
  data = analysis_df
)

# Define model specifications: MLE (unmarked/marked) + 4 Dispersion versions
n_versions <- list(
  list(n = "n_mle_unmarked", lambda = "lambda_mle_unmarked", label = "MLE (unmarked)"),
  list(n = "n_mle_marked", lambda = "lambda_mle_marked", label = "MLE (marked)"),
  list(n = "n_raw", lambda = "lambda_count", label = "Dispersion (raw)"),
  list(n = "n_adj", lambda = "lambda_count", label = "Dispersion (adj)"),
  list(n = "n_marks", lambda = "lambda_marks", label = "Dispersion (marks)"),
  list(n = "n_adj_marks", lambda = "lambda_marks", label = "Dispersion (adj+marks)")
)

basic_covars <- "+ Age + Gender + Race + WT"
full_covars <- "+ Age + Gender + Race + WT + Peak30 + BMI_cat + SmokeCigs + DrinkStatus + EducationAdult + MobilityProblem + Diabetes + CHF + CHD + Stroke + Cancer"

# Function to extract model results
extract_results <- function(model, var_name, data) {
  coefs <- summary(model)$coefficients
  
  beta <- coefs[var_name, 1]
  se <- coefs[var_name, 4]  # robust SE
  p <- coefs[var_name, 6]
  
  # Get SD for proper scaling
  sd_var <- sd(data[[var_name]], na.rm = TRUE)
  
  # HR per SD (not per 1 unit!)
  hr <- exp(beta * sd_var)
  ci_low <- exp((beta - 1.96 * se) * sd_var)
  ci_high <- exp((beta + 1.96 * se) * sd_var)
  
  list(hr = hr, ci_low = ci_low, ci_high = ci_high, p = p)
}

# Helper function to format HR with CI and significance stars
format_hr <- function(hr, ci_low, ci_high, p) {
  stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  sprintf("%.3f (%.3f-%.3f)%s", hr, ci_low, ci_high, stars)
}
format_hr <- function(hr, ci_low, ci_high, p) {
  stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.10, ".", ""))))
  sprintf("%.3f (%.3f-%.3f)%s p=%.4f", hr, ci_low, ci_high, stars, p)
}

# Function to calculate C-index from a fitted Cox model
# cvwrapr (Cross-validated, out-of-sample)
calculate_c_index <- function(model, model_formula, data) {
  model_vars <- all.vars(model_formula)
  model_data <- data[complete.cases(data[, model_vars]), , drop = FALSE]
  lp <- as.numeric(predict(model, newdata = model_data, type = "lp"))
  concordance_fit <- survival::concordance(
    Surv(surv_time, event) ~ lp,
    data = model_data,
    reverse = TRUE
  )
  
  c_index <- concordance_fit$concordance
  se <- sqrt(concordance_fit$var)
  
  list(
    c_index = c_index,
    se = se,
    ci_low = c_index - 1.96 * se,
    ci_high = c_index + 1.96 * se,
    n = concordance_fit$n
  )
}

format_c_index <- function(c_index, ci_low, ci_high) {
  sprintf("%.3f (%.3f-%.3f)", c_index, ci_low, ci_high)
}

compare_c_index <- function(model_y, formula_y, model_z, formula_z, data) {
  if (identical(model_y, model_z)) {
    model_vars <- all.vars(formula_y)
    model_data <- data[complete.cases(data[, model_vars]), , drop = FALSE]
    return(list(diff = 0, se = 0, z = NA_real_, p = NA_real_, n = nrow(model_data)))
  }
  
  model_vars <- union(all.vars(formula_y), all.vars(formula_z))
  model_data <- data[complete.cases(data[, model_vars]), , drop = FALSE]
  lp_y <- as.numeric(predict(model_y, newdata = model_data, type = "lp"))
  lp_z <- as.numeric(predict(model_z, newdata = model_data, type = "lp"))
  compare_fit <- compareC(
    timeX = model_data$surv_time,
    statusX = model_data$event,
    scoreY = -lp_y,
    scoreZ = -lp_z
  )
  
  list(
    diff = compare_fit$est.diff_c,
    se = sqrt(compare_fit$est.vardiff_c),
    z = compare_fit$zscore,
    p = compare_fit$pval,
    n = nrow(model_data)
  )
}

format_compare_c <- function(compare_result) {
  data.frame(
    diff = compare_result$diff,
    z = compare_result$z,
    p = compare_result$p,
    n = compare_result$n
  )
}

# Run all models and collect results
model_results <- lapply(n_versions, function(spec) {
  # Model 1
  f1 <- as.formula(paste0("Surv(surv_time, event) ~ ", spec$n, " + ", spec$lambda, basic_covars))
  m1 <- svycoxph(f1, design = survey_design)
  
  # Pass data to extract_results for SD calculation
  r1_n <- extract_results(m1, spec$n, analysis_df)
  r1_lambda <- extract_results(m1, spec$lambda, analysis_df)
  
  # Model 2
  f2 <- as.formula(paste0("Surv(surv_time, event) ~ ", spec$n, " + ", spec$lambda, full_covars))
  m2 <- svycoxph(f2, design = survey_design)
  
  r2_n <- extract_results(m2, spec$n, analysis_df)
  r2_lambda <- extract_results(m2, spec$lambda, analysis_df)
  r2_peak30 <- extract_results(m2, "Peak30", analysis_df)
  
  c1 <- calculate_c_index(m1, f1, analysis_df)
  c2 <- calculate_c_index(m2, f2, analysis_df)
  c_m2_vs_m1 <- compare_c_index(m2, f2, m1, f1, analysis_df)
  
  list(
    result = data.frame(
      Method = spec$label,
      n_M1 = format_hr(r1_n$hr, r1_n$ci_low, r1_n$ci_high, r1_n$p),
      lambda_M1 = format_hr(r1_lambda$hr, r1_lambda$ci_low, r1_lambda$ci_high, r1_lambda$p),
      n_M2 = format_hr(r2_n$hr, r2_n$ci_low, r2_n$ci_high, r2_n$p),
      lambda_M2 = format_hr(r2_lambda$hr, r2_lambda$ci_low, r2_lambda$ci_high, r2_lambda$p),
      Peak30_M2 = format_hr(r2_peak30$hr, r2_peak30$ci_low, r2_peak30$ci_high, r2_peak30$p),
      C_index_M1 = c1$c_index,
      C_index_M1_se = c1$se,
      C_index_M1_CI = format_c_index(c1$c_index, c1$ci_low, c1$ci_high),
      C_index_M1_n = c1$n,
      C_index_M2 = c2$c_index,
      C_index_M2_se = c2$se,
      C_index_M2_CI = format_c_index(c2$c_index, c2$ci_low, c2$ci_high),
      C_index_M2_n = c2$n,
      C_index_M2_minus_M1_p = c_m2_vs_m1$p,
      C_index_M2_minus_M1_z = c_m2_vs_m1$z,
      C_index_M2_minus_M1_compareC_n = c_m2_vs_m1$n
    ),
    m1 = m1,
    m2 = m2,
    f1 = f1,
    f2 = f2
  )
})

results_df <- bind_rows(map(model_results, "result"))
cindex_df <- results_df %>%
  select(Method, starts_with("C_index")) %>%
  mutate(
    C_index_M2_minus_M1 = C_index_M2 - C_index_M1
  )

paired_method_comparisons <- list(
  list(unmarked = "MLE (unmarked)", marked = "MLE (marked)", pair = "MLE: marked - unmarked"),
  list(unmarked = "Dispersion (raw)", marked = "Dispersion (marks)", pair = "Dispersion: marks - raw"),
  list(unmarked = "Dispersion (adj)", marked = "Dispersion (adj+marks)", pair = "Dispersion: adj+marks - adj")
)

paired_cindex_df <- bind_rows(lapply(paired_method_comparisons, function(pair_spec) {
  unmarked_index <- which(results_df$Method == pair_spec$unmarked)
  marked_index <- which(results_df$Method == pair_spec$marked)
  
  m1_compare <- compare_c_index(
    model_results[[marked_index]]$m1, model_results[[marked_index]]$f1,
    model_results[[unmarked_index]]$m1, model_results[[unmarked_index]]$f1,
    analysis_df
  )
  m2_compare <- compare_c_index(
    model_results[[marked_index]]$m2, model_results[[marked_index]]$f2,
    model_results[[unmarked_index]]$m2, model_results[[unmarked_index]]$f2,
    analysis_df
  )
  
  data.frame(
    Pair = pair_spec$pair,
    Unmarked_Method = pair_spec$unmarked,
    Marked_Method = pair_spec$marked,
    M1_C_index_unmarked = results_df$C_index_M1[unmarked_index],
    M1_C_index_marked = results_df$C_index_M1[marked_index],
    M1_marked_minus_unmarked = m1_compare$diff,
    M1_z = m1_compare$z,
    M1_p = m1_compare$p,
    M1_compareC_n = m1_compare$n,
    M2_C_index_unmarked = results_df$C_index_M2[unmarked_index],
    M2_C_index_marked = results_df$C_index_M2[marked_index],
    M2_marked_minus_unmarked = m2_compare$diff,
    M2_z = m2_compare$z,
    M2_p = m2_compare$p,
    M2_compareC_n = m2_compare$n
  )
}))

results_df <- results_df %>%
  select(-starts_with("C_index"))

# Print combined table
cat("\n=== Combined Table: Activity Metrics and Mortality ===\n")
cat("Model 1: Adjusted for age, gender, race, wear time\n")
cat("Model 2: + Peak30, BMI, smoking, alcohol, education, mobility, comorbidities\n")
cat("* p<0.05, ** p<0.01, *** p<0.001\n\n")
print(results_df, row.names = FALSE)

cat("\n=== C-index and Differences ===\n")
cat("C-index calculated from each Cox model linear predictor; M2-M1 p-values are from compareC.\n\n")
print(cindex_df, row.names = FALSE)

cat("\n=== Paired Marked vs Unmarked C-index Comparisons ===\n")
cat("Positive differences mean the marked model has a higher C-index than its paired unmarked model.\n\n")
print(paired_cindex_df, row.names = FALSE)

write.csv(results_df, "Output/mortality/hip_results_df.csv", row.names = FALSE)
write.csv(results_df, "Output/mortality/wrist_results_df.csv", row.names = FALSE)
write.csv(cindex_df, "Output/mortality/hip_cindex_df.csv", row.names = FALSE)
write.csv(cindex_df, "Output/mortality/wrist_cindex_df.csv", row.names = FALSE)
write.csv(paired_cindex_df, "Output/mortality/hip_paired_cindex_df.csv", row.names = FALSE)
write.csv(paired_cindex_df, "Output/mortality/wrist_paired_cindex_df.csv", row.names = FALSE)

#### 4 Visualizations ####
par(mfrow = c(1, 3))                           
                                                
  # n: MLE vs Dispersion                         
  plot(analysis_df$n_mle_marked, analysis_df$n_raw, main
   = paste0("n (r=", round(cor(analysis_df$n_mle_marked,
   analysis_df$n_raw, use = "complete.obs"), 3), 
  ")"), xlab = "n_mle", ylab = "n_star")          
  abline(lm(n_raw ~ n_mle_marked, data = analysis_df),  
  col = "blue")                                  
                                                 
  # lambda: MLE vs Dispersion                    
  plot(analysis_df$lambda_mle_marked,                   
  analysis_df$lambda_count, main = paste0("lambda
   (r=", round(cor(analysis_df$lambda_mle_marked,       
  analysis_df$lambda_count, use =                
  "complete.obs"), 3), ")"), xlab = "lambda_mle",
   ylab = "lambda_count")                        
  abline(lm(lambda_count ~ lambda_mle_marked, data =    
  analysis_df), col = "blue")                    
                                                 
  # mu: MLE vs Dispersion                        
  plot(analysis_df$mu_mle_marked,                       
  analysis_df$mu_star_raw, main = paste0("mu     
  (r=", round(cor(analysis_df$mu_mle_marked,            
  analysis_df$mu_star_raw, use = "complete.obs"),
   3), ")"), xlab = "mu_mle", ylab =             
  "mu_star")                                 
  abline(lm(mu_star_raw ~ mu_mle_marked, data =         
  analysis_df), col = "blue")                    
                                                 
  par(mfrow = c(1, 1))

  pdf("Output/dispersion/variables/hip_map.pdf", width = 12, height = 5)
  pdf("Output/dispersion/variables/wrist_map.pdf", width = 12, height = 5)
  par(mfrow = c(1, 3))
  plot(analysis_df$n_mle_marked, analysis_df$n_raw, main
   = paste0("n (r=", round(cor(analysis_df$n_mle_marked,
   analysis_df$n_raw, use = "complete.obs"), 3),
  ")"), xlab = "n_mle", ylab = "n_star")
  abline(lm(n_raw ~ n_mle_marked, data = analysis_df),
  col = "blue")
  plot(analysis_df$lambda_mle_marked,
  analysis_df$lambda_count, main = paste0("lambda
   (r=", round(cor(analysis_df$lambda_mle_marked,
  analysis_df$lambda_count, use =
  "complete.obs"), 3), ")"), xlab = "lambda_mle",
   ylab = "lambda_count")
  abline(lm(lambda_count ~ lambda_mle_marked, data =
  analysis_df), col = "blue")
  plot(analysis_df$mu_mle_marked,
  analysis_df$mu_star_raw, main = paste0("mu
  (r=", round(cor(analysis_df$mu_mle_marked,
  analysis_df$mu_star_raw, use = "complete.obs"),
   3), ")"), xlab = "mu_mle", ylab =
  "mu_star")
  abline(lm(mu_star_raw ~ mu_mle_marked, data =
  analysis_df), col = "blue")
  par(mfrow = c(1, 1))
  dev.off()

#### 4b Boxplots: Parameter Distributions Across Methods ####

# Reshape data for n estimates
n_long <- analysis_df %>%
  select(SEQN, n_raw, n_adj, n_marks, n_adj_marks, n_mle_unmarked, n_mle_marked) %>%
  pivot_longer(cols = -SEQN, names_to = "method", values_to = "n") %>%
  mutate(method = factor(method,
    levels = c("n_raw", "n_adj", "n_marks", "n_adj_marks", "n_mle_unmarked", "n_mle_marked"),
    labels = c("Disp\n(raw)", "Disp\n(adj)", "Disp\n(marks)", "Disp\n(adj+marks)", "MLE\n(unmarked)", "MLE\n(marked)")))

# Reshape data for mu estimates
mu_long <- analysis_df %>%
  select(SEQN, mu_star_raw, mu_star_adj, mu_star_marks, mu_star_adj_marks, mu_mle_unmarked, mu_mle_marked) %>%
  pivot_longer(cols = -SEQN, names_to = "method", values_to = "mu") %>%
  mutate(method = factor(method,
    levels = c("mu_star_raw", "mu_star_adj", "mu_star_marks", "mu_star_adj_marks", "mu_mle_unmarked", "mu_mle_marked"),
    labels = c("Disp\n(raw)", "Disp\n(adj)", "Disp\n(marks)", "Disp\n(adj+marks)", "MLE\n(unmarked)", "MLE\n(marked)")))

# Reshape data for lambda estimates
lambda_long <- analysis_df %>%
  select(SEQN, lambda_count, lambda_marks, lambda_mle_unmarked, lambda_mle_marked) %>%
  pivot_longer(cols = -SEQN, names_to = "method", values_to = "lambda") %>%
  mutate(method = factor(method,
    levels = c("lambda_count", "lambda_marks", "lambda_mle_unmarked", "lambda_mle_marked"),
    labels = c("Count-based", "Mark-based", "MLE\n(unmarked)", "MLE\n(marked)")))

# Boxplot for n estimates
p_n <- ggplot(n_long, aes(x = method, y = n, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Branching Ratio (n) Across Methods", x = "", y = "n") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))

# Boxplot for mu estimates
p_mu <- ggplot(mu_long, aes(x = method, y = mu, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Immigration Rate (mu) Across Methods", x = "", y = "mu (events/hour)") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))

# Boxplot for lambda estimates
p_lambda <- ggplot(lambda_long, aes(x = method, y = lambda, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Event Rate (lambda) Across Methods", x = "", y = "lambda (events/hour)") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))

# Display plots
print(p_n)
print(p_mu)
print(p_lambda)
ggsave("Output/dispersion/variables/hip_boxplot_n.pdf", p_n, width = 8, height = 6)
ggsave("Output/dispersion/variables/hip_boxplot_mu.pdf", p_mu, width = 8, height = 6)
ggsave("Output/dispersion/variables/hip_boxplot_lambda.pdf", p_lambda, width = 8, height = 6)
ggsave("Output/dispersion/variables/wrist_boxplot_n.pdf", p_n, width = 8, height = 6)
ggsave("Output/dispersion/variables/wrist_boxplot_mu.pdf", p_mu, width = 8, height = 6)
ggsave("Output/dispersion/variables/wrist_boxplot_lambda.pdf", p_lambda, width = 8, height = 6)

#### 5 Night Regularity Exploratory Analysis ####
model_night <- svycoxph((Surv(surv_time, event) ~
                         mean_sleep_cv +
                         mean_sleep_events +
                         days_early_activer +
                         days_late_activer +
                         Age + Gender + Race),
                       design = survey_design)        

print(summary(model_night))

# Extract night model results for LaTeX table
night_coefs <- summary(model_night)$coefficients
night_vars <- c("mean_sleep_cv",
                "mean_sleep_events",
                "days_early_activer", "days_late_activer",
                "Age", "GenderFemale",
                "RaceMexican American", "RaceOther Hispanic", "RaceBlack", "RaceOther")
night_labels <- c("Sleep consolidation (CV)",
                  "Mean nocturnal activity",
                  "Days with early activity", "Days with late activity",
                  "Age (per year)", "Female",
                  "Race: Mexican American", "Race: Other Hispanic", "Race: Black", "Race: Other")

night_table <- data.frame(
  Variable = night_labels,
  HR = sprintf("%.4f", exp(night_coefs[night_vars, 1])),
  CI = sprintf("%.4f--%.4f",
    exp(night_coefs[night_vars, 1] - 1.96 * night_coefs[night_vars, 4]),
    exp(night_coefs[night_vars, 1] + 1.96 * night_coefs[night_vars, 4])),
  p = ifelse(night_coefs[night_vars, 6] < 0.001, "<0.001",
    sprintf("%.4f", night_coefs[night_vars, 6])),
  sig = ifelse(night_coefs[night_vars, 6] < 0.001, "***",
    ifelse(night_coefs[night_vars, 6] < 0.01, "**",
    ifelse(night_coefs[night_vars, 6] < 0.05, "*",
      ifelse(night_coefs[night_vars, 6] < 0.10, ".", ""))))
)

print(night_table, row.names = FALSE)
write.csv(night_table, "Output/mortality/hip_night_table.csv", row.names = FALSE)
write.csv(night_table, "Output/mortality/wrist_night_table.csv", row.names = FALSE)
