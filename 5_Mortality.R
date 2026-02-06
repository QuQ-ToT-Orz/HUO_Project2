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
    MVPA = mean(MVPA, na.rm = TRUE),
    Peak30 = mean(Peak30, na.rm = TRUE),
    ASTP = mean(ASTP, na.rm = TRUE),
    SBout = mean(SBout, na.rm = TRUE),
    ABout = mean(ABout, na.rm = TRUE),
    ST = mean(ST, na.rm = TRUE),
    WT = mean(WT, na.rm = TRUE),
    .groups = "drop"
  )

# Extract night regularity metrics
night_metrics <- event_analysis %>%
  select(SEQN, days_early_activer, days_late_activer, days_consolidated_sleep) %>%
  distinct()

sleep_event_summary <- event_analysis %>%                     
    mutate(SEQN = as.numeric(as.character(SEQN))) %>%  # convert factor numeric                                            
    group_by(SEQN, WEEKDAY) %>%                                 
    summarise(                                                  
      sleep_events = sum(start >= 60 & start < 300),            
      .groups = "drop"                                          
    ) %>%                                                       
    group_by(SEQN) %>%                                          
    summarise(                                                  
      mean_sleep_events = mean(sleep_events, na.rm = TRUE),     
      .groups = "drop"                                          
    )     

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

  "TAC", "MVPA", "Peak30", "ASTP", "ST", "WT")
  library(corrplot)
  cor_matrix <- cor(analysis_df[, cor_vars], use = "pairwise.complete.obs")
  corrplot(cor_matrix, method = "color", type = "upper", addCoef.col = "black",
           number.cex = 0.5, tl.cex = 0.6, tl.srt = 45)

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
  weights = ~WTMEC2YR,
  # weights = ~I(as.numeric(WTMEC2YR) / 2),
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

# Run all models and collect results
results_list <- lapply(n_versions, function(spec) {
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
  
  data.frame(
    Method = spec$label,
    n_M1 = format_hr(r1_n$hr, r1_n$ci_low, r1_n$ci_high, r1_n$p),
    lambda_M1 = format_hr(r1_lambda$hr, r1_lambda$ci_low, r1_lambda$ci_high, r1_lambda$p),
    n_M2 = format_hr(r2_n$hr, r2_n$ci_low, r2_n$ci_high, r2_n$p),
    lambda_M2 = format_hr(r2_lambda$hr, r2_lambda$ci_low, r2_lambda$ci_high, r2_lambda$p),
    Peak30_M2 = format_hr(r2_peak30$hr, r2_peak30$ci_low, r2_peak30$ci_high, r2_peak30$p)
  )
})

results_df <- bind_rows(results_list)

# Print combined table
cat("\n=== Combined Table: Activity Metrics and Mortality ===\n")
cat("Model 1: Adjusted for age, gender, race, wear time\n")
cat("Model 2: + Peak30, BMI, smoking, alcohol, education, mobility, comorbidities\n")
cat("* p<0.05, ** p<0.01, *** p<0.001\n\n")

print(results_df, row.names = FALSE)

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


#### 5 Night Regularity Exploratory Analysis ####
 par(mfrow = c(1, 3))                                          
    hist(analysis_df$days_consolidated_sleep, main =            
  "Consolidated Sleep", xlab = "")                                                               
    hist(analysis_df$days_early_activer, main = "Early          
  Activation", xlab = "")                                       
    hist(analysis_df$days_late_activer, main = "Late            
  Activation", xlab = "")                                       
  par(mfrow = c(1, 1))     

model_night <- svycoxph((Surv(surv_time, event) ~
                         days_consolidated_sleep +
                         mean_sleep_events +
                         days_early_activer +
                         days_late_activer +
                         Age + Gender + Race),
                       design = survey_design)        

print(summary(model_night))
