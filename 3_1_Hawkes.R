#### Univariate marked hawkes ####
# capturing clustering/burstiness #
#### 0 Packages ####
gc()
rm(list = ls())
R.version

library(emhawkes)
library(hawkesbow)
library(hawkes)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)

# devtools::install_local("./HUO", force = TRUE, dependencies = TRUE)
# devtools::install_github("cmjt/stelfi")
# detach("package:HUO", unload = TRUE)
library(HUO)
source("./function_process.R")

#### 1 Data simulation ####
set.seed(42)
fn_mark <- function(...) {
  sample(1:5, 1)
}
h <- new("hspec",
  mu = 0.1, alpha = 0.5, beta = 1,
  rmark = fn_mark
)
res <- hsim(h, size = 1000)
res_simulation <- tibble(
  simulation_time = as.vector(res$N),
  mark = res$mark,
  arrival_simulation = res$arrival,
  inter_arrival_simulation = res$inter_arrival
)
res_simulation <- res_simulation %>% mutate(diff = c(0, diff(arrival_simulation)))
res_simulation <- res_simulation[-1, ]

iv <- list(mu = 0.1, alpha = 0.5, beta = 1)
fit <- fit_hawkes(
  times = res_simulation$arrival_simulation,
  marks = res_simulation$mark / mean(res_simulation$mark),
  parameters = iv
)
fit <- fit_hawkes(
  times = res_simulation$arrival_simulation,
  marks = res_simulation$mark / mean(res_simulation$mark)
)
fit <- fit_hawkes(
  times = res_simulation$arrival_simulation,
  marks = res_simulation$mark / mean(res_simulation$mark), model = 2
)
get_coefs(fit)
show_hawkes(fit)
show_hawkes_GOF(fit)
check_conv(fit)$converged

# !!!! douting !!!!
# # Simulate events for 24 hours (1440 minutes)
# # But force first 8 hours (480 minutes) to be empty
# times_active <- simulate_hawkes(mu, alpha, beta, 960) # Simulate for 16 hours
# times_full <- times_active + 480 # Shift events by 8 hours
# # Estimate parameters for both cases
# fit_full <- fit_hawkes(times_full, parameters = sv)
# fit_active <- fit_hawkes(times_active, parameters = sv)

#### 2 hip ####
load("../2025/data/runlength/event_analysis_old.RData")
load("../2025/data/count/Act_Analysis_old.RData")
load("../2025/data/count/data_analysis_old.RData")

daily_data <- event_analysis
daily_data <- event_analysis %>% mutate(original_value = 1) 
seqn <- unique(daily_data$SEQN)[1]
seqn = 40583
weekday <- unique(daily_data$WEEKDAY)[1]
hip_data <- transform_hawkes_data(daily_data, seqn, weekday, active_only = TRUE, data_type = "hip")
fit1_hip <- fit_hawkes(
  times = hip_data$times,
  marks = hip_data$marks, model = 1
)
get_coefs(fit1_hip)
show_hawkes(fit1_hip)
show_hawkes_GOF(fit1_hip)
show_interarrival_pattern(fit1_hip)

# Fit all 7 days together
all_days_list <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
all_days_times_list <- list()
all_days_marks_list <- list()
for (day in all_days_list) {
  if (day %in% unique(daily_data$WEEKDAY)) {
    day_data <- transform_hawkes_data(daily_data, seqn, day, data_type = "hip")
    all_days_times_list[[day]] <- day_data$times
    all_days_marks_list[[day]] <- day_data$marks
  }
}

# Grid search using helper function
best_fit_hip <- perform_hawkes_grid_search(
  times_list = all_days_times_list,
  marks_list = all_days_marks_list,
  fit_function = "fit_hawkes_multi_series",
  model = 1
)

get_coefs(best_fit_hip)
show_hawkes(best_fit_hip, per_series = T)
show_hawkes_GOF(best_fit_hip, per_series = T)
show_interarrival_pattern(best_fit_hip, per_series = T)
'''
marked
                 Estimate  Std. Error
log_mu        -3.53923117 0.142457277
logit_abratio  1.11442941 0.247603727
log_beta      -2.87911489 0.184769209
mu             0.02903564 0.004136338
alpha          0.04230432 0.006912638
beta           0.05618447 0.010381160

unmarked
                 Estimate  Std. Error
log_mu        -3.63032989 0.150684833
logit_abratio  1.25082364 0.263095588
log_beta      -2.88064019 0.172224747
mu             0.02650744 0.003994269
alpha          0.04361362 0.006729868
beta           0.05609884 0.009661608
'''

# Manual grid search for experimentation
a_par_values <- c(-2, -1, 0, 1, 2)
mu_values <- c(1, 0.5)
beta_values <- c(0.1, 1)

best_fit_manual <- NULL
best_objective <- Inf
best_params <- list(a_par = NA, mu_mult = NA, beta_mult = NA)

for (i in seq_along(a_par_values)) {
  for (j in seq_along(mu_values)) {
    for (k in seq_along(beta_values)) {
      try({
        mu_val <- mu_values[j] * mean(sapply(all_days_times_list, function(x) length(x) / max(x)))
        beta_val <- beta_values[k] / median(unlist(lapply(all_days_times_list, function(x) diff(x))))

        fit_candidate <- fit_hawkes_multi_series(
          times_list = all_days_times_list,
          marks_list = all_days_marks_list,
          model = 2,
          parameters = list(
            mu = mu_val,
            alpha = a_par_values[i],
            beta = beta_val
          )
        )

        # Check convergence and update best fit
        if (!is.null(fit_candidate) && TMB::sdreport(fit_candidate)$pdHess &&
          fit_candidate$objective < best_objective) {
          best_objective <- fit_candidate$objective
          best_fit_manual <- fit_candidate
          best_params <- list(
            a_par = a_par_values[i],
            mu_mult = mu_values[j],
            beta_mult = beta_values[k]
          )
        }
      })
    }
  }
}

# Compare results
if (!is.null(best_fit_manual)) {
  cat("Manual grid search best parameters:\n")
  print(best_params)
  get_coefs(best_fit_manual)
}

#### 3 wrist ####
load("../2025/data/runlength/event_analysis_new.RData")
load("../2025/data/mims/data_analysis_new.RData")
load("../2025/data/mims/Act_Analysis_new.RData")

daily_data <- event_analysis
daily_data <- event_analysis %>% mutate(original_value = 1) 
seqn <- unique(daily_data$SEQN)[1]
seqn = 64161
weekday <- unique(daily_data$WEEKDAY)[1]
wrist_data <- transform_hawkes_data(daily_data, seqn, weekday, active_only = TRUE, data_type = "wrist")
fit1_wrist <- fit_hawkes(
  times = wrist_data$times,
  marks = wrist_data$marks, model = 1
)
get_coefs(fit1_wrist)
show_hawkes(fit1_wrist)
show_hawkes_GOF(fit1_wrist)
show_interarrival_pattern(fit1_wrist)

# Fit all 7 days together
all_days_list <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
all_days_times_list <- list()
all_days_marks_list <- list()
for (day in all_days_list) {
  if (day %in% unique(daily_data$WEEKDAY)) {
    day_data <- transform_hawkes_data(daily_data, seqn, day, data_type = "wrist")
    all_days_times_list[[day]] <- day_data$times
    all_days_marks_list[[day]] <- day_data$marks
  }
}

# Grid search using helper function
best_fit_wrist <- perform_hawkes_grid_search(
  times_list = all_days_times_list,
  marks_list = all_days_marks_list,
  fit_function = "fit_hawkes_multi_series",
  model = 1
)

get_coefs(best_fit_wrist)
show_hawkes(best_fit_wrist, per_series = T)
show_hawkes_GOF(best_fit_wrist, per_series = T)
show_interarrival_pattern(best_fit_wrist, per_series = T)

'''
marked
                 Estimate  Std. Error
log_mu        -3.51911042 0.106014236
logit_abratio  1.54708514 0.203297858
log_beta      -1.91632813 0.081682702
mu             0.02962578 0.003140754
alpha          0.12132097 0.009777489
beta           0.14714627 0.012019305

unmarked
                 Estimate  Std. Error
log_mu        -3.60581534 0.110728118
logit_abratio  1.65498515 0.217135757
log_beta      -1.93452643 0.079469233
mu             0.02716529 0.003007961
alpha          0.12131081 0.009559128
beta           0.14449268 0.011482722
'''

# Manual grid search for experimentation
a_par_values <- c(-5, -2, -1, 0, 1, 2)
mu_values <- c(5, 1, 0.5)
beta_values <- c(0.1, 1)

best_fit_manual_wrist <- NULL
best_objective <- Inf
best_params <- list(a_par = NA, mu_mult = NA, beta_mult = NA)

for (i in seq_along(a_par_values)) {
  for (j in seq_along(mu_values)) {
    for (k in seq_along(beta_values)) {
      try({
        mu_val <- mu_values[j] * mean(sapply(all_days_times_list, function(x) length(x) / max(x)))
        beta_val <- beta_values[k] / median(unlist(lapply(all_days_times_list, function(x) diff(x))))

        fit_candidate <- fit_hawkes_multi_series(
          times_list = all_days_times_list,
          marks_list = all_days_marks_list,
          model = 2,
          parameters = list(
            mu = mu_val,
            alpha = a_par_values[i],
            beta = beta_val
          )
        )

        # Check convergence and update best fit
        if (!is.null(fit_candidate) && TMB::sdreport(fit_candidate)$pdHess &&
          fit_candidate$objective < best_objective) {
          best_objective <- fit_candidate$objective
          best_fit_manual_wrist <- fit_candidate
          best_params <- list(
            a_par = a_par_values[i],
            mu_mult = mu_values[j],
            beta_mult = beta_values[k]
          )
        }
      })
    }
  }
}

# Compare results
if (!is.null(best_fit_manual_wrist)) {
  cat("Manual grid search best parameters (wrist):\n")
  print(best_params)
  get_coefs(best_fit_manual_wrist)
}

