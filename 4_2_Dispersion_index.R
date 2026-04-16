#### 6 Dispersion Index Analysis ####
# drop + inclusive + all-bin normalization + event-time mu + 
# bin_size=120 + inclusive lambda + adj_marks uses mu_window #
gc()
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)

#### 1 Load Data ####
load(file = paste("./data/runlength/", "event_analysis_old.RData", sep = ""))
load(file = paste("./data/count/", "data_analysis_old.RData", sep = ""))

load(file = paste("./data/runlength/", "event_analysis_new.RData", sep = ""))
load(file = paste("./data/mims/", "data_analysis_new.RData", sep = ""))

circadian_bin_size <- 120
window_sizes <- c(15, 30, 60, 90, 120)

#### 2 Core Functions ####
sum_by_group <- function(values, group_idx, n_groups) {
  totals <- numeric(n_groups)
  grouped_sums <- tapply(values, group_idx, sum)
  totals[as.integer(names(grouped_sums))] <- as.numeric(grouped_sums)
  totals
}

# Single-day circadian baseline (normalized to mean=1) - count-based
estimate_circadian_baseline <- function(day_events, bin_size) {
  n_bins <- 1440 / bin_size
  bin_idx <- pmin(floor((day_events$start - 1) / bin_size) + 1, n_bins)
  bin_rates <- tabulate(bin_idx, nbins = n_bins)

  # Normalize across all bins so the daily expectation is properly distributed
  bin_rates <- bin_rates / mean(bin_rates)

  attr(bin_rates, "bin_size") <- bin_size
  return(bin_rates)
}

# Single-day circadian baseline (normalized to mean=1) - mark-weighted
estimate_circadian_baseline_marks <- function(day_events, bin_size) {
  n_bins <- 1440 / bin_size
  bin_idx <- pmin(floor((day_events$start - 1) / bin_size) + 1, n_bins)
  bin_rates <- sum_by_group(day_events$mark_norm, bin_idx, n_bins)

  bin_rates <- bin_rates / mean(bin_rates)

  attr(bin_rates, "bin_size") <- bin_size
  return(bin_rates)
}

# Dispersion for single day
compute_dispersion_single_day <- function(day_events, circadian_rates, circadian_rates_marks, window_sizes, bin_size) {
  event_starts <- day_events$start
  event_marks <- day_events$mark_norm
  day_start <- min(day_events$start)
  day_end <- max(day_events$start)
  obs_len <- day_end - day_start + 1
  n_bins <- length(circadian_rates)
  bin_idx_all <- pmin(floor((event_starts - 1) / bin_size) + 1, n_bins)

  results <- list()
  for (W in window_sizes) {
    if (obs_len < W) next
    n_windows <- floor(obs_len / W)
    if (n_windows < 3) next

    window_idx_all <- floor((event_starts - day_start) / W) + 1
    valid_events <- window_idx_all >= 1 & window_idx_all <= n_windows
    window_idx <- window_idx_all[valid_events]

    raw_count <- tabulate(window_idx, nbins = n_windows)
    raw_marks <- sum_by_group(event_marks[valid_events], window_idx, n_windows)

    mu_window <- rep(NA_real_, n_windows)
    mu_vals <- tapply(circadian_rates[bin_idx_all[valid_events]], window_idx, mean)
    mu_window[as.integer(names(mu_vals))] <- as.numeric(mu_vals)

    adj_count <- ifelse(!is.na(mu_window) & mu_window > 0, raw_count / mu_window, 0)
    adj_marks <- ifelse(!is.na(mu_window) & mu_window > 0, raw_marks / mu_window, 0)

    # D_raw: count-based (each event = 1)
    D_raw <- var(raw_count) / mean(raw_count)

    # D_marks: mark-based without circadian adjustment
    D_marks <- var(raw_marks) / mean(raw_marks)

    # D_adj
    D_adj <- var(adj_count) / mean(adj_count)

    # D_adj_marks
    D_adj_marks <- var(adj_marks) / mean(adj_marks)

    results[[as.character(W)]] <- data.frame(
      window_size = W, D_raw = D_raw, D_marks = D_marks, D_adj = D_adj, D_adj_marks = D_adj_marks
    )
  }
  bind_rows(results)
}

# Per day then average (circadian computed per day)
compute_dispersion <- function(day_events_list, window_sizes, bin_size) {

  day_results <- lapply(day_events_list, function(day_events) {
    circadian_rates <- estimate_circadian_baseline(day_events, bin_size)
    circadian_rates_marks <- estimate_circadian_baseline_marks(day_events, bin_size)
    compute_dispersion_single_day(day_events, circadian_rates, circadian_rates_marks, window_sizes, bin_size)
  })

  day_results <- Filter(function(x) nrow(x) > 0, day_results)
  if (length(day_results) == 0) {
    return(data.frame(
      window_size = numeric(0),
      D_raw = numeric(0),
      D_marks = numeric(0),
      D_adj = numeric(0),
      D_adj_marks = numeric(0)
    ))
  }

  bind_rows(day_results) %>%
    group_by(window_size) %>%
    summarise(D_raw = mean(D_raw, na.rm = TRUE),
              D_marks = mean(D_marks, na.rm = TRUE),
              D_adj = mean(D_adj, na.rm = TRUE),
              D_adj_marks = mean(D_adj_marks, na.rm = TRUE),
              .groups = "drop")
}

# Branching ratio: n = 1 - 1/sqrt(D)
# No capping: n < 0 for D < 1 (regularity), n > 0 for D > 1 (clustering)
branching_ratio <- function(D) {
  1 - 1/sqrt(D)
}

#### 3 Compute for All Participants ####

# Filter active events
active_events <- event_analysis %>%
  filter(categories != "sedentary") %>%
  group_by(SEQN, WEEKDAY) %>%
  mutate(
    # sqrt for hip
    mark_sqrt = (original_value),
    daily_mean_sqrt = mean(mark_sqrt),
    mark_norm = mark_sqrt / daily_mean_sqrt
  ) %>%
  ungroup()

daily_subject_metrics <- active_events %>%
  group_by(SEQN, WEEKDAY) %>%
  summarise(
    Q_day = mean(mark_norm^2) / mean(mark_norm),
    n_events = n(),
    sum_marks = sum(mark_norm),
    obs_period = max(start) - min(start) + 1,
    .groups = "drop"
  ) %>%
  mutate(
    obs_hours = obs_period / 60,
    lambda_count_day = n_events / obs_hours,
    lambda_marks_day = sum_marks / obs_hours
  )

person_metrics <- daily_subject_metrics %>%
  group_by(SEQN) %>%
  summarise(
    Q = mean(Q_day, na.rm = TRUE),
    lambda_count = mean(lambda_count_day, na.rm = TRUE),
    lambda_marks = mean(lambda_marks_day, na.rm = TRUE),
    .groups = "drop"
  )

person_events_list <- split(active_events, active_events$SEQN, drop = TRUE)
cat("n_participants:", length(person_events_list), "\n")

dispersion_results <- lapply(person_events_list, function(person_events) {
  person_disp <- compute_dispersion(split(person_events, person_events$WEEKDAY, drop = TRUE), window_sizes, circadian_bin_size)
  person_disp$SEQN <- person_events$SEQN[[1]]
  person_disp
})

dispersion_df <- bind_rows(dispersion_results) %>%
  left_join(select(person_metrics, SEQN, Q), by = "SEQN") %>%
  mutate(
    n_raw = branching_ratio(D_raw),
    n_adj = branching_ratio(D_adj),
    n_marks = 1 - sqrt(Q / D_marks),
    n_adj_marks = 1 - sqrt(Q / D_adj_marks)
  )

#### 4 Window Size Sensitivity ####
window_sensitivity <- dispersion_df %>%
  group_by(window_size) %>%
  summarise(
    mean_D_raw = mean(D_raw, na.rm = TRUE),
    sd_D_raw = sd(D_raw, na.rm = TRUE),
    mean_D_marks = mean(D_marks, na.rm = TRUE),
    sd_D_marks = sd(D_marks, na.rm = TRUE),
    mean_D_adj = mean(D_adj, na.rm = TRUE),
    sd_D_adj = sd(D_adj, na.rm = TRUE),
    mean_D_adj_marks = mean(D_adj_marks, na.rm = TRUE),
    sd_D_adj_marks = sd(D_adj_marks, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff_count = mean_D_raw - mean_D_adj,
    diff_marks = mean_D_marks - mean_D_adj_marks
  )

print(window_sensitivity)

p_windows <- ggplot(window_sensitivity, aes(x = window_size)) +
  # D_raw with error ribbon
  geom_ribbon(aes(ymin = mean_D_raw - sd_D_raw, ymax = mean_D_raw + sd_D_raw),
              fill = "red", alpha = 0.2) +
  geom_line(aes(y = mean_D_raw, color = "D_raw"), linewidth = 1) +
  geom_point(aes(y = mean_D_raw, color = "D_raw"), size = 3) +
  # D_marks with error ribbon
  geom_ribbon(aes(ymin = mean_D_marks - sd_D_marks, ymax = mean_D_marks + sd_D_marks),
              fill = "orange", alpha = 0.2) +
  geom_line(aes(y = mean_D_marks, color = "D_marks"), linewidth = 1) +
  geom_point(aes(y = mean_D_marks, color = "D_marks"), size = 3) +
  # D_adj with error ribbon
  geom_ribbon(aes(ymin = mean_D_adj - sd_D_adj, ymax = mean_D_adj + sd_D_adj),
              fill = "blue", alpha = 0.2) +
  geom_line(aes(y = mean_D_adj, color = "D_adj"), linewidth = 1) +
  geom_point(aes(y = mean_D_adj, color = "D_adj"), size = 3) +
  # D_adj_marks with error ribbon
  geom_ribbon(aes(ymin = mean_D_adj_marks - sd_D_adj_marks, ymax = mean_D_adj_marks + sd_D_adj_marks),
              fill = "green", alpha = 0.2) +
  geom_line(aes(y = mean_D_adj_marks, color = "D_adj_marks"), linewidth = 1) +
  geom_point(aes(y = mean_D_adj_marks, color = "D_adj_marks"), size = 3) +
  scale_color_manual(values = c("D_raw" = "red", "D_marks" = "orange", "D_adj" = "blue", "D_adj_marks" = "darkgreen")) +
  labs(x = "Window Size (min)", y = "Dispersion Index", color = "") +
  theme_minimal()
print(p_windows)
ggsave("Output/dispersion/hip_windows.pdf", p_windows, width = 10, height = 6)
ggsave("Output/dispersion/wrist_windows.pdf", p_windows, width = 10, height = 6)

#### 5 Create Summary and Save ####
primary_window <- 30

dispersion_summary <- dispersion_df %>%
  filter(window_size == primary_window) %>%
  left_join(select(person_metrics, SEQN, lambda_count, lambda_marks), by = "SEQN") %>%
  mutate(
    # Immigration rates using Hawkes relationship: μ* = λ(1-n)
    # 4 versions to match 4 versions of n
    mu_star_raw = lambda_count * (1 - n_raw),           # count, no circadian adj
    mu_star_adj = lambda_count * (1 - n_adj),           # count, circadian adj
    mu_star_marks = lambda_marks * (1 - n_marks),       # marks, no circadian adj
    mu_star_adj_marks = lambda_marks * (1 - n_adj_marks) # marks, circadian adj
  ) %>%
  select(SEQN, Q, D_raw, D_marks, D_adj, D_adj_marks,
         n_raw, n_marks, n_adj, n_adj_marks,
         lambda_count, lambda_marks,
         mu_star_raw, mu_star_adj, mu_star_marks, mu_star_adj_marks)

# Save results
dir_path <- "./data/dispersion/"
save(dispersion_df, file = paste0(dir_path, "dispersion_df_old.RData"))
save(dispersion_summary, file = paste0(dir_path, "dispersion_summary_old.RData"))

save(dispersion_df, file = paste0(dir_path, "dispersion_df_new.RData"))
save(dispersion_summary, file = paste0(dir_path, "dispersion_summary_new.RData"))
