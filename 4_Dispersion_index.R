#### 6 Dispersion Index Analysis ####
gc()
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)

#### 1 Load Data ####
load(file = paste("../2025/data/runlength/", "event_analysis_old.RData", sep = ""))
load(file = paste("../2025/data/count/", "data_analysis_old.RData", sep = ""))

load(file = paste("../2025/data/runlength/", "event_analysis_new.RData", sep = ""))
load(file = paste("../2025/data/mims/", "data_analysis_new.RData", sep = ""))

circadian_bin_size <- 30
window_sizes <- c(15, 30, 60, 90, 120)

#### 2 Core Functions ####
# Single-day circadian baseline (normalized to mean=1) - count-based
estimate_circadian_baseline <- function(day_events, bin_size) {
  n_bins <- 1440 / bin_size
  day_events <- day_events %>%
    mutate(bin_idx = pmin(floor((start - 1) / bin_size) + 1, n_bins))

  bin_counts <- day_events %>%
    group_by(bin_idx) %>%
    summarise(count = n(), .groups = "drop")

  bin_rates <- rep(0, n_bins)
  bin_rates[bin_counts$bin_idx] <- bin_counts$count

  observed_bins <- which(bin_rates > 0)
  bin_rates <- bin_rates / mean(bin_rates[observed_bins])

  attr(bin_rates, "bin_size") <- bin_size
  return(bin_rates)
}

# Single-day circadian baseline (normalized to mean=1) - mark-weighted
estimate_circadian_baseline_marks <- function(day_events, bin_size) {
  n_bins <- 1440 / bin_size
  day_events <- day_events %>%
    mutate(bin_idx = pmin(floor((start - 1) / bin_size) + 1, n_bins))

  bin_marks <- day_events %>%
    group_by(bin_idx) %>%
    summarise(total_marks = sum(mark_norm), .groups = "drop")

  bin_rates <- rep(0, n_bins)
  bin_rates[bin_marks$bin_idx] <- bin_marks$total_marks

  observed_bins <- which(bin_rates > 0)
  bin_rates <- bin_rates / mean(bin_rates[observed_bins])

  attr(bin_rates, "bin_size") <- bin_size
  return(bin_rates)
}

# Dispersion for single day
compute_dispersion_single_day <- function(day_events, circadian_rates, circadian_rates_marks, window_sizes, bin_size) {
  n_bins <- length(circadian_rates)
  day_start <- min(day_events$start)
  obs_period <- max(day_events$start) - day_start

  # Add bin_idx
  day_events <- day_events %>%
    mutate(
      bin_idx = pmin(floor((start - 1) / bin_size) + 1, n_bins)
    )

  results <- list()
  for (W in window_sizes) {
    if (obs_period < W) next
    n_windows <- floor(obs_period / W)
    if (n_windows < 3) next

    window_stats <- day_events %>%
      mutate(window_idx = pmin(floor((start - day_start) / W) + 1, n_windows)) %>%
      filter(window_idx <= n_windows) %>%
      group_by(window_idx) %>%
      summarise(
        raw_count = n(),
        raw_marks = sum(mark_norm),  # sum of normalized marks
        mu_window = mean(circadian_rates[bin_idx]),  # window's average circadian rate (count-based)
        mu_window_marks = mean(circadian_rates_marks[bin_idx]),  # window's average circadian rate (mark-weighted)
        adj_count = raw_count / mu_window,  
        adj_marks = raw_marks / mu_window_marks, 
        .groups = "drop"
      )

    all_windows <- data.frame(window_idx = 1:n_windows) %>%
      left_join(window_stats, by = "window_idx") %>%
      mutate(
        raw_count = ifelse(is.na(raw_count), 0, raw_count),
        raw_marks = ifelse(is.na(raw_marks), 0, raw_marks),
        adj_count = ifelse(is.na(adj_count), 0, adj_count),
        adj_marks = ifelse(is.na(adj_marks), 0, adj_marks)
      )

    # D_raw: count-based (each event = 1)
    D_raw <- var(all_windows$raw_count) / mean(all_windows$raw_count)

    # D_marks: mark-based without circadian adjustment
    D_marks <- var(all_windows$raw_marks) / mean(all_windows$raw_marks)

    # D_adj
    D_adj <- var(all_windows$adj_count) / mean(all_windows$adj_count)

    # D_adj_marks
    D_adj_marks <- var(all_windows$adj_marks) / mean(all_windows$adj_marks)

    results[[as.character(W)]] <- data.frame(
      window_size = W, D_raw = D_raw, D_marks = D_marks, D_adj = D_adj, D_adj_marks = D_adj_marks
    )
  }
  bind_rows(results)
}

# Per day then average (circadian computed per day)
compute_dispersion <- function(events_df, window_sizes, bin_size) {

  day_results <- lapply(unique(events_df$WEEKDAY), function(day) {
    day_events <- events_df %>% filter(WEEKDAY == day)
    # Compute circadian baselines for THIS day only
    circadian_rates <- estimate_circadian_baseline(day_events, bin_size)
    circadian_rates_marks <- estimate_circadian_baseline_marks(day_events, bin_size)
    compute_dispersion_single_day(day_events, circadian_rates, circadian_rates_marks, window_sizes, bin_size)
  })

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

seqn_list <- unique(active_events$SEQN)
cat("n_participants:", length(seqn_list), "\n")

dispersion_results <- lapply(seqn_list, function(seqn) {
  person_events <- active_events %>% filter(SEQN == seqn)

  person_disp <- compute_dispersion(person_events, window_sizes, circadian_bin_size)
  person_disp$SEQN <- seqn

  # Compute Q = E[m^2]/E[m] per day, then average across days
  daily_Q <- sapply(unique(person_events$WEEKDAY), function(day) {
    day_events <- person_events %>% filter(WEEKDAY == day)
    m <- day_events$mark_norm
    mean(m^2) / mean(m)
  })
  Q <- mean(daily_Q, na.rm = TRUE)
  person_disp$Q <- Q

  # Uncorrected branching ratios
  person_disp$n_raw <- branching_ratio(person_disp$D_raw)
  person_disp$n_adj <- branching_ratio(person_disp$D_adj)

  # Q-corrected branching ratios for marks: n = 1 - sqrt(Q/D)
  person_disp$n_marks <- 1 - sqrt(Q / person_disp$D_marks)
  person_disp$n_adj_marks <- 1 - sqrt(Q / person_disp$D_adj_marks)

  person_disp
})

dispersion_df <- bind_rows(dispersion_results)

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

ggplot(window_sensitivity, aes(x = window_size)) +
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

#### 5 Create Summary and Save ####
window_diagnostics <- active_events %>%
    group_by(SEQN, WEEKDAY) %>%
    summarise(
      obs_period = max(start) - min(start),
      n_events = n(),
      .groups = "drop"
    ) %>%
    mutate(
      windows_90 = floor(obs_period / 90),
      windows_60 = floor(obs_period / 60),
      events_per_window_90 = n_events / windows_90,
      events_per_window_60 = n_events / windows_60
    )

summary(window_diagnostics[, c("events_per_window_90", "events_per_window_60")])

primary_window <- 30

dispersion_summary <- dispersion_df %>%
  filter(window_size == primary_window) %>%
  mutate(
    # Compute raw event rate λ (events per hour)
    lambda_count = sapply(SEQN, function(id) {
      person_events <- active_events %>% filter(SEQN == id)

      daily_rates <- sapply(unique(person_events$WEEKDAY), function(day) {
        day_events <- person_events %>% filter(WEEKDAY == day)
        n_events <- nrow(day_events)
        obs_period <- max(day_events$start) - min(day_events$start)  # in minutes
        obs_hours <- obs_period / 60  # convert to hours
        n_events / obs_hours
      })

      mean(daily_rates, na.rm = TRUE)  # average across days
    }),
    # Compute raw mark rate λ_marks (sum of normalized marks per hour)
    lambda_marks = sapply(SEQN, function(id) {
      person_events <- active_events %>% filter(SEQN == id)

      daily_rates <- sapply(unique(person_events$WEEKDAY), function(day) {
        day_events <- person_events %>% filter(WEEKDAY == day)
        sum_marks <- sum(day_events$mark_norm)
        obs_period <- max(day_events$start) - min(day_events$start)  # in minutes
        obs_hours <- obs_period / 60  # convert to hours
        sum_marks / obs_hours
      })

      mean(daily_rates, na.rm = TRUE)  # average across days
    })
  ) %>%
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
dir_path <- "../2025/data/dispersion/"
save(dispersion_df, file = paste0(dir_path, "dispersion_df_old.RData"))
save(dispersion_summary, file = paste0(dir_path, "dispersion_summary_old.RData"))

save(dispersion_df, file = paste0(dir_path, "dispersion_df_new.RData"))
save(dispersion_summary, file = paste0(dir_path, "dispersion_summary_new.RData"))
