#### Simulation Test for Dispersion Index Functions ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(emhawkes)  # For proper Hawkes simulation

set.seed(42)

#### 1 Copy Core Functions from 4_Dispersion_index.R ####
# drop + inclusive + all-bin normalization + event-time mu + adj_marks uses mu_window
sum_by_group <- function(values, group_idx, n_groups) {
  totals <- numeric(n_groups)
  if (length(values) == 0) return(totals)

  grouped_sums <- tapply(values, group_idx, sum)
  totals[as.integer(names(grouped_sums))] <- as.numeric(grouped_sums)
  totals
}

estimate_circadian_baseline <- function(day_events, bin_size) {
  n_bins <- 1440 / bin_size
  bin_idx <- pmin(floor((day_events$start - 1) / bin_size) + 1, n_bins)
  bin_rates <- tabulate(bin_idx, nbins = n_bins)

  bin_rates <- bin_rates / mean(bin_rates)

  attr(bin_rates, "bin_size") <- bin_size
  return(bin_rates)
}

estimate_circadian_baseline_marks <- function(day_events, bin_size) {
  n_bins <- 1440 / bin_size
  bin_idx <- pmin(floor((day_events$start - 1) / bin_size) + 1, n_bins)
  bin_rates <- sum_by_group(day_events$mark_norm, bin_idx, n_bins)

  bin_rates <- bin_rates / mean(bin_rates)

  attr(bin_rates, "bin_size") <- bin_size
  return(bin_rates)
}

compute_dispersion_single_day <- function(day_events, circadian_rates, circadian_rates_marks,
                                          window_sizes, bin_size, obs_start, obs_end) {
  event_starts <- day_events$start
  event_marks <- day_events$mark_norm
  day_start <- obs_start
  day_end <- obs_end - 1
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

    D_raw <- var(raw_count) / mean(raw_count)
    D_marks <- var(raw_marks) / mean(raw_marks)
    D_adj <- var(adj_count) / mean(adj_count)
    D_adj_marks <- var(adj_marks) / mean(adj_marks)

    results[[as.character(W)]] <- data.frame(
      window_size = W, D_raw = D_raw, D_marks = D_marks, D_adj = D_adj, D_adj_marks = D_adj_marks
    )
  }
  bind_rows(results)
}

compute_dispersion <- function(events_df, window_sizes, bin_size, obs_start = 480, obs_end = 1320) {
  day_results <- lapply(split(events_df, events_df$WEEKDAY, drop = TRUE), function(day_events) {
    circadian_rates <- estimate_circadian_baseline(day_events, bin_size)
    circadian_rates_marks <- estimate_circadian_baseline_marks(day_events, bin_size)
    compute_dispersion_single_day(day_events, circadian_rates, circadian_rates_marks,
                                  window_sizes, bin_size, obs_start, obs_end)
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

branching_ratio <- function(D) {
  1 - 1/sqrt(D)
}

#### 2 Simulation Functions ####

# Generate Poisson (homogeneous) events - evenly distributed in time
generate_poisson_events <- function(n_participants, n_days, rate_per_hour,
                                     obs_start = 480, obs_end = 1320) {
  obs_period <- obs_end - obs_start
  expected_events <- rate_per_hour * (obs_period / 60)

  all_events <- list()

  for (seqn in 1:n_participants) {
    for (day in 1:n_days) {
      n_events <- rpois(1, expected_events)
      if (n_events == 0) n_events <- 1

      event_times <- sort(runif(n_events, min = obs_start, max = obs_end))
      marks <- rexp(n_events, rate = 1)

      day_events <- data.frame(
        SEQN = seqn,
        WEEKDAY = day,
        start = event_times,
        original_value = marks,
        categories = "active"
      ) %>%
        mutate(
          mark_sqrt = original_value,
          daily_mean_sqrt = mean(mark_sqrt),
          mark_norm = mark_sqrt / daily_mean_sqrt
        )

      all_events[[paste(seqn, day, sep = "_")]] <- day_events
    }
  }

  bind_rows(all_events)
}

# Generate clustered events using emhawkes package (proper Hawkes simulation)
generate_clustered_events <- function(n_participants, n_days, base_rate,
                                       branching_n = 0.6, obs_start = 480, obs_end = 1320) {
  obs_period <- obs_end - obs_start  # in minutes

  # Hawkes parameters (time in minutes)
  # mu = baseline rate per minute, alpha/beta = branching ratio
  # We want lambda = base_rate events/hour = base_rate/60 events/min
  # lambda = mu / (1-n), so mu = lambda * (1-n)
  lambda_per_min <- base_rate / 60
  mu <- lambda_per_min * (1 - branching_n)
  beta <- 0.15  # decay rate (per minute)
  alpha <- branching_n * beta  # so that n = alpha/beta

  # Create Hawkes specification
  h <- new("hspec", mu = mu, alpha = alpha, beta = beta)

  all_events <- list()

  for (seqn in 1:n_participants) {
    for (day in 1:n_days) {
      # Simulate Hawkes process
      # hsim generates events until size events or horizon time
      res <- hsim(h, size = 300, lambda_component0 = mu * 10)

      # Get arrival times and shift to observation window
      arrival_times <- res$arrival

      # Keep only events within observation period
      event_times <- arrival_times[arrival_times <= obs_period]
      event_times <- event_times + obs_start  # shift to observation window

      if (length(event_times) == 0) {
        event_times <- runif(1, obs_start, obs_end)
      }

      n_events <- length(event_times)
      marks <- rexp(n_events, rate = 1)

      day_events <- data.frame(
        SEQN = seqn,
        WEEKDAY = day,
        start = event_times,
        original_value = marks,
        categories = "active"
      ) %>%
        mutate(
          mark_sqrt = original_value,
          daily_mean_sqrt = mean(mark_sqrt),
          mark_norm = mark_sqrt / daily_mean_sqrt
        )

      all_events[[paste(seqn, day, sep = "_")]] <- day_events
    }
  }

  bind_rows(all_events)
}

# Generate regular (under-dispersed) events
generate_regular_events <- function(n_participants, n_days, n_events_per_day = 300,
                                     obs_start = 480, obs_end = 1320) {
  obs_period <- obs_end - obs_start

  all_events <- list()

  for (seqn in 1:n_participants) {
    for (day in 1:n_days) {
      spacing <- obs_period / n_events_per_day
      event_times <- seq(obs_start, obs_end - spacing, length.out = n_events_per_day)
      event_times <- event_times + runif(n_events_per_day, -spacing * 0.1, spacing * 0.1)
      event_times <- pmax(obs_start, pmin(obs_end, event_times))
      event_times <- sort(event_times)

      n_events <- length(event_times)
      marks <- rexp(n_events, rate = 1)

      day_events <- data.frame(
        SEQN = seqn,
        WEEKDAY = day,
        start = event_times,
        original_value = marks,
        categories = "active"
      ) %>%
        mutate(
          mark_sqrt = original_value,
          daily_mean_sqrt = mean(mark_sqrt),
          mark_norm = mark_sqrt / daily_mean_sqrt
        )

      all_events[[paste(seqn, day, sep = "_")]] <- day_events
    }
  }

  bind_rows(all_events)
}

#### 3 Run Simulations ####

cat("=== Simulation Tests for Dispersion Index Functions ===\n\n")

# Parameters
n_participants <- 50
n_days <- 7
rate_per_hour <- 20
window_sizes <- c(5, 10, 15, 30, 45, 60, 75, 90, 120)
bin_size <- 120

# Test 1: Poisson (uniform) process - expected D ≈ 1, n ≈ 0
cat("Test 1: Poisson Process (evenly distributed)\n")
cat("Expected: D ≈ 1, n ≈ 0\n")

poisson_events <- generate_poisson_events(n_participants, n_days, rate_per_hour)
cat("Generated", nrow(poisson_events), "events for", n_participants, "participants\n")
cat("Mean events per day:", mean(table(paste(poisson_events$SEQN, poisson_events$WEEKDAY))), "\n\n")

poisson_results <- lapply(1:n_participants, function(seqn) {
  person_events <- poisson_events %>% filter(SEQN == seqn)
  person_disp <- compute_dispersion(person_events, window_sizes, bin_size)
  person_disp$SEQN <- seqn
  person_disp$n_raw <- branching_ratio(person_disp$D_raw)
  person_disp
})

poisson_df <- bind_rows(poisson_results)

poisson_summary <- poisson_df %>%
  group_by(window_size) %>%
  summarise(
    mean_D_raw = mean(D_raw, na.rm = TRUE),
    sd_D_raw = sd(D_raw, na.rm = TRUE),
    mean_n_raw = mean(n_raw, na.rm = TRUE),
    .groups = "drop"
  )

# Test 2: Clustered process - expected D > 1, n > 0
cat("\nTest 2: Clustered Process (Hawkes-like)\n")
cat("Expected: D > 1, n > 0 (clustering)\n")

clustered_events <- generate_clustered_events(n_participants, n_days, rate_per_hour)
cat("Generated", nrow(clustered_events), "events\n")
cat("Mean events per day:", mean(table(paste(clustered_events$SEQN, clustered_events$WEEKDAY))), "\n\n")

clustered_results <- lapply(1:n_participants, function(seqn) {
  person_events <- clustered_events %>% filter(SEQN == seqn)
  person_disp <- compute_dispersion(person_events, window_sizes, bin_size)
  person_disp$SEQN <- seqn
  person_disp$n_raw <- branching_ratio(person_disp$D_raw)
  person_disp
})

clustered_df <- bind_rows(clustered_results)

clustered_summary <- clustered_df %>%
  group_by(window_size) %>%
  summarise(
    mean_D_raw = mean(D_raw, na.rm = TRUE),
    sd_D_raw = sd(D_raw, na.rm = TRUE),
    mean_n_raw = mean(n_raw, na.rm = TRUE),
    .groups = "drop"
  )

# Test 3: Regular (underdispersed) process - expected D < 1, n < 0
cat("\nTest 3: Regular Process (evenly spaced)\n")
cat("Expected: D < 1, n < 0 (regularity)\n")

regular_events <- generate_regular_events(n_participants, n_days, n_events_per_day = 280)

regular_results <- lapply(1:n_participants, function(seqn) {
  person_events <- regular_events %>% filter(SEQN == seqn)
  person_disp <- compute_dispersion(person_events, window_sizes, bin_size)
  person_disp$SEQN <- seqn
  person_disp$n_raw <- branching_ratio(person_disp$D_raw)
  person_disp
})

regular_df <- bind_rows(regular_results)

regular_summary <- regular_df %>%
  group_by(window_size) %>%
  summarise(
    mean_D_raw = mean(D_raw, na.rm = TRUE),
    sd_D_raw = sd(D_raw, na.rm = TRUE),
    mean_n_raw = mean(n_raw, na.rm = TRUE),
    .groups = "drop"
  )

#### 4 Visualization ####

all_results <- bind_rows(
  poisson_summary %>% mutate(process = "Poisson"),
  clustered_summary %>% mutate(process = "Clustered"),
  regular_summary %>% mutate(process = "Regular")
)

# Plot D_raw by window size
p1 <- ggplot(all_results, aes(x = window_size, y = mean_D_raw, color = process)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = mean_D_raw - sd_D_raw, ymax = mean_D_raw + sd_D_raw, fill = process),
              alpha = 0.2, color = NA) +
  labs(
    title = "Dispersion Index by Process Type",
    x = "Window Size (minutes)",
    y = "D_raw (Dispersion Index)",
    color = "Process",
    fill = "Process"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)
ggsave("Output/dispersion/simulation.pdf", p1, width = 12, height = 8)

poisson_summary
clustered_summary
regular_summary
'''
=== Simulation Tests for Dispersion Index Functions ===
Test 1: Poisson Process (evenly distributed)
Expected: D ≈ 1, n ≈ 0
Generated 98267 events for 50 participants
Mean events per day: 280.7629 

Test 2: Clustered Process (Hawkes-like)
Expected: D > 1, n > 0 (clustering)
Generated 98830 events
Mean events per day: 282.3714 

# A tibble: 9 × 4
  window_size mean_D_raw sd_D_raw mean_n_raw
        <dbl>      <dbl>    <dbl>      <dbl>
1           5       1.97    0.176      0.286
2          10       2.82    0.345      0.402
3          15       3.53    0.501      0.464
4          30       5.22    0.916      0.558
5          45       6.06    1.18       0.589
6          60       7.06    1.62       0.617
7          75       7.55    1.88       0.629
8          90       7.67    2.01       0.631
9         120       8.88    2.86       0.653
'''