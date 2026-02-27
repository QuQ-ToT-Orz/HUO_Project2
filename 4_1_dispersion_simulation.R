#### Simulation Test for Dispersion Index Functions ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(emhawkes)  # For proper Hawkes simulation

set.seed(42)

#### 1 Copy Core Functions from 4_Dispersion_index.R ####

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

compute_dispersion_single_day <- function(day_events, circadian_rates, circadian_rates_marks,
                                          window_sizes, bin_size, obs_start, obs_end) {
  n_bins <- length(circadian_rates)
  # Use fixed observation period instead of data-driven
  day_start <- obs_start
  obs_period <- obs_end - obs_start

  day_events <- day_events %>%
    mutate(bin_idx = pmin(floor((start - 1) / bin_size) + 1, n_bins))

  results <- list()
  for (W in window_sizes) {
    if (obs_period < W) next
    n_windows <- floor(obs_period / W)
    if (n_windows < 3) next

    window_stats <- day_events %>%
      mutate(window_idx = pmin(floor((start - day_start) / W) + 1, n_windows)) %>%
      filter(window_idx >= 1, window_idx <= n_windows) %>%
      group_by(window_idx) %>%
      summarise(
        raw_count = n(),
        raw_marks = sum(mark_norm),
        mu_window = mean(circadian_rates[bin_idx]),
        mu_window_marks = mean(circadian_rates_marks[bin_idx]),
        adj_count = raw_count / mu_window,
        adj_marks = raw_marks / mu_window,
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

    # Population variance (n denominator) for comparison against theoretical D
    pop_var <- function(x) mean((x - mean(x))^2)
    D_raw <- pop_var(all_windows$raw_count) / mean(all_windows$raw_count)
    D_marks <- pop_var(all_windows$raw_marks) / mean(all_windows$raw_marks)
    D_adj <- pop_var(all_windows$adj_count) / mean(all_windows$adj_count)
    D_adj_marks <- pop_var(all_windows$adj_marks) / mean(all_windows$adj_marks)

    results[[as.character(W)]] <- data.frame(
      window_size = W, D_raw = D_raw, D_marks = D_marks, D_adj = D_adj, D_adj_marks = D_adj_marks
    )
  }
  bind_rows(results)
}

compute_dispersion <- function(events_df, window_sizes, bin_size, obs_start = 480, obs_end = 1320) {
  day_results <- lapply(unique(events_df$WEEKDAY), function(day) {
    day_events <- events_df %>% filter(WEEKDAY == day)
    circadian_rates <- estimate_circadian_baseline(day_events, bin_size)
    circadian_rates_marks <- estimate_circadian_baseline_marks(day_events, bin_size)
    compute_dispersion_single_day(day_events, circadian_rates, circadian_rates_marks,
                                  window_sizes, bin_size, obs_start, obs_end)
  })

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
window_sizes <- c(5, 10, 15, 30, 45, 60, 90, 120)
bin_size <- 30

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

clustered_events <- generate_clustered_events(n_participants, n_days, rate_per_hour, branching_n = 0.5)
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
# A tibble: 8 × 4
  window_size mean_D_raw sd_D_raw mean_n_raw
        <dbl>      <dbl>    <dbl>      <dbl>
1           5       1.86    0.138      0.265
2          10       2.58    0.260      0.375
3          15       3.19    0.372      0.437
4          30       4.41    0.632      0.520
5          45       5.02    0.774      0.550
6          60       5.46    1.06       0.566
7          90       5.56    1.29       0.568
8         120       6.19    1.74       0.586
'''