#### Simulation for Dispersion Index ####
library(dplyr)
library(tidyr)
library(ggplot2)

#### 1 Simulation Functions ####
# Generate Poisson (homogeneous) events
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

# Periodic time-of-day multiplier, normalized to average 1 over the observed interval.
# amplitude = 0 recovers a constant baseline.
make_tod_multiplier <- function(amplitude = 0, peak_time = 14 * 60,
                                obs_start = 480, obs_end = 1320,
                                period = 1440, sharpness = 1) {

  peak_time <- as.numeric(peak_time)
  omega <- 2 * pi / period
  unscaled_multiplier <- function(time) {
    (1 + amplitude * cos(omega * (time - peak_time)))^sharpness
  }
  minute_midpoints <- seq(obs_start + 0.5, obs_end - 0.5, by = 1)
  normalizer <- mean(unscaled_multiplier(minute_midpoints))

  multiplier <- function(time) {
    unscaled_multiplier(time) / normalizer
  }
  attr(multiplier, "maximum") <-
    (1 + abs(amplitude))^sharpness / normalizer
  multiplier
}

# Expected total Hawkes rate under the periodic immigrant baseline and an exponential offspring kernel.
make_expected_hawkes_tod_multiplier <- function(
    branching_n, beta, amplitude = 0, peak_time = 14 * 60,
    obs_start = 480, obs_end = 1320, period = 1440, sharpness = 1) {
  immigrant_multiplier <- make_tod_multiplier(
    amplitude = amplitude,
    peak_time = peak_time,
    obs_start = obs_start,
    obs_end = obs_end,
    period = period,
    sharpness = sharpness
  )

  grid_step <- 1
  grid_times <- seq(grid_step / 2, period - grid_step / 2, by = grid_step)
  immigrant_grid <- immigrant_multiplier(grid_times)
  decay <- exp(-beta * (1 - branching_n) * grid_step)
  forcing_weight <- branching_n / (1 - branching_n) * (1 - decay)

  offspring_component <- 0
  for (immigrant_rate in immigrant_grid) {
    offspring_component <-
      decay * offspring_component + forcing_weight * immigrant_rate
  }
  offspring_component <- offspring_component / (1 - decay^length(grid_times))

  expected_grid <- numeric(length(grid_times))
  half_decay <- exp(-beta * (1 - branching_n) * grid_step / 2)
  half_forcing <- branching_n / (1 - branching_n) * (1 - half_decay)
  for (index in seq_along(grid_times)) {
    expected_grid[index] <- (1 - branching_n) * (
      immigrant_grid[index] +
        half_decay * offspring_component +
        half_forcing * immigrant_grid[index]
    )
    offspring_component <-
      decay * offspring_component + forcing_weight * immigrant_grid[index]
  }

  periodic_expected_rate <- function(time) {
    approx(
      x = c(grid_times - period, grid_times, grid_times + period),
      y = rep(expected_grid, 3),
      xout = time,
      rule = 2
    )$y
  }

  observed_midpoints <- seq(obs_start + 0.5, obs_end - 0.5, by = 1)
  normalizer <- mean(periodic_expected_rate(observed_midpoints))
  function(time) periodic_expected_rate(time) / normalizer
}

simulate_inhomogeneous_poisson <- function(rate_function, max_rate,
                                           interval_start, interval_end) {
  n_candidates <- rpois(1, max_rate * (interval_end - interval_start))

  candidates <- runif(n_candidates, interval_start, interval_end)
  accepted <- runif(n_candidates) <= rate_function(candidates) / max_rate
  sort(candidates[accepted])
}

# Cluster representation of an exponential Hawkes process.
# simulated marks E[M] = 1
simulate_hawkes_cluster_day <- function(
    base_rate, branching_n = 0.6, beta = 0.15,
    obs_start = 480, obs_end = 1320,
    tod_amplitude = 0, tod_peak = 14 * 60, tod_sharpness = 1,
    marked_excitation = FALSE, mark_shape = 1,
    burn_in = max(120, 15 / beta), max_events = 1000) {
  if (branching_n < 0 || branching_n >= 1) {
    stop("branching_n must be in [0, 1)")
  }
  if (base_rate <= 0 || beta <= 0) stop("base_rate and beta must be positive")

  immigrant_multiplier <- make_tod_multiplier(
    amplitude = tod_amplitude,
    peak_time = tod_peak,
    obs_start = obs_start,
    obs_end = obs_end,
    sharpness = tod_sharpness
  )
  mu0 <- (base_rate / 60) * (1 - branching_n)
  immigrant_rate <- function(time) mu0 * immigrant_multiplier(time)

  # Marks have E[M] = 1 and Var(M) = 1 / mark_shape. 
  draw_marks <- function(n) {
    if (is.infinite(mark_shape)) rep(1, n)
    else rgamma(n, shape = mark_shape, rate = mark_shape)
  }

  immigrant_times <- simulate_inhomogeneous_poisson(
    rate_function = immigrant_rate,
    max_rate = mu0 * attr(immigrant_multiplier, "maximum"),
    interval_start = obs_start - burn_in,
    interval_end = obs_end
  )
  if (length(immigrant_times) == 0) {
    return(data.frame(
      start = numeric(0), original_value = numeric(0),
      generation = integer(0), cluster_id = integer(0)
    ))
  }

  event_times <- immigrant_times
  event_marks <- draw_marks(length(immigrant_times))
  generations <- integer(length(immigrant_times))
  cluster_ids <- seq_along(immigrant_times)

  parent <- 1L
  while (parent <= length(event_times)) {
    offspring_mean <- branching_n
    if (marked_excitation) offspring_mean <- offspring_mean * event_marks[parent]
    n_offspring <- rpois(1, offspring_mean)

    if (n_offspring > 0) {
      child_times <- event_times[parent] + rexp(n_offspring, rate = beta)
      child_times <- child_times[child_times < obs_end]

      if (length(child_times) > 0) {
        event_times <- c(event_times, child_times)
        event_marks <- c(
          event_marks,
          draw_marks(length(child_times))
        )
        generations <- c(
          generations,
          rep.int(generations[parent] + 1L, length(child_times))
        )
        cluster_ids <- c(
          cluster_ids,
          rep.int(cluster_ids[parent], length(child_times))
        )
      }
    }
    parent <- parent + 1L
  }

  keep <- event_times >= obs_start & event_times < obs_end
  order_idx <- order(event_times[keep])
  data.frame(
    start = event_times[keep][order_idx],
    original_value = event_marks[keep][order_idx],
    generation = generations[keep][order_idx],
    cluster_id = cluster_ids[keep][order_idx]
  )
}

generate_hawkes_events <- function(
    n_participants, n_days, base_rate, branching_n = 0.6, beta = 0.15,
    obs_start = 480, obs_end = 1320,
    tod_amplitude = 0, tod_peak = 14 * 60, tod_sharpness = 1,
    marked_excitation = FALSE, mark_shape = 1,
    start_jitter = 0,
    person_amplitude = NULL, person_peak = NULL, person_rate = NULL) {

  all_events <- list()

  for (seqn in seq_len(n_participants)) {
    amp_s <- if (!is.null(person_amplitude)) person_amplitude[seqn] else tod_amplitude
    peak_s <- if (!is.null(person_peak)) person_peak[seqn] else tod_peak
    rate_s <- if (!is.null(person_rate)) person_rate[seqn] else base_rate
    for (day in seq_len(n_days)) {
      jit <- if (start_jitter > 0) runif(1, 0, start_jitter) else 0
      jit <- jit - start_jitter / 2
      day_events <- simulate_hawkes_cluster_day(
        base_rate = rate_s,
        branching_n = branching_n,
        beta = beta,
        obs_start = obs_start + jit,
        obs_end = obs_end   + jit,
        tod_amplitude = amp_s,
        tod_peak = peak_s,
        tod_sharpness = tod_sharpness,
        marked_excitation = marked_excitation,
        mark_shape = mark_shape
      )
      if (nrow(day_events) == 0) next

      day_events <- day_events %>%
        mutate(
          SEQN = seqn,
          WEEKDAY = day,
          categories = "active",
          mark_sqrt = original_value,
          daily_mean_sqrt = mean(mark_sqrt),
          mark_norm = mark_sqrt / daily_mean_sqrt
        ) %>%
        select(
          SEQN, WEEKDAY, start, original_value, categories,
          mark_sqrt, daily_mean_sqrt, mark_norm, generation, cluster_id
        )

      all_events[[paste(seqn, day, sep = "_")]] <- day_events
    }
  }

  bind_rows(all_events)
}

generate_clustered_events <- function(n_participants, n_days, base_rate,
                                      branching_n = 0.6,
                                      obs_start = 480, obs_end = 1320) {
  generate_hawkes_events(
    n_participants = n_participants,
    n_days = n_days,
    base_rate = base_rate,
    branching_n = branching_n,
    obs_start = obs_start,
    obs_end = obs_end,
    marked_excitation = FALSE,
    mark_shape = 1
  )
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

known_circadian_rates <- function(
    day_events, expected_rate_multiplier, bin_size) {
  n_bins <- 1440 / bin_size
  bin_rates <- vapply(seq_len(n_bins), function(bin_idx) {
    bin_start <- (bin_idx - 1) * bin_size
    minute_midpoints <- bin_start + seq(0.5, bin_size - 0.5, by = 1)
    mean(expected_rate_multiplier(minute_midpoints))
  }, numeric(1))

  covered_bins <- covered_bin_indices(
    min(day_events$start), max(day_events$start), bin_size, n_bins
  )
  bin_rates <- bin_rates / mean(bin_rates[covered_bins])
  attr(bin_rates, "bin_size") <- bin_size
  bin_rates
}

compute_dispersion_known_tod <- function(
    events_df, expected_rate_multiplier, window_sizes, bin_size) {
  day_results <- lapply(
    split(events_df, events_df$WEEKDAY, drop = TRUE),
    function(day_events) {
      circadian_rates <- known_circadian_rates(
        day_events, expected_rate_multiplier, bin_size
      )
      normalized_dispersion <- compute_dispersion_single_day(
        day_events = day_events,
        circadian_rates = circadian_rates,
        circadian_rates_marks = circadian_rates,
        window_sizes = window_sizes,
        bin_size = bin_size
      )
      normalized_dispersion
    }
  )

  day_results <- Filter(function(x) nrow(x) > 0, day_results)
  if (length(day_results) == 0) {
    return(data.frame(
      window_size = numeric(0),
      D_adj_known_tod = numeric(0),
      D_adj_marks_known_tod = numeric(0)
    ))
  }

  bind_rows(day_results) %>%
    group_by(window_size) %>%
    summarise(
      D_adj_known_tod = mean(D_adj, na.rm = TRUE),
      D_adj_marks_known_tod = mean(D_adj_marks, na.rm = TRUE),
      .groups = "drop"
    )
}

person_mark_factor <- function(person_events) {
  person_events %>%
    group_by(WEEKDAY) %>%
    summarise(
      Q_day = mean(mark_norm^2) / mean(mark_norm),
      .groups = "drop"
    ) %>%
    summarise(Q = mean(Q_day, na.rm = TRUE)) %>%
    pull(Q)
}

analyse_simulated_events <- function(
    events_df, window_sizes, bin_size,
    expected_rate_multiplier = NULL) {
  events_by_person <- split(events_df, events_df$SEQN, drop = TRUE)

  bind_rows(lapply(names(events_by_person), function(seqn_chr) {
    person_events <- events_by_person[[seqn_chr]]
    seqn <- person_events$SEQN[1]
    person_disp <- compute_dispersion(
      split(person_events, person_events$WEEKDAY, drop = TRUE),
      window_sizes = window_sizes,
      bin_size = bin_size
    )
    Q <- person_mark_factor(person_events)

    person_disp <- person_disp %>%
      mutate(
        SEQN = seqn,
        Q = Q,
        n_raw = branching_ratio(D_raw),
        n_adj = branching_ratio(D_adj),
        n_marks = 1 - sqrt(Q / D_marks),
        n_adj_marks = 1 - sqrt(Q / D_adj_marks)
      )

    if (!is.null(expected_rate_multiplier)) {
      known_tod_disp <- compute_dispersion_known_tod(
        person_events,
        expected_rate_multiplier = expected_rate_multiplier,
        window_sizes = window_sizes,
        bin_size = bin_size
      )
      person_disp <- person_disp %>%
        left_join(known_tod_disp, by = "window_size") %>%
        mutate(
          n_adj_known_tod = branching_ratio(D_adj_known_tod),
          n_adj_marks_known_tod = 1 - sqrt(Q / D_adj_marks_known_tod)
        )
    }

    person_disp
  }))
}

summarise_branching_estimates <- function(results_df, true_n) {
  estimate_columns <- intersect(
    c(
      "n_raw", "n_marks",
      "n_adj_known_tod", "n_adj_marks_known_tod"
    ),
    names(results_df)
  )

  results_df %>%
    group_by(window_size) %>%
    summarise(
      across(
        all_of(estimate_columns),
        list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) %>%
    mutate(true_n = true_n, .after = window_size)
}

#### 2 Run Simulations ####
simulation_seed <- 42
set_scenario_seed <- function(offset) {
  set.seed(simulation_seed + as.integer(offset))
}
source("function_dispersion.R")

# Parameters
n_participants <- 300
n_days <- 7
rate_per_hour <- 20
true_n <- 0.6
beta <- 0.15
# time-of-day baseline
tod_amplitude <- 0.67
tod_sharpness <- 1
tod_peak <- 14 * 60
tod_start_jitter <- 30
mark_variances <- c(0.25, 1, 2)
marked_tod_mark_variance <- 0.25
marked_tod_scenario <- paste0(
  "Marked Hawkes with TOD Baseline\n(Var(M) = ",
  format(marked_tod_mark_variance, trim = TRUE),
  ")"
)
window_sizes <- c(5, 10, 15, 30, 45, 60, 75, 90, 120)
bin_size <- 120

# Test 1: Poisson (uniform) process - expected D ≈ 1, n ≈ 0
cat("Test 1: Poisson Process \n")
cat("Expected: D ≈ 1, n ≈ 0\n")
set_scenario_seed(0)
poisson_events <- generate_poisson_events(n_participants, n_days, rate_per_hour)
cat("Generated", nrow(poisson_events), "events for", n_participants, "participants\n")
cat("Mean events per day:", mean(table(paste(poisson_events$SEQN, poisson_events$WEEKDAY))), "\n\n")
poisson_df <- analyse_simulated_events(
  poisson_events, window_sizes = window_sizes, bin_size = bin_size
)
poisson_summary <- summarise_branching_estimates(poisson_df, true_n = 0)

# Test 2: Clustered process - expected D > 1, n > 0
cat("\nTest 2: Stationary Unmarked Hawkes Process\n")
cat("Expected: n_raw approaches", true_n, "as the window size increases\n")
set_scenario_seed(0)
clustered_events <- generate_hawkes_events(
  n_participants = n_participants,
  n_days = n_days,
  base_rate = rate_per_hour,
  branching_n = true_n,
  beta = beta,
  marked_excitation = FALSE,
  mark_shape = 1e6,
  start_jitter = tod_start_jitter
)
cat("Generated", nrow(clustered_events), "events\n")
cat("Mean events per day:", mean(table(paste(clustered_events$SEQN, clustered_events$WEEKDAY))), "\n\n")

clustered_df <- analyse_simulated_events(
  clustered_events, window_sizes = window_sizes, bin_size = bin_size
)
clustered_summary <- summarise_branching_estimates(
  clustered_df, true_n = true_n
)

# Test 3: Regular (underdispersed) process - expected D < 1, n < 0
cat("\nTest 3: Regular Process (evenly spaced)\n")
cat("Expected: D < 1, n < 0 (regularity)\n")
set_scenario_seed(0)
regular_events <- generate_regular_events(n_participants, n_days, n_events_per_day = 280)
regular_df <- analyse_simulated_events(
  regular_events, window_sizes = window_sizes, bin_size = bin_size
)
regular_summary <- summarise_branching_estimates(regular_df, true_n = NA_real_)

#
# Sub-Test: different branching ratios across participants
#
cat("\nTest: Different Branching Ratios Across Participants\n")
cat("Expected: n_raw is biased low but tracks true n across participants\n")
set_scenario_seed(0)

person_n      <- runif(n_participants, 0.45, 0.85)
person_lambda <- pmax(8, rnorm(n_participants, rate_per_hour, 7.6))

diff_events <- bind_rows(lapply(seq_len(n_participants), function(i) {
  person_events <- generate_hawkes_events(
    n_participants    = 1,
    n_days            = n_days,
    base_rate         = person_lambda[i],
    branching_n       = person_n[i],
    beta              = beta,
    marked_excitation = FALSE,
    mark_shape        = 1e6
  )
  person_events$SEQN <- i
  person_events
}))

cat("Generated", nrow(diff_events), "events\n")
cat("Mean events per day:",
    mean(table(paste(diff_events$SEQN, diff_events$WEEKDAY))), "\n")

diff_df <- analyse_simulated_events(
  diff_events,
  window_sizes = window_sizes,
  bin_size     = bin_size
) %>%
  mutate(
    true_n      = person_n[SEQN],
    true_lambda = person_lambda[SEQN]
  )

# bias = level accuracy; 
# spearman = between-person ordering;
# cor_lambda should be ~0, i.e. the estimator is not just tracking activity rate.
diff_summary <- diff_df %>%
  group_by(window_size) %>%
  summarise(
    mean_true_n = mean(true_n),
    mean_n_raw  = mean(n_raw, na.rm = TRUE),
    bias        = mean(n_raw - true_n, na.rm = TRUE),
    spearman    = cor(n_raw, true_n, method = "spearman", use = "complete.obs"),
    slope       = coef(lm(n_raw ~ true_n))[2],
    cor_lambda  = cor(n_raw, true_lambda, use = "complete.obs"),
    .groups = "drop"
  )

# Test 4: stationary marked Hawkes process. 
cat("\nTest 4: Stationary Marked Hawkes Process\n")
cat(
  "Expected: n_marks approaches", true_n,
  "as Var(M) increases with E[M] fixed at 1\n"
)
marked_variance_df <- bind_rows(lapply(
  seq_along(mark_variances),
  function(index) {
    mark_variance <- mark_variances[index]
    set_scenario_seed(0)
    marked_events <- generate_hawkes_events(
      n_participants = n_participants,
      n_days = n_days,
      base_rate = rate_per_hour,
      branching_n = true_n,
      beta = beta,
      marked_excitation = TRUE,
      mark_shape = 1 / mark_variance,
      start_jitter = tod_start_jitter
    )
    analyse_simulated_events(
      marked_events,
      window_sizes = window_sizes,
      bin_size = bin_size
    ) %>%
      mutate(mark_variance = mark_variance)
  }
))

marked_variance_summary <- marked_variance_df %>%
  group_by(mark_variance, window_size) %>%
  summarise(
    true_n = true_n,
    n_marks_mean = mean(n_marks, na.rm = TRUE),
    n_marks_sd = sd(n_marks, na.rm = TRUE),
    .groups = "drop"
  )

# Test 5: time-varying baseline with unmarked excitation. Immigrants are generated by thinning mu(t), 
# while the offspring mechanism and true n stay unchanged.
cat("\nTest 5: Hawkes Process with Time-of-Day Baseline\n")
set_scenario_seed(0)
tod_expected_rate <- make_expected_hawkes_tod_multiplier(
  branching_n = true_n,
  beta = beta,
  amplitude = tod_amplitude,
  peak_time = tod_peak,
  sharpness = tod_sharpness
)
tod_events <- generate_hawkes_events(
  n_participants = n_participants,
  n_days = n_days,
  base_rate = rate_per_hour,
  branching_n = true_n,
  beta = beta,
  tod_amplitude = tod_amplitude,
  tod_peak = tod_peak,
  tod_sharpness = tod_sharpness,
  marked_excitation = FALSE,
  mark_shape = 1e6,
  start_jitter = tod_start_jitter
)
tod_df <- analyse_simulated_events(
  tod_events,
  window_sizes = window_sizes,
  bin_size = bin_size,
  expected_rate_multiplier = tod_expected_rate
)
tod_summary <- summarise_branching_estimates(tod_df, true_n = true_n)

# Test 6 combines both complications: positive marks modulate offspring and the immigrant rate follows the same time-of-day curve.
cat("\nTest 6: Marked Hawkes Process with Time-of-Day Baseline\n")
set_scenario_seed(0)
marked_tod_events <- generate_hawkes_events(
  n_participants = n_participants,
  n_days = n_days,
  base_rate = rate_per_hour,
  branching_n = true_n,
  beta = beta,
  tod_amplitude = tod_amplitude,
  tod_peak = tod_peak,
  tod_sharpness = tod_sharpness,
  marked_excitation = TRUE,
  mark_shape = 1 / marked_tod_mark_variance,
  start_jitter = tod_start_jitter
)
marked_tod_df <- analyse_simulated_events(
  marked_tod_events,
  window_sizes = window_sizes,
  bin_size = bin_size,
  expected_rate_multiplier = tod_expected_rate
)
marked_tod_summary <- summarise_branching_estimates(
  marked_tod_df, true_n = true_n
)

#### 3 Visualization ####

raw_dispersion_summary <- function(results_df, process_name) {
  results_df %>%
    group_by(window_size) %>%
    summarise(
      mean_D_raw = mean(D_raw, na.rm = TRUE),
      sd_D_raw = sd(D_raw, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(process = process_name)
}

all_results <- bind_rows(
  raw_dispersion_summary(poisson_df, "Poisson"),
  raw_dispersion_summary(clustered_df, "Clustered"),
  raw_dispersion_summary(regular_df, "Regular")
) %>%
  mutate(
    process = factor(
      process,
      levels = c(
        "Poisson",
        "Clustered",
        "Regular"
      )
    )
  )

simulation_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

# Plot D_raw for the Poisson, clustered, and regular process checks.
p1 <- ggplot(all_results, aes(x = window_size, y = mean_D_raw, color = process)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(
    ymin = mean_D_raw - sd_D_raw,
    ymax = mean_D_raw + sd_D_raw,
    fill = process
  ),
              alpha = 0.2, color = NA) +
  labs(
    title = "Dispersion Index Across Simulated Process Types",
    x = "Window size (minutes)",
    y = expression("Mean dispersion index, " * D),
    color = "Simulated process",
    fill = "Simulated process"
  ) +
  simulation_theme

print(p1)
ggsave("Output/dispersion/simulation.pdf", p1, width = 12, height = 8)

diff_plot_data <- diff_df %>% filter(window_size == 30)

diff_labels <- diff_summary %>%
  filter(window_size == 30) %>%
  mutate(
    label = sprintf("rho == %.2f ~~~ bias == '%+.3f'", spearman, bias)
  )

p2 <- ggplot(diff_plot_data, aes(x = true_n, y = n_raw)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(alpha = 0.45, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.9) +
  geom_text(
    data = diff_labels,
    aes(x = -Inf, y = Inf, label = label), parse = TRUE,
    inherit.aes = FALSE, hjust = -0.1, vjust = 1.6, size = 3.6
  ) +
  labs(
    title    = "Between-Participant Differences in the Branching Ratio",
    subtitle = expression("Window " * W * " = 30 min; dashed line is " * hat(n) * " = " * n),
    x        = expression("True branching ratio, " * n),
    y        = expression("Estimated branching ratio, " * hat(n))
  ) +
  simulation_theme

print(p2)
ggsave(
  "Output/dispersion/simulation_different_n.pdf",
  p2, width = 7, height = 6
)

validation_plot_data <- bind_rows(
  clustered_df %>%
    select(window_size, n_raw) %>%
    rename(`Event-count-based` = n_raw) %>%
    pivot_longer(-window_size, names_to = "estimator", values_to = "n_hat") %>%
    mutate(
      scenario = "Stationary Unmarked Hawkes",
      mark_variance = "Not varied"
    ),
  marked_variance_df %>%
    select(window_size, n_marks, mark_variance) %>%
    rename(`Intensity-mark-weighted` = n_marks) %>%
    pivot_longer(
      cols = `Intensity-mark-weighted`,
      names_to = "estimator",
      values_to = "n_hat"
    ) %>%
    mutate(
      scenario = "Stationary Marked Hawkes (E[M] = 1)",
      mark_variance = as.character(mark_variance)
    ),
  tod_df %>%
    select(
      window_size,
      n_raw,
      n_adj_known_tod
    ) %>%
    rename(
      `Event-count-based` = n_raw,
      `Event-count-based, TOD normalized` = n_adj_known_tod
    ) %>%
    pivot_longer(-window_size, names_to = "estimator", values_to = "n_hat") %>%
    mutate(
      scenario = "Unmarked Hawkes with TOD Baseline",
      mark_variance = "Not varied"
    ),
  marked_tod_df %>%
    select(
      window_size,
      n_marks,
      n_adj_marks_known_tod
    ) %>%
    rename(
      `Intensity-mark-weighted` = n_marks,
      `Intensity-mark-weighted, TOD normalized` = n_adj_marks_known_tod
    ) %>%
    pivot_longer(-window_size, names_to = "estimator", values_to = "n_hat") %>%
    mutate(
      scenario = marked_tod_scenario,
      mark_variance = "Not varied"
    )
) %>%
  group_by(scenario, estimator, mark_variance, window_size) %>%
  summarise(
    mean_n = mean(n_hat, na.rm = TRUE),
    sd_n = sd(n_hat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "Stationary Unmarked Hawkes",
        "Stationary Marked Hawkes (E[M] = 1)",
        "Unmarked Hawkes with TOD Baseline",
        marked_tod_scenario
      )
    ),
    estimator = factor(
      estimator,
      levels = c(
        "Event-count-based",
        "Event-count-based, TOD normalized",
        "Intensity-mark-weighted",
        "Intensity-mark-weighted, TOD normalized"
      )
    ),
    mark_variance = factor(
      mark_variance,
      levels = c("Not varied", "0.25", "1", "2")
    )
  )

# ===============================================================
# p3: branching-ratio recovery across scenarios
#
estimator_colors <- c(
  "Event-count-based"                       = "#E76F51",
  "Event-count-based, TOD normalized"       = "#7CAE00",
  "Intensity-mark-weighted"                 = "#00BFC4",
  "Intensity-mark-weighted, TOD normalized" = "#C77CFF"
)
marked_variance_plot_data <- validation_plot_data %>%
  filter(scenario == "Stationary Marked Hawkes (E[M] = 1)")

non_variance_plot_data <- validation_plot_data %>%
  filter(scenario != "Stationary Marked Hawkes (E[M] = 1)")

reference_label_unmarked <- "Stationary reference, unmarked"
reference_label_marked   <- paste0(
  "Stationary reference, marked (Var(M) = ",
  format(marked_tod_mark_variance, trim = TRUE), ")"
)

stationary_unmarked_ref <- validation_plot_data %>%
  filter(scenario == "Stationary Unmarked Hawkes") %>%
  transmute(
    window_size,
    ref_n    = mean_n,
    scenario = factor(
      "Unmarked Hawkes with TOD Baseline",
      levels = levels(validation_plot_data$scenario)
    )
  )

stationary_marked_ref <- validation_plot_data %>%
  filter(
    scenario == "Stationary Marked Hawkes (E[M] = 1)",
    as.character(mark_variance) == as.character(marked_tod_mark_variance)
  ) %>%
  transmute(
    window_size,
    ref_n    = mean_n,
    scenario = factor(
      marked_tod_scenario,
      levels = levels(validation_plot_data$scenario)
    )
  )


estimator_colors <- c(
  estimator_colors,
  setNames("gray35", reference_label_unmarked),
  setNames("gray62", reference_label_marked)
)

# ---- clipped view ---------------------------------------------------
y_clip_min    <- 0.40
y_visible_max <- max(
  validation_plot_data$mean_n + validation_plot_data$sd_n,
  na.rm = TRUE
)
y_span <- y_visible_max - y_clip_min

variance_legend_data <- data.frame(
  scenario = factor(
    "Stationary Marked Hawkes (E[M] = 1)",
    levels = levels(validation_plot_data$scenario)
  ),
  mark_variance = factor(
    c("0.25", "1", "2"),
    levels = levels(validation_plot_data$mark_variance)
  ),
  label   = c("0.25", "1", "2"),
  x       = 84,   # segment start
  xend    = 92,   # segment end
  x_point = 88,   # marker, mid-segment
  x_label = 95,   # value label
  y = y_clip_min + y_span * c(0.20, 0.13, 0.06)
)

variance_legend_title <- data.frame(
  scenario = factor(
    "Stationary Marked Hawkes (E[M] = 1)",
    levels = levels(validation_plot_data$scenario)
  ),
  x     = 84,
  y     = y_clip_min + y_span * 0.29,
  label = "Var(M)"
)

variance_legend_layers <- unlist(
  lapply(names(variance_colors), function(variance_value) {
    legend_row <- variance_legend_data %>%
      filter(as.character(mark_variance) == variance_value)
    variance_color <- variance_colors[[variance_value]]

    list(
      geom_segment(
        data = legend_row,
        aes(x = x, xend = xend, y = y, yend = y, linetype = mark_variance),
        inherit.aes = FALSE, color = variance_color, linewidth = 0.8
      ),
      geom_point(
        data = legend_row,
        aes(x = x_point, y = y, shape = mark_variance),
        inherit.aes = FALSE, color = variance_color, size = 2.0
      )
    )
  }),
  recursive = FALSE
)

stopifnot(all(variance_legend_data$y > y_clip_min))  # legend inside view

p3 <- ggplot(
  validation_plot_data,
  aes(
    x = window_size, y = mean_n,
    color = estimator, fill = estimator,
    linetype = mark_variance, shape = mark_variance,
    group = interaction(estimator, mark_variance)
  )
) +
  geom_hline(yintercept = true_n, linetype = "dashed", color = "gray40") +

  geom_line(
    data = stationary_unmarked_ref,
    aes(x = window_size, y = ref_n, color = reference_label_unmarked),
    inherit.aes = FALSE,
    linewidth = 0.7, linetype = "longdash"
  ) +
  geom_line(
    data = stationary_marked_ref,
    aes(x = window_size, y = ref_n, color = reference_label_marked),
    inherit.aes = FALSE,
    linewidth = 0.7, linetype = "dotdash"
  ) +

  geom_ribbon(
    data = non_variance_plot_data,
    aes(ymin = mean_n - sd_n, ymax = mean_n + sd_n),
    alpha = 0.15, color = NA
  ) +
  geom_line(data = non_variance_plot_data, linewidth = 1) +
  geom_point(data = non_variance_plot_data, size = 2) +
  marked_variance_layers +

  geom_text(
    data = variance_legend_title,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, hjust = 0, fontface = "bold", size = 3.0
  ) +
  variance_legend_layers +
  geom_text(
    data = variance_legend_data,
    aes(x = x_label, y = y, label = label),
    inherit.aes = FALSE, hjust = 0, size = 2.7
  ) +

  facet_wrap(~scenario, ncol = 2) +
  coord_cartesian(ylim = c(y_clip_min, NA)) +

  scale_linetype_manual(
    values = c("Not varied" = "solid", "0.25" = "solid",
               "1" = "dashed", "2" = "dotdash")
  ) +
  scale_shape_manual(
    values = c("Not varied" = 16, "0.25" = 16, "1" = 17, "2" = 15)
  ) +
  scale_color_manual(values = estimator_colors,
    breaks = c(
      setdiff(names(estimator_colors),
              c(reference_label_unmarked, reference_label_marked)),
      reference_label_unmarked,
      reference_label_marked
    )
  ) +
  scale_fill_manual(values = estimator_colors) +
  labs(
    title = "Branching Ratio Recovery Across Hawkes Simulations",
    x     = "Window size (minutes)",
    y     = expression("Mean estimated branching ratio, " * hat(n)),
    color = "Series"
  ) +
  guides(
    color    = guide_legend(nrow = 3, byrow = TRUE, order = 1),
    fill     = "none",   # colour has one more level; keeps one legend
    linetype = "none",   # drawn manually inside the panel
    shape    = "none"    # drawn manually inside the panel
  ) +
  simulation_theme

print(p3)
ggsave("Output/dispersion/simulation_marked_tod.pdf", p3, width = 10, height = 6)

print(all_results%>%filter(process == "Clustered"), n = Inf)
'''
# A tibble: 9 × 4
  window_size mean_D_raw sd_D_raw process  
        <dbl>      <dbl>    <dbl> <fct>    
1           5       1.67    0.103 Clustered
2          10       2.23    0.191 Clustered
3          15       2.69    0.265 Clustered
4          30       3.64    0.468 Clustered
5          45       4.25    0.621 Clustered
6          60       4.59    0.780 Clustered
7          75       4.90    0.934 Clustered
8          90       5.08    1.02  Clustered
9         120       5.26    1.32  Clustered
'''