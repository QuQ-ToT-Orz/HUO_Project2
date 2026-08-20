#### Shared Dispersion Index Functions ####

sum_by_group <- function(values, group_idx, n_groups) {
  totals <- numeric(n_groups)
  if (length(values) == 0) return(totals)

  grouped_sums <- tapply(values, group_idx, sum)
  totals[as.integer(names(grouped_sums))] <- as.numeric(grouped_sums)
  totals
}

covered_bin_indices <- function(segment_start, segment_end, bin_size, n_bins) {
  first_bin <- floor((segment_start - 1) / bin_size) + 1
  last_bin <- floor((segment_end - 1) / bin_size) + 1

  first_bin <- pmin(pmax(first_bin, 1), n_bins)
  last_bin <- pmin(pmax(last_bin, 1), n_bins)

  first_bin:last_bin
}

estimate_circadian_baseline <- function(day_events, bin_size,
                                        segment_start = min(day_events$start),
                                        segment_end = max(day_events$start)) {
  n_bins <- 1440 / bin_size
  bin_idx <- pmin(floor((day_events$start - 1) / bin_size) + 1, n_bins)
  bin_rates <- tabulate(bin_idx, nbins = n_bins)

  covered_bins <- covered_bin_indices(
    segment_start, segment_end, bin_size, n_bins
  )
  bin_rates <- bin_rates / mean(bin_rates[covered_bins])

  attr(bin_rates, "bin_size") <- bin_size
  bin_rates
}

estimate_circadian_baseline_marks <- function(
    day_events, bin_size,
    segment_start = min(day_events$start),
    segment_end = max(day_events$start)) {
  n_bins <- 1440 / bin_size
  bin_idx <- pmin(floor((day_events$start - 1) / bin_size) + 1, n_bins)
  bin_rates <- sum_by_group(day_events$mark_norm, bin_idx, n_bins)

  covered_bins <- covered_bin_indices(
    segment_start, segment_end, bin_size, n_bins
  )
  bin_rates <- bin_rates / mean(bin_rates[covered_bins])

  attr(bin_rates, "bin_size") <- bin_size
  bin_rates
}

compute_dispersion_single_day <- function(
    day_events, circadian_rates, circadian_rates_marks,
    window_sizes, bin_size) {
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
    mu_vals <- tapply(
      circadian_rates[bin_idx_all[valid_events]], window_idx, mean
    )
    mu_window[as.integer(names(mu_vals))] <- as.numeric(mu_vals)

    adj_count <- ifelse(
      !is.na(mu_window) & mu_window > 0,
      raw_count / mu_window,
      0
    )
    adj_marks <- ifelse(
      !is.na(mu_window) & mu_window > 0,
      raw_marks / mu_window,
      0
    )

    D_raw <- var(raw_count) / mean(raw_count)
    D_marks <- var(raw_marks) / mean(raw_marks)
    D_adj <- var(adj_count) / mean(adj_count)
    D_adj_marks <- var(adj_marks) / mean(adj_marks)

    results[[as.character(W)]] <- data.frame(
      window_size = W,
      D_raw = D_raw,
      D_marks = D_marks,
      D_adj = D_adj,
      D_adj_marks = D_adj_marks
    )
  }
  dplyr::bind_rows(results)
}

compute_dispersion <- function(day_events_list, window_sizes, bin_size) {
  day_results <- lapply(day_events_list, function(day_events) {
    segment_start <- min(day_events$start)
    segment_end <- max(day_events$start)
    circadian_rates <- estimate_circadian_baseline(
      day_events, bin_size, segment_start, segment_end
    )
    circadian_rates_marks <- estimate_circadian_baseline_marks(
      day_events, bin_size, segment_start, segment_end
    )
    compute_dispersion_single_day(
      day_events, circadian_rates, circadian_rates_marks,
      window_sizes, bin_size
    )
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

  combined_results <- dplyr::bind_rows(day_results)
  grouped_results <- dplyr::group_by(combined_results, window_size)
  dplyr::summarise(
    grouped_results,
    D_raw = mean(D_raw, na.rm = TRUE),
    D_marks = mean(D_marks, na.rm = TRUE),
    D_adj = mean(D_adj, na.rm = TRUE),
    D_adj_marks = mean(D_adj_marks, na.rm = TRUE),
    .groups = "drop"
  )
}

branching_ratio <- function(D) {
  1 - 1 / sqrt(D)
}
