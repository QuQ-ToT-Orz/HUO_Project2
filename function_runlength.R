#### 1  Generate RLE data with quantile binning ####
#' Runlength encoding with optimal structure for Hawkes modeling
#' @param n_bins No defined cutoffs
#' @param Act_Analysis Activity analysis data frame
#' 
runlength_single <- function(Act_Analysis, n_bins = 10, bin_method = "quantile", global_bins = FALSE) {
  act_mat <- as.matrix(Act_Analysis[, paste0("MIN", 1:1440)])
  act_mat[is.na(act_mat)] <- 0
  Peak30 <- Act_Analysis$Peak30

  # Cap all activity values at their respective Peak30
  act_mat_capped <- pmin(act_mat, Peak30)

  # Calculate global bin edges if requested
  if (global_bins && bin_method == "quantile") {
    all_activity_values <- as.vector(act_mat_capped)[as.vector(act_mat_capped) > 0]
    quantile_probs <- seq(0, 1, length.out = n_bins + 1)
    global_bin_edges <- quantile(all_activity_values, probs = quantile_probs, na.rm = TRUE)
    global_bin_edges <- unique(global_bin_edges)
    cat("Using global quantile bins. Number of unique bins:", length(global_bin_edges) - 1, "\n")
  }

  bout_mat <- map(seq_len(nrow(act_mat)), function(i) {
    activity_values <- act_mat_capped[i, ]
    current_peak30 <- Peak30[i]

    # Create bins using quantile or equal-width method
    if (bin_method == "quantile") {
      if (global_bins) {
        bin_edges <- global_bin_edges
      } else {
        non_zero_values <- activity_values[activity_values > 0]
        quantile_probs <- seq(0, 1, length.out = n_bins + 1)
        bin_edges <- quantile(non_zero_values, probs = quantile_probs, na.rm = TRUE)
      }
    } else {
      # Equal-width bins
      bin_edges <- seq(1, current_peak30, length.out = n_bins + 1)
    }
    bin_edges <- unique(sort(bin_edges))

    # sed_act_list <- map(seq(from = 1, to = floor(Peak30[i]), by = 5), function(j) {
    map_dfr(seq_along(bin_edges), function(bin_idx) {
      j <- bin_edges[bin_idx]
      x <- activity_values >= j
      mat <- rlenc(x)

      start_positions <- cumsum(c(1, mat$lengths[-length(mat$lengths)]))
      original_values <- activity_values[start_positions]

      tibble(
        SEQN = Act_Analysis$SEQN[i],
        WEEKDAY = Act_Analysis$WEEKDAY[i],
        run_id = bin_idx,
        lengths = mat$lengths,
        values = mat$values,
        start = start_positions,
        grid_value = j,
        original_value = original_values,
        end = start_positions + mat$lengths - 1,
        bin_edges = list(bin_edges)
      )
    })
  })

  bind_rows(bout_mat) %>%
    mutate(
      SEQN = as.factor(SEQN),
      WEEKDAY = as.factor(WEEKDAY),
      run_id = factor(run_id, levels = rev(levels(as.factor(run_id))))
    )
}

runlength <- function(Act_Analysis, n_bins = 10, bin_method = "quantile", global_bins = FALSE) {
  act_mat <- as.matrix(Act_Analysis[, paste0("MIN", 1:1440)])
  Peak30 <- Act_Analysis$Peak30

  # Cap all activity values at their respective Peak30
  act_mat_capped <- pmin(act_mat, Peak30)

  # Calculate global bin edges if requested
  if (global_bins && bin_method == "quantile") {
    all_activity_values <- as.vector(act_mat_capped)[as.vector(act_mat_capped) > 0]
    quantile_probs <- seq(0, 1, length.out = n_bins + 1)
    global_bin_edges <- quantile(all_activity_values, probs = quantile_probs, na.rm = TRUE)
    global_bin_edges <- unique(global_bin_edges)
    cat("Using global quantile bins. Number of unique bins:", length(global_bin_edges) - 1, "\n")
  }

  map_dfr(seq_len(nrow(act_mat)), function(i) {
    activity_values <- act_mat_capped[i, ]
    current_peak30 <- Peak30[i]

    # Create bins using quantile or equal-width method
    if (bin_method == "quantile") {
      if (global_bins) {
        bin_edges <- global_bin_edges
      } else {
        non_zero_values <- activity_values[activity_values > 0]
        quantile_probs <- seq(0, 1, length.out = n_bins + 1)
        bin_edges <- quantile(non_zero_values, probs = quantile_probs, na.rm = TRUE)
      }
    } else {
      # Equal-width bins
      bin_edges <- seq(1, current_peak30, length.out = n_bins + 1)
    }

    # Remove duplicate edges
    bin_edges <- unique(sort(bin_edges))

    # Assign each minute to a bin (grid value)
    # bin_centers <- (bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2
    bin_centers <- bin_edges
    bin_assignments <- cut(activity_values, breaks = bin_edges, labels = FALSE, include.lowest = TRUE)
    # Convert to grid values
    grid_value <- bin_centers[bin_assignments]
    # Perform RLE on grid values
    grid_rle <- rlenc(grid_value) # rle

    # Map grid values back to bin numbers
    bin_numbers <- sapply(grid_rle$values, function(val) {
      which(bin_centers == val)[1]
    })

    start_positions <- cumsum(c(1, grid_rle$lengths[-length(grid_rle$lengths)]))
    original_values <- activity_values[start_positions]

    tibble(
      SEQN = Act_Analysis$SEQN[i],
      WEEKDAY = Act_Analysis$WEEKDAY[i],
      run_id = bin_numbers, # bin number (1 to n_bins - 1)
      lengths = grid_rle$lengths,
      grid_value = grid_rle$values, # bin edge values
      original_value = original_values,
      start = start_positions,
      end = start_positions + grid_rle$lengths - 1,
      bin_edges = list(bin_edges)
    )
  })
}

#### 2 Generate RLE data with manual binning ####
#' Runlength encoding with optimal structure for Hawkes modeling
#'
#' @param Act_Analysis Activity analysis data frame
#' @param activity_cutoffs Device type ("hip" or "wrist")
#' @param force_optimal_structure Use 3-level optimal structure (default: TRUE)
runlength_manual <- function(Act_Analysis,
                          activity_cutoffs = c("hip", "wrist"),
                          force_optimal_structure = TRUE) {
  # Set cutoffs based on device type
  if (activity_cutoffs == "hip") {
    cutoff_values <- c(100, 2020)
  } else if (activity_cutoffs == "wrist") {
    cutoff_values <- c(10.558, 19.614)
  } else {
    stop("activity_cutoffs must be 'hip' or 'wrist'")
  }

  # Calculate dataset-level rest threshold
  act_mat <- as.matrix(Act_Analysis[, paste0("MIN", 1:1440)])
  act_mat[is.na(act_mat)] <- 0
  Peak30 <- Act_Analysis$Peak30

  all_activity <- as.vector(act_mat)
  positive_activity <- all_activity[all_activity > 0]
  rest_threshold <- quantile(positive_activity, probs = 0.1, na.rm = TRUE)

  sedentary_threshold <- cutoff_values[1]
  moderate_threshold <- cutoff_values[2]

  # Choose structure based on optimization flag
  if (force_optimal_structure) {
    # OPTIMAL: Use 3-level structure for both hip and wrist
    # sed_q1 <- rest_threshold + (sedentary_threshold - rest_threshold) * 0.3
    # light_threshold <- sedentary_threshold + (moderate_threshold - sedentary_threshold) * 0.3
    structure_used <- paste0("3level_", activity_cutoffs, "_optimized")
  } else {
    # DEFAULT: Use original 5-level for both (backward compatibility)
    light_threshold <- mean(c(sedentary_threshold, moderate_threshold))
    sed_q1 <- (rest_threshold + sedentary_threshold) * 0.5
    structure_used <- "5level_original"
  }

  # Process data with appropriate structure
  map_dfr(seq_len(nrow(act_mat)), function(i) {
    activity <- as.numeric(act_mat[i, ])
    peak30_value <- Peak30[i]

    # Cap activity values at peak30 for high intensity classification
    activity <- pmin(activity, peak30_value)

    if (force_optimal_structure) {
      # OPTIMAL: 3-level classification for both hip and wrist
      levels <- case_when(
        activity < rest_threshold ~ NA_real_,
        activity < sedentary_threshold ~ 1, # Sedentary
        activity < moderate_threshold ~ 2, # Active moderate
        activity >= moderate_threshold ~ 3, # Active high
        TRUE ~ NA_real_
      )

      # Categories for 3-level (balanced for Hawkes)
      category_mapping <- c("sedentary", "light_moderate", "high")
      grid_mapping <- c(sedentary_threshold, moderate_threshold, peak30_value)
      bin_edges <- c(rest_threshold, sedentary_threshold, moderate_threshold, peak30_value)
    } else {
      # DEFAULT: 5-level classification (original)
      levels <- case_when(
        activity < rest_threshold ~ NA_real_,
        activity < sed_q1 ~ 1,
        activity < sedentary_threshold ~ 2,
        activity < light_threshold ~ 3,
        activity < moderate_threshold ~ 4,
        activity >= moderate_threshold ~ 5,
        TRUE ~ NA_real_
      )

      # Categories for 5-level (original)
      category_mapping <- c("sedentary", "sedentary", "light_moderate", "light_moderate", "high")
      grid_mapping <- c(sed_q1, sedentary_threshold, light_threshold, moderate_threshold, peak30_value)
      bin_edges <- c(rest_threshold, sed_q1, sedentary_threshold, light_threshold, moderate_threshold, peak30_value)
    }

    # RLE and filter
    rle_result <- rle(levels)
    # long_bouts <- which(rle_result$lengths >= min_bout_length)
    long_bouts <- which(!is.na(rle_result$values))

    starts <- cumsum(c(1, rle_result$lengths[-length(rle_result$lengths)]))[long_bouts]

    tibble(
      SEQN = Act_Analysis$SEQN[i],
      WEEKDAY = Act_Analysis$WEEKDAY[i],
      start = starts,
      lengths = rle_result$lengths[long_bouts],
      run_id = rle_result$values[long_bouts],
      grid_value = grid_mapping[rle_result$values[long_bouts]],
      end = starts + rle_result$lengths[long_bouts] - 1,
      categories = category_mapping[rle_result$values[long_bouts]],
      original_value = activity[starts],
      bin_edges = list(bin_edges)
    )
  })
}

runlength_single_manual <- function(Act_Analysis,
                                 activity_cutoffs = c("hip", "wrist")) {
  # Set cutoffs based on device type
  if (activity_cutoffs == "hip") {
    cutoff_values <- c(100, 2020)
  } else if (activity_cutoffs == "wrist") {
    cutoff_values <- c(10.558, 19.614)
  } else {
    stop("activity_cutoffs must be 'hip' or 'wrist'")
  }

  # Use only sedentary threshold for binary transitions
  sedentary_threshold <- cutoff_values[1]

  act_mat <- as.matrix(Act_Analysis[, paste0("MIN", 1:1440)])
  act_mat[is.na(act_mat)] <- 0
  Peak30 <- Act_Analysis$Peak30

  map_dfr(seq_len(nrow(act_mat)), function(i) {
    activity <- as.numeric(act_mat[i, ])
    current_peak30 <- Peak30[i]
    # Cap activity values at peak30 for high intensity classification
    activity <- pmin(activity, current_peak30)

    # Filter to only active minutes (>0) first
    active_minutes <- which(activity > 0)
    # Get activity values for active minutes only
    active_values <- activity[active_minutes]

    # Binary classification on active minutes only: TRUE = active (>=threshold), FALSE = sedentary (<threshold)
    x <- active_values >= sedentary_threshold
    mat <- rlenc(x)

    # Map back to original minute positions
    start_positions <- active_minutes[cumsum(c(1, mat$lengths[-length(mat$lengths)]))]
    original_values <- active_values[cumsum(c(1, mat$lengths[-length(mat$lengths)]))]

    tibble(
      SEQN = Act_Analysis$SEQN[i],
      WEEKDAY = Act_Analysis$WEEKDAY[i],
      run_id = 1, # Only one threshold
      lengths = mat$lengths,
      values = mat$values, # TRUE = sedentary->active transition, FALSE = active->sedentary transition
      start = start_positions,
      grid_value = sedentary_threshold,
      original_value = original_values,
      end = start_positions + mat$lengths - 1,
      categories = ifelse(mat$values, "active", "sedentary") # TRUE = active bout, FALSE = sedentary bout
    )
  })
}

#### 3 Simple minute-level events ####
simple_events <- function(Act_Analysis, activity_cutoffs = c("hip", "wrist")) {
  # Set cutoffs based on device type
  if (activity_cutoffs == "hip") {
    cutoff_values <- c(100, 2020)
  } else if (activity_cutoffs == "wrist") {
    cutoff_values <- c(10.558, 19.614)
  } else {
    stop("activity_cutoffs must be 'hip' or 'wrist'")
  }

  sedentary_threshold <- cutoff_values[1]
  moderate_threshold <- cutoff_values[2]

  map_dfr(seq_len(nrow(Act_Analysis)), function(i) {
    activity <- as.numeric(Act_Analysis[i, paste0(
      "MIN",
      1:1440
    )])
    peak30_value <- Act_Analysis$Peak30[i]

    # Create categories: sedentary, light_moderate, high
    categories <- case_when(
      activity > 0 & activity < sedentary_threshold ~ "sedentary",
      activity >= sedentary_threshold & activity < moderate_threshold ~ "light_moderate",
      activity >= moderate_threshold ~ "high"
    )

    # For original_value: use peak30 for high activity values
    original_value <- ifelse(activity >= peak30_value,
      peak30_value, activity
    )

    # Filter to only include minutes with activity > 0
    valid_minutes <- which(activity > 0)

    # Create minute-level events
    tibble(
      SEQN = Act_Analysis$SEQN[i],
      WEEKDAY = Act_Analysis$WEEKDAY[i],
      start = valid_minutes,
      categories = categories[valid_minutes],
      original_value = original_value[valid_minutes]
    )
  })
}

#### Other helper functions ####
# Fast post-processing function to filter short bouts
filter_short_bouts <- function(events_data, min_bout_length = 2) {
  events_data %>%
    group_by(SEQN, WEEKDAY) %>%
    arrange(start) %>%
    mutate(
      # Detect category changes OR time gaps to identify bout boundaries
      category_change = categories != lag(categories,
        default =
          "different"
      ) |
        (start - lag(start, default = 0)) > 1,
      bout_id = cumsum(category_change)
    ) %>%
    group_by(SEQN, WEEKDAY, bout_id) %>%
    # Keep sedentary regardless of length, only filter active bouts with < min_bout_length minutes
    filter(n() >= min_bout_length) %>%
    ungroup() %>%
    # Clean up temporary columns
    select(-category_change, -bout_id) %>%
    arrange(SEQN, WEEKDAY, start)
}

# Remove outlier minutes using forecast::tsoutliers
clean_events_spikes <- function(events_data, activity_cutoffs = c("hip", "wrist")) {
  activity_cutoffs <- match.arg(activity_cutoffs)

  cutoff_values <- if (activity_cutoffs == "hip") {
    c(100, 2020)
  } else if (activity_cutoffs == "wrist") {
    c(10.558, 19.614)
  } else {
    stop("activity_cutoffs must be 'hip' or 'wrist'")
  }

  detect_outliers <- function(values) {
    n <- length(values)
    freq <- floor(n / 2)
    ts_values <- ts(values, frequency = freq)

    tryCatch({
      fit <- forecast::tsoutliers(ts_values)
      flags <- rep(FALSE, n)
      if (!is.null(fit$index) && length(fit$index) > 0) {
        idx <- fit$index[!is.na(fit$index)]
        idx <- idx[idx >= 1 & idx <= n]
        flags[idx] <- TRUE
      }
      flags
    }, error = function(e) rep(FALSE, n))
  }

  events_data %>%
    group_by(SEQN, WEEKDAY) %>%
    arrange(start) %>%
    mutate(
      outlier_flag = detect_outliers(original_value)
    ) %>%
    filter(!outlier_flag) %>%
    select(-outlier_flag) %>%
    ungroup()
}

combined_df <- function(rle_df, seqn = NULL, weekday = NULL) {
  # Apply filters if seqn or weekday are provided
  filtered_df <- rle_df
  if (!is.null(seqn)) {
    filtered_df <- filtered_df %>% filter(SEQN == seqn)
  }
  if (!is.null(weekday)) {
    filtered_df <- filtered_df %>% filter(WEEKDAY == weekday)
  }

  # Process the data
  filtered_df %>%
    filter(values == TRUE) %>%
    group_by(SEQN, WEEKDAY, start) %>%
    filter(run_id == max(as.character(run_id))) %>%
    ungroup() %>%
    arrange(SEQN, WEEKDAY, start) %>%
    mutate(run_id = as.numeric(as.character(run_id)))
}


# Intensity Categorization
categorize_activity_bouts <- function(
    rle_data,
    activity_cutoffs = c("hip", "wrist")) {
  if (activity_cutoffs == "hip") {
    activity_cutoff <- c(100, 2020)
  } else if (activity_cutoffs == "wrist") {
    activity_cutoff <- c(10.558, 19.614)
  }

  rle_data %>%
    mutate(
      categories = case_when(
        is.na(original_value) ~ NA,
        original_value < activity_cutoff[1] ~ "sedentary",
        original_value <= activity_cutoff[2] ~ "light_moderate",
        TRUE ~ "high"
      )
    )
}