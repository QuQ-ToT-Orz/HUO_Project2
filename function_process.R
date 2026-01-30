# Grid search helper function for Hawkes processes
perform_hawkes_grid_search <- function(times_data = NULL, marks_data = NULL, events_data = NULL,
                                      times_list = NULL, marks_list = NULL, events_list = NULL,
                                      a_par_values = c(-2, -1, 0, 1, 2),
                                      mu_values = c(1, 0.5),
                                      beta_values = c(0.1, 1),
                                      model = 1,
                                      fit_function = c("fit_hawkes", "fit_hawkes_multi_series")) {
  
  fit_function <- match.arg(fit_function)
  best_fit <- NULL
  best_objective <- Inf
  best_params <- list(a_par = NA, mu_mult = NA, beta_mult = NA)
  
  # Determine if this is multi-series or single-series
  is_multi_series <- !is.null(times_list)
  is_multivariate <- grepl("mhawkes", fit_function)
  
  # Loop through all parameter combinations
  for (i in seq_along(a_par_values)) {
    for (j in seq_along(mu_values)) {
      for (k in seq_along(beta_values)) {
        try({
          if (fit_function == "fit_hawkes") {
            # Single-series marked Hawkes
            mu_val <- mu_values[j] * mean(length(times_data) / max(times_data))
            beta_val <- beta_values[k] / median(diff(times_data))
            
            fit_candidate <- fit_hawkes(
              times = times_data,
              marks = marks_data,
              parameters = list(
                mu = mu_val,
                alpha = a_par_values[i],
                beta = beta_val
              ),
              model = model
            )
            
          } else if (fit_function == "fit_hawkes_multi_series") {
            # Multi-series marked Hawkes
            mu_val <- mu_values[j] * mean(sapply(times_list, function(x) length(x) / max(x)))
            beta_val <- beta_values[k] / median(unlist(lapply(times_list, function(x) diff(x))))
            
            fit_candidate <- fit_hawkes_multi_series(
              times_list = times_list,
              marks_list = marks_list,
              model = model,
              parameters = list(
                mu = mu_val,
                alpha = a_par_values[i],
                beta = beta_val
              )
            )
            
          }
          
          # Check convergence and update best fit
          if (!is.null(fit_candidate) && TMB::sdreport(fit_candidate)$pdHess && 
              fit_candidate$objective < best_objective) {
            best_objective <- fit_candidate$objective
            best_fit <- fit_candidate
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
  
  # Add best parameters to fit object
  if (!is.null(best_fit)) {
    best_fit$best_a_par_init <- best_params$a_par
    best_fit$best_mu_mult_init <- best_params$mu_mult
    best_fit$best_beta_mult_init <- best_params$beta_mult
  }
  
  return(best_fit)
}

transform_hawkes_data <- function(data, seqn, weekday, time_divisor = 1,
                                 jitter_factor = 0, active_only = TRUE, data_type) {
  # Filter data for specific seqn and weekday
  filtered_data <- data %>% filter(SEQN == seqn, WEEKDAY == weekday)

    # Filter for active events only if specified
  if (active_only) {
    filtered_data <- filtered_data %>%
      filter(categories != "sedentary") %>%
      mutate(start_other = start_normalized - min(start_normalized))
  }

  times <- filtered_data$start_other
  if (jitter_factor > 0 && length(times) > 1) {
    # Keep first event at 0, jitter only subsequent events
    jittered_times <- jitter(times, factor = jitter_factor)
    times[2:length(times)] <- jittered_times[2:length(times)]
    times[1] <- 0 # Ensure first event stays at 0
  }

  # Rescale time so typical inter-event time is O(1)
  times <- times / time_divisor
  marks <- filtered_data$original_value

  # Apply sqrt transformation for hip data
  if (data_type == "hip") {
    marks <- sqrt(marks)
    transformed_marks <- marks / mean(marks)
  } else if (data_type == "wrist") {
    transformed_marks <- marks / mean(marks)
  }

  return(list(
    times = times,
    marks = transformed_marks
  ))
}


check_conv <- function(fit) {
  sd_report <- TMB::sdreport(fit)
  coef_table <- get_coefs(fit)
  params <- coef_table[, 1]
  std_errors <- coef_table[, 2]

  # Check if pdHess is TRUE
  pdHess_ok <- sd_report$pdHess

  # Check standard errors of untransformed parameters (log_mu, logit_ratio, a_par)
  if (pdHess_ok) {
    # Get standard errors from sd_report for untransformed parameters
    untransformed_se <- summary(sd_report, "fixed")[, 2]

    # Check if any standard errors are too large (> 10)
    if (any(untransformed_se > 10, na.rm = TRUE)) {
      pdHess_ok <- FALSE
    }
  }
  # Overall convergence: check pdHess and standard error criteria
  converged <- pdHess_ok

  list(
    parameters = params,
    std_errors = std_errors,
    converged = converged
  )
}

# marked Hawkes process
process_run_marked <- function(
    seqn, weekday = NULL, runs_df, penalty_coef = 0,
    single_day = FALSE, data_type) {
  if (single_day && !is.null(weekday)) {
    # Single day fitting
    transformed_data <- transform_hawkes_data(runs_df, seqn, weekday, active_only = TRUE, data_type = data_type)

    # Initialization search
    best_fit <- perform_hawkes_grid_search(
      times_data = transformed_data$times,
      marks_data = transformed_data$marks,
      fit_function = "fit_hawkes"
    )

    if (is.null(best_fit)) {
      cat(sprintf("Warning: No converged fits found for SEQN %s, weekday %s. Returning NULL result.\n", seqn, weekday))
      return(list(
        id = as.character(seqn),
        weekday = as.character(weekday),
        parameters = NULL,
        std_errors = NULL,
        pdHess = FALSE,
        penalty_coef = penalty_coef,
        best_a_par_init = NA
      ))
    }

    # best_a_par_init is already set by the helper function
    results <- check_conv(best_fit)

    # Get KS statistics
    gof_results <- show_hawkes_GOF(best_fit, plot = FALSE, tests = TRUE, return_values = TRUE)
    ks_D <- gof_results$ks_D

    return(list(
      id = as.character(seqn),
      weekday = as.character(weekday),
      parameters = results$parameters,
      std_errors = results$std_errors,
      pdHess = results$converged,
      penalty_coef = penalty_coef,
      best_a_par_init = best_fit$best_a_par_init,
      ks_D = ks_D
    ))
  } else {
    # Multi-series fitting for all 7 days together
    all_days_list <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")

    # Create times and marks lists for all 7 days
    all_days_times_list <- list()
    all_days_marks_list <- list()
    for (day in all_days_list) {
      if (day %in% unique(runs_df$WEEKDAY)) {
        day_data <- transform_hawkes_data(runs_df, seqn, day, data_type = data_type)
        all_days_times_list[[day]] <- day_data$times
        all_days_marks_list[[day]] <- day_data$marks
      }
    }

    # Initialization search 
    best_fit <- perform_hawkes_grid_search(
      times_list = all_days_times_list,
      marks_list = all_days_marks_list,
      fit_function = "fit_hawkes_multi_series"
    )

    if (is.null(best_fit)) {
      cat(sprintf("Warning: No converged fits found for SEQN %s. Returning NULL result.\n", seqn))
      return(list(
        id = as.character(seqn),
        period = "all_days",
        parameters = NULL,
        std_errors = NULL,
        pdHess = FALSE,
        penalty_coef = penalty_coef,
        best_a_par_init = NA
      ))
    }

    # best_a_par_init is already set by the helper function
    results_all_days <- check_conv(best_fit)

    # Get KS statistics (per_series = TRUE for multi-series)
    gof_results <- show_hawkes_GOF(best_fit, plot = FALSE, tests = TRUE, return_values = TRUE, per_series = TRUE)
    ks_D <- gof_results$ks_D
    ks_per_series <- gof_results$ks_per_series

    return(list(
      id = as.character(seqn),
      period = "all_days",
      parameters = results_all_days$parameters,
      std_errors = results_all_days$std_errors,
      pdHess = results_all_days$converged,
      penalty_coef = penalty_coef,
      best_a_par_init = best_fit$best_a_par_init,
      ks_D = ks_D,
      ks_per_series = ks_per_series
    ))
  }
}


check_convergence <- function(fits_list) {
  pdHess_values <- sapply(fits_list, function(x) x$pdHess)
  names(pdHess_values) <- names(fits_list)

  successful <- sum(pdHess_values, na.rm = TRUE)
  success_rate <- round(mean(pdHess_values, na.rm = TRUE) * 100, 1)
  failed_fits <- names(pdHess_values)[!pdHess_values & !is.na(pdHess_values)]

  cat("Success rate:", success_rate, "% (", successful, "/", length(fits_list), ")\n")

  # Failure analysis
  if (length(failed_fits) > 0) {
    failure_analysis <- data.frame(
      fit_name = failed_fits,
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        seqn = str_extract(fit_name, "(?<=fit_)[0-9]+"),
        period = str_extract(fit_name, "(weekdays|weekends|all_days)")
      )

    print(table(failure_analysis$period))
    print(head(sort(table(failure_analysis$seqn), decreasing = TRUE), 10))
  } else {
    cat("No failed fits to analyze.\n")
  }

  return(list(
    convergence_status = pdHess_values,
    failed_fits = failed_fits,
    success_rate = success_rate
  ))
}

#' Process Hawkes fits data from list format to dataframe
#' @param fits_list List of Hawkes fit results
#' @return Processed dataframe with parameters
process_fits_data <- function(fits_list) {
  # Helper function to extract parameter by name
  get_param_by_name <- function(params, name, default = NA) {
    if (is.null(params) || is.null(names(params))) return(default)
    if (name %in% names(params)) return(params[name])
    return(default)
  }

  all_fits_data <- do.call(rbind, lapply(fits_list, function(fit) {
    params <- fit$parameters
    ses <- fit$std_errors

    # Extract transformed parameters
    mu <- get_param_by_name(params, "mu")
    alpha <- get_param_by_name(params, "alpha")
    beta <- get_param_by_name(params, "beta")

    # Extract standard errors
    mu_se <- get_param_by_name(ses, "mu")
    alpha_se <- get_param_by_name(ses, "alpha")
    beta_se <- get_param_by_name(ses, "beta")

    # Extract KS statistic
    ks_D <- if (!is.null(fit$ks_D)) fit$ks_D else NA

    # Extract raw parameters
    log_mu <- get_param_by_name(params, "log_mu")
    logit_abratio <- get_param_by_name(params, "logit_abratio")
    log_beta <- get_param_by_name(params, "log_beta")

    # Extract raw standard errors
    log_mu_se <- get_param_by_name(ses, "log_mu")
    logit_abratio_se <- get_param_by_name(ses, "logit_abratio")
    log_beta_se <- get_param_by_name(ses, "log_beta")

    data.frame(
      id = fit$id,
      mu = mu, alpha = alpha, beta = beta,
      mu_se = mu_se, alpha_se = alpha_se, beta_se = beta_se,
      pdHess = fit$pdHess,
      ks_D = ks_D,
      log_mu = log_mu, logit_abratio = logit_abratio, log_beta = log_beta,
      log_mu_se = log_mu_se, logit_abratio_se = logit_abratio_se, log_beta_se = log_beta_se,
      stringsAsFactors = FALSE
    )
  }))

  # Clean up and return
  all_fits_data <- all_fits_data[!is.na(all_fits_data$id), ]
  all_fits_data$id <- as.factor(all_fits_data$id)
  rownames(all_fits_data) <- NULL

  return(all_fits_data)
}