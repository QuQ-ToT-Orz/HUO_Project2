# Declare global variables for R CMD check (ggplot2 aes variables)
utils::globalVariables(c("time", "series"))

#' Creates plots for univariate Hawkes processes
#' Supports both single-series and multi-series models.
#' For multi-series models, intensity is computed for each series separately.
#'
#' @param obj TMB object from Hawkes process fitting
#' @param type Character string: "fitted" (intensity only),
#' "data" (histogram only), or "both" (both plots). Default: "fitted"
#' @param per_series Logical, for multi-series models:
#' if TRUE, show faceted plots for each series;
#' if FALSE (default), show combined/overlaid plot.
#' @return ggplot object or grid arrangement depending on type
#' @importFrom stats sd reshape ppoints quantile
#' @export
show_hawkes <- function(obj, type = c("fitted", "data", "both"), per_series = FALSE) {
  type <- match.arg(type)

  # Extract parameters and data
  times <- obj$env$data$times
  marks <- obj$env$data$marks
  series_id <- obj$env$data$series_id

  # Extract shared parameters
  params <- extract_univariate_params(obj)
  mu <- params$mu
  alpha <- params$alpha
  beta <- params$beta

  # Check if multi-series
  is_multi_series <- !is.null(series_id)

  if (is_multi_series) {
    # Multi-series: process each series separately
    n_series <- obj$env$data$n_series

    intensity_data_list <- list()
    event_data_list <- list()

    for (s in 0:(n_series - 1)) {
      # Get indices for this series
      idx <- which(series_id == s)
      if (length(idx) == 0) next

      times_s <- times[idx]
      marks_s <- if (!is.null(marks)) marks[idx] else rep(1, length(times_s))

      # Sort by time within series
      ord <- order(times_s)
      times_s <- times_s[ord]
      marks_s <- marks_s[ord]

      # Calculate intensity for this series
      max_time_s <- max(times_s)
      p_s <- seq(0, max_time_s, length.out = 200)
      intensity_s <- hawkes_intensity(
        times = times_s,
        mu = mu,
        alpha = alpha,
        beta = beta,
        p = p_s,
        marks = marks_s
      )

      intensity_data_list[[s + 1]] <- data.frame(
        time = p_s,
        intensity = intensity_s,
        series = paste("Series", s + 1)
      )

      event_data_list[[s + 1]] <- data.frame(
        times = times_s,
        series = paste("Series", s + 1)
      )
    }

    intensity_data <- do.call(rbind, intensity_data_list)
    event_data <- do.call(rbind, event_data_list)

    if (per_series) {
      # Faceted plots for each series
      intensity_plot <- ggplot2::ggplot(
        data = intensity_data,
        ggplot2::aes(x = time, y = intensity)
      ) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~series, scales = "free") +
        ggplot2::labs(x = "Time", y = expression(lambda(t)), title = "Intensity by Series") +
        ggplot2::theme_minimal()

      n_bins <- floor(min(30, nrow(event_data) / (5 * n_series)))
      n_bins <- max(1, n_bins)

      history_plot <- ggplot2::ggplot(
        data = event_data,
        ggplot2::aes(x = times)
      ) +
        ggplot2::geom_histogram(bins = n_bins) +
        ggplot2::facet_wrap(~series, scales = "free") +
        ggplot2::labs(x = "Time", y = "Count", title = "Event Distribution by Series") +
        ggplot2::theme_minimal()

    } else {
      # Combined/overlaid plots (default)
      intensity_plot <- ggplot2::ggplot(
        data = intensity_data,
        ggplot2::aes(x = time, y = intensity, color = series)
      ) +
        ggplot2::geom_line(alpha = 0.7) +
        ggplot2::labs(
          x = "Time",
          y = expression(lambda(t)),
          title = paste("Multi-Series Intensity (", n_series, " series)", sep = ""),
          color = "Series"
        ) +
        ggplot2::theme_minimal()

      n_bins <- floor(min(30, nrow(event_data) / 5))
      n_bins <- max(1, n_bins)

      history_plot <- ggplot2::ggplot(
        data = event_data,
        ggplot2::aes(x = times, fill = series)
      ) +
        ggplot2::geom_histogram(bins = n_bins, alpha = 0.5, position = "identity") +
        ggplot2::labs(
          x = "Time",
          y = "Count",
          title = paste("Event Distribution (", n_series, " series)", sep = ""),
          fill = "Series"
        ) +
        ggplot2::theme_minimal()
    }

  } else {
    # Single series: original behavior
    if (is.null(marks)) {
      marks <- rep(1, length(times))
    }

    # Calculate intensity over time range
    max_time <- max(times)
    p <- seq(0, max_time, length.out = 500)
    intensity <- hawkes_intensity(
      times = times,
      mu = mu,
      alpha = alpha,
      beta = beta,
      p = p,
      marks = marks
    )

    # Create intensity plot
    intensity_plot <- ggplot2::ggplot(
      data = data.frame(time = p, intensity = intensity),
      ggplot2::aes(x = time, y = intensity)
    ) +
      ggplot2::geom_line() +
      ggplot2::labs(x = "Time", y = expression(lambda(t))) +
      ggplot2::theme_minimal()

    # Create histogram of events
    n_bins <- floor(min(30, length(times) / 5))
    n_bins <- max(1, n_bins)

    history_plot <- ggplot2::ggplot(
      data = data.frame(times = times),
      ggplot2::aes(x = times)
    ) +
      ggplot2::geom_histogram(bins = n_bins) +
      ggplot2::labs(x = "Time", y = "Count") +
      ggplot2::theme_minimal()
  }

  # Return appropriate plot(s)
  switch(type,
    "fitted" = intensity_plot,
    "data" = history_plot,
    "both" = gridExtra::grid.arrange(intensity_plot, history_plot, ncol = 1)
  )
}

#' Goodness-of-fit diagnostics for univariate Hawkes process
#'
#' Uses the compensator transformation and checks for exponential interarrival times.
#' Supports both single-series and multi-series models. 
#' For multi-series models, the compensator is computed for each series independently,
#' then interarrivals are pooled for GOF testing (default) or tested per-series.
#'
#' @param obj TMB object from Hawkes process fitting
#' @param plot Logical, whether to produce diagnostic plots. Default: TRUE
#' @param tests Logical, whether to perform statistical tests. Default: TRUE
#' @param return_values Logical, whether to return computed values. Default: FALSE
#' @param per_series Logical, for multi-series models: 
#' if TRUE, run separate GOF tests for each series; 
#' if FALSE (default), pool interarrivals across series.
#' @return If return_values=TRUE, returns list with interarrival times (and per-series results for multi-series)
#' @export
show_hawkes_GOF <- function(obj, plot = TRUE, tests = TRUE, return_values = FALSE, per_series = FALSE) {
  # Extract parameters and data
  times <- obj$env$data$times
  marks <- obj$env$data$marks
  series_id <- obj$env$data$series_id

  # Extract shared parameters
  params <- extract_univariate_params(obj)
  mu <- params$mu
  alpha <- params$alpha
  beta <- params$beta

  # Check if multi-series
  is_multi_series <- !is.null(series_id)

  if (is_multi_series) {
    # Multi-series: process each series independently
    n_series <- obj$env$data$n_series
    if (is.null(n_series)) n_series <- length(unique(series_id))

    all_interarrivals <- list()
    all_U <- list()
    series_results <- list()

    for (s in 0:(n_series - 1)) {
      # Get indices for this series
      idx <- which(series_id == s)
      if (length(idx) < 2) next

      times_s <- times[idx]
      marks_s <- if (!is.null(marks)) marks[idx] else rep(1, length(times_s))

      # Sort by time within series
      ord <- order(times_s)
      times_s <- times_s[ord]
      marks_s <- marks_s[ord]

      # Calculate A vector for this series (Ozaki's algorithm)
      A_s <- numeric(length(times_s))
      for (i in 2:length(times_s)) {
        A_s[i] <- exp(-beta * (times_s[i] - times_s[i - 1])) * (marks_s[i - 1] + A_s[i - 1])
      }

      # Calculate compensator for this series
      compensator_s <- numeric(length(times_s))
      for (i in seq_along(times_s)) {
        compensator_s[i] <- (mu * times_s[i]) -
          ((alpha / beta) * A_s[i]) +
          ((alpha / beta) * (sum(marks_s[1:i]) - marks_s[i]))
      }

      # Calculate interarrivals for this series
      interarrivals_s <- compensator_s[2:length(compensator_s)] - compensator_s[1:(length(compensator_s) - 1)]
      all_interarrivals[[s + 1]] <- interarrivals_s

      # Calculate U for scatter plot
      U_s <- 1 - exp(-interarrivals_s)
      all_U[[s + 1]] <- U_s

      series_results[[paste0("series_", s + 1)]] <- list(
        times = times_s,
        compensator = compensator_s,
        interarrivals = interarrivals_s,
        n_events = length(times_s)
      )
    }

    # Pool all interarrivals
    interarrivals <- unlist(all_interarrivals)
    U_pooled <- unlist(all_U)

    if (per_series) {
      # Per-series GOF tests
      if (tests) {
        cat("Per-series GOF tests (", length(series_results), " series):\n")
        cat(paste(rep("=", 60), collapse = ""), "\n\n")

        test_results <- list()
        for (name in names(series_results)) {
          sr <- series_results[[name]]
          if (length(sr$interarrivals) < 5) {
            cat(name, ": insufficient data (", sr$n_events, " events)\n\n")
            test_results[[name]] <- list(
              n_events = sr$n_events,
              ks_pvalue = NA,
              ljungbox_pvalue = NA,
              status = "insufficient_data"
            )
            next
          }

          ks_test <- stats::ks.test(sr$interarrivals, "pexp", rate = 1)
          lb_test <- stats::Box.test(sr$interarrivals, type = "Ljung")

          # Diagnostic statistics for Exp(1): mean should be ~1, sd should be ~1
          ia_mean <- mean(sr$interarrivals)
          ia_sd <- sd(sr$interarrivals)

          cat(name, " (n=", sr$n_events, "):\n", sep = "")
          cat("  Interarrival mean:      ", format(ia_mean, digits = 4), " (expect ~1)\n")
          cat("  Interarrival SD:        ", format(ia_sd, digits = 4), " (expect ~1)\n")
          cat("  KS D statistic:         ", format(ks_test$statistic, digits = 4), "\n")
          cat("  Ljung-Box p-value:      ", format(lb_test$p.value, digits = 4), "\n\n")

          test_results[[name]] <- list(
            n_events = sr$n_events,
            ia_mean = ia_mean,
            ia_sd = ia_sd,
            ks_D = as.numeric(ks_test$statistic),
            ks_pvalue = ks_test$p.value,
            ljungbox_pvalue = lb_test$p.value,
            status = "success"
          )
        }

        # Summary table
        cat(paste(rep("=", 60), collapse = ""), "\n")
        cat("Summary:\n")
        ia_means <- sapply(test_results, function(x) x$ia_mean)
        ia_sds <- sapply(test_results, function(x) x$ia_sd)
        ks_Ds <- sapply(test_results, function(x) x$ks_D)
        lb_pvals <- sapply(test_results, function(x) x$ljungbox_pvalue)

        cat("  Mean of interarrivals (pooled):         ", format(mean(ia_means, na.rm = TRUE), digits = 4), " (expect ~1)\n")
        cat("  SD of interarrivals (pooled):           ", format(mean(ia_sds, na.rm = TRUE), digits = 4), " (expect ~1)\n")
        cat("  Average KS D statistic:                 ", format(mean(ks_Ds, na.rm = TRUE), digits = 4), "\n")
        cat("  Max KS D statistic:                     ", format(max(ks_Ds, na.rm = TRUE), digits = 4), "\n")
        n_failed_lb <- sum(lb_pvals < 0.05, na.rm = TRUE)
        cat("  Series failing Ljung-Box (p<0.05):      ", n_failed_lb, "/", length(test_results), "\n")

        # Store test results in series_results
        for (name in names(test_results)) {
          series_results[[name]]$ia_mean <- test_results[[name]]$ia_mean
          series_results[[name]]$ia_sd <- test_results[[name]]$ia_sd
          series_results[[name]]$ks_D <- test_results[[name]]$ks_D
          series_results[[name]]$ks_pvalue <- test_results[[name]]$ks_pvalue
          series_results[[name]]$ljungbox_pvalue <- test_results[[name]]$ljungbox_pvalue
        }
      }

      if (plot) {
        # 1. KS D statistic by series (D < 0.1 is good), with n on x-axis labels
        ks_data <- data.frame(
          series = sapply(names(series_results), function(name) {
            paste0(name, "\n(n=", series_results[[name]]$n_events, ")")
          }),
          D = sapply(series_results, function(x) {
            if (length(x$interarrivals) < 5) return(NA)
            stats::ks.test(x$interarrivals, "pexp", rate = 1)$statistic
          })
        )
        ks_data$quality <- ifelse(is.na(ks_data$D), "NA",
          ifelse(ks_data$D < 0.1, "D<0.1 (good)", "D>=0.1"))

        p1 <- ggplot2::ggplot(
          ks_data,
          ggplot2::aes(x = .data$series, y = .data$D, fill = .data$quality)
        ) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::scale_fill_manual(values = c("D<0.1 (good)" = "steelblue", "D>=0.1" = "coral", "NA" = "gray")) +
          ggplot2::labs(x = "Series", y = "KS D Statistic", title = "KS D Statistic by Series (lower is better)") +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")

        # 2. Faceted compensator vs observed events plots for each series
        comp_data_list <- list()
        for (name in names(series_results)) {
          sr <- series_results[[name]]
          if (length(sr$times) < 2) next
          label <- paste0(name, " (n=", sr$n_events, ")")
          n_events <- length(sr$times)
          comp_data_list[[name]] <- data.frame(
            xs = rep(sr$times, 2),
            val = c(seq_along(sr$times), sr$compensator),
            type = rep(c("observed", "compensator"), each = n_events),
            series = label
          )
        }
        comp_data_all <- do.call(rbind, comp_data_list)

        p2 <- ggplot2::ggplot(
          comp_data_all,
          ggplot2::aes(x = .data$xs, y = .data$val, colour = .data$type)
        ) +
          ggplot2::geom_line() +
          ggplot2::facet_wrap(~series, scales = "free") +
          ggplot2::labs(
            x = "Time",
            y = "Events",
            title = expression("Actual Events vs Compensator (" * Lambda * ") by Series")
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            strip.text = ggplot2::element_text(size = 9),
            legend.position = "bottom",
            axis.text.x = ggplot2::element_text(size = 6),
            axis.text.y = ggplot2::element_text(size = 6)
          )

        # 3. Faceted Q-Q plots for each series (with n in label)
        qq_data_list <- list()
        for (name in names(series_results)) {
          sr <- series_results[[name]]
          if (length(sr$interarrivals) < 5) next
          p_pts <- ppoints(min(100, length(sr$interarrivals)))
          q_obs <- quantile(sr$interarrivals, p = p_pts)
          label <- paste0(name, " (n=", sr$n_events, ")")
          qq_data_list[[name]] <- data.frame(
            theoretical = stats::qexp(p_pts),
            observed = q_obs,
            series = label
          )
        }
        qq_data_all <- do.call(rbind, qq_data_list)

        p3 <- ggplot2::ggplot(qq_data_all, ggplot2::aes(x = .data$theoretical, y = .data$observed)) +
          ggplot2::geom_point(size = 1, alpha = 0.6) +
          ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
          ggplot2::facet_wrap(~series, scales = "free") +
          ggplot2::labs(x = "Theoretical Quantiles (Exp(1))", y = "Observed Quantiles", title = "Q-Q Plots by Series") +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            strip.text = ggplot2::element_text(size = 9),
            axis.text.x = ggplot2::element_text(size = 6),
            axis.text.y = ggplot2::element_text(size = 6)
          )

        # Plot 1: KS D bar + Compensator vs observed
        print(gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(1, 2)))
        # Plot 2: KS D bar + Q-Q plots
        print(gridExtra::grid.arrange(p1, p3, ncol = 1, heights = c(1, 2)))
      }

      if (return_values) {
        # Compute KS D statistics per series
        ks_per_series <- sapply(series_results, function(x) {
          if (length(x$interarrivals) < 5) return(NA)
          as.numeric(stats::ks.test(x$interarrivals, "pexp", rate = 1)$statistic)
        })
        # Average KS D across series
        ks_D <- mean(ks_per_series, na.rm = TRUE)

        return(list(
          interarrivals = interarrivals,
          series_results = series_results,
          n_series = n_series,
          per_series = TRUE,
          ks_per_series = ks_per_series,
          ks_D = ks_D
        ))
      }

    } else {
      # Pooled GOF tests (default)
      if (tests) {
        ks_test <- stats::ks.test(interarrivals, "pexp", rate = 1)
        lb_test <- stats::Box.test(interarrivals, type = "Ljung")

        cat("Multi-series GOF (pooled interarrivals from", n_series, "series):\n")
        cat(paste(rep("=", 60), collapse = ""), "\n\n")
        cat("Sample size:              ", length(interarrivals), "\n")
        cat("Interarrival mean:        ", format(mean(interarrivals), digits = 4), " (expect ~1)\n")
        cat("Interarrival SD:          ", format(sd(interarrivals), digits = 4), " (expect ~1)\n")
        cat("KS D statistic:           ", format(ks_test$statistic, digits = 4), "\n")
        cat("Ljung-Box p-value:        ", format(lb_test$p.value, digits = 4), "\n")
      }

      if (plot) {
        # 1. Pooled Q-Q plot
        p <- ppoints(100)
        q <- quantile(interarrivals, p = p)
        qq_data <- data.frame(
          theoretical = stats::qexp(p),
          observed = q
        )

        p1 <- ggplot2::ggplot(
          qq_data,
          ggplot2::aes(x = .data$theoretical, y = .data$observed)
        ) +
          ggplot2::geom_point() +
          ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
          ggplot2::labs(
            x = "Theoretical Quantiles",
            y = "Observed Quantiles",
            title = "Q-Q Plot (Pooled)"
          ) +
          ggplot2::theme_minimal()

        # 2. Histogram of pooled interarrivals
        binwidth <- if (length(interarrivals) > 1500) 0.05 else 0.1
        hist_data <- data.frame(
          data = interarrivals[interarrivals < 4]
        )

        p2 <- ggplot2::ggplot(
          hist_data,
          ggplot2::aes(x = .data$data)
        ) +
          ggplot2::geom_histogram(binwidth = binwidth) +
          ggplot2::labs(
            x = "Interarrival Times",
            y = "Count",
            title = "Interarrival Times (Pooled)"
          ) +
          ggplot2::theme_minimal()

        # 3. Scatter plot of consecutive interarrivals (within each series, then pooled)
        scatter_data <- data.frame(
          x = U_pooled[-length(U_pooled)],
          y = U_pooled[-1]
        )

        p3 <- ggplot2::ggplot(
          scatter_data,
          ggplot2::aes(x = .data$x, y = .data$y)
        ) +
          ggplot2::geom_point(alpha = 0.5) +
          ggplot2::labs(
            x = expression("F(" * Lambda[k] ~ -Lambda[k - 1] * ")"),
            y = expression("F(" * Lambda[k + 1] ~ -Lambda[k] * ")"),
            title = "Consecutive Interarrivals (Pooled)"
          ) +
          ggplot2::theme_minimal()

        # 4. Per-series event count summary
        series_counts <- sapply(series_results, function(x) x$n_events)
        count_data <- data.frame(
          series = names(series_counts),
          n_events = as.numeric(series_counts)
        )

        p4 <- ggplot2::ggplot(
          count_data,
          ggplot2::aes(x = .data$series, y = .data$n_events)
        ) +
          ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
          ggplot2::labs(
            x = "Series",
            y = "Number of Events",
            title = "Events per Series"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

        gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
      }

      if (return_values) {
        return(list(
          interarrivals = interarrivals,
          series_results = series_results,
          n_series = n_series,
          per_series = FALSE,
          ks_D = as.numeric(ks_test$statistic)
        ))
      }
    }

  } else {
    # Single series: original behavior
    if (is.null(marks)) {
      marks <- rep(1, length(times))
    }

    # Calculate A vector using Ozaki's efficient algorithm (matches C++)
    A <- numeric(length(times))
    for (i in 2:length(times)) {
      A[i] <- exp(-beta * (times[i] - times[i - 1])) * (marks[i - 1] + A[i - 1])
    }

    # Calculate compensator: Λ(t_i) = μ*t_i - (α/β)*A_i + (α/β)*(∑marks[1:i] - marks[i])
    compensator <- numeric(length(times))
    for (i in seq_along(times)) {
      compensator[i] <- (mu * times[i]) -
        ((alpha / beta) * A[i]) +
        ((alpha / beta) * (sum(marks[1:i]) - marks[i]))
    }

    # Calculate interarrivals using explicit indexing
    interarrivals <- compensator[2:length(compensator)] - compensator[1:(length(compensator) - 1)]

    # Statistical tests
    if (tests) {
      ks_test <- stats::ks.test(interarrivals, "pexp", rate = 1)
      lb_test <- stats::Box.test(interarrivals, type = "Ljung")

      cat("Goodness-of-Fit Diagnostics:\n")
      cat(paste(rep("=", 60), collapse = ""), "\n\n")
      cat("Sample size:              ", length(interarrivals), "\n")
      cat("Interarrival mean:        ", format(mean(interarrivals), digits = 4), " (expect ~1)\n")
      cat("Interarrival SD:          ", format(sd(interarrivals), digits = 4), " (expect ~1)\n")
      cat("KS D statistic:           ", format(ks_test$statistic, digits = 4), "\n")
      cat("Ljung-Box p-value:        ", format(lb_test$p.value, digits = 4), "\n")
    }

    # Diagnostic plots
    if (plot) {
      # 1. Compensator vs observed events
      plot_data <- data.frame(xs = times, observed = seq_along(times), compensator = compensator)
      plot_data <- reshape(plot_data,
        direction = "long", idvar = "xs",
        varying = c("observed", "compensator"), v.names = "val",
        times = c("observed", "compensator"),
        new.row.names = 1:(2 * length(times))
      )

      p1 <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
          x = .data$xs, y = .data$val,
          colour = .data$time
        )
      ) +
        ggplot2::geom_line() +
        ggplot2::labs(
          x = "Time",
          y = "Events",
          title = expression("Actual Events and Compensator(" * Lambda * ")")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = c(0.8, 0.2))

      # 2. Q-Q plot for exponential distribution
      p <- ppoints(100)
      q <- quantile(interarrivals, p = p)
      qq_data <- data.frame(
        theoretical = stats::qexp(p),
        observed = q
      )

      p2 <- ggplot2::ggplot(
        qq_data,
        ggplot2::aes(x = .data$theoretical, y = .data$observed)
      ) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
        ggplot2::labs(
          x = "Theoretical Quantiles",
          y = "Observed Quantiles",
          title = "Q-Q Plot"
        ) +
        ggplot2::theme_minimal()

      # 3. Histogram with adaptive binwidth
      binwidth <- if (length(interarrivals) > 1500) 0.05 else 0.1
      hist_data <- data.frame(
        data = interarrivals[interarrivals < 4]
      )

      p3 <- ggplot2::ggplot(
        hist_data,
        ggplot2::aes(x = .data$data)
      ) +
        ggplot2::geom_histogram(binwidth = binwidth) +
        ggplot2::labs(
          x = "Interarrival Times",
          y = "Count",
          title = "Interarrival Times Distribution"
        ) +
        ggplot2::theme_minimal()

      # 4. Scatter plot of consecutive interarrivals
      U <- 1 - exp(-compensator[2:length(compensator)] + compensator[1:(length(compensator) - 1)])
      scatter_data <- data.frame(
        x = U[-length(U)],
        y = U[-1]
      )

      p4 <- ggplot2::ggplot(
        scatter_data,
        ggplot2::aes(x = .data$x, y = .data$y)
      ) +
        ggplot2::geom_point() +
        ggplot2::labs(
          x = expression("F(" * Lambda[k] ~ -Lambda[k - 1] * ")"),
          y = expression("F(" * Lambda[k + 1] ~ -Lambda[k] * ")"),
          title = "Consecutive Interarrivals"
        ) +
        ggplot2::theme_minimal()

      gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
    }

    if (return_values) {
      return(list(
        interarrivals = interarrivals,
        ks_D = as.numeric(ks_test$statistic)
      ))
    }
  }
}

#' Plots the time series of interarrival times for Hawkes processes,
#' showing patterns of regularity or clustering over the event sequence.
#'
#' Supports both single-series and multi-series models. 
#' For multi-series models, interarrivals are computed for each series independently.
#'
#' @param obj TMB object from Hawkes process fitting
#' @param return_values Logical, whether to return computed values. Default: FALSE
#' @param per_series Logical, for multi-series models: 
#' if TRUE, show separate plots for each series; 
#' if FALSE (default), pool and show combined plot.
#' @return If return_values=TRUE, returns list with interarrival times and process type
#' @export
show_interarrival_pattern <- function(obj, return_values = FALSE, per_series = FALSE) {
  # Extract parameters and data
  times <- obj$env$data$times
  marks <- obj$env$data$marks
  series_id <- obj$env$data$series_id

  # Extract shared parameters
  params <- extract_univariate_params(obj)
  mu <- params$mu
  alpha <- params$alpha
  beta <- params$beta

  # Check if multi-series
  is_multi_series <- !is.null(series_id)

  if (is_multi_series) {
    # Multi-series: process each series independently
    n_series <- obj$env$data$n_series
    if (is.null(n_series)) n_series <- length(unique(series_id))

    all_interarrivals <- list()
    series_results <- list()

    for (s in 0:(n_series - 1)) {
      # Get indices for this series
      idx <- which(series_id == s)
      if (length(idx) < 2) next

      times_s <- times[idx]
      marks_s <- if (!is.null(marks)) marks[idx] else rep(1, length(times_s))

      # Sort by time within series
      ord <- order(times_s)
      times_s <- times_s[ord]
      marks_s <- marks_s[ord]

      # Calculate A vector for this series (Ozaki's algorithm)
      A_s <- numeric(length(times_s))
      for (i in 2:length(times_s)) {
        A_s[i] <- exp(-beta * (times_s[i] - times_s[i - 1])) * (marks_s[i - 1] + A_s[i - 1])
      }

      # Calculate compensator for this series
      compensator_s <- numeric(length(times_s))
      for (i in seq_along(times_s)) {
        compensator_s[i] <- (mu * times_s[i]) -
          ((alpha / beta) * A_s[i]) +
          ((alpha / beta) * (sum(marks_s[1:i]) - marks_s[i]))
      }

      # Calculate interarrivals for this series
      interarrivals_s <- compensator_s[2:length(compensator_s)] - compensator_s[1:(length(compensator_s) - 1)]
      all_interarrivals[[s + 1]] <- interarrivals_s

      series_results[[paste0("series_", s + 1)]] <- list(
        times = times_s,
        compensator = compensator_s,
        interarrivals = interarrivals_s,
        n_events = length(times_s)
      )
    }

    if (per_series) {
      # Faceted plot for each series
      ts_data_list <- list()
      for (name in names(series_results)) {
        sr <- series_results[[name]]
        if (length(sr$interarrivals) < 2) next
        ts_data_list[[name]] <- data.frame(
          index = seq_along(sr$interarrivals),
          interarrival = sr$interarrivals,
          series = name
        )
      }
      ts_data <- do.call(rbind, ts_data_list)

      p <- ggplot2::ggplot(
        ts_data,
        ggplot2::aes(x = .data$index, y = .data$interarrival)
      ) +
        ggplot2::geom_line(alpha = 0.6) +
        ggplot2::geom_smooth(se = TRUE, color = if (alpha >= 0) "blue" else "red") +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        ggplot2::facet_wrap(~series, scales = "free") +
        ggplot2::labs(
          x = "Event Index",
          y = "Interarrival Time",
          title = "Interarrival Times by Series",
          subtitle = if (alpha >= 0) "Self-exciting" else "Self-inhibiting"
        ) +
        ggplot2::theme_minimal()

    } else {
      # Pooled plot
      interarrivals <- unlist(all_interarrivals)

      ts_data <- data.frame(
        index = seq_along(interarrivals),
        interarrival = interarrivals
      )

      p <- ggplot2::ggplot(
        ts_data,
        ggplot2::aes(x = .data$index, y = .data$interarrival)
      ) +
        ggplot2::geom_line(alpha = 0.6) +
        ggplot2::geom_smooth(se = TRUE, color = if (alpha >= 0) "blue" else "red") +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        ggplot2::labs(
          x = "Event Index (Pooled)",
          y = "Interarrival Time",
          title = paste("Interarrival Times (Pooled from", n_series, "series)"),
          subtitle = if (alpha >= 0) "Self-exciting" else "Self-inhibiting"
        ) +
        ggplot2::theme_minimal()
    }

    print(p)

    if (return_values) {
      return(list(
        interarrivals = unlist(all_interarrivals),
        series_results = series_results,
        process_type = if (alpha >= 0) "self-exciting" else "self-inhibiting",
        alpha = alpha,
        mu = mu,
        beta = beta,
        n_series = n_series
      ))
    }

  } else {
    # Single series: original behavior
    if (is.null(marks)) {
      marks <- rep(1, length(times))
    }

    # Calculate A vector using Ozaki's efficient algorithm (matches C++)
    A <- numeric(length(times))
    for (i in 2:length(times)) {
      A[i] <- exp(-beta * (times[i] - times[i - 1])) * (marks[i - 1] + A[i - 1])
    }

    # Calculate compensator
    compensator <- numeric(length(times))
    for (i in seq_along(times)) {
      compensator[i] <- (mu * times[i]) -
        ((alpha / beta) * A[i]) +
        ((alpha / beta) * (sum(marks[1:i]) - marks[i]))
    }

    # Calculate interarrivals
    interarrivals <- compensator[2:length(compensator)] - compensator[1:(length(compensator) - 1)]

    # Create time series plot of interarrivals
    ts_data <- data.frame(
      index = seq_along(interarrivals),
      interarrival = interarrivals
    )

    p <- ggplot2::ggplot(
      ts_data,
      ggplot2::aes(x = .data$index, y = .data$interarrival)
    ) +
      ggplot2::geom_line(alpha = 0.6) +
      ggplot2::geom_smooth(se = TRUE, color = if (alpha >= 0) "blue" else "red") +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      ggplot2::labs(
        x = "Event Index",
        y = "Interarrival Time",
        title = "Interarrival Times Over Event Sequence",
        subtitle = if (alpha >= 0) "Self-exciting: clustering may cause variability" else "Self-inhibiting: expect more regularity"
      ) +
      ggplot2::theme_minimal()

    print(p)

    if (return_values) {
      return(list(
        interarrivals = interarrivals,
        compensator = compensator,
        process_type = if (alpha >= 0) "self-exciting" else "self-inhibiting",
        alpha = alpha,
        mu = mu,
        beta = beta
      ))
    }
  }
}