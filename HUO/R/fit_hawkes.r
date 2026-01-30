#' @importFrom stats median

# Common bound settings with optional minimums and maximums
# @noRd
set_bounds <- function(n_mu, n_alpha, n_beta,
                       min_mu = NULL, max_mu = NULL,
                       min_ratio = NULL, max_ratio = NULL,
                       min_beta = NULL, max_beta = NULL) {
  lower_bounds <- c()
  upper_bounds <- c()

  # mu bounds (reasonable defaults: 1e-3 to 1e3)
  if (!is.null(min_mu)) {
    lower_bounds <- c(lower_bounds, rep(log(min_mu), n_mu))
  } else {
    lower_bounds <- c(lower_bounds, rep(log(1e-3), n_mu))
  }

  if (!is.null(max_mu)) {
    upper_bounds <- c(upper_bounds, rep(log(max_mu), n_mu))
  } else {
    upper_bounds <- c(upper_bounds, rep(log(1e3), n_mu))
  }

  # logit_abratio bounds (reasonable defaults: -10 to 10)
  if (!is.null(min_ratio)) {
    lower_bounds <- c(lower_bounds, rep(stats::qlogis(min_ratio), n_alpha))
  } else {
    lower_bounds <- c(lower_bounds, rep(-1e3, n_alpha))
  }

  if (!is.null(max_ratio)) {
    upper_bounds <- c(upper_bounds, rep(stats::qlogis(max_ratio), n_alpha))
  } else {
    upper_bounds <- c(upper_bounds, rep(1e3, n_alpha))
  }

  # beta bounds (reasonable defaults: 1e-6 to 1e6)
  if (!is.null(min_beta)) {
    lower_bounds <- c(lower_bounds, rep(log(min_beta), n_beta))
  } else {
    lower_bounds <- c(lower_bounds, rep(log(1e-6), n_beta))
  }

  if (!is.null(max_beta)) {
    upper_bounds <- c(upper_bounds, rep(log(max_beta), n_beta))
  } else {
    upper_bounds <- c(upper_bounds, rep(log(1e6), n_beta))
  }

  list(lower = lower_bounds, upper = upper_bounds)
}

# Bounds for model 2 (negative alpha allowed) 
# - only mu and beta bounds, no alpha bounds
set_bounds_model2 <- function(n_mu, n_alpha, n_beta,
                              min_mu = NULL, max_mu = NULL,
                              min_alpha = NULL, max_alpha = NULL,
                              min_beta = NULL, max_beta = NULL) {
  lower_bounds <- c()
  upper_bounds <- c()

  # mu bounds (reasonable defaults: 1e-3 to 1e3)
  if (!is.null(min_mu)) {
    lower_bounds <- c(lower_bounds, rep(log(min_mu), n_mu))
  } else {
    lower_bounds <- c(lower_bounds, rep(log(1e-3), n_mu))
  }

  if (!is.null(max_mu)) {
    upper_bounds <- c(upper_bounds, rep(log(max_mu), n_mu))
  } else {
    upper_bounds <- c(upper_bounds, rep(log(1e3), n_mu))
  }

  # a_par bounds (for tanh transformation, reasonable defaults: -10 to 10)
  if (!is.null(min_alpha)) {
    lower_bounds <- c(lower_bounds, rep(min_alpha, n_alpha))
  } else {
    lower_bounds <- c(lower_bounds, rep(-10, n_alpha))
  }

  if (!is.null(max_alpha)) {
    upper_bounds <- c(upper_bounds, rep(max_alpha, n_alpha))
  } else {
    upper_bounds <- c(upper_bounds, rep(10, n_alpha))
  }

  # beta bounds (reasonable defaults: 1e-6 to 1e6)
  if (!is.null(min_beta)) {
    lower_bounds <- c(lower_bounds, rep(log(min_beta), n_beta))
  } else {
    lower_bounds <- c(lower_bounds, rep(log(1e-6), n_beta))
  }

  if (!is.null(max_beta)) {
    upper_bounds <- c(upper_bounds, rep(log(max_beta), n_beta))
  } else {
    upper_bounds <- c(upper_bounds, rep(log(1e6), n_beta))
  }

  list(lower = lower_bounds, upper = upper_bounds)
}

# Helper function for univariate parameter validation and initialization
validate_params <- function(
    times, marks = NULL, model = 1,
    parameters = list(), min_events = 5) {
  # Validate times input
  if (!is.numeric(times)) {
    stop("times must be numeric")
  }
  if (length(times) < min_events) {
    stop("times must contain sufficient events")
  }
  if (any(diff(times) < 1.e-10)) {
    stop("times must be in ascending order with no simultaneous events")
  }
  if (any(is.na(times))) {
    stop("times cannot contain NA values")
  }

  # Validate marks if provided
  if (!is.null(marks)) {
    if (!is.numeric(marks)) {
      stop("marks must be numeric")
    }
    if (length(marks) != length(times)) {
      stop("marks must have same length as times")
    }
    if (any(is.na(marks))) {
      stop("marks cannot contain NA values")
    }
    if (any(marks < 0)) {
      stop("marks cannot be negative")
    }
    mean_marks <- mean(marks)
  }

  # Initialize parameters if not provided
  alpha <- parameters[["alpha"]]
  beta <- parameters[["beta"]]
  mu <- parameters[["mu"]]
  # Set default values if not provided
  if (is.null(mu)) {
    mu <- 0.5 * length(times) / max(times)
  }
  if (any(mu <= 0)) {
    stop("mu must be positive")
  }

  if (is.null(beta)) {
    # Calculate inter-arrival times
    inter_arrival <- diff(times)
    median_inter_arrival <- median(inter_arrival)
    beta <- 1 / median_inter_arrival
  }
  if (any(beta <= 0)) {
    stop("beta must be positive")
  }

  if (model == 1) {
    if (is.null(alpha)) {
      # Initialize alpha based on marks if provided
      if (!is.null(marks)) {
        alpha <- 0.3 * beta / mean(marks)
      } else {
        alpha <- 0.3 * beta
      }
    }

    if (any(alpha < 0)) {
      stop("alpha must be non-negative")
    }

    # Validate alpha/beta constraints with marks
    if (!is.null(marks)) {
      max_alpha <- beta / mean(marks)
      if (any(alpha >= max_alpha)) {
        stop(sprintf(
          "alpha = %.2f must be < beta/mean_marks = %.2f",
          alpha, max_alpha
        ))
      }
    } else {
      if (any(alpha >= beta)) {
        stop(sprintf(
          "alpha = %.2f must be < beta = %.2f",
          alpha, beta
        ))
      }
    }
  } else if (model == 2) {
    if (is.null(alpha)) {
      if (!is.null(marks)) {
        alpha <- 0.1 * beta / mean(marks)
      } else {
        alpha <- 0.1 * beta
      }
    }
  }

  # Validate parameter types
  if (!is.numeric(alpha) || !is.numeric(beta) || !is.numeric(mu)) {
    stop("alpha, beta, and mu must be numeric")
  }

  result <- list(
    alpha = alpha,
    beta = beta,
    mu = mu
  )

  # Add mean_marks if marks were provided
  if (!is.null(marks)) {
    result$mean_marks <- mean_marks
  }

  return(result)
}

#' Fit univariate Hawkes process
#' @param times Numeric vector of event times
#' @param marks Numeric vector of event marks (same length as times)
#' @param parameters List of initial parameter values:
#'   For model = 1: (mu, alpha, beta)
#'   For model = 2: (mu, a_par, beta)
#' @param model Numeric indicator specifying which model to fit:
#'   \itemize{
#'   \item \code{model = 1}, fits a Hawkes process with positive alpha (default);
#'   \item \code{model = 2}, fits a Hawkes process with alpha that can be negative.
#'   }
#' @param penalty_coef Numeric penalty coefficient for regularization (default: 0)
#' @param tmb_silent Logical, whether to suppress TMB output (default: TRUE)
#' @param ... Additional arguments passed to optimization
#' @return TMB object with fitted parameters
#' @export
fit_hawkes <- function(times,
                              marks,
                              parameters = list(),
                              model = 1,
                              penalty_coef = 0,
                              tmb_silent = TRUE,
                              ...) {
  validated_params <- validate_params(times, marks, 
    model = model, parameters = parameters)

  alpha <- validated_params$alpha
  beta <- validated_params$beta
  mu <- validated_params$mu
  mean_marks <- validated_params$mean_marks
  
  if (model == 1) {
    # Model 1: Positive alpha constraint
    # Setup TMB for positive alpha model
    obj <- TMB::MakeADFun(
      data = list(
        times = times,
        marks = marks,
        penalty_coef = penalty_coef,
        model_type = "hawkes"
      ),
      parameters = list(
        log_mu = log(mu),
        logit_abratio = stats::qlogis(alpha / (beta / mean_marks)),
        log_beta = log(beta)
      ),
      hessian = TRUE,
      DLL = "HUO",
      silent = tmb_silent
    )

    # Set bounds for positive alpha optimization
    bounds <- set_bounds(
      n_mu = 1,
      n_alpha = 1,
      n_beta = 1
    )

    opt <- stats::optim(
      par = obj$par,
      fn = obj$fn,
      gr = obj$gr,
      lower = bounds$lower,
      upper = bounds$upper,
      method = "L-BFGS-B",
      control = list(maxit = 20000, factr = 1e6, pgtol = 1e-12, lmm = 50),
      ...
    )
  } else if (model == 2) {
    # Model 2: Negative alpha allowed with a_par parameterization
    # Setup TMB for negative alpha model
    obj <- TMB::MakeADFun(
      data = list(
        times = times,
        marks = marks,
        penalty_coef = penalty_coef,
        model_type = "hawkes_neg"
      ),
      parameters = list(
        log_mu = log(mu),
        a_par = alpha,
        log_beta = log(beta)
      ),
      hessian = TRUE,
      DLL = "HUO",
      silent = tmb_silent
    )

    # Set bounds for negative alpha model
    bounds <- set_bounds_model2(
      n_mu = 1,
      n_alpha = 1,
      n_beta = 1
    )

    opt <- stats::optim(
      par = obj$par,
      fn = obj$fn,
      gr = obj$gr,
      # lower = bounds$lower,
      # upper = bounds$upper,
      method = "BFGS",
      control = list(maxit = 50000, factr = 1e7, pgtol = 1e-8, lmm = 20),
      ...
    )
  } else {
    stop("model must be 1 (positive alpha) or 2 (negative alpha allowed)")
  }

  # Store results
  obj$objective <- opt$value
  obj$model <- model

  return(obj)
}

# Validation function specifically for 
# multi-series univariate Hawkes with shared parameters
validate_multi_series_shared <- function(times_list, marks_list, 
  model = 1, parameters = list(), min_events = 5) {
  # Basic structure validation
  if (!is.list(times_list) || !is.list(marks_list)) {
    stop("times_list and marks_list must be lists")
  }

  n_series <- length(times_list)
  if (length(marks_list) != n_series) {
    stop("times_list and marks_list must have the same length")
  }

  if (n_series < 2) {
    stop("Need at least 2 series for multi-series fitting")
  }

  # Validate each series using existing validate_params
  validated_series <- list()

  for (i in 1:n_series) {
    if (length(times_list[[i]]) < min_events) {
      stop(paste("Series", i, ": need sufficient events per series"))
    }

    validated_series[[i]] <- validate_params(times_list[[i]], marks_list[[i]], model = model, parameters = list())
  }

  # Initialize shared parameters if not provided
  alpha <- parameters[["alpha"]]
  beta <- parameters[["beta"]]
  mu <- parameters[["mu"]]

  if (is.null(beta)) {
    beta <- mean(sapply(validated_series, function(x) x$beta))
  }

  if (is.null(mu)) {
    mu <- mean(sapply(validated_series, function(x) x$mu))
  }

  # Calculate overall mean marks across all series for shared alpha
  overall_mean_marks <- mean(unlist(marks_list))
  alpha <- 0.3 * beta / overall_mean_marks

  # Validate parameter types
  if (!is.numeric(alpha) || !is.numeric(beta) || !is.numeric(mu)) {
    stop("alpha, beta, and mu must be numeric")
  }

  # Additional parameter validations
  if (beta <= 0) {
    stop("beta must be positive")
  }
  if (mu <= 0) {
    stop("mu must be positive")
  }

  if (model == 1) {
    if (alpha < 0) {
      stop("alpha must be non-negative")
    }

    # Validate alpha/beta constraint
    max_alpha <- beta / overall_mean_marks
    if (alpha >= max_alpha) {
      stop(sprintf(
        "alpha = %.2f must be < beta/mean_marks = %.2f",
        alpha, max_alpha
      ))
    }
  }

  return(list(
    alpha = alpha,
    beta = beta,
    mu = mu,
    overall_mean_marks = overall_mean_marks
  ))
}

#' Fits multiple univariate Hawkes processes using shared parameters
#' across all series by combining their log-likelihoods.
#'
#' @param times_list List of numeric vectors, each containing event times for one series
#' @param marks_list List of numeric vectors, each containing event marks for one series
#' @param parameters List of initial parameter values - single values used for all series:
#'   For model = 1: (mu, alpha, beta)
#'   For model = 2: (mu, a_par, beta)
#' @param model Numeric indicator specifying which model to fit:
#'   \itemize{
#'   \item \code{model = 1}, fits a Hawkes process with positive alpha (default);
#'   \item \code{model = 2}, fits a Hawkes process with alpha that can be negative.
#'   }
#' @param penalty_coef Numeric penalty coefficient for regularization (default: 0)
#' @param tmb_silent Logical, whether to suppress TMB output (default: TRUE)
#' @param ... Additional arguments passed to optimization
#' @return TMB object with fitted shared parameters for all series
#' @export
fit_hawkes_multi_series <- function(times_list,
                                    marks_list,
                                    parameters = list(),
                                    model = 1,
                                    penalty_coef = 0,
                                    tmb_silent = TRUE,
                                    ...) {
  validated_params <- validate_multi_series_shared(times_list, marks_list, 
    model = model, parameters = parameters)

  # Extract validated shared parameters
  alpha <- validated_params$alpha
  beta <- validated_params$beta
  mu <- validated_params$mu
  overall_mean_marks <- validated_params$overall_mean_marks

  if (model == 1) {
    # Model 1: Positive alpha constraint
    n_series <- length(times_list)

    # Concatenate all series data
    times_concat <- unlist(times_list)
    marks_concat <- unlist(marks_list)

    # Create series identifiers (0-indexed for C++)
    series_id <- rep(0:(n_series - 1), sapply(times_list, length))

    # Setup TMB object for positive alpha
    obj <- TMB::MakeADFun(
      data = list(
        times = times_concat,
        marks = marks_concat,
        series_id = series_id, # 0-indexed
        penalty_coef = penalty_coef,
        n_series = n_series,
        model_type = "hawkes_multi"
      ),
      parameters = list(
        log_mu = log(mu),
        logit_abratio = stats::qlogis(alpha / (beta / overall_mean_marks)),
        log_beta = log(beta)
      ),
      hessian = TRUE,
      DLL = "HUO",
      silent = tmb_silent
    )

    # Set bounds for optimization
    bounds <- set_bounds(
      n_mu = 1,
      n_alpha = 1,
      n_beta = 1
    )

    opt <- stats::optim(
      par = obj$par,
      fn = obj$fn,
      gr = obj$gr,
      lower = bounds$lower,
      upper = bounds$upper,
      method = "L-BFGS-B",
      control = list(maxit = 20000, factr = 1e6, pgtol = 1e-12, lmm = 50),
      ...
    )
  } else if (model == 2) {
    # Model 2: Negative alpha allowed
    # Basic validation without alpha constraints
    n_series <- length(times_list)

    # Concatenate all series data
    times_concat <- unlist(times_list)
    marks_concat <- unlist(marks_list)

    # Create series identifiers (0-indexed for C++)
    series_id <- rep(0:(n_series - 1), sapply(times_list, length))

    # Setup TMB object for negative alpha
    obj <- TMB::MakeADFun(
      data = list(
        times = times_concat,
        marks = marks_concat,
        series_id = series_id, # 0-indexed
        penalty_coef = penalty_coef,
        n_series = n_series,
        model_type = "hawkes_multi_neg"
      ),
      parameters = list(
        log_mu = log(mu),
        a_par = alpha,
        log_beta = log(beta)
      ),
      hessian = TRUE,
      DLL = "HUO",
      silent = tmb_silent
    )

    # Set bounds for negative alpha model
    bounds <- set_bounds_model2(
      n_mu = 1,
      n_alpha = 1,
      n_beta = 1
    )

    opt <- stats::optim(
      par = obj$par,
      fn = obj$fn,
      gr = obj$gr,
      # lower = bounds$lower,
      # upper = bounds$upper,
      method = "BFGS",
      control = list(maxit = 50000, factr = 1e7, pgtol = 1e-8, lmm = 50),
      ...
    )
  } else {
    stop("model must be 1 (positive alpha) or 2 (negative alpha allowed)")
  }

  # Store results
  obj$objective <- opt$value
  obj$model <- model

  return(obj)
}