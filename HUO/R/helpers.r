#' Extract parameter estimates from TMB sdreport
#' @param obj TMB model object with fitted parameters
#' @return Matrix of parameter estimates and standard errors
#' @export
get_coefs <- function(obj) {
  table <- summary(TMB::sdreport(obj))
  return(table)
}

#' Extract parameters from univariate Hawkes TMB object
#' @param obj TMB model object from fit_hawkes or fit_hawkes_multi_series
#' @return List with mu, alpha, and beta parameter estimates
#' @export
extract_univariate_params <- function(obj) {
  # Use TMB::sdreport to get the transformed parameters
  sdr <- TMB::sdreport(obj)
  params <- summary(sdr, "report")

  # Extract the reported parameters (already transformed by C++)
  mu <- params[rownames(params) == "mu", "Estimate"]
  alpha <- params[rownames(params) == "alpha", "Estimate"]
  beta <- params[rownames(params) == "beta", "Estimate"]

  return(list(mu = mu, alpha = alpha, beta = beta))
}


#' Calculate univariate Hawkes intensity
#' @param times Numeric vector of observed event times
#' @param mu Baseline intensity rate (scalar)
#' @param alpha Excitation parameter (scalar)
#' @param beta Decay parameter (scalar)
#' @param p Numeric vector of time points at which to evaluate intensity (default: times)
#' @param marks Numeric vector of marks (default: NULL for unmarked process)
#' @return Numeric vector of intensity values at specified time points
#' @export
hawkes_intensity <- function(times, mu, alpha, beta, p = times, marks = NULL) {
  if (is.null(marks)) marks <- rep(1, length(times))

  # Efficient vectorized calculation for each evaluation point
  lam <- function(p_val, mu_val) {
    # Vectorized decay calculation for all past events
    decay_factors <- exp(-beta * (p_val - times))[times < p_val]
    past_marks <- marks[times < p_val]
    mu_val + alpha * sum(past_marks * decay_factors)
  }

  # Calculate intensity: λ(t) = μ + α∑_{t_j<t}m_j*e^{-β(t-t_j)}
  mus <- rep(mu, length(p))
  intensity <- numeric(length(p))

  for (i in seq_along(p)) {
    intensity[i] <- lam(p[i], mus[i])
  }

  return(intensity)
}
