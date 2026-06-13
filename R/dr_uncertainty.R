#' Build a dose-response parameter generator from database record
#'
#' @param rec A single-row data frame from \code{pathogen_db}.
#' @param use_uncertainty Logical. If TRUE, parameters are sampled from their
#'   uncertainty distributions; if FALSE, point estimates are used.
#'
#' @return A function that takes an integer \code{n} and returns a list
#'   with elements \code{type} and the sampled/constant parameter vectors.
#' @keywords internal
build_dr_generator <- function(rec, use_uncertainty = FALSE) {
  model <- rec$model

  # 1. Exponential model -------------------------------------------------------
  if (model == "exponential") {
    if (use_uncertainty && !is.na(rec$k_log_mean) && !is.na(rec$k_log_sd)) {
      f <- function(n) {
        k <- rlnorm(n, meanlog = rec$k_log_mean, sdlog = rec$k_log_sd)
        list(type = "exponential", k = k)
      }
    } else {
      k_fixed <- rec$k
      f <- function(n) {
        list(type = "exponential", k = rep(k_fixed, n))
      }
    }
    return(f)
  }

  # 2. Approximate beta-Poisson (N50-based) – strictly by model label --------
  if (model == "approx_beta_poisson") {
    if (is.na(rec$alpha) || is.na(rec$N50))
      stop("approx_beta_poisson requires alpha and N50 in the database.")

    alpha_unc <- use_uncertainty &&
      !is.na(rec$alpha_log_mean) && !is.na(rec$alpha_log_sd)
    n50_unc   <- use_uncertainty &&
      !is.na(rec$N50_log_mean) && !is.na(rec$N50_log_sd)

    if (use_uncertainty && (alpha_unc || n50_unc)) {
      f <- function(n) {
        alpha <- if (alpha_unc) rlnorm(n, rec$alpha_log_mean, rec$alpha_log_sd)
        else rep(rec$alpha, n)
        N50   <- if (n50_unc)   rlnorm(n, rec$N50_log_mean,   rec$N50_log_sd)
        else rep(rec$N50, n)
        list(type = "approx_beta_poisson", alpha = alpha, N50 = N50)
      }
    } else {
      f <- function(n) list(type = "approx_beta_poisson",
                            alpha = rep(rec$alpha, n),
                            N50   = rep(rec$N50, n))
    }
    return(f)
  }

  # 3. Exact beta-Poisson (beta-based) – strictly by model label --------------
  if (model %in% c("exact_beta_poisson", "beta_poisson")) {
    if (is.na(rec$alpha) || is.na(rec$beta))
      stop("exact_beta_poisson requires alpha and beta in the database.")

    alpha_unc <- use_uncertainty &&
      !is.na(rec$alpha_log_mean) && !is.na(rec$alpha_log_sd)
    beta_unc  <- use_uncertainty &&
      !is.na(rec$beta_log_mean) && !is.na(rec$beta_log_sd)

    if (use_uncertainty && (alpha_unc || beta_unc)) {
      f <- function(n) {
        alpha <- if (alpha_unc) rlnorm(n, rec$alpha_log_mean, rec$alpha_log_sd)
        else rep(rec$alpha, n)
        beta  <- if (beta_unc)  rlnorm(n, rec$beta_log_mean,  rec$beta_log_sd)
        else rep(rec$beta, n)
        list(type = "exact_beta_poisson", alpha = alpha, beta = beta)
      }
    } else {
      f <- function(n) list(type = "exact_beta_poisson",
                            alpha = rep(rec$alpha, n),
                            beta  = rep(rec$beta, n))
    }
    return(f)
  }

  stop("Unsupported dose-response model: ", model)
}
