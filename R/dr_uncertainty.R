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

  # 2. Beta-Poisson models (exact or approximate) ------------------------------
  if (model %in% c("exact_beta_poisson", "beta_poisson")) {

    # Helper: decide whether to use uncertainty for alpha and N50/beta
    alpha_unc <- use_uncertainty &&
      !is.na(rec$alpha_log_mean) && !is.na(rec$alpha_log_sd)
    n50_unc   <- use_uncertainty &&
      !is.na(rec$N50_log_mean) && !is.na(rec$N50_log_sd)
    beta_unc  <- use_uncertainty &&
      !is.na(rec$beta_log_mean) && !is.na(rec$beta_log_sd)

    # --- Case A: We have beta (directly or will compute from N50) ---
    # Check if a direct beta is available (point or uncertainty)
    have_beta <- !is.na(rec$beta) || (use_uncertainty && beta_unc)

    if (have_beta) {
      # Use exact beta-Poisson
      if (use_uncertainty && (alpha_unc || beta_unc)) {
        f <- function(n) {
          if (alpha_unc) {
            alpha <- rlnorm(n, meanlog = rec$alpha_log_mean, sdlog = rec$alpha_log_sd)
          } else {
            alpha <- rep(rec$alpha, n)
          }
          if (beta_unc) {
            beta <- rlnorm(n, meanlog = rec$beta_log_mean, sdlog = rec$beta_log_sd)
          } else {
            # If beta not uncertain, use fixed beta or compute from N50 if needed
            if (!is.na(rec$beta)) {
              beta <- rep(rec$beta, n)
            } else if (!is.na(rec$N50)) {
              # Compute beta from fixed alpha and N50 (should not happen if have_beta is TRUE)
              beta <- rec$N50 / (2^(1/alpha) - 1)
            } else {
              stop("No beta or N50 available.")
            }
          }
          list(type = "exact_beta_poisson", alpha = alpha, beta = beta)
        }
      } else {
        # Point estimates
        alpha_fixed <- rec$alpha
        beta_fixed  <- rec$beta
        if (is.na(beta_fixed) && !is.na(rec$N50)) {
          beta_fixed <- rec$N50 / (2^(1/alpha_fixed) - 1)
        }
        f <- function(n) {
          list(type = "exact_beta_poisson",
               alpha = rep(alpha_fixed, n),
               beta  = rep(beta_fixed, n))
        }
      }
      return(f)
    }

    # --- Case B: Only N50 is available, no beta -> approximate beta-Poisson ---
    if (!is.na(rec$N50)) {
      if (use_uncertainty && (alpha_unc || n50_unc)) {
        f <- function(n) {
          if (alpha_unc) {
            alpha <- rlnorm(n, meanlog = rec$alpha_log_mean, sdlog = rec$alpha_log_sd)
          } else {
            alpha <- rep(rec$alpha, n)
          }
          if (n50_unc) {
            N50 <- rlnorm(n, meanlog = rec$N50_log_mean, sdlog = rec$N50_log_sd)
          } else {
            N50 <- rep(rec$N50, n)
          }
          list(type = "approx_beta_poisson", alpha = alpha, N50 = N50)
        }
      } else {
        alpha_fixed <- rec$alpha
        N50_fixed <- rec$N50
        f <- function(n) {
          list(type = "approx_beta_poisson",
               alpha = rep(alpha_fixed, n),
               N50   = rep(N50_fixed, n))
        }
      }
      return(f)
    }

    # If neither beta nor N50 available
    stop("Beta-Poisson model requires either 'beta' or 'N50' in the database.")
  }

  stop("Unsupported dose-response model: ", model)
}
