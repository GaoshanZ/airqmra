# ===========================================================================
# 剂量‑反应模型拟合与比较 (支持二项似然、bootstrap、模型选择)
# ===========================================================================
#' Negative log-likelihood for exponential model
#' @keywords internal
exp_nll <- function(k, dose, n_total, n_dead) {
  p <- 1 - exp(-k * dose)
  p <- pmax(pmin(p, 1 - 1e-9), 1e-9)
  -sum(stats::dbinom(n_dead, size = n_total, prob = p, log = TRUE))
}

#' Negative log-likelihood for approximate beta-Poisson model
#' @keywords internal
apb_nll <- function(par, dose, n_total, n_dead) {
  alpha <- exp(par[1])
  N50   <- exp(par[2])
  p <- 1 - (1 + dose * ((2^(1/alpha) - 1) / N50))^(-alpha)
  p <- pmax(pmin(p, 1 - 1e-9), 1e-9)
  -sum(stats::dbinom(n_dead, size = n_total, prob = p, log = TRUE))
}

#' Fit an exponential dose-response model
#'
#' Fits \eqn{P(dose) = 1 - \exp(-k \cdot dose)} using binomial maximum
#' likelihood (when \code{total} is supplied) or least squares.
#'
#' @param dose Numeric vector of dose.
#' @param response Numeric vector of response fractions (0–1) OR dead counts
#'   when \code{total} is supplied.
#' @param total Optional integer vector of total subjects per dose.
#' @param start Starting value for \code{k}. Default \code{1e-6}.
#' @param ... Not used.
#'
#' @return An object of class \code{"dr_fit"}.
#' @export
fit_exponential <- function(dose, response, total = NULL,
                            start = 1e-6, ...) {
  if (!is.null(total)) {
    n_dead <- round(response * total)
    # Use optimize for robust 1D optimization
    opt <- stats::optimize(f = exp_nll, interval = c(1e-12, 1),
                           dose = dose, n_total = total, n_dead = n_dead)
    k_hat <- opt$minimum
    # Hessian approximated by finite difference
    hess <- numDeriv::hessian(exp_nll, k_hat, dose = dose,
                              n_total = total, n_dead = n_dead)
    method_used <- "binomial MLE"
  } else {
    warning("No 'total' supplied; using non‑linear least squares.")
    obj <- function(k) sum((response - (1 - exp(-k * dose)))^2)
    opt <- stats::optimize(f = obj, interval = c(1e-12, 1))
    k_hat <- opt$minimum
    hess <- NULL
    method_used <- "least squares"
  }

  se <- tryCatch(sqrt(solve(hess)[1,1]), error = function(e) NA_real_)
  aic <- if (!is.null(total)) 2 * 1 + 2 * exp_nll(k_hat, dose, total, n_dead) else NA

  structure(
    list(par = c(k = k_hat),
         value = ifelse(is.null(total), opt$objective, exp_nll(k_hat, dose, total, n_dead)),
         convergence = 0L,   # optimize always converges
         hessian = hess,
         se = se,
         dose = dose, response = response, total = total,
         model = "exponential",
         method = method_used,
         aic = aic,
         bootstrap = NULL),
    class = "dr_fit"
  )
}

#' Fit an approximate beta-Poisson dose-response model
#'
#' Fits \deqn{P = 1 - \left(1 + dose \cdot \frac{2^{1/\alpha} - 1}{N_{50}}\right)^{-\alpha}}
#' using binomial maximum likelihood (when \code{total} is supplied) or least squares.
#'
#' @param dose,response,total As in \code{\link{fit_exponential}}.
#' @param start Named list \code{list(alpha = 0.2, N50 = 1000)}.
#' @param method Optimization method passed to \code{optim}. Default
#'   \code{"Nelder-Mead"}.
#' @param ... Additional arguments passed to \code{optim}.
#'
#' @return An object of class \code{"dr_fit"}.
#' @export
fit_approx_beta_poisson <- function(dose, response, total = NULL,
                                    start = list(alpha = 0.2, N50 = 1000),
                                    method = "Nelder-Mead", ...) {
  if (!is.null(total)) {
    n_dead <- round(response * total)
    start_par <- log(c(start$alpha, start$N50))
    opt <- stats::optim(start_par, fn = apb_nll,
                        dose = dose, n_total = total, n_dead = n_dead,
                        method = method, hessian = TRUE, ...)
    method_used <- "binomial MLE"
  } else {
    warning("No 'total' supplied; using non‑linear least squares.")
    obj <- function(par) {
      alpha <- exp(par[1]); N50 <- exp(par[2])
      p <- 1 - (1 + dose * ((2^(1/alpha) - 1) / N50))^(-alpha)
      sum((response - p)^2)
    }
    start_par <- log(c(start$alpha, start$N50))
    opt <- stats::optim(start_par, fn = obj, method = method, hessian = TRUE, ...)
    method_used <- "least squares"
  }

  est <- exp(opt$par); names(est) <- c("alpha", "N50")
  se <- tryCatch(sqrt(diag(solve(opt$hessian))), error = function(e) rep(NA_real_, 2))
  aic <- if (!is.null(total)) 2 * 2 + 2 * opt$value else NA

  structure(
    list(par = est,
         value = opt$value,
         convergence = opt$convergence,
         hessian = opt$hessian,
         se = se,
         dose = dose, response = response, total = total,
         model = "approx_beta_poisson",
         method = method_used,
         aic = aic,
         bootstrap = NULL),
    class = "dr_fit"
  )
}


#' Unified dose‑response model fitting with optional bootstrap
#'
#' Fits either exponential or approximate beta‑Poisson model and, if requested,
#' conducts a parametric bootstrap to obtain parameter percentile confidence
#' intervals.
#'
#' @param dose,response,total,start,model,... As in the individual fit functions.
#' @param bootstrap Logical. If \code{TRUE} (and \code{total} is provided),
#'   parametric bootstrap is performed. Default \code{FALSE}.
#' @param nboot Number of bootstrap replicates. Default 1000.
#'
#' @return An object of class \code{"dr_fit"} with additional element
#'   \code{bootstrap} (a matrix of bootstrap samples) if bootstrapping was
#'   performed.
#' @export
fit_dr_model <- function(dose, response, total = NULL,
                         model = c("exponential", "approx_beta_poisson"),
                         start = NULL, bootstrap = FALSE, nboot = 1000, ...) {
  model <- match.arg(model)

  if (model == "exponential") {
    if (is.null(start)) start <- 1e-6
    fit <- fit_exponential(dose, response, total = total, start = start, ...)
  } else {
    if (is.null(start)) start <- list(alpha = 0.2, N50 = 1000)
    fit <- fit_approx_beta_poisson(dose, response, total = total,
                                   start = start, ...)
  }

  if (bootstrap && !is.null(total)) {
    boot_mat <- parametric_bootstrap(fit, nboot = nboot)
    fit$bootstrap <- boot_mat
    # add bootstrap percentiles
    fit$boot_percentiles <- apply(boot_mat, 2, quantile,
                                  probs = c(0.005, 0.025, 0.05, 0.95, 0.975, 0.995))
  }

  fit
}

#' Parametric bootstrap for a dose‑response model
#'
#' @param fit A \code{"dr_fit"} object.
#' @param nboot Number of bootstrap replicates.
#'
#' @return A matrix with rows corresponding to bootstrap samples and columns
#'   to parameters.
#' @keywords internal
parametric_bootstrap <- function(fit, nboot = 1000) {
  dose <- fit$dose
  total <- fit$total
  if (is.null(total)) stop("Bootstrap requires 'total' to be available.")

  boot_samples <- matrix(NA, nrow = nboot, ncol = length(fit$par))
  colnames(boot_samples) <- names(fit$par)

  for (i in seq_len(nboot)) {
    # Generate new response from fitted model
    if (fit$model == "exponential") {
      p <- 1 - exp(-fit$par["k"] * dose)
    } else {
      alpha <- fit$par["alpha"]; N50 <- fit$par["N50"]
      p <- 1 - (1 + dose * ((2^(1/alpha) - 1) / N50))^(-alpha)
    }
    dead_sim <- stats::rbinom(length(dose), size = total, prob = p)

    # Re‑fit model
    tryCatch({
      if (fit$model == "exponential") {
        boot_fit <- fit_exponential(dose, dead_sim / total, total = total,
                                    start = fit$par["k"], method = "Brent",
                                    lower = 1e-12, upper = 1)
        boot_samples[i, "k"] <- boot_fit$par["k"]
      } else {
        boot_fit <- fit_approx_beta_poisson(dose, dead_sim / total, total = total,
                                            start = list(alpha = fit$par["alpha"],
                                                         N50 = fit$par["N50"]),
                                            method = "Nelder-Mead")
        boot_samples[i, c("alpha","N50")] <- boot_fit$par[c("alpha","N50")]
      }
    }, error = function(e) {
      # leave NA if optimization fails
    })
  }

  boot_samples
}

#' @export
summary.dr_fit <- function(object, ...) {
  cat("Dose‑response model fit (", object$model, ")\n", sep = "")
  cat("Method: ", object$method, "\n")
  cat("-------------------------------------------\n")
  cat("Estimated parameters:\n")
  print(object$par)
  cat("\nStd. errors:\n")
  print(object$se)

  if (!is.null(object$bootstrap)) {
    cat("\nBootstrap percentiles:\n")
    print(object$boot_percentiles)
  }

  if (!is.null(object$aic)) {
    cat("\nAIC:", round(object$aic, 2), "\n")
  }
  cat("Convergence:", object$convergence, "(0 = successful)\n")
  invisible(object)
}

#' @export
plot.dr_fit <- function(x, ...) {
  dose_seq <- 10^seq(log10(min(x$dose) * 0.5), log10(max(x$dose) * 2), length.out = 100)

  if (x$model == "exponential") {
    pred <- 1 - exp(-x$par["k"] * dose_seq)
  } else if (x$model == "approx_beta_poisson") {
    alpha <- x$par["alpha"]; N50 <- x$par["N50"]
    pred <- 1 - (1 + dose_seq * ((2^(1/alpha) - 1) / N50))^(-alpha)
  } else stop("Unknown model type")

  plot(x$dose, x$response, log = "x", xlab = "Dose", ylab = "Response",
       main = paste0("Fitted ", x$model), col = "blue", pch = 16)
  lines(dose_seq, pred, col = "red", lwd = 2)

  # Add confidence bounds if bootstrap available
  if (!is.null(x$bootstrap)) {
    # Predict for each bootstrap sample and take quantiles
    nboot <- nrow(x$bootstrap)
    pred_mat <- matrix(NA, nrow = nboot, ncol = length(dose_seq))
    for (i in seq_len(nboot)) {
      if (x$model == "exponential") {
        k_b <- x$bootstrap[i, "k"]
        if (!is.na(k_b)) pred_mat[i, ] <- 1 - exp(-k_b * dose_seq)
      } else {
        a_b <- x$bootstrap[i, "alpha"]; n_b <- x$bootstrap[i, "N50"]
        if (!all(is.na(c(a_b, n_b))))
          pred_mat[i, ] <- 1 - (1 + dose_seq * ((2^(1/a_b) - 1) / n_b))^(-a_b)
      }
    }
    q_low <- apply(pred_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
    q_high <- apply(pred_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    lines(dose_seq, q_low, lty = 2, col = "red")
    lines(dose_seq, q_high, lty = 2, col = "red")
  }
  invisible(x)
}

#' Compare exponential and approximate beta-Poisson models
#'
#' Performs a likelihood ratio test (when \code{total} is supplied) or
#' compares sums of squares to help choose between the two models.
#'
#' @param dose Numeric vector of dose.
#' @param response Numeric vector of response fractions or dead counts.
#' @param total Optional integer vector of total subjects.
#'
#' @return A list with elements \code{exponential}, \code{beta_poisson},
#'   \code{delta_deviance}, \code{p_value} (if \code{total} supplied),
#'   and a model recommendation.
#' @export
compare_dr_models <- function(dose, response, total = NULL) {
  fit_exp <- fit_dr_model(dose, response, total = total, model = "exponential")
  fit_apb <- fit_dr_model(dose, response, total = total,
                          model = "approx_beta_poisson")

  if (!is.null(total)) {
    # Deviance = 2 * negative log likelihood
    dev_exp <- fit_exp$value * 2
    dev_apb <- fit_apb$value * 2
    delta <- dev_exp - dev_apb
    df <- 1   # approximate beta-Poisson has one extra parameter (N50)
    p_value <- stats::pchisq(delta, df, lower.tail = FALSE)

    cat("Model comparison (likelihood ratio test)\n")
    cat("  Exponential deviance:", round(dev_exp, 2), "\n")
    cat("  Approx beta-Poisson deviance:", round(dev_apb, 2), "\n")
    cat("  Delta deviance:", round(delta, 2), "df =", df,
        "p-value =", format.pval(p_value), "\n")
    if (p_value < 0.05)
      cat("  -> Approx beta-Poisson fits significantly better.\n")
    else
      cat("  -> No significant improvement; exponential is adequate.\n")

    invisible(list(exponential = fit_exp, beta_poisson = fit_apb,
                   delta_deviance = delta, p_value = p_value))
  } else {
    cat("Model comparison (sum of squares)\n")
    ss_exp <- fit_exp$value
    ss_apb <- fit_apb$value
    cat("  Exponential SS:", round(ss_exp, 4),
        "  Approx beta-Poisson SS:", round(ss_apb, 4), "\n")
    invisible(list(exponential = fit_exp, beta_poisson = fit_apb))
  }
}
