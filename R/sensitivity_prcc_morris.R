# =============================================================================
# 模块五：PRCC 与 Morris 敏感性分析
# =============================================================================

# 内部函数：计算 PRCC
.compute_prcc <- function(X, Y) {
  p <- ncol(X)
  X_rank <- apply(X, 2, rank)
  Y_rank <- rank(Y)
  X_std <- scale(X_rank)
  Y_std <- scale(Y_rank)
  cor_mat <- cor(cbind(Y_std, X_std))
  inv_cor <- tryCatch(solve(cor_mat), error = function(e) NULL)
  if (is.null(inv_cor)) return(rep(NA, p))
  prcc <- -inv_cor[1, -1] / sqrt(inv_cor[1,1] * diag(inv_cor[-1, -1]))
  names(prcc) <- colnames(X)
  prcc
}

#' Partial Rank Correlation Coefficient (PRCC) sensitivity analysis
#'
#' Computes PRCC between each input variable and the infection risk using
#' a 1D Monte Carlo sample. PRCC measures the strength of the monotonic
#' relationship while controlling for the effects of other variables.
#' Confidence intervals are obtained via bootstrap.
#'
#' @param scenario An `airqmra_scenario` object (point-estimate mode).
#' @param n_iter Integer. Number of Monte Carlo iterations.
#' @param nboot Integer. Number of bootstrap replicates for confidence
#'   intervals. Set to 0 to skip bootstrap.
#' @param conf Numeric. Confidence level (0–1). Default 0.95.
#'
#' @return A list with elements \code{prcc} (named vector of PRCC values),
#'   \code{ranking} (data frame sorted by absolute PRCC), and optionally
#'   \code{bootstrap} (matrix of bootstrap replicates) and \code{ci}
#'   (confidence intervals).
#' @export
#'
#' @examples
#' scen <- new_scenario(pathogen_id = "legionella pneumophila",
#'                      exposure_group = "Adult", microenvironment = "Outdoor",
#'                      use_uncertainty = FALSE)
#' prcc_res <- sensitivity_prcc(scen, n_iter = 3000, nboot = 200)
#' print(prcc_res$ranking)
#' tornado_plot(prcc_res, title = "PRCC Sensitivity")
sensitivity_prcc <- function(scenario, n_iter = 5000, nboot = 0, conf = 0.95) {
  # Generate Monte Carlo samples
  C   <- scenario$C_air_dist(n_iter)
  IR  <- scenario$IR_dist(n_iter)
  t   <- scenario$t_exp_dist(n_iter)
  surv <- scenario$survival_dist(n_iter)
  depo <- scenario$deposition_dist(n_iter)
  AG   <- scenario$AG_dist(n_iter)

  dose <- C * IR * t * surv * depo * AG

  dr <- scenario$dr_pars
  type <- dr$type
  if (type == "exponential") {
    risk <- 1 - exp(-dr$k * dose)
  } else if (type == "exact_beta_poisson") {
    risk <- exact_beta_poisson(dose, dr$alpha, dr$beta)
  } else if (type == "approx_beta_poisson") {
    risk <- approx_beta_poisson(dose, dr$alpha, dr$N50)
  } else {
    stop("Unknown dose-response type: ", type)
  }

  # Input matrix (exclude AG to avoid constant column)
  X <- cbind(C, IR, t, surv, depo)
  colnames(X) <- c("C", "IR", "t", "survival", "deposition")
  Y <- risk

  prcc_vals <- .compute_prcc(X, Y)

  result <- list(prcc = prcc_vals)

  if (nboot > 0) {
    boot_mat <- matrix(NA, nrow = nboot, ncol = ncol(X))
    colnames(boot_mat) <- colnames(X)
    for (b in seq_len(nboot)) {
      idx <- sample.int(n_iter, replace = TRUE)
      boot_mat[b, ] <- .compute_prcc(X[idx, , drop = FALSE], Y[idx])
    }
    result$bootstrap <- boot_mat
    alpha <- 1 - conf
    ci_low <- apply(boot_mat, 2, quantile, probs = alpha / 2, na.rm = TRUE)
    ci_high <- apply(boot_mat, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
    result$ci <- rbind(ci_low, ci_high)
  }

  # Ranking data frame
  ranking <- data.frame(
    Variable = names(prcc_vals),
    PRCC = prcc_vals,
    AbsPRCC = abs(prcc_vals),
    stringsAsFactors = FALSE
  )
  ranking <- ranking[order(-ranking$AbsPRCC), ]
  rownames(ranking) <- NULL
  result$ranking <- ranking

  class(result) <- c("prcc_sensitivity", "list")
  result
}

#' @export
print.prcc_sensitivity <- function(x, ...) {
  cat("PRCC Sensitivity Analysis\n")
  cat("-------------------------\n")
  print(x$ranking)
  if (!is.null(x$ci)) {
    cat("\nBootstrap confidence intervals:\n")
    print(x$ci)
  }
  invisible(x)
}

#' Morris elementary effects screening (built‑in implementation)
#'
#' Performs Morris sensitivity analysis using a self‑contained implementation
#' of the classic Morris method (Saltelli et al., 2004).  The algorithm does
#' not require the \pkg{sensitivity} package, avoiding potential compatibility
#' issues and ensuring full control over the analysis.
#'
#' @param scenario An `airqmra_scenario` object (point‑estimate mode).
#' @param n_levels Integer. Number of levels in the Morris design. Default 8.
#' @param n_trajectories Integer. Number of trajectories. Default 50.
#'
#' @return An object of class \code{"morris_airqmra"} containing:
#'   \item{ee}{Data frame with rows = factors, columns = \code{mu}, \code{mu.star}, \code{sigma}.}
#'   \item{X}{Design matrix used for evaluation.}
#'   \item{y}{Model output at design points.}
#'   \item{r}{Factor range matrix.}
#' @export
#'
#' @examples
#' \dontrun{
#' scen <- new_scenario(pathogen_id = "legionella pneumophila",
#'                      exposure_group = "Adult", microenvironment = "Outdoor")
#' mor <- sensitivity_morris(scen, n_levels = 8, n_trajectories = 30)
#' plot(mor)
#' tornado_plot(mor, title = "Morris Sensitivity")
#' }
sensitivity_morris <- function(scenario, n_levels = 8, n_trajectories = 50) {
  # ------------------------------------------------------------------
  # 1. Model function
  # ------------------------------------------------------------------
  dr <- scenario$dr_pars
  model_func <- function(X) {
    dose <- X[, 1] * X[, 2] * X[, 3] * X[, 4] * X[, 5]
    if (dr$type == "exponential") {
      1 - exp(-dr$k * dose)
    } else if (dr$type == "exact_beta_poisson") {
      exact_beta_poisson(dose, dr$alpha, dr$beta)
    } else if (dr$type == "approx_beta_poisson") {
      approx_beta_poisson(dose, dr$alpha, dr$N50)
    } else {
      stop("Unknown dose-response type.")
    }
  }

  # ------------------------------------------------------------------
  # 2. Factor ranges (1st and 99th percentiles)
  # ------------------------------------------------------------------
  set.seed(1)
  n_sample <- 10000
  C_range   <- quantile(scenario$C_air_dist(n_sample), c(0.01, 0.99))
  IR_range  <- quantile(scenario$IR_dist(n_sample), c(0.01, 0.99))
  t_range   <- quantile(scenario$t_exp_dist(n_sample), c(0.01, 0.99))
  surv_range <- quantile(scenario$survival_dist(n_sample), c(0.01, 0.99))
  depo_range <- quantile(scenario$deposition_dist(n_sample), c(0.01, 0.99))

  factors <- c("C", "IR", "t", "survival", "deposition")
  p <- length(factors)
  r <- rbind(C_range, IR_range, t_range, surv_range, depo_range)
  rownames(r) <- factors
  colnames(r) <- c("min", "max")

  # ------------------------------------------------------------------
  # 3. Morris sampling matrix construction (classic algorithm)
  # ------------------------------------------------------------------
  delta <- (r[, "max"] - r[, "min"]) / (n_levels - 1)

  ee_mat <- matrix(NA, nrow = n_trajectories, ncol = p)
  colnames(ee_mat) <- factors

  for (t in seq_len(n_trajectories)) {
    start <- runif(p, min = 0, max = 1)
    perm <- sample(p)
    z_current <- start

    for (j in perm) {
      z_pert <- z_current
      z_pert[j] <- z_current[j] + 1 / (n_levels - 1)

      if (z_pert[j] > 1) {
        z_pert[j] <- z_current[j] - 1 / (n_levels - 1)
      }

      base_orig <- z_current * (r[, "max"] - r[, "min"]) + r[, "min"]
      pert_orig <- z_pert   * (r[, "max"] - r[, "min"]) + r[, "min"]

      y_base <- model_func(matrix(base_orig, nrow = 1))
      y_pert <- model_func(matrix(pert_orig, nrow = 1))

      ee_mat[t, j] <- (y_pert - y_base) / delta[j]

      z_current <- z_pert
    }
  }

  # ------------------------------------------------------------------
  # 4. Compute summary statistics
  # ------------------------------------------------------------------
  mu      <- colMeans(ee_mat)
  mu_star <- colMeans(abs(ee_mat))
  sigma   <- apply(ee_mat, 2, sd)

  ee <- data.frame(mu = mu, mu.star = mu_star, sigma = sigma)
  rownames(ee) <- factors

  # ------------------------------------------------------------------
  # 5. Return object
  # ------------------------------------------------------------------
  result <- list(ee = ee, r = r, call = match.call())
  class(result) <- c("morris_airqmra", "list")
  result
}

#' @export
plot.morris_airqmra <- function(x, ...) {
  ee <- x$ee
  mu_star <- ee[, "mu.star"]
  sigma   <- ee[, "sigma"]
  vars    <- rownames(ee)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- data.frame(Variable = vars, mu_star = mu_star, sigma = sigma)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = mu_star, y = sigma, label = Variable)) +
      ggplot2::geom_point() +
      ggplot2::geom_text(ggplot2::aes(vjust = -1)) +
      ggplot2::labs(x = expression(mu*"*"), y = expression(sigma),
                    title = "Morris Elementary Effects") +
      ggplot2::theme_minimal()
    print(p)
  } else {
    plot(mu_star, sigma, xlab = expression(mu*"*"), ylab = expression(sigma),
         main = "Morris Elementary Effects")
    text(mu_star, sigma, labels = vars, pos = 3)
  }
  invisible(x)
}

#' Tornado plot for Morris sensitivity results (with log-scale option)
#'
#' Displays \eqn{\mu^*} (the mean of absolute elementary effects) as a
#' bar chart, providing a ranking of variable importance. Because Morris
#' \eqn{\mu^*} values can span many orders of magnitude, the plot uses a
#' log10 x‑axis by default. Set \code{log_scale = FALSE} to use the original
#' linear scale.
#'
#' @param sens A \code{morris_airqmra} object.
#' @param title Plot title.
#' @param xlab X-axis label (defaults to \eqn{\mu^*}).
#' @param log_scale Logical. If \code{TRUE} (default), the x‑axis is
#'   log10‑transformed.
#'
#' @return A ggplot2 object.
#' @export
tornado_plot.morris_airqmra <- function(sens,
                                        title = "Morris Sensitivity (mu*)",
                                        xlab = expression(mu*"*"),
                                        log_scale = TRUE) {
  ee <- sens$ee
  df <- data.frame(Variable = rownames(ee), mu_star = ee[, "mu.star"])

  p <- ggplot2::ggplot(df, ggplot2::aes(x = reorder(Variable, mu_star), y = mu_star)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = "", y = xlab) +
    ggplot2::theme_minimal()

  if (log_scale) {
    p <- p + ggplot2::scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  }

  return(p)
}

#' @export
print.morris_airqmra <- function(x, ...) {
  cat("Morris Elementary Effects\n")
  cat("-------------------------\n")
  print(x$ee)
  invisible(x)
}

#' Tornado plot for PRCC sensitivity results
#'
#' @param sens A \code{prcc_sensitivity} object.
#' @param title Plot title.
#' @param xlab X-axis label.
#'
#' @return A ggplot2 object.
#' @export
tornado_plot.prcc_sensitivity <- function(sens, title = "PRCC Sensitivity",
                                          xlab = "Partial Rank Correlation") {
  ranking <- sens$ranking

  ci_df <- data.frame(
    Variable = ranking$Variable,
    ci_low   = sens$ci["ci_low", ranking$Variable],
    ci_high  = sens$ci["ci_high", ranking$Variable]
  )
  ci_df$PRCC <- ranking$PRCC[match(ci_df$Variable, ranking$Variable)]

  p <- ggplot2::ggplot(ranking, ggplot2::aes(x = reorder(Variable, PRCC), y = PRCC)) +
    ggplot2::geom_col(ggplot2::aes(fill = PRCC > 0)) +
    ggplot2::geom_errorbar(
      data = ci_df,
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      width = 0.2,
      inherit.aes = TRUE
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = "", y = xlab) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                               guide = "none")
  return(p)
}

# =============================================================================
# 新增：全局敏感性分析（PRCC 和 Morris 的全局版本）
# =============================================================================

#' Global PRCC sensitivity analysis (variability + uncertainty)
#'
#' Calculates partial rank correlation coefficients for all input variables
#' (exposure and dose‑response parameters) using the recorded samples from a
#' 2D Monte Carlo simulation.
#'
#' @param res2d An `airqmra_2d` object created with `record_inputs = TRUE`.
#' @param nboot Integer. Number of bootstrap replicates for confidence
#'   intervals (0 to skip). Default 0.
#' @param conf Numeric. Confidence level (0–1). Default 0.95.
#'
#' @return An object of class \code{"prcc_sensitivity"} with elements
#'   \code{prcc}, \code{ranking}, and (if bootstrap) \code{ci}.
#' @export
sensitivity_global_prcc <- function(res2d, nboot = 0, conf = 0.95) {
  if (is.null(res2d$C_matrix) || is.null(res2d$param_list)) {
    stop("res2d must be created with record_inputs = TRUE.")
  }
  n_var <- nrow(res2d$C_matrix)
  n_unc <- ncol(res2d$C_matrix)

  C_vec <- as.vector(res2d$C_matrix)
  IR_vec <- as.vector(res2d$IR_matrix)
  t_vec <- as.vector(res2d$t_matrix)
  surv_vec <- as.vector(res2d$survival_matrix)
  depo_vec <- as.vector(res2d$deposition_matrix)

  param_df <- do.call(rbind, lapply(res2d$param_list, function(p) data.frame(t(p))))
  param_vec_list <- lapply(param_df, function(x) rep(x, each = n_var))

  X <- cbind(C_vec, IR_vec, t_vec, surv_vec, depo_vec)
  colnames(X) <- c("C", "IR", "t", "survival", "deposition")
  for (nm in names(param_vec_list)) {
    X <- cbind(X, param_vec_list[[nm]])
    colnames(X)[ncol(X)] <- nm
  }
  Y <- res2d$overall_risk

  prcc_vals <- .compute_prcc(X, Y)

  result <- list(prcc = prcc_vals)

  if (nboot > 0) {
    boot_mat <- matrix(NA, nrow = nboot, ncol = ncol(X))
    colnames(boot_mat) <- colnames(X)
    n_total <- nrow(X)
    for (b in seq_len(nboot)) {
      idx <- sample.int(n_total, replace = TRUE)
      boot_mat[b, ] <- .compute_prcc(X[idx, , drop = FALSE], Y[idx])
    }
    alpha <- 1 - conf
    ci_low <- apply(boot_mat, 2, quantile, probs = alpha / 2, na.rm = TRUE)
    ci_high <- apply(boot_mat, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
    result$ci <- rbind(ci_low, ci_high)
  }

  ranking <- data.frame(
    Variable = names(prcc_vals),
    PRCC = prcc_vals,
    AbsPRCC = abs(prcc_vals),
    stringsAsFactors = FALSE
  )
  ranking <- ranking[order(-ranking$AbsPRCC), ]
  rownames(ranking) <- NULL
  result$ranking <- ranking

  class(result) <- c("prcc_sensitivity", "list")
  result
}

#' Global Morris elementary effects screening
#'
#' Performs Morris sensitivity analysis on the full set of input variables
#' (exposure factors and, if the scenario includes uncertainty, dose‑response
#' parameters). Uses a self‑contained implementation.
#'
#' @param scenario An `airqmra_scenario` object. If `use_uncertainty = TRUE`,
#'   dose‑response parameters are added as additional input factors with their
#'   uncertainty ranges.
#' @param n_levels Integer. Number of levels for the Morris design. Default 8.
#' @param n_trajectories Integer. Number of trajectories. Default 50.
#'
#' @return An object of class \code{"morris_airqmra"}.
#' @export
sensitivity_global_morris <- function(scenario, n_levels = 8, n_trajectories = 50) {
  # ---- 1. Identify factors ----
  base_factors <- c("C", "IR", "t", "survival", "deposition")
  param_factors <- character(0)
  param_ranges <- list()

  if (!is.null(scenario$dr_generator)) {
    set.seed(1)
    param_samples <- scenario$dr_generator(10000)
    dr_sample <- scenario$dr_generator(1)
    param_factors <- setdiff(names(dr_sample), "type")
    for (pname in param_factors) {
      vals <- param_samples[[pname]]
      param_ranges[[pname]] <- quantile(vals, c(0.01, 0.99))
    }
  }

  all_factors <- c(base_factors, param_factors)
  p <- length(all_factors)

  # ---- 2. Build factor range matrix ----
  set.seed(1)
  n_sample <- 10000
  C_range   <- quantile(scenario$C_air_dist(n_sample), c(0.01, 0.99))
  IR_range  <- quantile(scenario$IR_dist(n_sample), c(0.01, 0.99))
  t_range   <- quantile(scenario$t_exp_dist(n_sample), c(0.01, 0.99))
  surv_range <- quantile(scenario$survival_dist(n_sample), c(0.01, 0.99))
  depo_range <- quantile(scenario$deposition_dist(n_sample), c(0.01, 0.99))

  r <- rbind(C_range, IR_range, t_range, surv_range, depo_range)
  if (length(param_factors) > 0) {
    for (pname in param_factors) {
      r <- rbind(r, param_ranges[[pname]])
    }
  }
  rownames(r) <- all_factors
  colnames(r) <- c("min", "max")

  # ---- 3. Model function ----
  model_func <- function(X) {
    colnames(X) <- all_factors
    dose <- X[, "C"] * X[, "IR"] * X[, "t"] * X[, "survival"] * X[, "deposition"]
    if (length(param_factors) == 0) {
      dr <- scenario$dr_pars
      if (dr$type == "exponential") {
        1 - exp(-dr$k * dose)
      } else if (dr$type == "exact_beta_poisson") {
        exact_beta_poisson(dose, dr$alpha, dr$beta)
      } else if (dr$type == "approx_beta_poisson") {
        approx_beta_poisson(dose, dr$alpha, dr$N50)
      } else stop("Unknown model")
    } else {
      if ("k" %in% all_factors) {
        k_vec <- X[, "k"]
        1 - exp(-k_vec * dose)
      } else if ("alpha" %in% all_factors && "beta" %in% all_factors) {
        alpha_vec <- X[, "alpha"]
        beta_vec <- X[, "beta"]
        exact_beta_poisson(dose, alpha_vec, beta_vec)
      } else if ("alpha" %in% all_factors && "N50" %in% all_factors) {
        alpha_vec <- X[, "alpha"]
        N50_vec <- X[, "N50"]
        approx_beta_poisson(dose, alpha_vec, N50_vec)
      } else {
        stop("Unknown parameter combination")
      }
    }
  }

  # ---- 4. Morris core ----
  delta <- (r[, "max"] - r[, "min"]) / (n_levels - 1)
  ee_mat <- matrix(NA, nrow = n_trajectories, ncol = p)
  colnames(ee_mat) <- all_factors

  for (t in seq_len(n_trajectories)) {
    start <- runif(p, min = 0, max = 1)
    perm <- sample(p)
    z_current <- start

    for (j in perm) {
      z_pert <- z_current
      z_pert[j] <- z_current[j] + 1 / (n_levels - 1)
      if (z_pert[j] > 1) z_pert[j] <- z_current[j] - 1 / (n_levels - 1)

      base_orig <- z_current * (r[, "max"] - r[, "min"]) + r[, "min"]
      pert_orig <- z_pert   * (r[, "max"] - r[, "min"]) + r[, "min"]

      y_base <- model_func(matrix(base_orig, nrow = 1))
      y_pert <- model_func(matrix(pert_orig, nrow = 1))

      ee_mat[t, j] <- (y_pert - y_base) / delta[j]
      z_current <- z_pert
    }
  }

  mu      <- colMeans(ee_mat)
  mu_star <- colMeans(abs(ee_mat))
  sigma   <- apply(ee_mat, 2, sd)
  ee <- data.frame(mu = mu, mu.star = mu_star, sigma = sigma)
  rownames(ee) <- all_factors

  result <- list(ee = ee, r = r, call = match.call())
  class(result) <- c("morris_airqmra", "list")
  result
}
