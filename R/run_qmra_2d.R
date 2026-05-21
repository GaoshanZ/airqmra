#' Run 2D Monte Carlo QMRA (variability + uncertainty)
#'
#' Performs a two-dimensional Monte Carlo simulation:
#' \itemize{
#'   \item Outer loop (uncertainty iterations): samples dose-response parameters
#'         from their uncertainty distributions.
#'   \item Inner loop (variability iterations): for each parameter set, simulates
#'         exposure variability to produce a risk distribution.
#' }
#'
#' @param scenario An `airqmra_scenario` object. Must contain a `dr_generator`
#'   function (created by `new_scenario(..., use_uncertainty = TRUE)`).
#' @param n_var Number of variability iterations (inner loop). Default 5000.
#' @param n_unc Number of uncertainty iterations (outer loop). Default 1000.
#' @param probs Numeric vector of probabilities for quantile calculation.
#'   Default `c(0.025, 0.5, 0.975)`.
#' @param frequency Optional. If provided, annual risk is also computed using
#'   this frequency (passed to `annual_risk()`). Can be a fixed integer or a function.
#'
#' @return An object of class `airqmra_2d`, which is a list containing:
#'   \item{risk_matrix}{Matrix of single-exposure risks (rows = variability, cols = uncertainty).}
#'   \item{risk_quantiles}{Matrix of risk quantiles (rows = probs, columns = uncertainty iterations).}
#'   \item{mean_risk}{Numeric vector of mean risk per uncertainty iteration.}
#'   \item{summary}{Data frame with overall mean, median, and confidence limits of the mean risk.}
#'   \item{annual_quantiles}{If `frequency` is provided, annual risk quantiles matrix.}
#'   \item{annual_mean}{If `frequency` is provided, vector of annual mean risk per iteration.}
#'   \item{annual_summary}{If `frequency` is provided, data frame with annual mean, median, CI.}
#' @export
#'
#' @examples
#' # Legionella with uncertainty
#' scen <- new_scenario(pathogen_id = "legionella pneumophila",
#'                      use_uncertainty = TRUE)
#' res2d <- run_qmra_2d(scen, n_var = 2000, n_unc = 500, frequency = 240)
#' summary(res2d)
#' plot(res2d)
run_qmra_2d <- function(scenario, n_var = 5000, n_unc = 1000,
                        probs = c(0.025, 0.5, 0.975),
                        frequency = NULL) {

  if (is.null(scenario$dr_generator)) {
    stop("Scenario must be created with 'use_uncertainty = TRUE' to enable 2D simulation.")
  }

  # 预分配矩阵
  risk_matrix <- matrix(NA_real_, nrow = n_var, ncol = n_unc)
  mean_risk <- numeric(n_unc)

  if (!is.null(frequency)) {
    annual_matrix <- matrix(NA_real_, nrow = n_var, ncol = n_unc)
    annual_mean <- numeric(n_unc)
  }

  for (i in seq_len(n_unc)) {
    # 从生成器获取一组固定的剂量-反应参数（标量）
    dr_params <- scenario$dr_generator(1)
    type <- dr_params$type
    dr_args_fixed <- dr_params[names(dr_params) != "type"]

    # 构造一个临时的固定参数场景
    scen_fixed <- scenario
    scen_fixed$dr_pars <- c(list(type = type), dr_args_fixed)
    scen_fixed$dr_generator <- NULL

    # 变异性模拟
    single_risk <- run_qmra_1d(scen_fixed, n_iter = n_var)
    risk_matrix[, i] <- single_risk
    mean_risk[i] <- mean(single_risk)

    if (!is.null(frequency)) {
      annual_risk_vec <- annual_risk(single_risk, frequency = frequency)
      annual_matrix[, i] <- annual_risk_vec
      annual_mean[i] <- mean(annual_risk_vec)
    }
  }

  # 计算风险分位数（每列一个不确定性样本）
  risk_quantiles <- apply(risk_matrix, 2, quantile, probs = probs, na.rm = TRUE)

  overall_mean <- mean(mean_risk)
  overall_ci <- quantile(mean_risk, probs = c(probs[1], probs[length(probs)]), na.rm = TRUE)

  res <- list(
    risk_matrix    = risk_matrix,
    risk_quantiles = risk_quantiles,
    mean_risk      = mean_risk,
    summary        = data.frame(
      Mean   = overall_mean,
      Median = median(mean_risk),
      CI_low = overall_ci[1],
      CI_high = overall_ci[2]
    )
  )

  if (!is.null(frequency)) {
    annual_quantiles <- apply(annual_matrix, 2, quantile, probs = probs, na.rm = TRUE)
    annual_overall_mean <- mean(annual_mean)
    annual_overall_ci <- quantile(annual_mean, probs = c(probs[1], probs[length(probs)]), na.rm = TRUE)
    res$annual_quantiles <- annual_quantiles
    res$annual_mean <- annual_mean
    res$annual_summary <- data.frame(
      Mean   = annual_overall_mean,
      Median = median(annual_mean),
      CI_low = annual_overall_ci[1],
      CI_high = annual_overall_ci[2]
    )
  }

  class(res) <- c("airqmra_2d", "list")
  res
}

#' @export
summary.airqmra_2d <- function(object, ...) {
  cat("2D Monte Carlo QMRA Results\n")
  cat("---------------------------\n")
  cat("Uncertainty iterations:", ncol(object$risk_matrix), "\n")
  cat("Variability iterations per U sample:", nrow(object$risk_matrix), "\n\n")

  cat("Single-Exposure Risk (mean across variability):\n")
  print(object$summary, row.names = FALSE)

  if (!is.null(object$annual_summary)) {
    cat("\nAnnual Risk (mean across variability):\n")
    print(object$annual_summary, row.names = FALSE)
  }

  invisible(object)
}

#' @export
plot.airqmra_2d <- function(x, ...) {
  hist(x$mean_risk, main = "Uncertainty in Mean Risk",
       xlab = "Mean infection risk per exposure", col = "lightblue")
}
