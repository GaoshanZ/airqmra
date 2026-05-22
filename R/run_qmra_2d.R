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
#' @param record_inputs Logical. If `TRUE`, exposure variability matrices (C, IR, t,
#'   survival, deposition) are also recorded and returned (for global sensitivity
#'   analysis). Default `FALSE`.
#'
#' @return An object of class `airqmra_2d`, which is a list containing:
#'   \item{risk_matrix}{Matrix of single-exposure risks (rows = variability, cols = uncertainty).}
#'   \item{overall_risk}{Numeric vector combining all variability and uncertainty samples.}
#'   \item{risk_quantiles}{Matrix of risk quantiles (rows = probs, columns = uncertainty iterations).}
#'   \item{mean_risk}{Numeric vector of mean risk per uncertainty iteration.}
#'   \item{summary}{Data frame with overall mean, median, and confidence limits of the mean risk.}
#'   \item{param_list}{List (length = n_unc) of named vectors containing the dose‑response
#'         parameter values actually used in each uncertainty iteration. Always returned.}
#'   \item{annual_quantiles, annual_mean, annual_summary, overall_annual_risk}{If `frequency` is provided.}
#'   \item{C_matrix, IR_matrix, t_matrix, survival_matrix, deposition_matrix}{If `record_inputs = TRUE`,
#'         matrices of exposure inputs (rows = variability, cols = uncertainty).}
#'   \item{annual_matrix}{If \code{frequency} is provided, matrix of annual risks (rows = variability, cols = uncertainty).}
#' @export
#'
#' @examples
#' # Legionella with uncertainty
#' scen <- new_scenario(pathogen_id = "legionella pneumophila",
#'                      use_uncertainty = TRUE)
#' res2d <- run_qmra_2d(scen, n_var = 2000, n_unc = 500, frequency = 240)
#' summary(res2d)
#' hist(res2d$overall_risk, main="Combined variability & uncertainty")
#'
#' # With global sensitivity inputs recorded
#' res2d_global <- run_qmra_2d(scen, n_var = 500, n_unc = 200,
#'                             frequency = NULL, record_inputs = TRUE)
run_qmra_2d <- function(scenario, n_var = 5000, n_unc = 1000,
                        probs = c(0.025, 0.5, 0.975),
                        frequency = NULL,
                        record_inputs = FALSE) {

  if (is.null(scenario$dr_generator)) {
    stop("Scenario must be created with 'use_uncertainty = TRUE' to enable 2D simulation.")
  }

  # 预分配输出矩阵
  risk_matrix <- matrix(NA_real_, nrow = n_var, ncol = n_unc)
  mean_risk <- numeric(n_unc)

  if (!is.null(frequency)) {
    annual_matrix <- matrix(NA_real_, nrow = n_var, ncol = n_unc)
    annual_mean <- numeric(n_unc)
  }

  # 如果需要记录输入，预分配存储
  if (record_inputs) {
    C_matrix <- matrix(NA_real_, nrow = n_var, ncol = n_unc)
    IR_matrix <- matrix(NA_real_, nrow = n_var, ncol = n_unc)
    t_matrix <- matrix(NA_real_, nrow = n_var, ncol = n_unc)
    survival_matrix <- matrix(NA_real_, nrow = n_var, ncol = n_unc)
    deposition_matrix <- matrix(NA_real_, nrow = n_var, ncol = n_unc)
  }

  # 始终记录参数值（无论 record_inputs 如何）
  param_list <- vector("list", n_unc)

  for (i in seq_len(n_unc)) {
    # 获取一组固定的剂量‑反应参数（标量）
    dr_params <- scenario$dr_generator(1)
    type <- dr_params$type
    dr_args_fixed <- dr_params[names(dr_params) != "type"]

    # 保存参数值
    if (type == "exponential") {
      param_list[[i]] <- c(k = dr_args_fixed$k)
    } else if (type == "exact_beta_poisson") {
      param_list[[i]] <- c(alpha = dr_args_fixed$alpha, beta = dr_args_fixed$beta)
    } else if (type == "approx_beta_poisson") {
      param_list[[i]] <- c(alpha = dr_args_fixed$alpha, N50 = dr_args_fixed$N50)
    }

    # 将参数固定到场景副本中
    scen_fixed <- scenario
    scen_fixed$dr_pars <- c(list(type = type), dr_args_fixed)
    scen_fixed$dr_generator <- NULL

    # 生成当前迭代的所有暴露变异性样本
    C <- scen_fixed$C_air_dist(n_var)
    IR <- scen_fixed$IR_dist(n_var)
    t <- scen_fixed$t_exp_dist(n_var)
    surv <- scen_fixed$survival_dist(n_var)
    depo <- scen_fixed$deposition_dist(n_var)
    AG <- scen_fixed$AG_dist(n_var)

    # 计算有效剂量和风险
    dose <- C * IR * t * surv * depo * AG
    single_risk <- do.call(dose_response, c(list(dose = dose, type = type), dr_args_fixed))

    risk_matrix[, i] <- single_risk
    mean_risk[i] <- mean(single_risk)

    # 记录年风险（如需要）
    if (!is.null(frequency)) {
      annual_risk_vec <- annual_risk(single_risk, frequency = frequency)
      annual_matrix[, i] <- annual_risk_vec
      annual_mean[i] <- mean(annual_risk_vec)
    }

    # 记录输入变量（仅当 record_inputs = TRUE）
    if (record_inputs) {
      C_matrix[, i] <- C
      IR_matrix[, i] <- IR
      t_matrix[, i] <- t
      survival_matrix[, i] <- surv
      deposition_matrix[, i] <- depo
    }
  }

  # 合并变异性和不确定性样本
  overall_risk <- as.vector(risk_matrix)

  risk_quantiles <- apply(risk_matrix, 2, quantile, probs = probs, na.rm = TRUE)
  overall_mean <- mean(mean_risk)
  overall_ci <- quantile(mean_risk, probs = c(probs[1], probs[length(probs)]), na.rm = TRUE)

  res <- list(
    risk_matrix    = risk_matrix,
    overall_risk   = overall_risk,
    risk_quantiles = risk_quantiles,
    mean_risk      = mean_risk,
    summary        = data.frame(
      Mean   = overall_mean,
      Median = median(mean_risk),
      CI_low = overall_ci[1],
      CI_high = overall_ci[2]
    ),
    param_list     = param_list   # 始终返回
  )

  if (!is.null(frequency)) {
    annual_quantiles <- apply(annual_matrix, 2, quantile, probs = probs, na.rm = TRUE)
    annual_overall_mean <- mean(annual_mean)
    annual_overall_ci <- quantile(annual_mean, probs = c(probs[1], probs[length(probs)]), na.rm = TRUE)
    overall_annual_risk <- as.vector(annual_matrix)
    res$annual_matrix <- annual_matrix
    res$annual_quantiles <- annual_quantiles
    res$annual_mean <- annual_mean
    res$annual_summary <- data.frame(
      Mean   = annual_overall_mean,
      Median = median(annual_mean),
      CI_low = annual_overall_ci[1],
      CI_high = annual_overall_ci[2]
    )
    res$overall_annual_risk <- overall_annual_risk
  }

  # 附加记录数据（如果开启）
  if (record_inputs) {
    res$C_matrix <- C_matrix
    res$IR_matrix <- IR_matrix
    res$t_matrix <- t_matrix
    res$survival_matrix <- survival_matrix
    res$deposition_matrix <- deposition_matrix
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

  cat("Overall single-exposure risk (combined variability & uncertainty):\n")
  print(summary(object$overall_risk))
  cat("\n")

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
  hist(x$overall_risk, main = "Combined Variability & Uncertainty",
       xlab = "Single-exposure infection risk", col = "lightblue")
}
