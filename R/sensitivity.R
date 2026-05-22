#' Compute Spearman correlation sensitivity for a set of inputs and an output
#'
#' @param data A data frame containing input variables and the output variable.
#' @param output_col Character. Name of the output column in `data`.
#' @param method Character. Correlation method, default `"spearman"`.
#'
#' @return A list with elements:
#'   \item{correlations}{Named vector of correlation coefficients between each
#'   input and the output.}
#'   \item{ranking}{Data frame with inputs sorted by absolute correlation.}
#' @export
#'
#' @examples
#' df <- data.frame(x = rnorm(100), y = rnorm(100), z = rnorm(100))
#' df$result <- 2*df$x + 0.5*df$y
#' sens <- compute_sensitivity(df, "result")
#' print(sens$ranking)
compute_sensitivity <- function(data, output_col, method = "spearman") {
  if (!output_col %in% names(data)) {
    stop("Output column '", output_col, "' not found in data.")
  }
  output <- data[[output_col]]
  inputs <- data[, setdiff(names(data), output_col), drop = FALSE]

  cors <- sapply(inputs, function(x) cor(x, output, method = method, use = "complete.obs"))

  ranking <- data.frame(
    Variable = names(cors),
    Correlation = cors,
    AbsCorrelation = abs(cors)
  )
  ranking <- ranking[order(-ranking$AbsCorrelation), ]
  rownames(ranking) <- NULL

  list(correlations = cors, ranking = ranking)
}


#' Sensitivity analysis for a 1D (point‑estimate) QMRA scenario
#'
#' Runs a 1D Monte Carlo simulation, records all input random samples and the
#' resulting risk, and computes Spearman rank correlation coefficients.
#' This function is intended for point‑estimate scenarios (`use_uncertainty = FALSE`).
#' **Note:** AG is excluded from the sensitivity analysis because it is fixed (no variability).
#'
#' @param scenario An `airqmra_scenario` object (from `new_scenario`).
#' @param n_iter Number of Monte Carlo iterations.
#'
#' @return A list with:
#'   \item{sensitivity}{Output from \code{\link{compute_sensitivity}}.}
#'   \item{data}{Data frame containing all sampled inputs and the risk.}
#' @export
sensitivity_1d <- function(scenario, n_iter = 5000) {
  # 生成所有输入变量的样本
  C <- scenario$C_air_dist(n_iter)
  IR <- scenario$IR_dist(n_iter)
  t <- scenario$t_exp_dist(n_iter)
  survival <- scenario$survival_dist(n_iter)
  deposition <- scenario$deposition_dist(n_iter)
  # AG 固定为 1，虽然乘入剂量，但不作为敏感性变量
  AG <- scenario$AG_dist(n_iter)

  dose <- C * IR * t * survival * deposition * AG   # 剂量计算保持不变

  dr <- scenario$dr_pars
  type <- dr$type
  if (type == "exponential") {
    k <- dr$k
    risk <- 1 - exp(-k * dose)
  } else if (type == "exact_beta_poisson") {
    alpha <- dr$alpha
    beta <- dr$beta
    risk <- exact_beta_poisson(dose, alpha, beta)
  } else if (type == "approx_beta_poisson") {
    alpha <- dr$alpha
    N50 <- dr$N50
    risk <- approx_beta_poisson(dose, alpha, N50)
  } else {
    stop("Unknown dose-response type: ", type)
  }

  # 构建数据框：不包含 AG
  df <- data.frame(
    C = C,
    IR = IR,
    t = t,
    survival = survival,
    deposition = deposition,
    risk = risk
  )

  sens <- compute_sensitivity(df, output_col = "risk")
  list(sensitivity = sens, data = df)
}


#' Sensitivity analysis for 2D Monte Carlo results (parameter uncertainty)
#'
#' Calculates Spearman correlation between dose-response parameter values
#' (e.g., k) and the mean single-exposure risk obtained from a 2D simulation.
#'
#' @param param_values Numeric vector of parameter values across uncertainty iterations.
#' @param mean_risk Numeric vector of mean risk (or other statistic) from each uncertainty iteration.
#' @param param_name Character. Name of the parameter (e.g., "k").
#'
#' @return A list with elements \code{correlations} and \code{ranking}.
#' @export
sensitivity_2d <- function(param_values, mean_risk, param_name) {
  df <- data.frame(param = param_values, mean_risk = mean_risk)
  names(df)[1] <- param_name
  compute_sensitivity(df, output_col = "mean_risk")
}


#' Global sensitivity analysis combining variability and uncertainty
#'
#' Uses the recorded inputs and overall risk from a 2D simulation to compute
#' Spearman correlations for all input variables (exposure and dose-response parameters).
#'
#' @param res2d An object of class `airqmra_2d` that was created with
#'   `run_qmra_2d(..., record_inputs = TRUE)`.
#' @param method Correlation method, default `"spearman"`.
#'
#' @return A list as from \code{\link{compute_sensitivity}} with correlations and ranking.
#' @export
#'
#' @examples
#' res2d <- run_qmra_2d(scen_unc, n_var = 500, n_unc = 200,
#'                      frequency = NULL, record_inputs = TRUE)
#' sens_global <- sensitivity_global(res2d)
#' tornado_plot(sens_global, title = "Global Sensitivity")
sensitivity_global <- function(res2d, method = "spearman") {
  if (is.null(res2d$C_matrix) || is.null(res2d$param_list)) {
    stop("The 2D simulation must be run with 'record_inputs = TRUE' to use global sensitivity.")
  }

  n_var <- nrow(res2d$C_matrix)
  n_unc <- ncol(res2d$C_matrix)

  # 将变异性输入矩阵拉直为向量
  C_vec <- as.vector(res2d$C_matrix)
  IR_vec <- as.vector(res2d$IR_matrix)
  t_vec <- as.vector(res2d$t_matrix)
  survival_vec <- as.vector(res2d$survival_matrix)
  deposition_vec <- as.vector(res2d$deposition_matrix)

  # 参数向量：匹配长度 n_var * n_unc
  param_df <- do.call(rbind, lapply(res2d$param_list, function(p) {
    data.frame(t(p))
  }))
  # 将每一列复制 n_var 次，以匹配变异性样本
  param_vec_list <- lapply(param_df, function(x) {
    rep(x, each = n_var)
  })

  # 整体风险向量
  overall_risk <- res2d$overall_risk

  # 构造数据框
  df <- data.frame(
    C = C_vec,
    IR = IR_vec,
    t = t_vec,
    survival = survival_vec,
    deposition = deposition_vec
  )
  # 添加上参数列
  for (nm in names(param_vec_list)) {
    df[[nm]] <- param_vec_list[[nm]]
  }
  df$risk <- overall_risk

  compute_sensitivity(df, output_col = "risk", method = method)
}


#' Create a tornado (sensitivity) plot
#'
#' @param sens A sensitivity object from \code{\link{compute_sensitivity}},
#'   \code{\link{sensitivity_1d}}, \code{\link{sensitivity_2d}}, or
#'   \code{\link{sensitivity_global}}.
#' @param title Character. Plot title.
#' @param xlab Character. Label for x-axis.
#'
#' @return A ggplot object (if ggplot2 is installed) or a base R barplot.
#' @export
tornado_plot <- function(sens, title = "Sensitivity Analysis", xlab = "Spearman Correlation") {
  # 提取 ranking 数据框
  if ("sensitivity" %in% names(sens)) {
    ranking <- sens$sensitivity$ranking
  } else if ("ranking" %in% names(sens)) {
    ranking <- sens$ranking
  } else {
    stop("Input must contain a 'ranking' data frame.")
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(ranking, ggplot2::aes(x = reorder(Variable, Correlation),
                                               y = Correlation, fill = Correlation > 0)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(title = title, x = "", y = xlab) +
      ggplot2::theme_minimal() +
      ggplot2::scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"), guide = "none")
    return(p)
  } else {
    # 基础绘图
    barplot(ranking$Correlation, names.arg = ranking$Variable, horiz = TRUE,
            col = ifelse(ranking$Correlation > 0, "steelblue", "firebrick"),
            main = title, xlab = xlab)
  }
}



#' Sensitivity analysis for 2D Monte Carlo results (convenience wrapper)
#'
#' @param res2d An object of class `airqmra_2d`.
#' @return A list with correlations and ranking for the primary dose‑response parameter.
#' @export
sensitivity_2d_auto <- function(res2d) {
  if (is.null(res2d$param_list)) stop("res2d must contain param_list (always present).")
  # 提取第一个参数（对于指数模型是 k，贝塔模型可能是 alpha）
  param_name <- names(res2d$param_list[[1]])[1]
  param_values <- sapply(res2d$param_list, function(x) x[param_name])
  mean_risk <- res2d$mean_risk
  sensitivity_2d(param_values, mean_risk, param_name = param_name)
}
