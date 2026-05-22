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


#' Sensitivity analysis for a 1D (point-estimate) QMRA scenario
#'
#' Runs a 1D Monte Carlo simulation, records all input random samples and the
#' resulting risk, and computes Spearman rank correlation coefficients.
#' This function is intended for point-estimate scenarios (`use_uncertainty = FALSE`).
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
  # Generate all input variable samples
  C <- scenario$C_air_dist(n_iter)
  IR <- scenario$IR_dist(n_iter)
  t <- scenario$t_exp_dist(n_iter)
  survival <- scenario$survival_dist(n_iter)
  deposition <- scenario$deposition_dist(n_iter)
  # AG fixed to 1, multiplied in dose but excluded from sensitivity
  AG <- scenario$AG_dist(n_iter)

  dose <- C * IR * t * survival * deposition * AG

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

  # Build data frame without AG
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
#' scen_unc <- new_scenario(pathogen_id = "legionella pneumophila",
#'                          exposure_group = "Adult",
#'                          microenvironment = "Outdoor",
#'                          population_deposition = "Male",
#'                          use_uncertainty = TRUE)
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

  # Flatten variability input matrices into vectors
  C_vec <- as.vector(res2d$C_matrix)
  IR_vec <- as.vector(res2d$IR_matrix)
  t_vec <- as.vector(res2d$t_matrix)
  survival_vec <- as.vector(res2d$survival_matrix)
  deposition_vec <- as.vector(res2d$deposition_matrix)

  # Parameter vectors: match length n_var * n_unc
  param_df <- do.call(rbind, lapply(res2d$param_list, function(p) {
    data.frame(t(p))
  }))
  # Replicate each column n_var times
  param_vec_list <- lapply(param_df, function(x) {
    rep(x, each = n_var)
  })

  # Overall risk vector
  overall_risk <- res2d$overall_risk

  # Construct data frame
  df <- data.frame(
    C = C_vec,
    IR = IR_vec,
    t = t_vec,
    survival = survival_vec,
    deposition = deposition_vec
  )
  # Add parameter columns
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
  # Extract ranking data frame
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
    # Base R barplot
    barplot(ranking$Correlation, names.arg = ranking$Variable, horiz = TRUE,
            col = ifelse(ranking$Correlation > 0, "steelblue", "firebrick"),
            main = title, xlab = xlab)
  }
}


#' Sensitivity analysis for 2D Monte Carlo results (convenience wrapper)
#'
#' @param res2d An object of class `airqmra_2d`.
#' @return A list with correlations and ranking for the primary dose-response parameter.
#' @export
sensitivity_2d_auto <- function(res2d) {
  if (is.null(res2d$param_list)) stop("res2d must contain param_list (always present).")
  # Extract the first parameter (k for exponential, alpha for beta models)
  param_name <- names(res2d$param_list[[1]])[1]
  param_values <- sapply(res2d$param_list, function(x) x[param_name])
  mean_risk <- res2d$mean_risk
  sensitivity_2d(param_values, mean_risk, param_name = param_name)
}


#' Global sensitivity analysis for annual risk (combined variability & uncertainty)
#'
#' Uses recorded exposure inputs and annual risk from a 2D simulation to compute
#' Spearman correlations for all input variables (exposure and dose-response parameters).
#'
#' @param res2d An object of class `airqmra_2d` created with
#'   `run_qmra_2d(..., frequency = <your frequency>, record_inputs = TRUE)`.
#' @param method Correlation method, default `"spearman"`.
#'
#' @return A list as from \code{\link{compute_sensitivity}} with correlations and ranking.
#' @export
sensitivity_global_annual <- function(res2d, method = "spearman") {
  if (is.null(res2d$annual_matrix)) {
    stop("The 2D simulation must include annual risk (supply 'frequency' to run_qmra_2d).")
  }
  if (is.null(res2d$C_matrix) || is.null(res2d$param_list)) {
    stop("The 2D simulation must be run with 'record_inputs = TRUE'.")
  }

  n_var <- nrow(res2d$C_matrix)
  n_unc <- ncol(res2d$C_matrix)

  # Flatten variability inputs
  C_vec <- as.vector(res2d$C_matrix)
  IR_vec <- as.vector(res2d$IR_matrix)
  t_vec <- as.vector(res2d$t_matrix)
  survival_vec <- as.vector(res2d$survival_matrix)
  deposition_vec <- as.vector(res2d$deposition_matrix)

  # Expand parameter vectors
  param_df <- do.call(rbind, lapply(res2d$param_list, function(p) data.frame(t(p))))
  param_vec_list <- lapply(param_df, function(x) rep(x, each = n_var))

  # Overall annual risk vector
  overall_annual_risk <- as.vector(res2d$annual_matrix)

  # Construct data frame
  df <- data.frame(
    C = C_vec,
    IR = IR_vec,
    t = t_vec,
    survival = survival_vec,
    deposition = deposition_vec
  )
  for (nm in names(param_vec_list)) {
    df[[nm]] <- param_vec_list[[nm]]
  }
  df$annual_risk <- overall_annual_risk

  compute_sensitivity(df, output_col = "annual_risk", method = method)
}


#' Sensitivity analysis for annual risk (variability only, point-estimate mode)
#'
#' Runs a 1D Monte Carlo simulation, converts single-exposure risk to annual risk
#' using the given frequency, and computes Spearman correlations between exposure
#' variables and the annual risk.
#'
#' @param scenario An `airqmra_scenario` object (from `new_scenario` with `use_uncertainty = FALSE`).
#' @param n_iter Number of Monte Carlo iterations.
#' @param frequency Either a fixed integer or a function for annual exposure frequency.
#'
#' @return A list as from `compute_sensitivity`, with additional element `annual_risk`.
#' @export
sensitivity_1d_annual <- function(scenario, n_iter = 5000, frequency) {
  # Generate all input variable samples
  C <- scenario$C_air_dist(n_iter)
  IR <- scenario$IR_dist(n_iter)
  t <- scenario$t_exp_dist(n_iter)
  survival <- scenario$survival_dist(n_iter)
  deposition <- scenario$deposition_dist(n_iter)
  AG <- scenario$AG_dist(n_iter)

  dose <- C * IR * t * survival * deposition * AG

  dr <- scenario$dr_pars
  type <- dr$type
  if (type == "exponential") {
    k <- dr$k
    single_risk <- 1 - exp(-k * dose)
  } else if (type == "exact_beta_poisson") {
    alpha <- dr$alpha
    beta <- dr$beta
    single_risk <- exact_beta_poisson(dose, alpha, beta)
  } else if (type == "approx_beta_poisson") {
    alpha <- dr$alpha
    N50 <- dr$N50
    single_risk <- approx_beta_poisson(dose, alpha, N50)
  } else {
    stop("Unknown dose-response type: ", type)
  }

  # Calculate annual risk
  annual_risk_vec <- annual_risk(single_risk, frequency = frequency)

  # Build data frame (no AG)
  df <- data.frame(
    C = C,
    IR = IR,
    t = t,
    survival = survival,
    deposition = deposition,
    annual_risk = annual_risk_vec
  )

  sens <- compute_sensitivity(df, output_col = "annual_risk")
  list(sensitivity = sens, data = df, annual_risk = annual_risk_vec)
}


#' Sensitivity of annual risk to dose-response parameter uncertainty
#'
#' Uses the recorded parameter values and annual mean risk from a 2D simulation
#' to calculate Spearman correlation between the dose-response parameter (e.g., k)
#' and the mean annual risk per uncertainty iteration.
#'
#' @param res2d An object of class `airqmra_2d` created with `frequency` specified.
#' @return A list as from `compute_sensitivity`.
#' @export
sensitivity_2d_annual_auto <- function(res2d) {
  if (is.null(res2d$annual_mean)) {
    stop("The 2D simulation must include annual risk (supply 'frequency' to run_qmra_2d).")
  }
  if (is.null(res2d$param_list)) {
    stop("res2d must contain param_list (always present).")
  }
  # Extract first parameter name
  param_name <- names(res2d$param_list[[1]])[1]
  param_values <- sapply(res2d$param_list, function(x) x[param_name])
  df <- data.frame(param = param_values, mean_annual_risk = res2d$annual_mean)
  names(df)[1] <- param_name
  compute_sensitivity(df, output_col = "mean_annual_risk")
}
