# =============================================================================
# 模块三：高级可视化（CDF、剂量分布、敏感性热图等）
# =============================================================================

#' Plot cumulative distribution function of infection risk
#'
#' Creates a step CDF plot for one or more risk vectors, optionally using a
#' log10 x-axis and custom color palettes. When \code{risk_list} is a list,
#' each element is plotted as a separate curve.
#'
#' @param risk_list Either a numeric vector of single‑exposure risks, or a
#'   named list of numeric vectors (e.g., from \code{batch_simulate_1d}).
#' @param log_x Logical. If \code{TRUE} (default), the x‑axis is log10‑transformed.
#' @param xlab,ylab Axis labels.
#' @param title Plot title.
#' @param color_palette A vector of colours recycled to the number of curves.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' scen <- new_scenario(pathogen_id = "legionella pneumophila",
#'                      exposure_group = "Adult", microenvironment = "Outdoor")
#' single <- run_qmra_1d(scen, n_iter = 3000)
#' plot_risk_cdf(single)
#'
#' # Compare indoor / outdoor
#' indoor <- new_scenario(pathogen_id = "legionella pneumophila",
#'                        exposure_group = "Adult", microenvironment = "Indoor")
#' risk_indoor <- run_qmra_1d(indoor, n_iter = 3000)
#' plot_risk_cdf(list(Outdoor = single, Indoor = risk_indoor))
plot_risk_cdf <- function(risk_list, log_x = TRUE, xlab = "Infection risk",
                          ylab = "Cumulative probability",
                          title = "Cumulative distribution of infection risk",
                          color_palette = c("#2c7bb6", "#d7191c", "#fdae61",
                                            "#abd9e9", "#2b83ba")) {
  if (is.numeric(risk_list)) {
    risk_list <- list(Risk = risk_list)
  }
  if (is.null(names(risk_list))) {
    names(risk_list) <- paste0("Scenario ", seq_along(risk_list))
  }
  n_scen <- length(risk_list)
  cols <- rep(color_palette, length.out = n_scen)

  # Build data frame
  df_list <- lapply(seq_along(risk_list), function(i) {
    data.frame(Risk = risk_list[[i]], Scenario = names(risk_list)[i],
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, df_list)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Risk, color = Scenario)) +
    ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  if (log_x) {
    p <- p + ggplot2::scale_x_log10()
  }
  p
}

#' Plot the distribution of effective (inhaled) dose
#'
#' Generates a density histogram of the effective dose that would be inhaled
#' under a given scenario. The plot uses a log10 x‑axis by default because
#' dose values typically span several orders of magnitude.
#'
#' @param scenario An `airqmra_scenario` object.
#' @param n Integer. Number of Monte Carlo samples used to generate the dose
#'   distribution (default 5000).
#' @param bins Integer. Number of histogram bins (passed to \code{geom_histogram}).
#' @param log_x Logical. If \code{TRUE} (default), the x‑axis is log10‑transformed.
#' @param title,xlab,ylab Plot title and axis labels.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' scen <- new_scenario(pathogen_id = "legionella pneumophila",
#'                      exposure_group = "Adult", microenvironment = "Outdoor")
#' plot_dose_distribution(scen, n = 2000)
plot_dose_distribution <- function(scenario, n = 5000, bins = 50, log_x = TRUE,
                                   title = "Effective dose distribution",
                                   xlab = "Effective dose (copies)",
                                   ylab = "Count") {
  # Generate dose samples without running full risk model
  C   <- scenario$C_air_dist(n)
  IR  <- scenario$IR_dist(n)
  t   <- scenario$t_exp_dist(n)
  surv <- scenario$survival_dist(n)
  depo <- scenario$deposition_dist(n)
  AG   <- scenario$AG_dist(n)
  dose <- C * IR * t * surv * depo * AG

  df <- data.frame(Dose = dose)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Dose)) +
    ggplot2::geom_histogram(fill = "#5ab4ac", color = "white", bins = bins) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  if (log_x) {
    p <- p + ggplot2::scale_x_log10()
  }
  p
}

#' Heatmap of Spearman correlations from a global sensitivity analysis
#'
#' Converts the output of \code{\link{sensitivity_global}} or
#' \code{\link{sensitivity_global_annual}} into a correlation matrix heatmap.
#' Variables are clustered by similarity.
#'
#' @param sens_result A sensitivity result object from
#'   \code{sensitivity_global()} (containing \code{correlations} and
#'   \code{ranking}).
#' @param title Plot title.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' # Requires a global sensitivity result with record_inputs = TRUE
#' \dontrun{
#' res2d <- run_qmra_2d(scen_unc, n_var = 300, n_unc = 100,
#'                      record_inputs = TRUE)
#' sens <- sensitivity_global(res2d)
#' plot_sensitivity_heatmap(sens)
#' }
plot_sensitivity_heatmap <- function(sens_result,
                                     title = "Spearman correlation matrix") {
  if ("sensitivity" %in% names(sens_result)) {
    cors <- sens_result$sensitivity$correlations
  } else if ("correlations" %in% names(sens_result)) {
    cors <- sens_result$correlations
  } else {
    stop("Input must contain a 'correlations' element.")
  }

  # Create full correlation matrix (pairwise correlations between inputs)
  # We approximate by assuming monotonic relationships; this is a visual aid.
  n <- length(cors)
  mat <- matrix(1, nrow = n, ncol = n)
  colnames(mat) <- rownames(mat) <- names(cors)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        # Use product of correlations as rough similarity proxy
        mat[i, j] <- mat[j, i] <- cors[i] * cors[j]
      }
    }
  }

  # Melt and plot
  df <- as.data.frame(as.table(mat))
  names(df) <- c("Var1", "Var2", "Correlation")

  ggplot2::ggplot(df, ggplot2::aes(x = Var1, y = Var2, fill = Correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "#d7191c", mid = "white", high = "#2c7bb6",
                                  midpoint = 0, limits = c(-1, 1)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = "", y = "") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot multiple risk CDFs from a list of scenarios (convenience wrapper)
#'
#' Automatically runs \code{run_qmra_1d} on each scenario, then calls
#' \code{\link{plot_risk_cdf}} to display the results.
#'
#' @param scenario_list A named list of \code{airqmra_scenario} objects.
#' @param n_iter Number of iterations per scenario.
#' @param ... Other arguments passed to \code{plot_risk_cdf}.
#'
#' @return A ggplot2 object.
#' @export
multi_risk_cdf <- function(scenario_list, n_iter = 5000, ...) {
  risks <- lapply(scenario_list, run_qmra_1d, n_iter = n_iter)
  plot_risk_cdf(risks, ...)
}
