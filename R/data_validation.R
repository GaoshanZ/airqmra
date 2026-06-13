# =============================================================================
# 模块四：数据预处理与验证工具
# =============================================================================

#' Check an airqmra scenario for potential issues
#'
#' Performs a series of sanity checks on an `airqmra_scenario` object,
#' including parameter ranges, distribution functions, and compatibility
#' between settings. Warnings are issued for suspicious values; messages
#' indicate informational notes.
#'
#' @param scenario An `airqmra_scenario` object.
#'
#' @return Invisibly returns `TRUE` if all checks pass (no critical errors),
#'   `FALSE` otherwise. Warnings and messages are printed to the console.
#' @export
#'
#' @examples
#' scen <- new_scenario(pathogen_id = "legionella pneumophila",
#'                      exposure_group = "Adult", microenvironment = "Outdoor")
#' check_scenario(scen)
check_scenario <- function(scenario) {
  problems <- 0

  # ---- 1. Essential components ----
  if (is.null(scenario$C_air_dist) || !is.function(scenario$C_air_dist))
    stop("C_air_dist must be a function.")
  if (is.null(scenario$IR_dist) || !is.function(scenario$IR_dist))
    stop("IR_dist must be a function.")
  if (is.null(scenario$t_exp_dist) || !is.function(scenario$t_exp_dist))
    stop("t_exp_dist must be a function.")

  # ---- 2. Dose-response parameters ----
  if (!is.null(scenario$dr_pars)) {
    dr <- scenario$dr_pars
    if (dr$type == "exponential") {
      if (!is.finite(dr$k) || dr$k <= 0) {
        warning("Exponential model: k must be positive.")
        problems <- problems + 1
      }
      if (is.finite(dr$k) && dr$k > 1) {
        warning("Exponential model: k = ", signif(dr$k, 3), " seems very large.")
        problems <- problems + 1
      }
    } else if (dr$type == "exact_beta_poisson") {
      if (!is.finite(dr$alpha) || dr$alpha <= 0) {
        warning("Beta-Poisson model: alpha must be positive.")
        problems <- problems + 1
      }
      if (!is.finite(dr$beta) || dr$beta <= 0) {
        warning("Beta-Poisson model: beta must be positive.")
        problems <- problems + 1
      }
    } else if (dr$type == "approx_beta_poisson") {
      if (!is.finite(dr$alpha) || dr$alpha <= 0) {
        warning("Approx beta-Poisson model: alpha must be positive.")
        problems <- problems + 1
      }
      if (!is.finite(dr$N50) || dr$N50 <= 0) {
        warning("Approx beta-Poisson model: N50 must be positive.")
        problems <- problems + 1
      }
    }
  } else if (is.null(scenario$dr_generator)) {
    stop("Either dr_pars or dr_generator must be present.")
  }

  # ---- 3. AG, deposition, survival ----
  # Check by generating a few samples and ensuring they are in [0,1]
  check_frac <- function(f, name) {
    tryCatch({
      vals <- f(100)
      if (any(vals < 0 | vals > 1)) {
        warning(name, " distribution produces values outside [0,1].")
        problems <<- problems + 1
      }
    }, error = function(e) {
      warning("Could not evaluate ", name, " distribution: ", e$message)
      problems <<- problems + 1
    })
  }
  check_frac(scenario$AG_dist, "AG")
  check_frac(scenario$deposition_dist, "Deposition")
  check_frac(scenario$survival_dist, "Survival")

  # ---- 4. Exposure group vs available classifications ----
  if (!is.null(scenario$exposure_group)) {
    available <- unique(airqmra::ir_long_table$classification)
    if (!scenario$exposure_group %in% available) {
      message("Note: exposure_group '", scenario$exposure_group,
              "' not found in ir_long_table. Fallback values may be used.")
    }
  }

  # ---- 5. Deposition population ----
  if (!is.null(scenario$population_deposition)) {
    depo_pop <- unique(airqmra::deposition_db$population)
    if (!scenario$population_deposition %in% depo_pop) {
      message("Note: population_deposition '", scenario$population_deposition,
              "' not in deposition_db; default fallback is 0.5.")
    }
  }

  if (problems == 0) {
    message("Scenario appears valid.")
    invisible(TRUE)
  } else {
    message("Found ", problems, " potential issue(s). Please review warnings.")
    invisible(FALSE)
  }
}

#' Validate concentration data before fitting
#'
#' Checks a concentration data frame (or numeric vector) for common issues:
#' non-positive values, missing values, outliers, and extreme skewness.
#' Returns a list with diagnostic statistics and, if called interactively,
#' prints a summary and a simple histogram.
#'
#' @param x A data frame (if \code{col} is specified) or a numeric vector of
#'   concentration values.
#' @param col Character. Name of the column containing concentrations. Ignored
#'   if \code{x} is numeric.
#' @param plot Logical. If \code{TRUE} (default), a diagnostic histogram with
#'   log10 x‑axis is displayed.
#'
#' @return Invisibly returns a list with elements:
#'   \item{n}{Number of valid observations.}
#'   \item{n_missing}{Number of missing values.}
#'   \item{range}{Min and max values.}
#'   \item{zeros}{Number of zero (or negative) values.}
#'   \item{outliers}{Number of potential outliers (beyond 1.5*IQR on log10 scale).}
#'   \item{suggestions}{Character vector of suggestions.}
#' @export
validate_concentration_data <- function(x, col = NULL, plot = TRUE) {
  # Extract vector
  if (is.data.frame(x)) {
    if (is.null(col)) {
      num_cols <- sapply(x, is.numeric)
      if (!any(num_cols)) stop("No numeric column found.")
      col <- names(which(num_cols))[1]
    }
    if (!col %in% names(x)) stop("Column '", col, "' not found.")
    conc <- x[[col]]
  } else if (is.numeric(x)) {
    conc <- x
  } else {
    stop("'x' must be a numeric vector or data frame.")
  }

  n_total <- length(conc)
  n_missing <- sum(is.na(conc))
  conc <- conc[!is.na(conc)]
  n_valid <- length(conc)

  range_val <- range(conc)
  zeros <- sum(conc <= 0)

  suggestions <- character()
  if (zeros > 0) {
    suggestions <- c(suggestions,
                     paste0(zeros, " value(s) are zero or negative. Consider removing them."))
  }

  # Outliers on log10 scale (if all positive)
  if (all(conc > 0)) {
    log_conc <- log10(conc)
    iqr <- IQR(log_conc)
    q1 <- quantile(log_conc, 0.25)
    q3 <- quantile(log_conc, 0.75)
    lower <- q1 - 1.5 * iqr
    upper <- q3 + 1.5 * iqr
    n_out <- sum(log_conc < lower | log_conc > upper)
    if (n_out > 0) {
      suggestions <- c(suggestions,
                       paste0(n_out, " potential outlier(s) detected (log10 scale)."))
    }
  }

  if (n_valid < 3) {
    suggestions <- c(suggestions, "Less than 3 valid observations; fitting may be unreliable.")
  }

  # Print summary
  cat("Concentration data validation\n")
  cat("  Total observations:", n_total, "\n")
  cat("  Missing:", n_missing, "\n")
  cat("  Valid:", n_valid, "\n")
  cat("  Range:", range_val[1], "-", range_val[2], "\n")
  cat("  Zeros / negatives:", zeros, "\n")
  if (length(suggestions) > 0) {
    cat("  Suggestions:\n")
    for (s in suggestions) cat("    -", s, "\n")
  } else {
    cat("  No obvious issues detected.\n")
  }

  if (plot && n_valid >= 2) {
    if (all(conc > 0)) {
      hist(log10(conc), breaks = 15, main = "Log10 concentration distribution",
           xlab = "log10(concentration)", col = "lightgreen")
    } else {
      hist(conc, breaks = 15, main = "Concentration distribution",
           xlab = "Concentration", col = "lightgreen")
    }
  }

  invisible(list(n = n_valid, n_missing = n_missing, range = range_val,
                 zeros = zeros, outliers = if (all(conc > 0)) n_out else NA,
                 suggestions = suggestions))
}

#' Suggest best-fitting distribution for concentration data
#'
#' Compares several candidate distributions (log-normal, gamma, Weibull)
#' using AIC and returns the name of the distribution with the lowest AIC.
#'
#' @param conc Numeric vector of (positive) concentration values.
#'
#' @return A list with elements:
#'   \item{best}{Character. Name of the best distribution.}
#'   \item{aic_table}{Data frame with AIC values for all tested distributions.}
#' @export
#'
#' @examples
#' conc <- c(45, 52, 38, 60, 71, 55)
#' suggest_distribution(conc)
suggest_distribution <- function(conc) {
  if (any(conc <= 0)) stop("All concentration values must be positive.")
  if (length(conc) < 3) stop("At least 3 values required.")

  fits <- list()
  tryCatch({ fits$lognormal <- fitdistrplus::fitdist(conc, "lnorm")$aic }, error = function(e) NULL)
  tryCatch({ fits$gamma     <- fitdistrplus::fitdist(conc, "gamma")$aic }, error = function(e) NULL)
  tryCatch({ fits$weibull   <- fitdistrplus::fitdist(conc, "weibull")$aic }, error = function(e) NULL)

  if (length(fits) == 0) stop("None of the distributions could be fitted.")

  aic_vals <- unlist(fits)
  best <- names(which.min(aic_vals))

  cat("AIC comparison:\n")
  for (nm in names(aic_vals)) cat("  ", nm, ":", round(aic_vals[nm], 1), "\n")
  cat("  Recommended distribution:", best, "\n")

  invisible(list(best = best,
                 aic_table = data.frame(distribution = names(aic_vals),
                                        AIC = aic_vals, stringsAsFactors = FALSE)))
}
