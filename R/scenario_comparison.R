# =============================================================================
# 模块二：情景比较与批量评估
# =============================================================================

#' Create a grid of exposure scenarios from parameter combinations
#'
#' Given lists of parameter values, this function creates all combinations
#' and returns a list of `airqmra_scenario` objects. Parameters not specified
#' will use the default values from `new_scenario`.
#'
#' @param pathogen_id Character vector (length 1). Pathogen identifier.
#' @param C_air_dist A function, or a list of functions. If a list, each
#'   element will be used as the concentration distribution in one scenario.
#'   If a single function, it will be reused for all scenarios.
#' @param exposure_group Character vector. Population groups.
#' @param microenvironment Character vector. Microenvironments.
#' @param season Character vector. Seasons.
#' @param activity_template Character vector. Activity template IDs.
#' @param activity_intensity Character vector. Activity intensities.
#' @param ... Other named arguments that will be passed directly to
#'   `new_scenario`. If a value is a list, it will be expanded across
#'   scenarios; otherwise it is constant.
#'
#' @return A list of `airqmra_scenario` objects, with names describing the
#'   varying parameters.
#' @export
#'
#' @examples
#' # Compare two concentration distributions (indoor vs outdoor)
#' conc1 <- function(n) rlnorm(n, log(50), log(1.5))
#' conc2 <- function(n) rlnorm(n, log(100), log(1.5))
#' scenarios <- create_scenario_grid(
#'   pathogen_id = "legionella pneumophila",
#'   C_air_dist = list(conc1, conc2),
#'   microenvironment = c("Indoor", "Outdoor"),
#'   exposure_group = "Adult"
#' )
#' length(scenarios)  # 2 * 2 = 4
create_scenario_grid <- function(pathogen_id,
                                 C_air_dist = function(n) rlnorm(n, log(50), log(1.8)),
                                 exposure_group = "Adult",
                                 microenvironment = "Indoor",
                                 season = "Annual",
                                 activity_template = NULL,
                                 activity_intensity = NULL,
                                 ...) {

  # Collect varying arguments: those that are lists will be expanded
  var_args <- list(
    C_air_dist = if (is.function(C_air_dist)) list(C_air_dist) else C_air_dist,
    exposure_group = as.list(exposure_group),
    microenvironment = as.list(microenvironment),
    season = as.list(season),
    activity_template = if (is.null(activity_template)) list(NULL) else as.list(activity_template),
    activity_intensity = if (is.null(activity_intensity)) list(NULL) else as.list(activity_intensity)
  )

  # Additional arguments
  dots <- list(...)
  for (nm in names(dots)) {
    val <- dots[[nm]]
    if (is.list(val) && !is.function(val)) {
      var_args[[nm]] <- val
    } else {
      var_args[[nm]] <- list(val)
    }
  }

  # Expand all combinations
  grid <- expand.grid(var_args, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

  scenarios <- vector("list", nrow(grid))
  names_list <- character(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    args <- list(pathogen_id = pathogen_id)
    row <- grid[i, , drop = FALSE]

    # Build call arguments
    args$C_air_dist <- row$C_air_dist[[1]]
    args$exposure_group <- row$exposure_group[[1]]
    args$microenvironment <- row$microenvironment[[1]]
    args$season <- row$season[[1]]
    if (!is.null(row$activity_template[[1]])) args$activity_template <- row$activity_template[[1]]
    if (!is.null(row$activity_intensity[[1]])) args$activity_intensity <- row$activity_intensity[[1]]

    # Add extra arguments
    extra_names <- setdiff(names(row), c("C_air_dist","exposure_group","microenvironment",
                                         "season","activity_template","activity_intensity"))
    for (en in extra_names) {
      args[[en]] <- row[[en]][[1]]
    }

    scenarios[[i]] <- do.call(new_scenario, args)

    # Build descriptive name
    desc_parts <- character()
    if (length(unique(grid$exposure_group)) > 1) desc_parts <- c(desc_parts, paste0("Grp=", args$exposure_group))
    if (length(unique(grid$microenvironment)) > 1) desc_parts <- c(desc_parts, args$microenvironment)
    if (length(unique(grid$season)) > 1) desc_parts <- c(desc_parts, args$season)
    if (length(unique(grid$C_air_dist)) > 1) desc_parts <- c(desc_parts, paste0("Conc", i))
    if (length(unique(grid$activity_template)) > 1) desc_parts <- c(desc_parts, args$activity_template)
    names_list[i] <- if (length(desc_parts) > 0) paste(desc_parts, collapse = ", ") else paste0("Scenario_", i)
  }

  names(scenarios) <- names_list
  scenarios
}

#' Batch run 1D Monte Carlo simulations
#'
#' Runs `run_qmra_1d` on a list of scenarios and returns a list of raw risk
#' vectors (optionally plus annual risk) and a summary table.
#'
#' @param scenario_list A list of `airqmra_scenario` objects, e.g. from
#'   `create_scenario_grid`.
#' @param n_iter Number of iterations per scenario.
#' @param frequency Annual frequency, passed to `annual_risk`. If NULL,
#'   annual risk is not computed.
#'
#' @return A list with two elements:
#'   \item{raw}{List of risk vectors (and annual if frequency provided).}
#'   \item{summary}{Data frame with quantiles and means.}
#' @export
batch_simulate_1d <- function(scenario_list, n_iter = 5000, frequency = NULL) {
  raw <- vector("list", length(scenario_list))
  names(raw) <- names(scenario_list)

  summary_rows <- vector("list", length(scenario_list))

  for (i in seq_along(scenario_list)) {
    scen <- scenario_list[[i]]
    single <- run_qmra_1d(scen, n_iter = n_iter)

    if (!is.null(frequency)) {
      annual <- annual_risk(single, frequency = frequency)
      raw[[i]] <- list(single = single, annual = annual)
      summ <- data.frame(
        Scenario = names(scenario_list)[i],
        Median_single = median(single),
        Mean_single = mean(single),
        P95_single = quantile(single, 0.95),
        Median_annual = median(annual),
        Mean_annual = mean(annual),
        P95_annual = quantile(annual, 0.95)
      )
    } else {
      raw[[i]] <- list(single = single)
      summ <- data.frame(
        Scenario = names(scenario_list)[i],
        Median = median(single),
        Mean = mean(single),
        P95 = quantile(single, 0.95)
      )
    }
    summary_rows[[i]] <- summ
  }

  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL

  list(raw = raw, summary = summary_table)
}

#' Compare multiple QMRA simulation results
#'
#' Takes either a list of raw risk vectors (from `batch_simulate_1d`) or
#' a list of `airqmra_2d` objects, and creates a summary comparison table.
#'
#' @param result_list A list where each element is either a numeric vector
#'   (single risk) or a list containing `single`/`overall_risk`.
#' @param probs Numeric vector of quantiles to compute. Default
#'   `c(0.5, 0.95, 0.99)`.
#' @param annual Logical. If TRUE, extracts the annual risk component
#'   (the element must have a name `annual`).
#'
#' @return A data frame with rows = scenarios, columns = quantiles + mean.
#' @export
compare_scenarios <- function(result_list, probs = c(0.5, 0.95, 0.99),
                              annual = FALSE) {
  scenario_names <- names(result_list)
  if (is.null(scenario_names)) scenario_names <- paste0("Scenario", seq_along(result_list))

  out <- data.frame(Scenario = scenario_names, stringsAsFactors = FALSE)

  for (i in seq_along(result_list)) {
    entry <- result_list[[i]]
    if (annual) {
      if (is.list(entry) && !is.null(entry$annual)) {
        risk_vec <- entry$annual
      } else if (is.list(entry) && !is.null(entry$overall_annual_risk)) {
        risk_vec <- entry$overall_annual_risk
      } else {
        stop("Annual risk not found in element ", i)
      }
    } else {
      if (is.numeric(entry)) {
        risk_vec <- entry
      } else if (is.list(entry) && !is.null(entry$single)) {
        risk_vec <- entry$single
      } else if (is.list(entry) && !is.null(entry$overall_risk)) {
        risk_vec <- entry$overall_risk
      } else {
        stop("Cannot extract risk vector from element ", i)
      }
    }
    quants <- quantile(risk_vec, probs = probs, na.rm = TRUE)
    out[i, paste0("Q", probs*100)] <- quants
    out[i, "Mean"] <- mean(risk_vec, na.rm = TRUE)
  }
  out
}

#' Plot scenario comparison
#'
#' Creates a ggplot2 boxplot or density plot comparing risk distributions
#' across scenarios.
#'
#' @param result_list A list of risk vectors or objects containing risk vectors.
#' @param plot_type Character. "box" or "density". Default "box".
#' @param annual Logical. If TRUE, plot annual risk.
#' @param log_scale Logical. If TRUE, use log10 y-axis. Default TRUE.
#'
#' @return A ggplot2 object.
#' @export
plot_scenario_comparison <- function(result_list, plot_type = c("box", "density"),
                                     annual = FALSE, log_scale = TRUE) {
  plot_type <- match.arg(plot_type)

  # Extract risk vectors and scenario names
  scenarios <- names(result_list)
  if (is.null(scenarios)) scenarios <- paste0("Scenario", seq_along(result_list))

  risk_df <- data.frame()
  for (i in seq_along(result_list)) {
    entry <- result_list[[i]]
    if (annual) {
      if (is.list(entry) && !is.null(entry$annual)) {
        vec <- entry$annual
      } else if (is.list(entry) && !is.null(entry$overall_annual_risk)) {
        vec <- entry$overall_annual_risk
      } else stop("Annual risk not found in element ", i)
    } else {
      if (is.numeric(entry)) vec <- entry
      else if (is.list(entry) && !is.null(entry$single)) vec <- entry$single
      else if (is.list(entry) && !is.null(entry$overall_risk)) vec <- entry$overall_risk
      else stop("Cannot extract risk vector from element ", i)
    }
    risk_df <- rbind(risk_df,
                     data.frame(Scenario = scenarios[i], Risk = vec,
                                stringsAsFactors = FALSE))
  }

  if (plot_type == "box") {
    p <- ggplot2::ggplot(risk_df, ggplot2::aes(x = Scenario, y = Risk)) +
      ggplot2::geom_boxplot(fill = "steelblue", alpha = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    p <- ggplot2::ggplot(risk_df, ggplot2::aes(x = Risk, color = Scenario, fill = Scenario)) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::theme_minimal()
  }

  if (log_scale) {
    p <- p + ggplot2::scale_y_log10(labels = scales::scientific)
  }

  p + ggplot2::labs(title = ifelse(annual, "Annual Risk Comparison", "Single Exposure Risk Comparison"),
                    y = "Risk", x = "")
}

#' Annual risk summary table from multiple 2D simulation results
#'
#' Extracts the annual risk summary (mean, median, confidence limits) from
#' a list of `airqmra_2d` objects and combines them into a data frame.
#'
#' @param res2d_list A named list of `airqmra_2d` objects.
#'
#' @return A data frame with columns: Scenario, Mean_annual, Median_annual,
#'   CI_low, CI_high.
#' @export
annual_risk_table <- function(res2d_list) {
  if (is.null(names(res2d_list))) {
    names(res2d_list) <- paste0("Scenario", seq_along(res2d_list))
  }

  tbl <- data.frame(Scenario = names(res2d_list),
                    Mean_annual = numeric(length(res2d_list)),
                    Median_annual = numeric(length(res2d_list)),
                    CI_low = numeric(length(res2d_list)),
                    CI_high = numeric(length(res2d_list)),
                    stringsAsFactors = FALSE)

  for (i in seq_along(res2d_list)) {
    obj <- res2d_list[[i]]
    if (is.null(obj$annual_summary)) {
      stop("Element ", i, " does not contain annual_summary. Did you provide a 2D result with frequency?")
    }
    tbl$Mean_annual[i]   <- obj$annual_summary$Mean
    tbl$Median_annual[i] <- obj$annual_summary$Median
    tbl$CI_low[i]        <- obj$annual_summary$CI_low
    tbl$CI_high[i]       <- obj$annual_summary$CI_high
  }
  tbl
}
