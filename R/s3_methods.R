# ============================================================================
# S3 methods for airqmra objects
# ============================================================================

#' @export
print.airqmra_scenario <- function(x, ...) {
  cat("airqmra inhalation exposure scenario\n")
  cat("--------------------------------------\n")
  cat("Pathogen      :", x$pathogen, "\n")
  if (!is.null(x$pathogen_id)) cat("Pathogen ID   :", x$pathogen_id, "\n")
  cat("Exposure group:", ifelse(is.null(x$exposure_group), "User-defined", x$exposure_group), "\n")
  cat("Microenvironment:", x$microenvironment, "\n")
  cat("Season        :", x$season, "\n")
  if (!is.null(x$activity_template))
    cat("Activity template:", x$activity_template, " (", x$activity_intensity, ")\n", sep = "")
  if (!is.null(x$dr_pars)) {
    cat("Dose-response model:", x$dr_pars$type, "\n")
    if (x$dr_pars$type == "exponential") {
      cat("  k =", signif(x$dr_pars$k, 3), "\n")
    } else if (x$dr_pars$type == "exact_beta_poisson") {
      cat("  alpha =", x$dr_pars$alpha, " beta =", signif(x$dr_pars$beta, 3), "\n")
    } else if (x$dr_pars$type == "approx_beta_poisson") {
      cat("  alpha =", x$dr_pars$alpha, " N50 =", signif(x$dr_pars$N50, 3), "\n")
    }
  }
  if (!is.null(x$dr_generator)) cat("Uncertainty: ENABLED\n")
  invisible(x)
}

#' @export
summary.airqmra_scenario <- function(object, ...) {
  print(object)
}

#' Update parameters of an existing scenario
#'
#' Takes an existing scenario object and changes one or more exposure-related
#' parameters without re-specifying all settings. Only the parameters that are
#' explicitly provided are updated; all others are kept as in the original
#' scenario.  Dose-response parameters cannot be changed via this function;
#' to modify them, create a new scenario with `new_scenario()`.
#'
#' @param scenario An `airqmra_scenario` object.
#' @param pathogen_name Character. A descriptive name for the scenario.
#' @param exposure_group Character. Population classification for exposure
#'   parameters (e.g., "Adult", "Male", "18Y-44Y").
#' @param microenvironment Character. "Indoor" or "Outdoor".
#' @param season Character. "Annual", "Spring & Autumn", "Summer", "Winter".
#' @param activity_template Character. Activity template ID
#'   (e.g., "office_worker").
#' @param activity_intensity Character. Activity intensity
#'   (e.g., "Light_intensity").
#' @param C_air_dist Function. A function that generates air concentration
#'   samples (copies/m^3).
#' @param IR_dist Function. A function that generates inhalation rate
#'   samples (m^3/h).
#' @param t_exp Single numeric or function. Exposure duration in hours.
#' @param AG Numeric or function. Fraction of deposited pathogens that
#'   initiate infection.
#' @param deposition Numeric or function. Fraction of inhaled particles that
#'   deposit in the lung.
#' @param survival Numeric or function. Fraction of airborne pathogens that
#'   remain viable upon inhalation.
#' @param population_deposition Character. Population group for deposition
#'   data lookup (e.g., "Male", "Female", "Child", "Senior").
#'
#' @return An updated `airqmra_scenario` object.
#'
#' @export
update_scenario <- function(scenario,
                            pathogen_name = NULL,
                            exposure_group = NULL,
                            microenvironment = NULL,
                            season        = NULL,
                            activity_template = NULL,
                            activity_intensity = NULL,
                            C_air_dist    = NULL,
                            IR_dist       = NULL,
                            t_exp         = NULL,
                            AG            = NULL,
                            deposition    = NULL,
                            survival      = NULL,
                            population_deposition = NULL) {

  new_scen <- scenario  # copy original list

  # update scalar fields
  if (!is.null(pathogen_name))     new_scen$pathogen <- pathogen_name
  if (!is.null(exposure_group))    new_scen$exposure_group <- exposure_group
  if (!is.null(microenvironment))  new_scen$microenvironment <- microenvironment
  if (!is.null(season))            new_scen$season <- season
  if (!is.null(activity_template)) new_scen$activity_template <- activity_template
  if (!is.null(activity_intensity)) new_scen$activity_intensity <- activity_intensity
  if (!is.null(population_deposition)) new_scen$population_deposition <- population_deposition

  # update distribution functions (only when explicitly provided)
  if (!is.null(C_air_dist)) new_scen$C_air_dist <- C_air_dist
  if (!is.null(IR_dist))    new_scen$IR_dist <- IR_dist
  if (!is.null(t_exp))      new_scen$t_exp_dist <- if (is.function(t_exp)) t_exp else function(n) rep(t_exp, n)
  if (!is.null(AG)) {
    if (is.function(AG)) new_scen$AG_dist <- AG
    else new_scen$AG_dist <- function(n) rep(AG, n)
  }
  if (!is.null(deposition)) {
    if (is.function(deposition)) new_scen$deposition_dist <- deposition
    else new_scen$deposition_dist <- function(n) rep(deposition, n)
  }
  if (!is.null(survival)) {
    if (is.function(survival)) new_scen$survival_dist <- survival
    else new_scen$survival_dist <- function(n) rep(survival, n)
  }

  # If exposure-related parameters are changed, one would ideally re-load
  # IR/t_exp from the database. Here we keep it simple: if the user does not
  # provide new distribution functions, the originals are kept.
  # Users can override by supplying custom functions.

  class(new_scen) <- c("airqmra_scenario", "list")
  new_scen
}

# ---- Methods for 1D simulation results ----

#' @export
print.airqmra_result <- function(x, ...) {
  cat("airqmra 1D Monte Carlo risk vector (", length(x), " iterations)\n", sep = "")
  print(summary(x))
  invisible(x)
}

#' @export
summary.airqmra_result <- function(object, ...) {
  # Already returns numeric summary; no need to override
  NextMethod()
}

#' @export
plot.airqmra_result <- function(x, ...) {
  # simple density plot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- data.frame(Risk = x)
    ggplot2::ggplot(df, ggplot2::aes(x = Risk)) +
      ggplot2::geom_density(fill = "steelblue", alpha = 0.4) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "Single-exposure risk distribution", x = "Risk", y = "Density") +
      ggplot2::theme_light()
  } else {
    hist(x, main = "Single-exposure risk", xlab = "Risk", col = "steelblue")
  }
}

# ---- Methods for 2D simulation results ----

#' @export
print.airqmra_2d <- function(x, ...) {
  cat("airqmra 2D Monte Carlo results\n")
  cat("  Uncertainty iterations:", ncol(x$risk_matrix), "\n")
  cat("  Variability iterations per U sample:", nrow(x$risk_matrix), "\n")
  if (!is.null(x$annual_summary)) {
    cat("  Annual risk summary (mean across variability):\n")
    print(x$annual_summary, row.names = FALSE)
  }
  invisible(x)
}
