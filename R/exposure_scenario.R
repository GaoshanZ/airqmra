#' Create an inhalation QMRA exposure scenario
#'
#' @param pathogen_name Character. A descriptive name for the scenario.
#'   If not provided, `pathogen_id` or a generic name will be used.
#' @param pathogen_id Character. English ID of the pathogen as in the built-in database.
#'   If `dr_pars` is NULL, dose-response parameters will be loaded from the database.
#' @param C_air_dist Function. A function that generates air concentration samples (copies/m^3).
#'   Defaults to a lognormal distribution with GM=50, GSD=1.8.
#' @param IR_dist Function. A function that generates inhalation rate samples (m^3/h).
#'   If not provided manually and `exposure_group` is given, it will be loaded automatically
#'   from the built-in exposure database (converted from m^3/d to m^3/h).
#' @param t_exp Single numeric or function for exposure duration (hours).
#'   If not provided manually and `exposure_group` is given, it will be loaded automatically
#'   from the time-activity database (converted from min/d to hours) or from an activity template.
#' @param dr_pars List. Dose-response parameters (e.g., `list(type = "exponential", k = 0.0054)`).
#'   If NULL, parameters will be loaded from the pathogen database using `pathogen_id`.
#'   Ignored when `use_uncertainty = TRUE` (the parameter generator takes over).
#' @param dr_choice Integer or character. Choice of parameter set when multiple records exist.
#' @param exposure_group Character. Population classification for exposure parameters
#'   (e.g., "Adult", "Male", "18Y-44Y"). If provided, `IR_dist` and `t_exp` will be
#'   overwritten with database values unless manually specified.
#' @param microenvironment Character. "Indoor" or "Outdoor". Default "Indoor".
#' @param season Character. "Annual", "Spring & Autumn", "Summer", "Winter". Default "Annual".
#' @param activity_template Character. Activity template ID (e.g., "office_worker").
#'   If provided together with `activity_intensity`, the exposure duration will be set to
#'   the fixed minutes of that intensity from the template (converted to hours).
#' @param activity_intensity Character. Activity intensity (e.g., "Light_intensity").
#'   Used with `activity_template` to set a fixed exposure duration.
#' @param AG Numeric or function. Fraction of deposited pathogens that initiate infection.
#'   Default is 1 (no additional reduction). Accepts a fixed value or a distribution function.
#' @param deposition Numeric or function. Fraction of inhaled particles that deposit in the lung.
#'   If NULL (default), loaded automatically from the deposition database based on
#'   `population_deposition` (site "PUL").
#' @param survival Numeric or function. Fraction of airborne pathogens that remain viable upon inhalation.
#'   If NULL (default), loaded automatically from the survival database.
#' @param population_deposition Character. Population group for deposition data lookup.
#'   One of "Male", "Female", "Child", "Senior". Default "Male".
#' @param use_uncertainty Logical. If TRUE, dose-response parameters (alpha, beta, N50, k) are sampled
#'   from their uncertainty distributions (log-normal) stored in the pathogen database.
#'   If FALSE (default), point estimates are used.
#' @param fallback Character vector. Fallback classifications for exposure data,
#'   in priority order. Default `c("Adult")`.
#'
#' @return An object of class `airqmra_scenario`.
#' @export
new_scenario <- function(pathogen_name = NULL,
                         pathogen_id   = NULL,
                         C_air_dist    = function(n) stats::rlnorm(n, meanlog = log(50), sdlog = log(1.8)),
                         IR_dist       = NULL,
                         t_exp         = NULL,
                         dr_pars       = NULL,
                         dr_choice     = NULL,
                         exposure_group = NULL,
                         microenvironment = "Indoor",
                         season        = "Annual",
                         activity_template = NULL,
                         activity_intensity = NULL,
                         AG            = 1.0,
                         deposition    = NULL,
                         survival      = NULL,
                         population_deposition = "Male",
                         use_uncertainty = FALSE,
                         fallback      = c("Adult")) {

  # ---------- dose-response parameter handling ----------
  if (!is.null(dr_pars)) {
    if (!is.null(pathogen_id)) {
      message("Both 'pathogen_id' and 'dr_pars' provided. Using 'dr_pars' and ignoring 'pathogen_id'.")
    }
  } else {
    if (is.null(pathogen_id)) {
      stop("Either 'dr_pars' or 'pathogen_id' must be provided.")
    }
    rec <- pathogen_info(pathogen_id)
    if (nrow(rec) == 0) {
      stop("No pathogen found for '", pathogen_id, "'.")
    }
    # multi-parameter selection logic
    if (nrow(rec) == 1) {
      chosen <- 1L
    } else {
      if (is.null(dr_choice)) {
        if (interactive()) {
          labels <- paste0(
            "[", seq_len(nrow(rec)), "] ",
            rec$model,
            ifelse(rec$model == "exponential",
                   paste0(" (k=", format(rec$k, digits=3), ")"),
                   paste0(" (alpha=", rec$alpha, ", beta=", format(rec$beta, digits=3), ")")
            )
          )
          cat("\nMultiple parameter sets found for '", pathogen_id, "':\n", sep="")
          choice_idx <- menu(labels, title = "Select a parameter set by number:")
          if (choice_idx == 0) stop("Selection aborted by user.")
          chosen <- choice_idx
        } else {
          msg <- paste0("Multiple parameter sets found for '", pathogen_id, "'. ",
                        "Please specify 'dr_choice' as an integer or model type.\nAvailable options:\n",
                        paste(capture.output(print(rec[, c("model","alpha","beta","N50","k")])), collapse="\n"))
          stop(msg)
        }
      } else {
        if (is.numeric(dr_choice)) {
          idx <- as.integer(dr_choice)
          if (idx < 1 || idx > nrow(rec)) stop("'dr_choice' out of range.")
          chosen <- idx
        } else if (is.character(dr_choice)) {
          model_match <- which(tolower(rec$model) == tolower(dr_choice))
          if (length(model_match) == 0) stop("Model '", dr_choice, "' not found.")
          chosen <- model_match[1]
        } else {
          stop("'dr_choice' must be integer or character.")
        }
      }
    }
    chosen_row <- rec[chosen, ]
    model <- chosen_row$model

    # build parameters based on uncertainty setting
    if (use_uncertainty) {
      dr_generator <- build_dr_generator(chosen_row, use_uncertainty = TRUE)
      dr_pars <- NULL   # not used in uncertainty mode
    } else {
      dr_generator <- NULL
      dr_pars <- list(type = model)
      if (model == "exponential") {
        if (is.na(chosen_row$k)) stop("Parameter 'k' is missing.")
        dr_pars$k <- chosen_row$k
      } else if (model == "exact_beta_poisson") {
        if (!is.na(chosen_row$alpha) && !is.na(chosen_row$beta)) {
          dr_pars$alpha <- chosen_row$alpha
          dr_pars$beta  <- chosen_row$beta
        } else if (!is.na(chosen_row$alpha) && !is.na(chosen_row$N50)) {
          dr_pars$type <- "approx_beta_poisson"
          dr_pars$alpha <- chosen_row$alpha
          dr_pars$N50   <- chosen_row$N50
        } else {
          stop("Beta-Poisson model requires either 'beta' or 'N50' in the database.")
        }
      } else {
        stop("Unknown model type: ", model)
      }
    }
  }

  # ---------- exposure parameter auto-loading ----------
  user_IR <- !is.null(IR_dist)
  user_t  <- !is.null(t_exp)

  if (!is.null(exposure_group)) {
    # 1. inhalation rate (long-term, m^3/d -> m^3/h)
    if (!user_IR) {
      ir_info <- get_ir_long(exposure_group, fallback = fallback)
      if (!is.null(ir_info)) {
        IR_dist <- function(n) stats::rlnorm(n,
                                             meanlog = ir_info$log_mean - log(24),
                                             sdlog   = ir_info$log_sd)
        message("IR_dist set from database for group '", exposure_group,
                "' (converted to m^3/h).")
      }
    }

    # 2. exposure time
    if (!user_t) {
      if (!is.null(activity_template) && !is.null(activity_intensity)) {
        # use fixed duration from activity template
        budget <- get_activity_budget(activity_template)
        act_min <- budget$minutes[activity_intensity]
        if (is.na(act_min)) {
          stop("Activity intensity '", activity_intensity, "' not in template '", activity_template, "'.")
        }
        t_exp <- act_min / 60   # convert to hours
        message("Exposure duration set to ", t_exp, " hours (from template '",
                activity_template, "', intensity '", activity_intensity, "').")
      } else {
        # use time-activity database (indoor/outdoor, minute distribution)
        time_info <- get_time(exposure_group, microenvironment, season, fallback = fallback)
        if (!is.null(time_info)) {
          # minutes -> hours
          t_exp <- function(n) stats::rlnorm(n,
                                             meanlog = time_info$log_mean - log(60),
                                             sdlog   = time_info$log_sd)
          message("Exposure duration distribution set from time-activity database for group '",
                  exposure_group, "', microenvironment '", microenvironment,
                  "', season '", season, "'.")
        }
      }
    }
  }

  # ---------- survival rate loading ----------
  if (is.null(survival)) {
    surv_rec <- get_survival(pathogen_id)
    survival_dist <- switch(surv_rec$dist,
                            lognormal = function(n) stats::rlnorm(n,
                                                                  meanlog = surv_rec$log_mean,
                                                                  sdlog   = surv_rec$log_sd),
                            constant  = function(n) rep(surv_rec$value, n),
                            stop("Unsupported survival distribution: ", surv_rec$dist)
    )
    message("Survival distribution loaded from database.")
  } else if (is.function(survival)) {
    survival_dist <- survival
  } else {
    survival_dist <- function(n) rep(survival, n)
  }

  # ---------- deposition rate loading ----------
  if (is.null(deposition)) {
    dep_rec <- get_deposition(population_deposition, site = "PUL")
    if (!is.null(dep_rec)) {
      deposition_dist <- switch(dep_rec$dist,
                                lognormal = function(n) stats::rlnorm(n,
                                                                      meanlog = dep_rec$log_mean,
                                                                      sdlog   = dep_rec$log_sd),
                                stop("Unsupported deposition distribution: ", dep_rec$dist)
      )
      message("Deposition distribution loaded for population '", population_deposition, "'.")
    } else {
      deposition_dist <- function(n) rep(0.5, n)  # safe fallback
    }
  } else if (is.function(deposition)) {
    deposition_dist <- deposition
  } else {
    deposition_dist <- function(n) rep(deposition, n)
  }

  # ---------- AG fixed to 1 (no additional multiplication by default) ----------
  if (is.function(AG)) {
    AG_dist <- AG
  } else {
    AG_dist <- function(n) rep(AG, n)
  }

  # ---------- scenario naming ----------
  if (is.null(pathogen_name)) {
    if (!is.null(pathogen_id)) {
      pathogen_name <- pathogen_id
    } else {
      pathogen_name <- "Unnamed scenario"
    }
  }

  # ---------- set final defaults ----------
  if (is.null(IR_dist)) {
    IR_dist <- function(n) stats::rnorm(n, mean = 0.7, sd = 0.1)
  }
  if (is.null(t_exp)) {
    t_exp <- 1.0
  }

  # ---------- build scenario object ----------
  scenario <- list(
    pathogen         = pathogen_name,
    pathogen_id      = pathogen_id,
    C_air_dist       = C_air_dist,
    IR_dist          = IR_dist,
    t_exp_dist       = if (is.function(t_exp)) t_exp else function(n) rep(t_exp, n),
    AG_dist          = AG_dist,
    deposition_dist  = deposition_dist,
    survival_dist    = survival_dist,
    dr_pars          = dr_pars,
    dr_generator     = dr_generator   # uncertainty parameter generator (may be NULL)
  )
  class(scenario) <- c("airqmra_scenario", "list")
  scenario
}


#' Run 1D Monte Carlo QMRA simulation
#'
#' @param scenario An `airqmra_scenario` object.
#' @param n_iter Number of Monte Carlo iterations.
#' @return A numeric vector of infection risk with class `airqmra_result`.
#' @export
run_qmra_1d <- function(scenario, n_iter = 1e5) {
  # Generate random samples for all input variables
  C   <- scenario$C_air_dist(n_iter)
  IR  <- scenario$IR_dist(n_iter)
  t   <- scenario$t_exp_dist(n_iter)

  # Calculate effective dose (viable, deposited, and invasion-adjusted)
  dose <- C * IR * t *
    scenario$survival_dist(n_iter) *
    scenario$deposition_dist(n_iter) *
    scenario$AG_dist(n_iter)

  # Obtain dose-response parameters (compatible with point estimate and uncertainty generator)
  if (!is.null(scenario$dr_generator)) {
    # uncertainty mode: generator returns a list with type and parameter vectors
    dr_params <- scenario$dr_generator(n_iter)
    type <- dr_params$type
    dr_args <- dr_params[names(dr_params) != "type"]
    risk <- do.call(dose_response, c(list(dose = dose, type = type), dr_args))
  } else {
    # point estimate mode: use the fixed parameter list pre-configured in the scenario
    dr  <- scenario$dr_pars
    type <- dr$type
    dr_args <- dr[names(dr) != "type"]
    risk <- do.call(dose_response, c(list(dose = dose, type = type), dr_args))
  }

  class(risk) <- c("airqmra_result", "numeric")
  risk
}
