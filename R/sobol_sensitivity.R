#' Sobol' variance-based global sensitivity analysis (variability + uncertainty)
#'
#' Performs a Sobol' sensitivity analysis on single-exposure (or annual) infection risk.
#' The input variables include exposure factors (C, IR, t, survival, deposition) and,
#' when available, the dose-response parameters (k for exponential; alpha and N50 for
#' approximate beta-Poisson; alpha and beta for exact beta-Poisson) with their
#' uncertainty distributions taken from the pathogen database.
#'
#' @param scenario An `airqmra_scenario` object created with `use_uncertainty = FALSE`
#'   (point-estimate mode) or `use_uncertainty = TRUE` (uncertainty mode). The function
#'   will automatically extract the necessary information from the scenario.
#' @param n Integer. Base sample size for Sobol' design. Total model evaluations
#'   will be `n * (p + 2)` where `p` is the total number of input variables.
#' @param order Integer. Order of Sobol' indices to compute (1 for first-order,
#'   2 for first + total, 3 for first + second + total). Default 2.
#' @param nboot Integer. Number of bootstrap replicates for confidence intervals.
#' @param conf Numeric. Confidence level for intervals (0-1). Default 0.95.
#' @param frequency Optional. If provided, the model function returns annual risk
#'   instead of single-exposure risk (using `annual_risk()`). Can be a fixed integer
#'   or a function.
#'
#' @return An object of class `sobol` (from the `sensitivity` package) containing
#'   the estimated Sobol' indices and their bootstrap confidence intervals.
#' @export
#'
#' @examples
#' # Legionella outdoor scenario with parameter uncertainty
#' scen <- new_scenario(pathogen_id = "legionella pneumophila",
#'                      exposure_group = "Adult",
#'                      microenvironment = "Outdoor",
#'                      population_deposition = "Male",
#'                      use_uncertainty = TRUE)
#' sobol_res <- sobol_sensitivity(scen, n = 500, order = 2, nboot = 100)
#' print(sobol_res)
#' plot(sobol_res)
#'
#' # Annual risk Sobol' analysis
#' sobol_annual <- sobol_sensitivity(scen, n = 500, frequency = 20)
#' plot(sobol_annual)
sobol_sensitivity <- function(scenario, n = 1000, order = 2, nboot = 100, conf = 0.95,
                              frequency = NULL) {
  if (!requireNamespace("sensitivity", quietly = TRUE)) {
    stop("Package 'sensitivity' is required for Sobol' analysis. Please install it.")
  }

  # ---- 1. define exposure variables ----
  var_exp <- c("C", "IR", "t", "survival", "deposition")
  p_exp <- length(var_exp)

  # ---- 2. obtain pathogen info and model type ----
  pathogen_id <- scenario$pathogen_id
  if (is.null(pathogen_id)) {
    stop("Scenario must have a pathogen_id to identify dose-response parameters.")
  }
  rec <- pathogen_info(pathogen_id)
  if (nrow(rec) == 0) stop("Pathogen not found in database.")
  # if multiple records, use the first one (simplified here; can be aligned with scenario selection)
  chosen_row <- rec[1, ]
  model_type <- chosen_row$model   # directly use model type from database
  if (is.null(model_type) || is.na(model_type)) {
    stop("Model type not found in pathogen database.")
  }

  # ---- 3. extract parameter uncertainty information ----
  param_info <- list()   # each element: name, log_mean, log_sd
  if (model_type == "exponential") {
    if (!is.na(chosen_row$k_log_mean) && !is.na(chosen_row$k_log_sd)) {
      param_info[[1]] <- list(name = "k",
                              log_mean = chosen_row$k_log_mean,
                              log_sd   = chosen_row$k_log_sd)
    }
  } else if (model_type == "exact_beta_poisson") {
    if (!is.na(chosen_row$alpha_log_mean) && !is.na(chosen_row$alpha_log_sd)) {
      param_info[[length(param_info)+1]] <- list(name = "alpha",
                                                 log_mean = chosen_row$alpha_log_mean,
                                                 log_sd   = chosen_row$alpha_log_sd)
    }
    if (!is.na(chosen_row$beta_log_mean) && !is.na(chosen_row$beta_log_sd)) {
      param_info[[length(param_info)+1]] <- list(name = "beta",
                                                 log_mean = chosen_row$beta_log_mean,
                                                 log_sd   = chosen_row$beta_log_sd)
    }
  } else if (model_type == "approx_beta_poisson") {
    if (!is.na(chosen_row$alpha_log_mean) && !is.na(chosen_row$alpha_log_sd)) {
      param_info[[length(param_info)+1]] <- list(name = "alpha",
                                                 log_mean = chosen_row$alpha_log_mean,
                                                 log_sd   = chosen_row$alpha_log_sd)
    }
    if (!is.na(chosen_row$N50_log_mean) && !is.na(chosen_row$N50_log_sd)) {
      param_info[[length(param_info)+1]] <- list(name = "N50",
                                                 log_mean = chosen_row$N50_log_mean,
                                                 log_sd   = chosen_row$N50_log_sd)
    }
  } else {
    stop("Unsupported dose-response model type: ", model_type)
  }

  # all input variable names
  var_names <- c(var_exp, sapply(param_info, `[[`, "name"))
  p_total <- length(var_names)

  # ---- 4. generate Sobol' design matrices ----
  X1 <- matrix(NA_real_, nrow = n, ncol = p_total)
  X2 <- matrix(NA_real_, nrow = n, ncol = p_total)
  colnames(X1) <- colnames(X2) <- var_names

  # exposure variable sampling
  X1[, 1] <- scenario$C_air_dist(n)
  X1[, 2] <- scenario$IR_dist(n)
  X1[, 3] <- scenario$t_exp_dist(n)
  X1[, 4] <- scenario$survival_dist(n)
  X1[, 5] <- scenario$deposition_dist(n)

  X2[, 1] <- scenario$C_air_dist(n)
  X2[, 2] <- scenario$IR_dist(n)
  X2[, 3] <- scenario$t_exp_dist(n)
  X2[, 4] <- scenario$survival_dist(n)
  X2[, 5] <- scenario$deposition_dist(n)

  # parameter uncertainty sampling (log-normal)
  for (j in seq_along(param_info)) {
    col_idx <- p_exp + j
    X1[, col_idx] <- rlnorm(n, meanlog = param_info[[j]]$log_mean,
                            sdlog = param_info[[j]]$log_sd)
    X2[, col_idx] <- rlnorm(n, meanlog = param_info[[j]]$log_mean,
                            sdlog = param_info[[j]]$log_sd)
  }

  # ---- 5. model function ----
  model_func <- function(X) {
    dose <- X[, 1] * X[, 2] * X[, 3] * X[, 4] * X[, 5]   # AG = 1

    if (model_type == "exponential") {
      k_vec <- X[, "k"]
      single <- 1 - exp(-k_vec * dose)
    } else if (model_type == "exact_beta_poisson") {
      alpha_vec <- X[, "alpha"]
      beta_vec  <- X[, "beta"]
      single <- exact_beta_poisson(dose, alpha_vec, beta_vec)
    } else if (model_type == "approx_beta_poisson") {
      alpha_vec <- X[, "alpha"]
      N50_vec   <- X[, "N50"]
      single <- approx_beta_poisson(dose, alpha_vec, N50_vec)
    } else {
      stop("Unknown model type.")
    }

    if (!is.null(frequency)) {
      annual_risk(single, frequency = frequency)
    } else {
      single
    }
  }

  # ---- 6. run Sobol' analysis ----
  sensitivity::sobol(model = model_func, X1 = X1, X2 = X2,
                     order = order, nboot = nboot, conf = conf)
}
