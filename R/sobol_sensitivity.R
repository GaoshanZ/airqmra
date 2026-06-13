#' Sobol' variance-based global sensitivity analysis (variability + uncertainty)
#'
#' Performs a Sobol' sensitivity analysis on single-exposure (or annual) infection
#' risk. The input variables include exposure factors (C, IR, t, survival,
#' deposition) and, when available, the dose-response parameters (k for
#' exponential; alpha and N50 for approximate beta-Poisson; alpha and beta for
#' exact beta-Poisson) with their uncertainty distributions taken from the
#' pathogen database.
#'
#' @param scenario An `airqmra_scenario` object created with
#'   `use_uncertainty = FALSE` (point-estimate mode) or `use_uncertainty = TRUE`
#'   (uncertainty mode).
#' @param n Integer. Base sample size for the Sobol' design. Default 1000.
#' @param order Integer. Order of Sobol' indices to compute when `method = "sobol"`
#'   (1 for first-order + total, 2 for first + second, 3 for first + second +
#'   third + total). Ignored when `method = "sobol2007"`. Default 2.
#' @param nboot Integer. Number of bootstrap replicates for confidence intervals.
#'   Set to 0 to skip bootstrap and avoid potential instability with extremely
#'   low risk values.
#' @param conf Numeric. Confidence level for intervals (0-1). Default 0.95.
#' @param frequency Optional. If provided, the model function returns annual risk
#'   instead of single-exposure risk (using `annual_risk()`).
#' @param method Character. Sobol' algorithm to use. `"sobol2007"` (default)
#'   computes first-order and total-effect indices directly and is recommended
#'   for most applications. `"sobol"` uses the older algorithm and allows
#'   specifying higher-order interactions via `order`.
#'
#' @return An object of class `"sobol2007"` or `"sobol"` (from the `sensitivity`
#'   package) containing the estimated indices and, if bootstrap was successful,
#'   bootstrap confidence intervals.
#' @export
sobol_sensitivity <- function(scenario,
                              n = 1000,
                              order = 2,
                              nboot = 100,
                              conf = 0.95,
                              frequency = NULL,
                              method = c("sobol2007", "sobol")) {
  method <- match.arg(method)

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
  chosen_row <- rec[1, ]
  model_type <- chosen_row$model
  if (is.null(model_type) || is.na(model_type)) {
    stop("Model type not found in pathogen database.")
  }

  # ---- 3. extract parameter uncertainty information ----
  param_info <- list()
  if (model_type == "exponential") {
    if (!is.na(chosen_row$k_log_mean) && !is.na(chosen_row$k_log_sd)) {
      param_info[[1]] <- list(name = "k",
                              log_mean = chosen_row$k_log_mean,
                              log_sd   = chosen_row$k_log_sd)
    } else if (!is.na(chosen_row$k)) {
      param_info[[1]] <- list(name = "k",
                              log_mean = log(chosen_row$k),
                              log_sd   = 0.01)
    }
  } else if (model_type == "exact_beta_poisson") {
    if (!is.na(chosen_row$alpha_log_mean) && !is.na(chosen_row$alpha_log_sd)) {
      param_info[[length(param_info)+1]] <- list(name = "alpha",
                                                 log_mean = chosen_row$alpha_log_mean,
                                                 log_sd   = chosen_row$alpha_log_sd)
    } else if (!is.na(chosen_row$alpha)) {
      param_info[[length(param_info)+1]] <- list(name = "alpha",
                                                 log_mean = log(chosen_row$alpha),
                                                 log_sd   = 0.01)
    }
    if (!is.na(chosen_row$beta_log_mean) && !is.na(chosen_row$beta_log_sd)) {
      param_info[[length(param_info)+1]] <- list(name = "beta",
                                                 log_mean = chosen_row$beta_log_mean,
                                                 log_sd   = chosen_row$beta_log_sd)
    } else if (!is.na(chosen_row$beta)) {
      param_info[[length(param_info)+1]] <- list(name = "beta",
                                                 log_mean = log(chosen_row$beta),
                                                 log_sd   = 0.01)
    }
  } else if (model_type == "approx_beta_poisson") {
    if (!is.na(chosen_row$alpha_log_mean) && !is.na(chosen_row$alpha_log_sd)) {
      param_info[[length(param_info)+1]] <- list(name = "alpha",
                                                 log_mean = chosen_row$alpha_log_mean,
                                                 log_sd   = chosen_row$alpha_log_sd)
    } else if (!is.na(chosen_row$alpha)) {
      param_info[[length(param_info)+1]] <- list(name = "alpha",
                                                 log_mean = log(chosen_row$alpha),
                                                 log_sd   = 0.01)
    }
    if (!is.na(chosen_row$N50_log_mean) && !is.na(chosen_row$N50_log_sd)) {
      param_info[[length(param_info)+1]] <- list(name = "N50",
                                                 log_mean = chosen_row$N50_log_mean,
                                                 log_sd   = chosen_row$N50_log_sd)
    } else if (!is.na(chosen_row$N50)) {
      param_info[[length(param_info)+1]] <- list(name = "N50",
                                                 log_mean = log(chosen_row$N50),
                                                 log_sd   = 0.01)
    }
  } else {
    stop("Unsupported dose-response model type: ", model_type)
  }

  var_names <- c(var_exp, sapply(param_info, `[[`, "name"))
  p_total <- length(var_names)

  # ---- 4. generate Sobol' design matrices ----
  X1 <- matrix(NA_real_, nrow = n, ncol = p_total)
  X2 <- matrix(NA_real_, nrow = n, ncol = p_total)
  colnames(X1) <- colnames(X2) <- var_names

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

  for (j in seq_along(param_info)) {
    col_idx <- p_exp + j
    X1[, col_idx] <- rlnorm(n, meanlog = param_info[[j]]$log_mean,
                            sdlog = param_info[[j]]$log_sd)
    X2[, col_idx] <- rlnorm(n, meanlog = param_info[[j]]$log_mean,
                            sdlog = param_info[[j]]$log_sd)
  }

  # ---- 5. model function (uses colnames(X) for parameter detection) ----
  model_func <- function(X) {
    dose <- X[, "C"] * X[, "IR"] * X[, "t"] * X[, "survival"] * X[, "deposition"]
    cn <- colnames(X)   # available column names

    if (model_type == "exponential") {
      k_vec <- if ("k" %in% cn) X[, "k"] else rep(chosen_row$k, nrow(X))
      single <- 1 - exp(-k_vec * dose)
    } else if (model_type == "exact_beta_poisson") {
      alpha_vec <- if ("alpha" %in% cn) X[, "alpha"] else rep(chosen_row$alpha, nrow(X))
      beta_vec  <- if ("beta" %in% cn) X[, "beta"] else rep(chosen_row$beta, nrow(X))
      single <- exact_beta_poisson(dose, alpha_vec, beta_vec)
    } else if (model_type == "approx_beta_poisson") {
      alpha_vec <- if ("alpha" %in% cn) X[, "alpha"] else rep(chosen_row$alpha, nrow(X))
      N50_vec   <- if ("N50" %in% cn) X[, "N50"] else rep(chosen_row$N50, nrow(X))
      # Numerical safeguard: clip extreme parameters
      alpha_vec <- pmax(pmin(alpha_vec, 1e4), 1e-4)
      N50_vec   <- pmax(N50_vec, 1e-2)
      single <- approx_beta_poisson(dose, alpha_vec, N50_vec)
    } else {
      stop("Unknown model type.")
    }
    single <- pmax(pmin(single, 1), 0)
    if (!is.null(frequency)) {
      annual_risk(single, frequency = frequency)
    } else {
      single
    }
  }

  # ---- 6. run Sobol' analysis with bootstrap fallback ----
  run_sobol <- function(use_boot) {
    if (method == "sobol2007") {
      sensitivity::sobol2007(model = model_func, X1 = X1, X2 = X2,
                             nboot = if (use_boot) nboot else 0, conf = conf)
    } else {
      sensitivity::sobol(model = model_func, X1 = X1, X2 = X2,
                         order = order, nboot = if (use_boot) nboot else 0, conf = conf)
    }
  }

  if (nboot > 0) {
    res <- tryCatch(run_sobol(TRUE), error = function(e) NULL)
    if (is.null(res)) {
      warning("Bootstrap failed (possibly due to near-zero variance of model output). Falling back to no bootstrap.")
      res <- run_sobol(FALSE)
    }
  } else {
    res <- run_sobol(FALSE)
  }
  res
}

#' @export
tornado_plot.sobol2007 <- function(sens, title = "Sobol' Variance Decomposition",
                                   xlab = "Sobol' Index", ...) {
  # extract first and total effect data
  S <- as.data.frame(sens$S) %>% rownames_to_column("Variable") %>%
    rename(original = original)
  T <- as.data.frame(sens$T) %>% rownames_to_column("Variable") %>%
    rename(original = original)

  has_ci <- all(c("min. c.i.", "max. c.i.") %in% names(S))
  if (has_ci) {
    S <- S %>% rename(min_ci = `min. c.i.`, max_ci = `max. c.i.`)
    T <- T %>% rename(min_ci = `min. c.i.`, max_ci = `max. c.i.`)
  }

  combined <- dplyr::bind_rows(S %>% mutate(Type = "First-order"),
                               T %>% mutate(Type = "Total-effect"))

  p <- ggplot2::ggplot(combined, ggplot2::aes(x = reorder(Variable, original), y = original, fill = Type)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("First-order" = "steelblue", "Total-effect" = "coral")) +
    ggplot2::labs(title = title,
                  subtitle = if (has_ci) "First-order and Total-effect Indices (95% CI)" else "First-order and Total-effect Indices (no bootstrap CI)",
                  x = "", y = xlab) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )

  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, min_ci), ymax = max_ci),
                                    position = ggplot2::position_dodge(0.7), width = 0.15, color = "grey30")
  }
  p
}

#' @export
tornado_plot.sobol <- function(sens, title = "Sobol' Variance Decomposition",
                               xlab = "Sobol' Index", ...) {
  tornado_plot.sobol2007(sens, title, xlab, ...)
}
