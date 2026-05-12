#' Create an inhalation QMRA exposure scenario
#'
#' @param pathogen_name Character. Name of the pathogen (e.g., "SARS-CoV-2").
#' @param C_air_dist A function that generates air concentration samples (copies/m³).
#' @param IR_dist A function that generates inhalation rate samples (m³/h).
#' @param t_exp Single numeric or function for exposure duration (h).
#' @param dr_pars List of dose-response parameters.
#'
#' @return An object of class `airqmra_scenario`.
#' @export
new_scenario <- function(pathogen_name,
                         C_air_dist = function(n) stats::rlnorm(n, meanlog = log(50), sdlog = log(1.8)),
                         IR_dist    = function(n) stats::rnorm(n, mean = 0.7, sd = 0.1),
                         t_exp      = 1.0,
                         dr_pars    = list(type = "exponential", k = 0.0054)) {
  scenario <- list(
    pathogen     = pathogen_name,
    C_air_dist   = C_air_dist,
    IR_dist      = IR_dist,
    t_exp_dist   = if (is.function(t_exp)) t_exp else function(n) rep(t_exp, n),
    dr_pars      = dr_pars
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

  # Calculate inhaled dose
  dose <- C * IR * t

  # Extract dose-response type and parameters
  dr  <- scenario$dr_pars
  type <- dr$type
  dr_args <- dr[names(dr) != "type"]  # all parameters except 'type'

  # Dispatch to dose_response
  risk <- do.call(dose_response, c(list(dose = dose, type = type), dr_args))

  class(risk) <- c("airqmra_result", "numeric")
  risk
}
