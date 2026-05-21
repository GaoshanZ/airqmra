#' Calculate disease burden from annual infection risk
#'
#' Computes disease burden (DB) as the product of annual infection risk,
#' the probability of illness given infection, and a health burden factor.
#'
#' @param annual_risk Numeric vector of annual infection probabilities
#'   (e.g., output from \code{\link{annual_risk}}). Values should be between 0 and 1.
#' @param HB Numeric or function. Health burden factor (e.g., DALY per case).
#'   Default \code{0.0026} (point estimate). Can also be a function that generates
#'   random values, e.g., \code{function(n) rlnorm(n, log(0.0026), 0.1)}.
#' @param P_ill_inf Numeric or function. Probability of illness given infection.
#'   Default \code{0.2} (point estimate). Can also be a distribution function.
#'
#' @return A numeric vector of disease burden values (same length as \code{annual_risk}),
#'   with class \code{airqmra_result} for compatibility with summary and plot methods.
#' @export
#'
#' @examples
#' # First, obtain single exposure risk
#' scen <- new_scenario(pathogen_id = "legionella pneumophila",
#'                      exposure_group = "Adult",
#'                      population_deposition = "Male")
#' single <- run_qmra_1d(scen, n_iter = 2000)
#'
#' # Calculate annual risk (e.g., 240 exposures per year)
#' annual <- annual_risk(single, frequency = 240)
#'
#' # Disease burden with default point estimates
#' db <- disease_burden(annual)
#' summary(db)
#' hist(db, breaks = 50)
#'
#' # Using distributions for HB and P_ill_inf
#' db_dist <- disease_burden(annual,
#'                           HB = function(n) rlnorm(n, log(0.0026), 0.1),
#'                           P_ill_inf = function(n) rbeta(n, 2, 8))
#' summary(db_dist)
disease_burden <- function(annual_risk, HB = 0.0026, P_ill_inf = 0.2) {
  if (!is.numeric(annual_risk)) {
    stop("'annual_risk' must be a numeric vector.")
  }

  n <- length(annual_risk)

  # 处理 HB：可以是数值或函数
  if (is.function(HB)) {
    hb_val <- HB(n)
  } else if (is.numeric(HB) && length(HB) == 1) {
    hb_val <- rep(HB, n)
  } else {
    stop("'HB' must be a single numeric value or a function taking 'n' and returning n values.")
  }

  # 处理 P_ill_inf
  if (is.function(P_ill_inf)) {
    p_ill <- P_ill_inf(n)
  } else if (is.numeric(P_ill_inf) && length(P_ill_inf) == 1) {
    p_ill <- rep(P_ill_inf, n)
  } else {
    stop("'P_ill_inf' must be a single numeric value or a function taking 'n' and returning n values.")
  }

  # 计算疾病负担
  db <- annual_risk * p_ill * hb_val

  # 确保非负（理论上都非负）
  db <- pmax(0, db)

  class(db) <- c("airqmra_result", "numeric")
  db
}
