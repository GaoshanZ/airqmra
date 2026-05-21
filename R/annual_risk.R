#' Calculate annual infection risk from single-exposure risks
#'
#' @param single_risk A numeric vector of single-exposure infection probabilities
#'   (e.g., output from \code{\link{run_qmra_1d}}).
#' @param frequency Either a fixed integer (number of exposures per year) or
#'   a function that takes an integer \code{n} and returns \code{n} random
#'   frequency values (e.g., \code{function(n) rpois(n, 240)}).
#' @param method Character. Method to compute annual risk:
#'   \code{"exact"} uses \eqn{1 - (1 - p)^k} for each pair of probability and
#'   frequency; \code{"approx"} uses the Poisson approximation
#'   \eqn{1 - e^{-k \cdot p}} (valid for small \eqn{p} and large \eqn{k}).
#'
#' @return A numeric vector of annual infection risks, with class \code{airqmra_result}
#'   for compatibility with summary and plot methods.
#' @export
#'
#' @examples
#' # Single exposure risk from a scenario
#' scen <- new_scenario(pathogen_id = "legionella pneumophila",
#'                      exposure_group = "Adult",
#'                      population_deposition = "Male")
#' single_risk <- run_qmra_1d(scen, n_iter = 5000)
#'
#' # Annual risk with fixed 240 exposures per year
#' annual <- annual_risk(single_risk, frequency = 240)
#' summary(annual)
#' hist(annual, breaks = 50)
#'
#' # Annual risk with Poisson-distributed frequency
#' freq_func <- function(n) rpois(n, lambda = 240)
#' annual_pois <- annual_risk(single_risk, frequency = freq_func)
#' summary(annual_pois)
annual_risk <- function(single_risk, frequency, method = c("exact", "approx")) {
  method <- match.arg(method)

  if (!is.numeric(single_risk)) {
    stop("'single_risk' must be a numeric vector.")
  }

  n <- length(single_risk)

  # 生成频率向量
  if (is.function(frequency)) {
    k <- frequency(n)   # 每个模拟对应的暴露次数
  } else if (is.numeric(frequency) && length(frequency) == 1) {
    k <- rep(frequency, n)
  } else {
    stop("'frequency' must be a single integer or a function taking 'n' and returning n values.")
  }

  # 计算年风险
  if (method == "exact") {
    annual <- 1 - (1 - single_risk)^k
  } else if (method == "approx") {
    annual <- 1 - exp(-k * single_risk)
  }

  # 确保落在 [0,1]
  annual <- pmax(0, pmin(1, annual))

  class(annual) <- c("airqmra_result", "numeric")
  annual
}
