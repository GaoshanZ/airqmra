#' @importFrom stats cor median quantile rlnorm reorder IQR na.omit rgamma runif rweibull sd predict
#' @importFrom graphics barplot hist lines text
#' @importFrom utils capture.output menu read.csv read.table globalVariables
#' @importFrom dplyr %>% mutate rename bind_rows
#' @importFrom tibble rownames_to_column
#' @importFrom grid unit
#' @importFrom rlang .data
NULL

# Declare non-standard evaluation variables used in ggplot2 aesthetics
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Variable", "Correlation",
    ".x", "Dose", "PRCC", "Risk", "Scenario",
    "Var1", "Var2", "ci_high", "ci_low", "mu_star",
    "original", "min. c.i.", "max. c.i.", "Type", "min_ci", "max_ci"
  ))
}
