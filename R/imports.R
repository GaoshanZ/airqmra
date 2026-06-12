#' @importFrom stats cor median quantile rlnorm reorder IQR na.omit rgamma runif rweibull sd
#' @importFrom graphics barplot hist lines text
#' @importFrom utils capture.output menu read.csv read.table
#' @importFrom utils globalVariables
#' @importFrom stats predict
NULL

# Declare non-standard evaluation variables used in ggplot2 aesthetics
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Variable", "Correlation",
    ".x", "Dose", "PRCC", "Risk", "Scenario",
    "Var1", "Var2", "ci_high", "ci_low", "mu_star"
  ))
}
