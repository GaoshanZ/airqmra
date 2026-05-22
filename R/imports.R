#' @importFrom stats cor median quantile rlnorm reorder
#' @importFrom graphics barplot hist
#' @importFrom utils capture.output menu
#' @importFrom utils globalVariables
NULL

# Declare non-standard evaluation variables used in ggplot2 aesthetics
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Variable", "Correlation"))
}
