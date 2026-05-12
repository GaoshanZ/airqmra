#' View all dose-response parameter sets for a pathogen
#'
#' @param pathogen_id Character. English ID of the pathogen.
#' @return A data frame with available records (model, alpha, beta, N50, k).
#' @export
view_params <- function(pathogen_id) {
  rec <- pathogen_info(pathogen_id)
  if (nrow(rec) == 0) stop("Pathogen not found.")
  rec[, c("model", "alpha", "beta", "N50", "k")]
}
