#' Query the built-in pathogen dose-response parameter database
#'
#' @param name Character. English or Chinese name (or partial match) of the pathogen.
#'   Case-insensitive. Example: "Legionella", "军团菌", "Staphylococcus".
#' @return A data frame (tibble) with columns: pathogen_id, pathogen_cn, model, alpha, beta,
#'   N50, k, primer_F, primer_R. Returns an empty tibble if no match found.
#'
#' @examples
#' pathogen_info("escherichia")
#' pathogen_info("军团菌")
#' @export
pathogen_info <- function(name) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function. Please install it.")
  }
  db <- airqmra::pathogen_db

  # 同时匹配英文ID和中文名
  idx <- grepl(tolower(name), tolower(db$pathogen_id)) |
    grepl(name, db$pathogen_cn, ignore.case = TRUE)

  dplyr::filter(db, idx)
}
