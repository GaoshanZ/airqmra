#' Launch the airqmra Shiny web application
#'
#' @export
run_shiny <- function() {
  app_dir <- system.file("shiny", package = "airqmra")
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try re-installing the package.")
  }
  shiny::runApp(app_dir, display.mode = "normal")
}
