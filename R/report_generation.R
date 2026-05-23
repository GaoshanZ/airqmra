#' Generate a standardized QMRA report (multi-format, bilingual)
#'
#' This function generates an infection risk assessment report based on
#' the T/CSES XXXX-XXXX technical guideline, supporting HTML, Word, and PDF
#' output, with Chinese or English text.
#'
#' @param scenario An `airqmra_scenario` object.
#' @param risk_results A named list containing analysis results.
#' @param output_file Character. Name of the output file (extension will be
#'   adjusted automatically).
#' @param output_format Character. One of `"html"`, `"word"`, or `"pdf"`.
#' @param language Character. `"zh"` (Chinese) or `"en"` (English).
#' @param ... Additional arguments passed to `rmarkdown::render`.
#'
#' @return Invisibly returns the path to the generated report.
#' @export
generate_risk_report <- function(scenario, risk_results,
                                 output_file = "airqmra_report.html",
                                 output_format = c("html", "word", "pdf"),
                                 language = c("zh", "en"), ...) {
  output_format <- match.arg(output_format)
  language <- match.arg(language)

  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is required. Please install it.")

  # Choose template file
  template_name <- paste0("skeleton_", language, ".Rmd")
  template_path <- system.file(
    file.path("rmarkdown", "templates", "risk_report", "skeleton", template_name),
    package = "airqmra"
  )
  if (template_path == "") {
    stop("Report template not found for language '", language, "'. Please reinstall 'airqmra'.")
  }

  # Set output format
  fmt <- switch(output_format,
                html = rmarkdown::html_document(css = system.file("rmarkdown/templates/risk_report/style.css", package = "airqmra")),
                word = rmarkdown::word_document(),
                pdf  = rmarkdown::pdf_document(latex_engine = "xelatex", toc = FALSE)
  )

  # Adjust output file extension
  ext <- switch(output_format, html = ".html", word = ".docx", pdf = ".pdf")
  if (!grepl(paste0("\\", ext, "$"), output_file)) {
    output_file <- paste0(tools::file_path_sans_ext(output_file), ext)
  }

  rmarkdown::render(
    input = template_path,
    output_format = fmt,
    output_file = output_file,
    params = list(scenario = scenario, results = risk_results),
    envir = new.env(parent = globalenv()),
    ...
  )

  message("Report generated: ", normalizePath(output_file))
  invisible(output_file)
}
