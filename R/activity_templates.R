#' Get daily activity time budget from a template
#'
#' Each template provides typical minutes per day for different activity intensities.
#' The values are fixed recommended durations, not distribution parameters.
#'
#' @param id Character. Template ID, e.g. "office_worker", "student".
#'
#' @return A list with elements:
#'   \item{minutes}{Named numeric vector of activity intensities and their daily minutes.}
#'   \item{time_unit}{Character, unit of time (e.g. "min/d").}
#' @export
#'
#' @examples
#' get_activity_budget("office_worker")
get_activity_budget <- function(id = "office_worker") {
  db <- airqmra::activity_templates

  if (!id %in% db$template_id) {
    stop("Template '", id, "' not found. Available: ",
         paste(db$template_id, collapse = ", "))
  }

  row <- db[db$template_id == id, ]

  # 提取活动强度列（除了 template_id 和 time_unit）
  act_cols <- setdiff(names(row), c("template_id", "time_unit"))
  minutes   <- as.numeric(row[1, act_cols])
  names(minutes) <- act_cols

  list(
    minutes   = minutes,
    time_unit = row$time_unit
  )
}
