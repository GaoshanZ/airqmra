# 暴露参数智能查询函数（单维度 + 回退）

#' Internal: match a classification in a table with fallback
#'
#' @param db Data frame, must contain column 'classification'.
#' @param group Character. Desired classification, e.g. "Adult", "18Y-44Y", "Urban".
#' @param fallback Character vector. Fallback classifications in priority order.
#'   Defaults to c("Adult").
#'
#' @return The first matched row, or NULL if nothing matches.
#' @keywords internal
match_exposure <- function(db, group, fallback = c("Adult")) {
  # 直接匹配
  idx <- which(db$classification == group)
  if (length(idx) > 0) return(db[idx[1], , drop = FALSE])

  # 模糊匹配（包含关键词）
  idx <- grep(group, db$classification, ignore.case = TRUE)
  if (length(idx) > 0) return(db[idx[1], , drop = FALSE])

  # 按回退列表顺序尝试
  for (fb in fallback) {
    idx <- which(db$classification == fb)
    if (length(idx) > 0) {
      message("No exact match for '", group, "', using fallback '", fb, "'.")
      return(db[idx[1], , drop = FALSE])
    }
    idx <- grep(fb, db$classification, ignore.case = TRUE)
    if (length(idx) > 0) {
      message("No exact match for '", group, "', using fallback '", fb, "'.")
      return(db[idx[1], , drop = FALSE])
    }
  }

  NULL
}


#' Get long‑term inhalation rate (daily average, m3/d)
#'
#' @param group Character. Population classification (e.g., "Adult", "18Y-44Y",
#'   "Male", "Urban"). Single dimension only.
#' @param fallback Character vector. Fallback classifications in priority order.
#'   Default: c("Adult").
#'
#' @return A list with elements: dist, log_mean, log_sd, unit, ref.
#'   Returns NULL if no record found.
#' @export
#'
#' @examples
#' get_ir_long("18Y-44Y")
#' get_ir_long("Beijing", fallback = c("Urban", "Adult"))
get_ir_long <- function(group, fallback = c("Adult")) {
  rec <- match_exposure(airqmra::ir_long_table, group, fallback)
  if (is.null(rec)) return(NULL)
  list(
    dist     = rec$IR_dist,
    log_mean = rec$IR_log_mean,
    log_sd   = rec$IR_log_sd,
    unit     = rec$IR_unit,
    ref      = rec$ref
  )
}


#' Get short‑term inhalation rate for a specific activity intensity (L/min)
#'
#' @param group Character. Population classification.
#' @param activity Character. Activity intensity, e.g. "Sleeping",
#'   "Sedentary", "Light-intensity", "Moderate-intensity",
#'   "High-intensity", "Very_High-intensity".
#' @param fallback Character vector. Fallback classifications.
#'
#' @return A list with elements: dist, log_mean, log_sd, unit, ref.
#' @export
get_ir_short <- function(group, activity, fallback = c("Adult")) {
  rec <- match_exposure(airqmra::ir_short_table, group, fallback)
  if (is.null(rec)) return(NULL)
  # 再按活动强度筛选
  rec <- rec[rec$activity_intensity == activity, , drop = FALSE]
  if (nrow(rec) == 0) {
    warning("Activity '", activity, "' not found for group '", group, "'.")
    return(NULL)
  }
  list(
    dist     = rec$IR_dist,
    log_mean = rec$IR_log_mean,
    log_sd   = rec$IR_log_sd,
    unit     = rec$IR_unit,
    ref      = rec$ref
  )
}


#' Get time‑activity pattern (minutes per day) for a population group
#'
#' @param group Character. Population classification.
#' @param microenvironment Character. "Indoor" or "Outdoor".
#' @param season Character. One of "Annual", "Spring & Autumn", "Summer", "Winter".
#' @param fallback Character vector. Fallback classifications.
#'
#' @return A list with elements: dist, log_mean, log_sd, unit, ref.
#' @export
get_time <- function(group, microenvironment = "Indoor", season = "Annual",
                     fallback = c("Adult")) {
  rec <- match_exposure(airqmra::time_table, group, fallback)
  if (is.null(rec)) return(NULL)
  # 再按季节和微环境筛选
  rec <- rec[rec$season == season & rec$microenvironment == microenvironment, , drop = FALSE]
  if (nrow(rec) == 0) {
    # 如果指定季节没有，尝试 Annual
    if (season != "Annual") {
      rec <- match_exposure(airqmra::time_table, group, fallback)
      if (!is.null(rec)) {
        rec <- rec[rec$season == "Annual" & rec$microenvironment == microenvironment, , drop = FALSE]
        if (nrow(rec) > 0) {
          message("Season '", season, "' not found, using 'Annual'.")
        }
      }
    }
  }
  if (nrow(rec) == 0) {
    warning("No time data found for group '", group, "', microenvironment '",
            microenvironment, "', season '", season, "'.")
    return(NULL)
  }
  list(
    dist     = rec$time_dist,
    log_mean = rec$time_log_mean,
    log_sd   = rec$time_log_sd,
    unit     = rec$time_unit,
    ref      = rec$ref
  )
}
