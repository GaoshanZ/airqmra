# 暴露参数智能查询函数（单维度 + 回退）

#' Internal: match a classification label with fallback
#'
#' @param db Data frame containing a 'classification' column.
#' @param group Character. Desired classification, e.g. "Adult", "Urban".
#' @param fallback Character vector. Fallback classification labels in priority order.
#'
#' @return A character string (the matched classification label), or NULL if nothing found.
#' @keywords internal
match_exposure_label <- function(db, group, fallback = c("Adult")) {
  # 1. 精确匹配
  if (group %in% db$classification) return(group)

  # 2. 模糊匹配（包含关键词）
  idx <- grep(group, db$classification, ignore.case = TRUE)
  if (length(idx) > 0) return(db$classification[idx[1]])

  # 3. 按回退列表尝试
  for (fb in fallback) {
    if (fb %in% db$classification) {
      message("No match for '", group, "', using fallback '", fb, "'.")
      return(fb)
    }
    idx <- grep(fb, db$classification, ignore.case = TRUE)
    if (length(idx) > 0) {
      message("No match for '", group, "', using fallback '", fb, "'.")
      return(db$classification[idx[1]])
    }
  }

  NULL
}


#' Get long‑term inhalation rate (daily average, m3/d)
#'
#' @param group Character. Population classification (e.g., "Adult", "18Y-44Y").
#' @param fallback Character vector. Fallback classifications.
#'
#' @return A list with elements: dist, log_mean, log_sd, unit, ref.
#' @export
get_ir_long <- function(group, fallback = c("Adult")) {
  label <- match_exposure_label(airqmra::ir_long_table, group, fallback)
  if (is.null(label)) return(NULL)

  rec <- airqmra::ir_long_table[airqmra::ir_long_table$classification == label, ][1, ]
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
#' @param activity Character. Activity intensity, e.g. "Sleeping", "Light_intensity".
#' @param fallback Character vector. Fallback classifications.
#'
#' @return A list with elements: dist, log_mean, log_sd, unit, ref.
#' @export
get_ir_short <- function(group, activity, fallback = c("Adult")) {
  label <- match_exposure_label(airqmra::ir_short_table, group, fallback)
  if (is.null(label)) return(NULL)

  rec <- airqmra::ir_short_table[
    airqmra::ir_short_table$classification == label &
      airqmra::ir_short_table$activity_intensity == activity,
  ]

  if (nrow(rec) == 0) {
    warning("Activity '", activity, "' not found for group '", label, "'.")
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


#' Get time‑activity pattern (minutes per day)
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
  label <- match_exposure_label(airqmra::time_table, group, fallback)
  if (is.null(label)) return(NULL)

  # 在匹配的分类下筛选季节和微环境
  rec <- airqmra::time_table[
    airqmra::time_table$classification == label &
      airqmra::time_table$microenvironment == microenvironment &
      airqmra::time_table$season == season,
  ]

  # 如果指定季节没有，回退到 Annual
  if (nrow(rec) == 0 && season != "Annual") {
    rec <- airqmra::time_table[
      airqmra::time_table$classification == label &
        airqmra::time_table$microenvironment == microenvironment &
        airqmra::time_table$season == "Annual",
    ]
    if (nrow(rec) > 0) {
      message("Season '", season, "' not found for '", label, "', using 'Annual'.")
    }
  }

  if (nrow(rec) == 0) {
    warning("No time data found for group '", label, "', microenvironment '",
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

#' Get survival fraction distribution
#'
#' @param pathogen_id Character or NULL. Pathogen ID for specific survival data.
#'   If NULL, returns the default (all-pathogen) survival distribution.
#'
#' @return A list with dist, log_mean, log_sd, value, ref.
#' @export
get_survival <- function(pathogen_id = NULL) {
  db <- airqmra::survival_db

  # 优先匹配指定病原体
  if (!is.null(pathogen_id)) {
    rec <- db[!is.na(db$pathogen_id) & db$pathogen_id == pathogen_id, ]
    if (nrow(rec) > 0) return(as.list(rec[1, ]))
  }

  # 回退到通用值（pathogen_id 为 NA）
  rec <- db[is.na(db$pathogen_id), ]
  if (nrow(rec) == 0) stop("No default survival parameters found in database.")
  as.list(rec[1, ])
}

#' Get pulmonary deposition fraction for a population group
#'
#' @param population Character. "Male", "Female", "Child", "Senior".
#' @param site Character. Deposition site code, default "PUL".
#'
#' @return A list with dist, log_mean, log_sd, deposition_site.
#' @export
get_deposition <- function(population = "Male", site = "PUL") {
  db <- airqmra::deposition_db
  rec <- db[db$population == population & db$site == site, ]
  if (nrow(rec) == 0) {
    warning("No deposition data for population '", population, "', site '", site, "'.")
    return(NULL)
  }
  list(
    dist      = rec$dist,
    log_mean  = rec$log_mean,
    log_sd    = rec$log_sd,
    deposition_site = rec$deposition_site
  )
}
