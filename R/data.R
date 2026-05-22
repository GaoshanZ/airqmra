#' Pathogen Dose-Response Database
#'
#' A dataset containing dose-response model parameters and uncertainty
#' distributions for various pathogens.
#'
#' @format A data frame with columns for pathogen identifiers, model type,
#'   point estimates (alpha, beta, N50, k), uncertainty distribution parameters
#'   (log-mean and log-sd for each parameter), exposure pathways, response type,
#'   and literature reference.
#' @source QMRA Wiki (https://qmrawiki.msu.edu) and published literature.
"pathogen_db"

#' Long-term Inhalation Rate Table
#'
#' Daily average inhalation rates (m^3/d) by population classification.
#'
#' @format A data frame with classification, IR_dist, IR_log_mean, IR_log_sd,
#'   IR_unit, and ref columns.
"ir_long_table"

#' Short-term Inhalation Rate Table
#'
#' Short-term inhalation rates (L/min) by population and activity.
#'
#' @format A data frame with classification, activity_intensity, IR_dist,
#'   IR_log_mean, IR_log_sd, IR_unit, and ref columns.
"ir_short_table"

#' Time-Activity Table
#'
#' Time spent indoors or outdoors by population, season, and microenvironment.
#'
#' @format A data frame with classification, season, microenvironment,
#'   time_dist, time_log_mean, time_log_sd, time_unit, and ref columns.
"time_table"

#' Deposition Fraction Database
#'
#' Pulmonary deposition fractions by population and site.
#'
#' @format A data frame with population, site, dist, log_mean, log_sd, and
#'   deposition_site columns.
"deposition_db"

#' Survival Fraction Database
#'
#' Airborne pathogen survival parameters.
#'
#' @format A data frame with pathogen_id, dist, log_mean, log_sd, value, and
#'   ref columns.
"survival_db"

#' Activity Time Budget Templates
#'
#' Pre-defined daily activity time budgets (min/d) for typical population groups.
#'
#' @format A data frame with template_id, activity intensity columns (Sleeping,
#'   Sedentary, Light_intensity, Moderate_intensity, High_intensity,
#'   Very_High_intensity), and time_unit column.
"activity_templates"
