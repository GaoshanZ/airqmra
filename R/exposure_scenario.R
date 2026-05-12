#' Create an inhalation QMRA exposure scenario
#'
#' @param pathogen_name Character. A descriptive name for the scenario (e.g., "Office outbreak").
#'   If not provided but `pathogen_id` is given, `pathogen_id` will be used as the name.
#' @param pathogen_id Character. English ID of the pathogen as in the built-in database
#'   (e.g., "legionella pneumophila", "escherichia coli"). If provided and `dr_pars` is NULL,
#'   dose-response parameters will be loaded automatically from the database.
#' @param C_air_dist A function that generates air concentration samples (copies/m³).
#' @param IR_dist A function that generates inhalation rate samples (m³/h).
#' @param t_exp Single numeric or function for exposure duration (h).
#' @param dr_pars List of dose-response parameters (e.g., `list(type = "exponential", k = 0.0054)`).
#'   If `NULL`, the function will try to load parameters from the built-in database using `pathogen_id`.
#'
#' @return An object of class `airqmra_scenario`.
#' @export
new_scenario <- function(pathogen_name = NULL,
                         pathogen_id   = NULL,
                         C_air_dist    = function(n) stats::rlnorm(n, meanlog = log(50), sdlog = log(1.8)),
                         IR_dist       = function(n) stats::rnorm(n, mean = 0.7, sd = 0.1),
                         t_exp         = 1.0,
                         dr_pars       = NULL) {

  # --- 处理剂量-反应参数 ---
  if (!is.null(dr_pars)) {
    # 用户手动提供了 dr_pars，直接使用
    if (!is.null(pathogen_id)) {
      message("Both 'pathogen_id' and 'dr_pars' provided. Using 'dr_pars' and ignoring 'pathogen_id'.")
    }
  } else {
    # dr_pars 为空，尝试从数据库获取
    if (is.null(pathogen_id)) {
      stop("Either 'dr_pars' or 'pathogen_id' must be provided.")
    }

    # 查询数据库
    rec <- pathogen_info(pathogen_id)
    if (nrow(rec) == 0) {
      stop("No pathogen found for '", pathogen_id, "'. Check spelling or use 'pathogen_info()' to explore the database.")
    }

    # 如果有多条匹配，默认选用第一条并提示
    if (nrow(rec) > 1) {
      message("Multiple records found for '", pathogen_id, "'. Using the first one:\n",
              "  model: ", rec$model[1],
              ", alpha: ", rec$alpha[1],
              ", beta: ", rec$beta[1],
              ", N50: ", rec$N50[1],
              ", k: ", rec$k[1])
      rec <- rec[1, ]
    }

    # 根据模型类型构建 dr_pars
    model <- rec$model
    dr_pars <- list(type = model)

    if (model == "exponential") {
      if (is.na(rec$k)) stop("Parameter 'k' is missing in the database for this pathogen.")
      dr_pars$k <- rec$k
    } else if (model == "exact_beta_poisson") {
      if (is.na(rec$alpha) || is.na(rec$beta)) stop("Parameters 'alpha' or 'beta' are missing in the database.")
      dr_pars$alpha <- rec$alpha
      dr_pars$beta  <- rec$beta
    } else {
      stop("Unknown model type in database: ", model)
    }
  }

  # --- 处理场景名称 ---
  if (is.null(pathogen_name)) {
    if (!is.null(pathogen_id)) {
      pathogen_name <- pathogen_id
    } else {
      pathogen_name <- "Unnamed scenario"
    }
  }

  # --- 构建场景对象 ---
  scenario <- list(
    pathogen     = pathogen_name,
    pathogen_id  = pathogen_id,
    C_air_dist   = C_air_dist,
    IR_dist      = IR_dist,
    t_exp_dist   = if (is.function(t_exp)) t_exp else function(n) rep(t_exp, n),
    dr_pars      = dr_pars
  )
  class(scenario) <- c("airqmra_scenario", "list")
  scenario
}




#' Run 1D Monte Carlo QMRA simulation
#'
#' @param scenario An `airqmra_scenario` object.
#' @param n_iter Number of Monte Carlo iterations.
#' @return A numeric vector of infection risk with class `airqmra_result`.
#' @export
run_qmra_1d <- function(scenario, n_iter = 1e5) {
  # Generate random samples for all input variables
  C   <- scenario$C_air_dist(n_iter)
  IR  <- scenario$IR_dist(n_iter)
  t   <- scenario$t_exp_dist(n_iter)

  # Calculate inhaled dose
  dose <- C * IR * t

  # Extract dose-response type and parameters
  dr  <- scenario$dr_pars
  type <- dr$type
  dr_args <- dr[names(dr) != "type"]  # all parameters except 'type'

  # Dispatch to dose_response
  risk <- do.call(dose_response, c(list(dose = dose, type = type), dr_args))

  class(risk) <- c("airqmra_result", "numeric")
  risk
}
