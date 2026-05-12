#' Create an inhalation QMRA exposure scenario
#'
#' @param pathogen_name Character. A descriptive name for the scenario. If not provided but
#'   `pathogen_id` is given, `pathogen_id` will be used as the name.
#' @param pathogen_id Character. English ID of the pathogen as in the built-in database
#'   (e.g., "legionella pneumophila", "escherichia coli"). If provided and `dr_pars` is NULL,
#'   dose-response parameters will be loaded automatically from the database.
#' @param C_air_dist Function. A function that generates air concentration samples (copies/m³).
#' @param IR_dist Function. A function that generates inhalation rate samples (m³/h).
#' @param t_exp Single numeric or function for exposure duration (hours).
#' @param dr_pars List. Dose-response parameters (e.g., `list(type = "exponential", k = 0.0054)`).
#'   If `NULL`, the function will load parameters from the built-in database using `pathogen_id`.
#' @param dr_choice Integer or character. Choice of which parameter set to use when the pathogen
#'   has multiple records in the database. An integer selects the row number (from the
#'   `pathogen_info()` result); a character string can be `"exponential"` or
#'   `"exact_beta_poisson"` to select the first record of that model type. If `NULL` and
#'   multiple records exist, the function will ask interactively in an R session, or raise an
#'   error in non-interactive scripts.
#'
#' @return An object of class `airqmra_scenario`.
#' @export
new_scenario <- function(pathogen_name = NULL,
                         pathogen_id   = NULL,
                         C_air_dist    = function(n) stats::rlnorm(n, meanlog = log(50), sdlog = log(1.8)),
                         IR_dist       = function(n) stats::rnorm(n, mean = 0.7, sd = 0.1),
                         t_exp         = 1.0,
                         dr_pars       = NULL,
                         dr_choice     = NULL) {

  # --- 如果手动提供了 dr_pars，优先使用 ---
  if (!is.null(dr_pars)) {
    if (!is.null(pathogen_id)) {
      message("Both 'pathogen_id' and 'dr_pars' provided. Using 'dr_pars' and ignoring 'pathogen_id'.")
    }
  } else {
    # --- dr_pars 为空，必须提供 pathogen_id ---
    if (is.null(pathogen_id)) {
      stop("Either 'dr_pars' or 'pathogen_id' must be provided.")
    }

    # 查询数据库
    rec <- pathogen_info(pathogen_id)
    if (nrow(rec) == 0) {
      stop("No pathogen found for '", pathogen_id, "'. Check spelling or use 'pathogen_info()' to explore the database.")
    }

    # --- 处理多套参数的选择逻辑 ---
    if (nrow(rec) == 1) {
      chosen <- 1L
      if (!is.null(dr_choice) && dr_choice != 1) {
        warning("Only one parameter set exists; ignoring 'dr_choice'.")
      }
    } else {
      # 有多条记录
      if (is.null(dr_choice)) {
        if (interactive()) {
          # 创建清晰的菜单项
          labels <- paste0(
            "[", seq_len(nrow(rec)), "] ",
            rec$model,
            ifelse(rec$model == "exponential",
                   paste0(" (k=", format(rec$k, digits=3), ")"),
                   paste0(" (alpha=", rec$alpha, ", beta=", format(rec$beta, digits=3), ")")
            )
          )
          cat("\nMultiple parameter sets found for '", pathogen_id, "':\n", sep="")
          cat("Choose a number and press Enter (0 to abort):\n")
          choice_idx <- menu(labels, title = "Select a parameter set by number:")
          if (choice_idx == 0) {
            stop("Selection aborted by user.")
          }
          chosen <- choice_idx
        } else {
          # 非交互环境：报错并列出选项
          msg <- paste0(
            "Multiple parameter sets found for '", pathogen_id, "'. ",
            "Please specify 'dr_choice' as an integer (1-", nrow(rec), ") or model type.\n",
            "Available options:\n",
            paste(capture.output(print(rec[, c("model","alpha","beta","N50","k")])), collapse="\n")
          )
          stop(msg)
        }
      } else {
        # dr_choice 已提供
        if (is.numeric(dr_choice)) {
          idx <- as.integer(dr_choice)
          if (idx < 1 || idx > nrow(rec)) {
            stop("'dr_choice' must be between 1 and ", nrow(rec), ".")
          }
          chosen <- idx
        } else if (is.character(dr_choice)) {
          # 按模型名匹配，取第一个匹配的
          model_match <- which(tolower(rec$model) == tolower(dr_choice))
          if (length(model_match) == 0) {
            stop("No record with model '", dr_choice, "'. Available models: ",
                 paste(unique(rec$model), collapse=", "))
          }
          chosen <- model_match[1]
        } else {
          stop("'dr_choice' must be an integer index or a model name (e.g., 'exponential').")
        }
      }
    }

    # --- 从选定的行构建 dr_pars ---
    chosen_row <- rec[chosen, ]
    model <- chosen_row$model
    dr_pars <- list(type = model)

    if (model == "exponential") {
      if (is.na(chosen_row$k)) stop("Parameter 'k' is missing in the database for this record.")
      dr_pars$k <- chosen_row$k
    } else if (model == "exact_beta_poisson") {
      if (is.na(chosen_row$alpha) || is.na(chosen_row$beta))
        stop("Parameters 'alpha' or 'beta' are missing in the database for this record.")
      dr_pars$alpha <- chosen_row$alpha
      dr_pars$beta  <- chosen_row$beta
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
