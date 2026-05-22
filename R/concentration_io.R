#' Import and fit concentration data
#'
#' Reads concentration data from a CSV file or data frame, fits a distribution,
#' and returns a random number generator function that can be directly passed to
#' \code{C_air_dist} in \code{\link{new_scenario}}.
#'
#' @param x A file path (character) or a data frame containing concentration data.
#' @param col Column name or index (if data frame/file with header). If \code{NULL}
#'   (default), the first numeric column is used automatically.
#' @param method Fitting method:
#'   \itemize{
#'     \item \code{"empirical"}: bootstrap sampling from the observed values.
#'     \item \code{"lognormal"}: fit a lognormal distribution.
#'     \item \code{"auto"}: automatically select the best among lognormal,
#'           Gamma, and Weibull based on AIC.
#'   }
#' @param ... Additional arguments passed to \code{read.csv} (when \code{x} is a
#'   file path).
#'
#' @return A function with signature \code{function(n)} that generates \code{n}
#'   random concentration values (copies/m³). The function has an attribute
#'   \code{"distribution"} containing fitted distribution details.
#'
#' @examples
#' # From a CSV file (assumes column "conc")
#' conc_gen <- import_concentration("data/field_samples.csv", col = "conc")
#' scen <- new_scenario(..., C_air_dist = conc_gen)
#'
#' # From a data frame
#' mydata <- data.frame(conc = c(50, 42, 67, 31))
#' conc_gen <- import_concentration(mydata, method = "lognormal")
#' conc_gen(10)
#'
#' @export
import_concentration <- function(x, col = NULL, method = c("empirical", "lognormal", "auto"), ...) {
  method <- match.arg(method)

  # ---- 1. 读入数据 ----
  if (is.character(x)) {
    data <- read.csv(x, stringsAsFactors = FALSE, ...)
  } else if (is.data.frame(x)) {
    data <- x
  } else {
    stop("'x' must be a file path or a data frame.")
  }

  # ---- 2. 提取浓度列 ----
  if (is.null(col)) {
    # 自动寻找第一个数值列
    num_cols <- sapply(data, is.numeric)
    if (!any(num_cols)) stop("No numeric column found in the data.")
    col <- names(which(num_cols))[1]
    message("Using column '", col, "' as concentration.")
  } else {
    if (is.numeric(col)) {
      col <- names(data)[col]
    }
    if (!col %in% names(data)) stop("Column '", col, "' not found in the data.")
  }

  conc <- data[[col]]
  conc <- na.omit(conc)
  if (any(conc < 0)) stop("Concentration values must be non-negative.")
  if (length(conc) < 3) stop("Need at least 3 concentration values for fitting.")

  # ---- 3. 分布拟合与生成器构建 ----
  dist_info <- list()

  if (method == "empirical") {
    f <- function(n) sample(conc, size = n, replace = TRUE)
    dist_info$type <- "empirical"
    dist_info$data <- conc
  } else if (method == "lognormal") {
    # 要求全部为正数
    if (any(conc == 0)) stop("Lognormal fitting requires all concentrations > 0.")
    fit <- fitdistrplus::fitdist(conc, "lnorm")
    mlog <- fit$estimate["meanlog"]
    slog <- fit$estimate["sdlog"]
    f <- function(n) rlnorm(n, meanlog = mlog, sdlog = slog)
    dist_info$type <- "lognormal"
    dist_info$params <- list(meanlog = mlog, sdlog = slog)
    dist_info$AIC <- fit$aic
  } else if (method == "auto") {
    if (any(conc == 0)) stop("Auto fitting requires all concentrations > 0 for parametric distributions.")
    fits <- list()
    tryCatch({ fits$lnorm <- fitdistrplus::fitdist(conc, "lnorm") }, error = function(e) NULL)
    tryCatch({ fits$gamma <- fitdistrplus::fitdist(conc, "gamma") }, error = function(e) NULL)
    tryCatch({ fits$weibull <- fitdistrplus::fitdist(conc, "weibull") }, error = function(e) NULL)
    if (length(fits) == 0) stop("None of the distributions (lognormal, gamma, weibull) could be fitted.")

    aics <- sapply(fits, function(x) x$aic)
    best <- which.min(aics)
    best_fit <- fits[[best]]
    best_name <- names(fits)[best]

    f <- switch(best_name,
                lnorm = function(n) rlnorm(n, meanlog = best_fit$estimate["meanlog"], sdlog = best_fit$estimate["sdlog"]),
                gamma = function(n) rgamma(n, shape = best_fit$estimate["shape"], rate = best_fit$estimate["rate"]),
                weibull = function(n) rweibull(n, shape = best_fit$estimate["shape"], scale = best_fit$estimate["scale"])
    )
    dist_info$type <- best_name
    dist_info$params <- as.list(best_fit$estimate)
    dist_info$AIC <- best_fit$aic
  }

  # 将分布信息附加为函数属性
  attr(f, "distribution") <- dist_info
  class(f) <- c("airqmra_concentration", class(f))
  return(f)
}

#' Print method for concentration generator
#' @export
print.airqmra_concentration <- function(x, ...) {
  info <- attr(x, "distribution")
  cat("Concentration random generator\n")
  cat("  Fitting method: ", info$type, "\n", sep = "")
  if (info$type == "empirical") {
    cat("  Based on ", length(info$data), " observed values.\n", sep = "")
  } else {
    cat("  Parameters: ", paste(names(info$params), "=", signif(unlist(info$params), 3)), "\n", sep = "")
    if (!is.null(info$AIC)) cat("  AIC: ", round(info$AIC, 1), "\n", sep = "")
  }
  invisible(x)
}
