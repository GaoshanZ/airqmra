#' Import and fit concentration data
#'
#' Reads concentration data from a CSV or TXT file or a data frame, fits a
#' distribution, and returns a random number generator function that can be
#' directly passed to \code{C_air_dist} in \code{\link{new_scenario}}.
#'
#' @param x A file path (character) or a data frame containing concentration data.
#' @param col Column name or index. If \code{NULL} (default), the first numeric
#'   column is used automatically.
#' @param method Fitting method:
#'   \itemize{
#'     \item \code{"empirical"}: bootstrap sampling from the observed values.
#'     \item \code{"lognormal"}: fit a lognormal distribution.
#'     \item \code{"auto"}: automatically select the best among lognormal,
#'           Gamma, and Weibull based on AIC.
#'   }
#' @param sep Field separator (only used when \code{x} is a file path). Defaults
#'   to \code{","} for CSV files; for \code{.txt} files it defaults to \code{"\\t"}.
#'   You can override this manually.
#' @param ... Additional arguments passed to \code{read.table}.
#'
#' @return A function with signature \code{function(n)} that generates \code{n}
#'   random concentration values (copies/m^3). The function has an attribute
#'   \code{"distribution"} containing fitted distribution details.
#' @export
#'
#' @examples
#' # From a data frame (in-memory)
#' mydata <- data.frame(conc = c(50, 42, 67, 31))
#' conc_gen <- import_concentration(mydata, col = "conc", method = "lognormal")
#' conc_gen(5)
#'
#' # Empirical bootstrap
#' conc_gen_emp <- import_concentration(mydata, method = "empirical")
#' print(conc_gen_emp)
import_concentration <- function(x, col = NULL, method = c("empirical", "lognormal", "auto"),
                                 sep = NULL, ...) {
  method <- match.arg(method)

  # 1. Read data
  if (is.character(x)) {
    if (is.null(sep)) {
      if (grepl("\\.txt$", x, ignore.case = TRUE)) {
        sep <- "\t"
        message("Detected .txt file, using tab as separator.")
      } else {
        sep <- ","
      }
    }
    data <- read.table(x, header = TRUE, sep = sep, stringsAsFactors = FALSE, ...)
  } else if (is.data.frame(x)) {
    data <- x
  } else {
    stop("'x' must be a file path or a data frame.")
  }

  # 2. Extract concentration column
  if (is.null(col)) {
    num_cols <- sapply(data, is.numeric)
    if (!any(num_cols)) stop("No numeric column found in the data.")
    col <- names(which(num_cols))[1]
    message("Using column '", col, "' as concentration.")
  } else {
    if (is.numeric(col)) col <- names(data)[col]
    if (!col %in% names(data)) stop("Column '", col, "' not found in the data.")
  }

  conc <- data[[col]]
  conc <- na.omit(conc)
  if (any(conc < 0)) stop("Concentration values must be non-negative.")
  if (length(conc) < 3) stop("Need at least 3 concentration values for fitting.")

  # 3. Distribution fitting and generator construction
  dist_info <- list()

  if (method == "empirical") {
    f <- function(n) sample(conc, size = n, replace = TRUE)
    dist_info$type <- "empirical"
    dist_info$data <- conc
  } else if (method == "lognormal") {
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

  attr(f, "distribution") <- dist_info
  class(f) <- c("airqmra_concentration", class(f))
  return(f)
}

#' Print method for concentration generator
#' @param x An object of class `airqmra_concentration`.
#' @param ... Further arguments passed to or from other methods.
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
