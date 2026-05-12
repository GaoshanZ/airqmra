#' Exact beta-Poisson dose-response model (safe version)
#'
#' Computes infection probability using the exact beta-Poisson model
#' via the confluent hypergeometric function, with numerical safeguards.
#'
#' @param dose Numeric vector of inhaled doses (copies or CFU). Must be >= 0.
#' @param alpha Non-negative shape parameter.
#' @param beta Non-negative shape parameter.
#'
#' @return Numeric vector of infection probabilities (0–1).
#' @examples
#' exact_beta_poisson(dose = c(0, 1, 10, 1000), alpha = 0.45, beta = 3.2)
#' @export
exact_beta_poisson <- function(dose, alpha, beta) {
  if (!requireNamespace("hypergeo", quietly = TRUE)) {
    stop("Package 'hypergeo' is required for exact beta-Poisson. Please install it.")
  }

  # 初始化概率向量
  prob <- numeric(length(dose))

  # 负剂量或零剂量直接返回 0
  prob[dose <= 0] <- 0

  # 极大剂量直接返回 1（避免数值爆炸）
  # 阈值设定：当 dose > 1e10 时，1F1 几乎为 0，感染概率 ≈ 1
  huge <- dose > 1e10
  prob[huge] <- 1

  # 对正常范围内的剂量计算精确值
  normal <- dose > 0 & !huge
  if (any(normal)) {
    # 计算 1F1，并自动处理可能的错误
    f1 <- hypergeo::genhypergeo(alpha, alpha + beta, -dose[normal])
    # 如果计算结果为 NA/NaN/Inf，或因误差导致 >1 或 <0，强制夹紧
    prob[normal] <- 1 - pmax(0, pmin(1, f1, na.rm = TRUE), na.rm = TRUE)
  }

  # 最后全局夹紧，确保所有值在 [0,1]
  prob <- pmax(0, pmin(1, prob))
  return(prob)
}


#' Dose-response model dispatcher
#'
#' Selects the appropriate dose-response model based on the type specified.
#' Currently supports:
#'   - `"exponential"`: P = 1 - exp(-k * dose)
#'   - `"exact_beta_poisson"`: P = 1 - 1F1(alpha, alpha + beta, -dose)
#'
#' @param dose Numeric vector of dose.
#' @param type Character. One of "exponential" or "exact_beta_poisson".
#' @param ... Additional parameters passed to the underlying model
#'   (e.g., `k` for exponential, `alpha` and `beta` for exact_beta_poisson).
#'
#' @return Numeric vector of infection probabilities.
#' @export
dose_response <- function(dose, type = c("exponential", "exact_beta_poisson"), ...) {
  type <- match.arg(type)
  switch(type,
         exponential = {
           args <- list(...)
           k <- args$k
           if (is.null(k)) stop("Parameter 'k' is required for exponential model.")
           1 - exp(-k * dose)
         },
         exact_beta_poisson = {
           args <- list(...)
           alpha <- args$alpha
           beta  <- args$beta
           if (is.null(alpha) || is.null(beta))
             stop("Parameters 'alpha' and 'beta' are required for exact beta-Poisson.")
           exact_beta_poisson(dose, alpha, beta)
         }
  )
}
