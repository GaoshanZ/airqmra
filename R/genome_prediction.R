# =============================================================================
# Genome-nformed dose-esponse parameter prediction
# =============================================================================

#' Compute adhesion gene density from VFDB annotation
#'
#' Reads a DIAMOND/ABRicate VFDB annotation table and a genome FASTA file,
#' counts unique adhesion-elated genes, and returns the density per million
#' base pairs (genes / Mb).  The VFDB file may have a header line or not.
#'
#' @param vfdb_file Path to the VFDB annotation file.
#' @param genome_fasta Path to the genome FASTA file.
#' @param keyword_pattern Regular expression for adhesion genes.
#' @param min_identity Minimum percent identity.
#' @param max_evalue Maximum E-value.
#' @param has_header Logical; does the VFDB file contain a header line?
#'   If \code{TRUE} (default), column names are taken from the first line and
#'   must match the standard DIAMOND output names.  If \code{FALSE}, columns
#'   are assumed to be in the standard DIAMOND order and named automatically.
#' @param col_sseqid,col_stitle,col_pident,col_evalue Column names to use when
#'   \code{has_header = TRUE}.  Ignored when \code{has_header = FALSE}.
#'
#' @return A single numeric value: adhesion gene density (genes / Mb).
#' @export
compute_adhesion_density <- function(
    vfdb_file,
    genome_fasta,
    keyword_pattern = "adhesin|fimbri|pilus|intimin|afimbrial",
    min_identity    = 40,
    max_evalue      = 1e-5,
    has_header      = TRUE,
    col_sseqid      = "sseqid",
    col_stitle      = "stitle",
    col_pident      = "pident",
    col_evalue      = "evalue") {

  # ---- read VFDB table ----
  if (has_header) {
    vfdb <- read.table(vfdb_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")
    needed <- c(col_sseqid, col_stitle, col_pident, col_evalue)
    missing_cols <- setdiff(needed, names(vfdb))
    if (length(missing_cols) > 0) {
      stop("The following required columns are missing: ",
           paste(missing_cols, collapse = ", "),
           ".  Try setting has_header = FALSE or specify the correct column names.")
    }
  } else {
    # Standard DIAMOND output order (no header)
    vfdb <- read.table(vfdb_file, header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "",
                       col.names = c("qseqid","sseqid","pident","length",
                                     "mismatch","gapopen","qstart","qend",
                                     "sstart","send","evalue","bitscore","stitle"))
    col_sseqid <- "sseqid"
    col_stitle <- "stitle"
    col_pident <- "pident"
    col_evalue <- "evalue"
  }

  # ---- filter by quality ----
  vfdb_filt <- subset(vfdb,
                      get(col_pident) >= min_identity &
                        get(col_evalue) <= max_evalue)

  # ---- search adhesion keywords ----
  is_adhesion <- grepl(keyword_pattern, vfdb_filt[[col_stitle]], ignore.case = TRUE)
  adhesion_hits <- vfdb_filt[is_adhesion, , drop = FALSE]

  # ---- count unique adhesion genes ----
  adhesion_genes <- unique(adhesion_hits[[col_sseqid]])
  n_adhesion <- length(adhesion_genes)

  # ---- calculate genome length in Mb ----
  # simple FASTA length counter (no external dependency)
  genome_length_bp <- 0
  con <- file(genome_fasta, "r")
  on.exit(close(con))
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (!startsWith(line, ">")) {
      genome_length_bp <- genome_length_bp + nchar(line)
    }
  }
  genome_mb <- genome_length_bp / 1e6

  # ---- density ----
  adhesion_per_mb <- n_adhesion / genome_mb
  return(adhesion_per_mb)
}


#' Predict alpha parameter from adhesion density
#'
#' @param adhesion_per_mb Numeric. Adhesion gene density (genes / Mb).
#' @param use_db_N50 Logical. If \code{TRUE} (default), the N50 distribution
#'   parameters are taken from the built-n \code{pathogen_db} for
#'   \emph{Escherichia coli}.  If \code{FALSE}, a hard-coded default N50 is
#'   used (see details).
#' @param N50_record Integer. Which \emph{E. coli} record to use when
#'   multiple records exist in the database (default 1).
#'
#' @return A list with \code{alpha_est}, prediction interval, log-normal
#'   approximation parameters, model information, and \code{N50_default}
#'   (containing \code{value}, \code{log_mean}, \code{log_sd}).
#' @export
predict_alpha_from_genome <- function(adhesion_per_mb,
                                      use_db_N50 = TRUE,
                                      N50_record = 1L) {
  # ---- load pre-trained model ----
  model_path <- system.file("extdata", "model_adhesion_10strains.rds",
                            package = "airqmra")
  if (!file.exists(model_path)) {
    stop("Model file not found. Please reinstall the 'airqmra' package.")
  }
  fit <- readRDS(model_path)

  # ---- predict alpha ----
  newdata <- data.frame(adhesion_per_mb = adhesion_per_mb)
  pred <- predict(fit, newdata, interval = "prediction", level = 0.95)
  alpha_est <- as.numeric(pred[1, "fit"])
  alpha_lwr <- as.numeric(pred[1, "lwr"])
  alpha_upr <- as.numeric(pred[1, "upr"])

  sigma <- summary(fit)$sigma
  log_mean <- log(max(alpha_est, 0.0001))
  log_sd   <- sigma / alpha_est

  # ---- model info ----
  sm <- summary(fit)
  r_squared <- sm$r.squared
  rmse <- sigma
  spearman_rho <- 0.759   # from your training

  # ---- N50 default (from database or hard-coded) ----
  if (use_db_N50) {
    rec <- pathogen_info("escherichia coli")
    # keep only records that have N50 (approx or exact beta-Poisson)
    rec_n50 <- rec[!is.na(rec$N50), , drop = FALSE]
    if (nrow(rec_n50) == 0) {
      stop("No Escherichia coli record with N50 found in the database.")
    }
    if (N50_record > nrow(rec_n50)) {
      stop("N50_record exceeds available records (", nrow(rec_n50), ").")
    }
    chosen <- rec_n50[N50_record, ]
    N50_value <- chosen$N50
    # Use uncertainty parameters if available, otherwise estimate a log_sd
    if (!is.na(chosen$N50_log_mean) && !is.na(chosen$N50_log_sd)) {
      N50_log_mean_ln <- chosen$N50_log_mean   # already natural log
      N50_log_sd_ln   <- chosen$N50_log_sd
    } else {
      N50_log_mean_ln <- log(N50_value)
      # rough estimate: assume 50% coefficient of variation on log scale
      N50_log_sd_ln   <- 0.5
    }
  } else {
    # fallback hard-coded values (original training set)
    N50_value       <- 494083.1          # 10^5.6938
    N50_log_mean_ln <- log(N50_value)
    N50_log_sd_ln   <- 1.8476 * log(10)  # converted from log10 sd
  }

  list(
    alpha_est   = alpha_est,
    alpha_lwr   = alpha_lwr,
    alpha_upr   = alpha_upr,
    log_mean    = log_mean,
    log_sd      = log_sd,
    model_info  = list(
      name         = "Adhesion-density linear model",
      R_squared    = r_squared,
      RMSE         = rmse,
      Spearman_rho = spearman_rho
    ),
    N50_default = list(
      value    = N50_value,
      log_mean = N50_log_mean_ln,
      log_sd   = N50_log_sd_ln
    )
  )
}

#' Create a dose-response generator from genome-predicted parameters
#'
#' Converts the output of \code{\link{predict_alpha_from_genome}} into a
#' function that can be passed as \code{dr_generator} to
#' \code{\link{new_scenario}}.  The generator samples alpha and N50 from
#' log-normal distributions using the predicted log-mean and log-sd.
#'
#' @param pred A list returned by \code{\link{predict_alpha_from_genome}}.
#'
#' @return A function with signature \code{function(n)} that returns a list
#'   with elements \code{type = "approx_beta_poisson"}, \code{alpha} (vector
#'   of length \code{n}), and \code{N50} (vector of length \code{n}).
#' @export
#'
#' @examples
#' \dontrun{
#' pred <- predict_alpha_from_genome(adhesion_dens)
#' gen <- make_prediction_generator(pred)
#' scen <- new_scenario(pathogen_name = "E. coli predicted",
#'                      exposure_group = "Adult",
#'                      dr_generator = gen, ...)
#' }
make_prediction_generator <- function(pred) {
  force(pred)
  function(n) {
    alpha_samp <- rlnorm(n, meanlog = pred$log_mean, sdlog = pred$log_sd)
    N50_samp   <- rlnorm(n, meanlog = pred$N50_default$log_mean,
                         sdlog   = pred$N50_default$log_sd)
    list(
      type  = "approx_beta_poisson",
      alpha = alpha_samp,
      N50   = N50_samp
    )
  }
}
