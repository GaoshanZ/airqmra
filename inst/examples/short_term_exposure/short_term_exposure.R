# =============================================================================
# airqmra Short-Term Exposure Risk Assessment Example
# This script reads concentration data from the built-in example dataset,
# builds outdoor short-term exposure scenarios with user-defined activity
# durations, runs 1D Monte Carlo simulations (with parameter uncertainty),
# and produces publication-quality comparison plots.
#
# Output figures are saved in the current working directory.
# =============================================================================

# ------------------------- User-adjustable parameters -----------------------
# Pathogen ID (must match the airqmra database)
pathogen_id       <- "staphylococcus aureus"   # also "escherichia coli", "pseudomonas aeruginosa"
# Column name in the concentration file
conc_col_name     <- "S_aureus"                # also "E_coli", "P_aeruginosa"
# Display name for figure titles
display_name      <- "Staphylococcus aureus"
# Prefix for output image files
file_prefix       <- "Saureus_short"
# Target population (must exist in ir_long_table)
population_group  <- "Adult"
# Outdoor activity durations (hours) for each intensity level
# Keys = activity intensity (must match ir_short_table activity_intensity)
# Values = vectors of exposure durations (hours)
outdoor_durations <- list(
  "Sedentary"          = c(0.5, 1.0, 2.0),
  "Light_intensity"    = c(0.5, 1.0, 2.0),
  "Moderate_intensity" = c(0.5, 1.0, 2.0)
)
# Number of times each activity is repeated per year (for annual risk)
ANNUAL_REPEATS    <- 30
# Number of 1D Monte Carlo iterations
N_ITER_1D         <- 500
# Risk threshold for reference line (comment out to remove)
THRESHOLD         <- 1e-4
# ----------------------------------------------------------------------------

# Load required packages
library(airqmra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(fitdistrplus)
library(scales)

# ---------- Data loading ----------------------------------------------------
# Locate the example concentration file shipped with the package
data_file <- system.file("examples/short_term_exposure/data/Concentration.txt",
                         package = "airqmra")
if (data_file == "") stop("Example data not found. Please reinstall the airqmra package.")

conc_data <- read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
conc_vec <- conc_data[[conc_col_name]]
conc_vec <- conc_vec[conc_vec > 0]
if (length(conc_vec) < 3) stop("Not enough concentration values to fit a distribution.")

# Fit lognormal distribution to the concentration data
fit <- fitdistrplus::fitdist(conc_vec, "lnorm")
mu_conc    <- fit$estimate["meanlog"]
sigma_conc <- fit$estimate["sdlog"]
cat(conc_col_name, "concentration lognormal fit: meanlog =", mu_conc,
    "sdlog =", sigma_conc, "\n")

# Concentration sampling function
C_fun <- function(n) rlnorm(n, meanlog = mu_conc, sdlog = sigma_conc)

# ---------- Build scenarios and run 1D simulations for each duration ---------
risk_df <- data.frame()

for (act in names(outdoor_durations)) {
  durations <- outdoor_durations[[act]]
  for (dur in durations) {
    scenario_label <- paste0(act, " (", dur, "h)")

    # Obtain short-term inhalation rate for the chosen population and activity
    ir_short <- get_ir_short(population_group, act)
    if (!is.null(ir_short)) {
      # Convert L/min -> m^3/h (1 L/min = 0.06 m^3/h)
      IR_dist <- function(n) rlnorm(n,
                                    meanlog = ir_short$log_mean + log(0.06),
                                    sdlog   = ir_short$log_sd)
    } else {
      # Fallback to long-term average inhalation rate
      warning("No short-term IR for '", population_group, "' doing '", act,
              "', falling back to long-term average.")
      ir_long <- get_ir_long(population_group)
      IR_dist <- function(n) rlnorm(n, meanlog = ir_long$log_mean - log(24),
                                    sdlog   = ir_long$log_sd)
    }

    # Build point-estimate scenario (parameter uncertainty enabled)
    scen <- new_scenario(
      pathogen_id      = pathogen_id,
      exposure_group   = population_group,
      microenvironment = "Outdoor",
      C_air_dist       = C_fun,
      IR_dist          = IR_dist,
      t_exp            = dur,
      use_uncertainty  = TRUE          # include parameter uncertainty
    )

    set.seed(123)
    single <- run_qmra_1d(scen, n_iter = N_ITER_1D)
    annual <- annual_risk(single, frequency = ANNUAL_REPEATS)

    temp <- data.frame(
      Activity      = act,
      Duration_h    = dur,
      Scenario      = scenario_label,
      SingleRisk    = single,
      AnnualRisk    = annual
    )
    risk_df <- rbind(risk_df, temp)

    cat(sprintf("%-30s: single median = %.2e, annual median = %.2e\n",
                scenario_label, median(single), median(annual)))
  }
}

# ---------- Summary statistics (median & 95% CI) ----------------------------
summary_df <- risk_df %>%
  group_by(Activity, Duration_h, Scenario) %>%
  summarise(
    Median_single  = median(SingleRisk),
    CI_low_single  = quantile(SingleRisk, 0.025),
    CI_high_single = quantile(SingleRisk, 0.975),
    Median_annual  = median(AnnualRisk),
    CI_low_annual  = quantile(AnnualRisk, 0.025),
    CI_high_annual = quantile(AnnualRisk, 0.975),
    .groups = "drop"
  )

# =============================================================================
# Publication-quality figures
# =============================================================================

# Custom journal theme
theme_journal <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "grey30", linewidth = 0.8),
      axis.ticks         = element_line(color = "grey30"),
      axis.text          = element_text(color = "black"),
      axis.title         = element_text(face = "bold", size = 13),
      plot.title         = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle      = element_text(size = 11, color = "grey40", hjust = 0.5),
      legend.position    = "right",
      legend.title       = element_blank()
    )
}

# ---- Single-exposure risk comparison (pointrange by activity and duration) ----
p_single <- ggplot(summary_df, aes(x = factor(Duration_h), y = Median_single,
                                   color = Activity, group = Activity)) +
  geom_pointrange(aes(ymin = CI_low_single, ymax = CI_high_single),
                  size = 1.2, linewidth = 1.0,
                  position = position_dodge(width = 0.5)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title = paste(display_name, "Short-term Exposure Risk\n(", population_group, ")"),
       subtitle = "Single exposure risk vs. duration",
       x = "Exposure duration (hours)",
       y = "Single infection probability (log scale)") +
  theme_journal()
ggsave(paste0(file_prefix, "_single_risk_by_duration.png"), p_single,
       width = 10, height = 6, dpi = 600)
print(p_single)

# ---- Annual risk comparison (pointrange by activity and duration) ----------
p_annual <- ggplot(summary_df, aes(x = factor(Duration_h), y = Median_annual,
                                   color = Activity, group = Activity)) +
  geom_pointrange(aes(ymin = CI_low_annual, ymax = CI_high_annual),
                  size = 1.2, linewidth = 1.0,
                  position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = THRESHOLD, linetype = "dashed",
             color = "#d73027", linewidth = 0.7, alpha = 0.8) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title = paste(display_name, "Annual Risk by Duration\n(", population_group, ")"),
       subtitle = paste("Annual repeats:", ANNUAL_REPEATS),
       x = "Exposure duration (hours)",
       y = "Annual infection probability (log scale)") +
  theme_journal()
ggsave(paste0(file_prefix, "_annual_risk_by_duration.png"), p_annual,
       width = 10, height = 6, dpi = 600)
print(p_annual)

cat("\nShort-term exposure analysis complete. Figures saved in the current working directory.\n")
