# =============================================================================
# airqmra Long-Term Exposure Risk Assessment Example
# This script reads concentration data from the built-in example dataset,
# performs 2D Monte Carlo simulation for multiple population groups,
# generates publication-quality figures, and conducts global sensitivity
# analysis with uncertainty propagation visualization.
#
# Output figures are saved in the current working directory.
# =============================================================================

# ------------------------- User-adjustable parameters -----------------------
# Pathogen ID (must match the airqmra database)
pathogen_id     <- "staphylococcus aureus"   # also "escherichia coli", "pseudomonas aeruginosa"
# Column name in the concentration file
conc_col_name   <- "S_aureus"                # also "E_coli", "P_aeruginosa"
# Prefix for output image files
file_prefix     <- "Saureus_long"
# Display name for figure titles
display_name    <- "Staphylococcus aureus"
# Population groups to compare
groups          <- c("9Y-12Y", "18Y-44Y", "60Y-79Y")
# Representative group for sensitivity analysis
sensitivity_group <- "18Y-44Y"
# 2D Monte Carlo simulation sizes (adjust according to computational resources)
N_VAR           <- 100   # variability iterations
N_UNC           <- 50    # uncertainty iterations
# Global sensitivity simulation sizes (recorded inputs, relatively smaller)
N_VAR_GLOBAL    <- 300
N_UNC_GLOBAL    <- 100
# Annual exposure frequency (internal use only, not displayed on plots)
ANNUAL_FREQ     <- 20
# Risk threshold for reference line on comparison plot
THRESHOLD       <- 1e-4
# ----------------------------------------------------------------------------

# Load required packages
library(airqmra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(fitdistrplus)
library(scales)
library(cowplot)   # for combining subplots

# ---------- Data loading ----------------------------------------------------
# Locate the example concentration file shipped with the package
data_file <- system.file("examples/long_term_exposure/data/Concentration.txt",
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

# ---------- Validate population groups --------------------------------------
valid_groups <- intersect(groups, unique(airqmra::ir_long_table$classification))
if (length(valid_groups) < length(groups)) {
  warning("The following groups were not found in the database and have been filtered out: ",
          paste(setdiff(groups, valid_groups), collapse = ", "))
}
if (length(valid_groups) == 0) stop("No valid population groups. Check the 'groups' parameter.")
cat("Assessing groups:", paste(valid_groups, collapse = ", "), "\n")

# ---------- 2D Monte Carlo simulation ---------------------------------------
res2d_list       <- list()
annual_means_list <- list()
for (grp in valid_groups) {
  cat("Running 2D simulation:", pathogen_id, "|", grp, "\n")
  scen_unc <- new_scenario(
    pathogen_id      = pathogen_id,
    exposure_group   = grp,
    microenvironment = "Outdoor",
    season           = "Annual",
    C_air_dist       = C_fun,
    use_uncertainty  = TRUE,          # enable parameter uncertainty
    dr_generator     = NULL
  )
  res2d <- run_qmra_2d(scen_unc, n_var = N_VAR, n_unc = N_UNC,
                       frequency = ANNUAL_FREQ)
  res2d_list[[grp]]      <- res2d
  annual_means_list[[grp]] <- res2d$annual_mean
}

# ---------- Annual risk summary (median & 95% CI) ---------------------------
summary_annual <- data.frame(
  Group   = valid_groups,
  Median  = sapply(annual_means_list, median),
  CI_low  = sapply(annual_means_list, quantile, probs = 0.025),
  CI_high = sapply(annual_means_list, quantile, probs = 0.975),
  row.names = NULL
)
cat("\n========== Annual risk summary (mean across variability) ==========\n")
print(summary_annual)

# ---------- Density data ----------------------------------------------------
density_data <- bind_rows(lapply(valid_groups, function(g) {
  data.frame(Group = g, AnnualRisk = res2d_list[[g]]$overall_annual_risk)
}))

# =============================================================================
# Publication-quality figures: density + comparison
# =============================================================================

# Dynamic colour palette
n_colors <- length(valid_groups)
if (n_colors <= 9) {
  pal <- RColorBrewer::brewer.pal(n_colors, "Set1")
} else {
  pal <- scales::hue_pal()(n_colors)
}
names(pal) <- valid_groups

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
      legend.position    = "bottom",
      legend.title       = element_blank()
    )
}

# ---- Density plot ----
p_density <- ggplot(density_data, aes(x = AnnualRisk, color = Group, fill = Group)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(title = paste(display_name, "Annual Infection Risk Density"),
       subtitle = "Outdoor long-term exposure",
       x = "Annual infection probability (log scale)", y = "Density") +
  theme_journal()
ggsave(paste0(file_prefix, "_annual_risk_density.png"), p_density,
       width = 10, height = 6, dpi = 600)
print(p_density)

# ---- Comparison plot (point range + reference line) ----
p_compare <- ggplot(summary_annual, aes(x = Group, y = Median, color = Group)) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high),
                  size = 4, linewidth = 1.8) +
  geom_hline(yintercept = THRESHOLD, linetype = "dashed",
             color = "#d73027", linewidth = 0.7, alpha = 0.8) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = pal) +
  labs(title = paste(display_name, "Annual Infection Risk Comparison"),
       subtitle = "Median and 95% uncertainty interval",
       x = "Population group", y = "Annual infection probability (log scale)") +
  theme_journal() +
  theme(legend.position = "none")
ggsave(paste0(file_prefix, "_risk_comparison.png"), p_compare,
       width = 8, height = 6, dpi = 600)
print(p_compare)

# =============================================================================
# Global sensitivity analysis (on the sensitivity_group)
# =============================================================================
cat("\n========== Global sensitivity analysis (", sensitivity_group, " group) ==========\n")
scen_global <- new_scenario(
  pathogen_id      = pathogen_id,
  exposure_group   = sensitivity_group,
  microenvironment = "Outdoor",
  season           = "Annual",
  C_air_dist       = C_fun,
  use_uncertainty  = TRUE,
  dr_generator     = NULL
)
set.seed(123)
res2d_global <- run_qmra_2d(scen_global, n_var = N_VAR_GLOBAL, n_unc = N_UNC_GLOBAL,
                            frequency = ANNUAL_FREQ, record_inputs = TRUE)

# ---- Global Spearman ----
sens_global_spearman <- sensitivity_global(res2d_global, method = "spearman")
cat("Global Spearman ranking:\n")
print(sens_global_spearman$ranking)
p_spearman <- tornado_plot(sens_global_spearman,
                           title = paste("Global Spearman (", display_name, ",", sensitivity_group, ")"))
ggsave(paste0(file_prefix, "_global_spearman.png"), p_spearman,
       width = 10, height = 6, dpi = 300)
print(p_spearman)

# ---- Global PRCC ----
sens_global_prcc <- sensitivity_global_prcc(res2d_global, nboot = 0)
cat("Global PRCC ranking:\n")
print(sens_global_prcc$ranking)
p_prcc <- tornado_plot(sens_global_prcc,
                       title = paste("Global PRCC (", display_name, ",", sensitivity_group, ")"))
ggsave(paste0(file_prefix, "_global_prcc.png"), p_prcc,
       width = 10, height = 6, dpi = 300)
print(p_prcc)

# ---- Global Morris ----
scen_morris <- new_scenario(
  pathogen_id      = pathogen_id,
  exposure_group   = sensitivity_group,
  microenvironment = "Outdoor",
  season           = "Annual",
  C_air_dist       = C_fun,
  use_uncertainty  = TRUE,
  dr_generator     = NULL
)
sens_morris_global <- sensitivity_global_morris(scen_morris, n_levels = 8, n_trajectories = 50)
cat("Morris mu* values:\n")
print(sens_morris_global$ee)
p_morris <- tornado_plot(sens_morris_global,
                         title = paste("Global Morris (", display_name, ",", sensitivity_group, ")"))
ggsave(paste0(file_prefix, "_global_morris.png"), p_morris,
       width = 10, height = 6, dpi = 300)
print(p_morris)

# ---- Sobol' variance decomposition ----
sobol_res <- sobol_sensitivity(scen_global, n = 100, nboot = 0)
cat("First-order indices:\n"); print(sobol_res$S)
cat("Total-effect indices:\n"); print(sobol_res$T)
p_sobol <- tornado_plot(sobol_res,
                        title = paste("Sobol' variance decomposition (", display_name, ",", sensitivity_group, ")"))
ggsave(paste0(file_prefix, "_sobol_global.png"), p_sobol,
       width = 10, height = 6, dpi = 300)
print(p_sobol)

# =============================================================================
# Parameter uncertainty propagation visualization
# =============================================================================
cat("\n========== Plotting parameter uncertainty propagation ==========\n")

res2d_demo <- res2d_list[[sensitivity_group]]
mean_annual_vec <- res2d_demo$annual_mean
param_list <- res2d_demo$param_list
param_name <- names(param_list[[1]])[1]
param_vals <- sapply(param_list, function(x) x[[param_name]])

if (length(names(param_list[[1]])) > 1) {
  param_name2 <- names(param_list[[1]])[2]
  param_vals2 <- sapply(param_list, function(x) x[[param_name2]])
} else {
  param_name2 <- NULL
  param_vals2 <- NULL
}

df_unc <- data.frame(mean_annual = mean_annual_vec, param1 = param_vals)
if (!is.null(param_name2)) df_unc$param2 <- param_vals2

# Panel A: Uncertainty distribution of mean annual risk
p_A <- ggplot(df_unc, aes(x = mean_annual)) +
  geom_density(fill = "#5e3c99", alpha = 0.4, color = "#5e3c99", linewidth = 1) +
  scale_x_log10() +
  labs(title = "Uncertainty distribution of mean annual risk",
       subtitle = paste("Pathogen:", display_name, "| Group:", sensitivity_group),
       x = "Mean annual infection probability (log scale)", y = "Density") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Panel B: Main parameter vs. mean annual risk
p_B <- ggplot(df_unc, aes(x = param1, y = mean_annual)) +
  geom_point(alpha = 0.4, color = "#2c7bb6", size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#d7191c", fill = "#d7191c", alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = paste(param_name, "(log scale)"),
       y = "Mean annual infection probability (log scale)") +
  theme_classic(base_size = 12)

if (!is.null(param_name2)) {
  # Panel C: Second parameter vs. mean annual risk (if present)
  p_C <- ggplot(df_unc, aes(x = param2, y = mean_annual)) +
    geom_point(alpha = 0.4, color = "#fdae61", size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "#d7191c", fill = "#d7191c", alpha = 0.2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = paste(param_name2, "(log scale)"),
         y = "Mean annual infection probability (log scale)") +
    theme_classic(base_size = 12)
  p_combined <- cowplot::plot_grid(p_A, p_B, p_C, ncol = 1, align = "v",
                                   rel_heights = c(1, 1, 1))
} else {
  p_combined <- cowplot::plot_grid(p_A, p_B, ncol = 1, align = "v",
                                   rel_heights = c(1, 1))
}

ggsave(paste0(file_prefix, "_uncertainty_propagation.png"), p_combined,
       width = 8, height = 6 + ifelse(!is.null(param_name2), 3, 0), dpi = 600)
print(p_combined)

cat("\nAll analyses complete. Figures saved in the current working directory.\n")
