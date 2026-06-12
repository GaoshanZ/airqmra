# airqmra: Inhalation Quantitative Microbial Risk Assessment in R

`airqmra` is a professional R package designed for **inhalation infection risk assessment of bioaerosols in the atmospheric environment**.
It integrates a Chinese population exposure factor database, multi‑pathogen dose‑response models, one‑ and two‑dimensional Monte Carlo simulation, annual risk and disease burden calculation, and a comprehensive sensitivity analysis toolkit (Spearman, PRCC, Morris, Sobol’).
The package also provides an interactive Shiny web application and automated report generation following the national technical guideline (T/CSES XXXX‑XXXX).

## Installation

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("GaoshanZ/airqmra")
```

## Quick start

```r
library(airqmra)

# Create an exposure scenario: adult outdoor inhalation of Legionella pneumophila
scen <- new_scenario(
  pathogen_id         = "legionella pneumophila",
  exposure_group      = "Adult",
  microenvironment    = "Outdoor",
  population_deposition = "Male",
  use_uncertainty     = FALSE
)

# 1D Monte Carlo simulation (variability only)
single_risk <- run_qmra_1d(scen, n_iter = 5000)
annual      <- annual_risk(single_risk, frequency = 20)

# 2D Monte Carlo (variability + parameter uncertainty)
scen_unc <- new_scenario(
  pathogen_id         = "legionella pneumophila",
  exposure_group      = "Adult",
  microenvironment    = "Outdoor",
  population_deposition = "Male",
  use_uncertainty     = TRUE
)
res2d <- run_qmra_2d(scen_unc, n_var = 2000, n_unc = 500, frequency = 20)

# Sensitivity analyses
sens_spearman <- sensitivity_1d(scen, n_iter = 3000)
tornado_plot(sens_spearman, title = "Spearman sensitivity")

prcc_res <- sensitivity_prcc(scen, n_iter = 3000, nboot = 200)
tornado_plot(prcc_res, title = "PRCC sensitivity")

mor <- sensitivity_morris(scen, n_levels = 8, n_trajectories = 30)
tornado_plot(mor, title = "Morris sensitivity (mu*)", log_scale = TRUE)

sobol_res <- sobol_sensitivity(scen_unc, n = 1000)
plot(sobol_res)

# Generate a standard‑compliant risk assessment report
generate_risk_report(scen, list(single_risk = single_risk,
                                annual_risk  = annual,
                                sensitivity_1d = sens_spearman),
                     language = "en")
```

## Key features

### 📦 Built‑in databases
- Dose‑response parameters for over 15 common pathogens (exponential, exact beta‑Poisson, approximate beta‑Poisson), including log‑normal uncertainty distributions.
- Chinese population exposure factors: long‑term / short‑term inhalation rates, indoor/outdoor time‑activity patterns, pulmonary deposition fractions, airborne survival rates, and activity‑template scenarios.
- All parameters are CSV‑driven and can be extended by the user.

### 🔬 Dose‑response models
- Exponential model, exact beta‑Poisson model (with numerical safeguards), and approximate beta‑Poisson model (N₅₀‑based).
- **Model‑fitting module**: maximum likelihood estimation from raw dose‑response data, bootstrap percentile confidence intervals, and likelihood‑ratio model comparison.

### 🧬 Genome‑informed dose‑response prediction (proof‑of‑concept)
`airqmra` incorporates a first‑of‑its‑kind module that predicts the α parameter of the approximate beta‑Poisson model from bacterial genome annotations. The prediction relies on a linear regression between adhesion gene density (number of adhesion genes per million base pairs) and α, trained on 10 representative pathogenic strains.

**Key functions**:
- `compute_adhesion_density(vfdb_file, genome_fasta, ...)` – compute adhesion gene density from a VFDB annotation file and a genome FASTA.
- `predict_alpha_from_genome(adhesion_per_mb)` – predict α and its uncertainty using the pre‑trained model.
- `make_prediction_generator(pred)` – convert the prediction into a `dr_generator` suitable for `new_scenario()`.

**Example**:
```r
# compute adhesion density
adhesion_dens <- compute_adhesion_density("strain_vfdb.tsv", "strain.fna",
                                         has_header = FALSE)
# predict alpha
pred <- predict_alpha_from_genome(adhesion_dens)
print(pred)

# create a custom dose‑response generator
gen <- make_prediction_generator(pred)
scen <- new_scenario(pathogen_name = "Predicted strain",
                     exposure_group = "Adult",
                     microenvironment = "Outdoor",
                     dr_generator = gen)
```

> **Note**: The current model is a proof‑of‑concept with limited training data. It is applicable to high‑adhesion, low‑invasiveness pathogens (e.g., *E. coli*, *K. aerogenes*). Predictions should be interpreted with caution. This functionality is available only from the R command line and is not yet integrated into the Shiny application.

### 🎲 Monte Carlo simulation
- **1D simulation**: propagates exposure variability (concentration, inhalation rate, exposure time, deposition fraction, survival fraction).
- **2D simulation**: separates parameter uncertainty (k, α, β, N₅₀) from population variability; outputs confidence intervals for the mean risk and the combined risk distribution.
- Annual infection risk calculation (supports fixed frequencies or frequency distributions).

### 📊 Sensitivity analysis
- **Spearman rank correlation** (1D variability / 2D parameter uncertainty / global mixed).
- **PRCC (Partial Rank Correlation Coefficient)**: monotonic influence after controlling for other variables, with bootstrap confidence intervals.
- **Morris elementary effects**: self‑contained implementation for screening important variables and interactions.
- **Sobol’ variance decomposition** (based on the sobol2007 algorithm): first‑order and total‑effect indices with full interaction quantification.
- Tornado charts (`tornado_plot`) and standard scatter plots.

### 🌐 Data import and validation
- Concentration data import: supports CSV/TXT files, automatic distribution fitting (lognormal, gamma, Weibull), and data diagnostics.
- Scenario validation: `check_scenario()` automatically evaluates parameter plausibility.
- Flexible `update_scenario()` method for batch scenario analyses.

### 📄 Report generation
- One‑click generation of bilingual (Chinese/English) risk assessment reports in HTML, Word, or PDF format, strictly following the framework of “Technical Guideline for Health Risk Assessment of Bioaerosols in the Atmospheric Environment” (T/CSES XXXX‑XXXX).
- Adaptive content: the report automatically includes Spearman, PRCC, Morris, Sobol’, and 2D uncertainty analyses when the corresponding results are provided.

### 🖥️ Shiny interactive interface
- A full‑featured web application (`run_shiny()`) supporting scenario configuration, 1D/2D simulation, multiple sensitivity analyses, scenario comparison, and report download.

## Selected pathogens

| Pathogen | Model | Point‑estimate parameters |
|----------|-------|---------------------------|
| *Legionella pneumophila* | exponential | k = 4.99×10⁻⁵ |
| *Staphylococcus aureus* | exponential | k = 7.64×10⁻⁸ |
| *Pseudomonas aeruginosa* | approx. beta‑Poisson | α = 0.19, N₅₀ = 1.85×10⁴ |
| *Aspergillus fumigatus* | exact beta‑Poisson | α = 0.09, β = 14 |
| *Influenza virus* | approx. beta‑Poisson | α = 0.581, N₅₀ = 9.45×10⁵ |

(For the full list see the built‑in `pathogen_db` dataset.)

## Dependencies

- `dplyr`, `tidyr`, `readr` – data manipulation
- `hypergeo` – exact beta‑Poisson model
- `sensitivity` – Sobol’ analysis
- `shiny`, `shinythemes` – interactive interface
- `rmarkdown`, `knitr` – report generation
- `ggplot2`, `scales` – visualization
- R ≥ 4.0

## Citation

If you use `airqmra` in your research, please cite:

```bibtex
@software{airqmra,
  author = {Gaoshan Zhang},
  title  = {airqmra: Inhalation Quantitative Microbial Risk Assessment in R},
  year   = {2025},
  url    = {https://github.com/GaoshanZ/airqmra},
  note   = {R package version 0.0.0.9000}
}
```

Main data sources for exposure factors and pathogen parameters:
- Duan, X.L. et al. (2014, 2016) *Highlights of the Chinese Exposure Factors Handbook*.
- [QMRA Wiki](https://qmrawiki.msu.edu) and published literature.

## Contributions and feedback

Bug reports and suggestions are welcome via [GitHub Issues](https://github.com/GaoshanZ/airqmra/issues).
Author email: 2024029025@chd.edu.cn

---

**Making risk analysis transparent – from black box to scientific reproducibility!**
