# Long-Term Outdoor Exposure QMRA Example

This example demonstrates how to use the `airqmra` package to perform a
long‑term (annual) inhalation infection risk assessment for outdoor
exposure to a specific airborne pathogen.  It uses the built‑in example
concentration data (`Concentration.txt`) and covers multiple population
groups (children 9–12, adults 18–44, elderly 60–79).

## What the script does

- Reads pathogen concentration data and fits a lognormal distribution.
- Creates outdoor exposure scenarios for the selected population groups
  with **parameter uncertainty** enabled (2D Monte Carlo).
- Computes annual infection risks and compares them across groups
  (median + 95% confidence interval).
- Performs **global sensitivity analysis** on the adult group:
  Spearman, PRCC, Morris, and Sobol' variance decomposition.
- Visualises the propagation of dose‑response parameter uncertainty
  on the mean annual risk estimate.

## How to run

1. Install `airqmra` and its dependencies (see the package README).
2. Copy the `long_term_exposure.R` script and the `data` folder to a
   writable directory.
3. Open the script in R or RStudio.
4. Adjust the user‑defined parameters at the top of the script if needed
   (pathogen, population groups, simulation sizes, etc.).
5. Run the entire script.  All output figures will be saved in your
   current working directory as PNG files.

## Input data

The file `data/Concentration.txt` contains tab‑separated columns:
`Sample_ID`, `E_coli`, `P_aeruginosa`, and `S_aureus`.
Each row provides the air concentration (copies/m^3) of the respective
pathogen in one outdoor sample.
By default, the script analyses *Staphylococcus aureus* (column `S_aureus`).

## Output figures

| File | Content |
|------|---------|
| `*_annual_risk_density.png` | Density plot of annual infection risks for the three populations |
| `*_risk_comparison.png` | Point‑range plot showing median and 95% CI, with 10⁻⁴ benchmark line |
| `*_global_spearman.png` | Global Spearman tornado chart |
| `*_global_prcc.png` | Global PRCC tornado chart |
| `*_global_morris.png` | Global Morris sensitivity (mu*) |
| `*_sobol_global.png` | Sobol' variance decomposition (first‑order and total effects) |
| `*_uncertainty_propagation.png` | Uncertainty propagation: distribution of mean risk and parameter‑risk relationship |

## Customisation

All key parameters are clearly labelled at the top of the script.
You can easily change:
- The pathogen (set `pathogen_id` and `conc_col_name`)
- The population groups (edit the `groups` vector)
- The simulation resolution (`N_VAR`, `N_UNC`, `N_VAR_GLOBAL`, `N_UNC_GLOBAL`)
- The representative group for sensitivity analysis (`sensitivity_group`)

## Expected run time

With the default settings (10,000 variability iterations and 5,000
uncertainty iterations for three groups), the 2D simulation may take
several minutes to tens of minutes depending on your hardware.
The global sensitivity analysis (using smaller recorded matrices) runs
more quickly.

## References

Zhang, G. (2025). `airqmra`: Inhalation Quantitative Microbial Risk
Assessment in R.  R package version 0.0.0.9000.
https://github.com/GaoshanZ/airqmra
