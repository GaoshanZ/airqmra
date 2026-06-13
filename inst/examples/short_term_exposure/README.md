# Short-Term Outdoor Exposure QMRA Example

This example demonstrates how to use the `airqmra` package to assess
outdoor short‑term inhalation infection risks for different activity
intensities and exposure durations.  It uses the built‑in example
concentration data (`Concentration.txt`) and allows flexible
customisation of activities, durations, and the target population.

## What the script does

- Reads pathogen concentration data and fits a lognormal distribution.
- Builds outdoor exposure scenarios for user‑specified activity
  intensities (e.g., Sedentary, Light, Moderate) and multiple
  exposure durations (hours).
- Uses short‑term inhalation rates (L/min → m³/h) for the chosen
  population, with parameter uncertainty enabled (1D Monte Carlo).
- Calculates both single‑exposure and annual infection risks
  (annual repeats are set by the user).
- Generates publication‑quality comparison plots showing risk
  as a function of activity intensity and duration.

## How to run

1. Install `airqmra` and its dependencies (see the package README).
2. Copy the `short_term_exposure.R` script and the `data` folder to a
   writable directory.
3. Open the script in R or RStudio.
4. Adjust the user‑defined parameters at the top of the script
   (pathogen, population, activity durations, annual repeats, etc.).
5. Run the entire script.  Output figures are saved as PNG files
   in your current working directory.

## Input data

The file `data/Concentration.txt` contains tab‑separated columns:
`Sample_ID`, `E_coli`, `P_aeruginosa`, and `S_aureus`.
Each row provides the air concentration (copies/m³) of the respective
pathogen in one outdoor sample.
By default, the script analyses *Staphylococcus aureus* (column `S_aureus`).

## Output figures

| File | Content |
|------|---------|
| `*_single_risk_by_duration.png` | Single‑exposure risk (median + 95% CI) across activities and durations |
| `*_annual_risk_by_duration.png` | Annual risk (median + 95% CI) with a 10⁻⁴ benchmark line |

## Customisation

The following parameters at the top of the script can be freely changed:
- `pathogen_id` and `conc_col_name` – select a different pathogen
- `population_group` – any valid classification from `ir_long_table`
- `outdoor_durations` – list each activity intensity and its duration vector
- `ANNUAL_REPEATS` – how many times the activity is performed per year
- `N_ITER_1D` – Monte Carlo iterations

## Expected run time

With the default settings (5000 iterations for about 8 activity‑duration
combinations), the script usually completes within a minute on a
standard laptop.

## References

Zhang, G. (2025). `airqmra`: Inhalation Quantitative Microbial Risk
Assessment in R.  R package version 0.0.0.9000.
https://github.com/GaoshanZ/airqmra
