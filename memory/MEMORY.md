# RAPTOR Hotspot — Project Memory

## Project
Shiny app for Bayesian traffic collision hotspot prediction (RAPTOR).
Working directory: `/Users/joematthews/Desktop/RAPTORHotspot/`

## Key Files
- `app.R` — complete single-file Shiny app
- `modelscript.R` — original JAGS model reference (not used directly by app)
- `ExampleData.csv` — example upload file (~424KB, many correlated predictors)
- `README.Rmd` — description of app functionality

## App Structure (6 tabs)
1. Introduction — static description
2. Data Upload — fileInput CSV, column mapping (ID/Year/Count), validation
3. Site Selection — checkboxGroupInput, collision trends plot
4. Simulation Settings — MCMC params (nAdapt/nBurnin/nIter/nThin), Run button, preflight checks, log
5. Results — DT table (clickable row → time-series plot), download button
6. Site Warnings — threshold numericInput, colour-coded exceedance probability table

## Model (JAGS)
- Global SPF: `glm.nb(count ~ all_numeric_predictors, data = all_data)` via MASS
- Per-site JAGS model with spike-and-slab trend (`alpha = alpha_n * alpha_z`)
- `sigma` = site-specific intercept offset; `tau` = overdispersion decay
- `pred` = future year count prediction (posterior predictive)
- `lambda[1..n]` = Poisson mean at each observed time point
- `predmu` = SPF value projected one year forward via year coefficient

## R Packages Required
shiny, MASS, dplyr, ggplot2, DT, rjags
System dependency: JAGS (mcmc-jags.sourceforge.net)

## Data Format
Tidy CSV: one row per site per year.
Columns: Site ID, Year (numeric), Collision count (non-negative int), + any numeric predictors.
Example data has many correlated columns (Volume + LogVolume etc.) — GLM may warn about collinearity but runs.

## Time-Series Plot
- Black points = observed counts
- Blue line + ribbon = posterior mean lambda ± 95% CI
- Orange dashed = SPF (global model) fitted values
- Red diamond + error bar = predicted count for next year ± 95% CI
