# Supplmentary material for "Computing Injury Probability from Microscopic Traffic Data: An Energy-Oriented Approach"

## ⚠️ Disclaimer

**This code is provided for research purposes only.** It is intended solely for academic research, experimentation, and educational use. The code is provided "as is" without warranty of any kind, express or implied. The authors and contributors assume no liability for any consequences arising from the use of this code.

**Please note:**
- This version of the code is not intended for production or commercial use.
- Results should be validated independently before being used in any safety-critical applications.
- Users are responsible for ensuring the code is appropriate for their specific use case.

---

## Project Description
   
This repository contains the R code for reproducing the results in several papers related to the usage of SMoS to compute injury probability. The project implements bivariate extreme value theory (EVT) methods for estimating crash and injury probabilities.

### Key Features

- **Bivariate Peak-Over-Threshold (BPOT) methods** for modeling joint extreme values of conflict proximity and conflict severity, covering the logistic, asymmetric logistic, negative logistic, asymmetric negative logistic, bilogistic, negative bilogistic, Coles-Tawn and Husler-Reiss dependence models
- **Non-stationary BPOT (BPOT-NS)** via `fbvpot_ns()` with per-observation scale driven by covariates
- **Conditional Bivariate POT (CBPOT) methods** using various copula families (Gumbel, Clayton, Gaussian, BB1, Tawn T1)
- **Injury probability estimation** integrating extreme value models with injury severity functions
- **Model diagnostics** including goodness-of-fit tests and tail dependence analysis

---

## Requirements

### R Packages

R >= 4.5.0 required. Packages are managed with [uvr](https://github.com/nbafrank/uvr), **not** renv. The `.Rprofile` auto-links `.uvr/library/` into `.libPaths()`, but packages are not installed until you run:

```r
uvr::sync()
```

---

## Workspace Structure

*How to use the repo*

The project is organized as follows:

- **results/**: Manuscript notebooks (R Markdown, with knitted `.html`/`.nb.html` outputs). Each notebook demonstrates a specific aspect of the analysis:
    - `execution.Rmd`: Demonstration of the main BPOT and CBPOT functions.
    - `estimation.Rmd`: Comparison of estimation methods (censored / Poisson / marginal) for BPOT.
    - `dependenceModels.Rmd`: Comparison of BPOT dependence models and CBPOT copulas for the same sampling approach.
    - `sensitivity.Rmd`: Sensitivity of crash severity / injury probability to the dependence strength.
    - `supplimentary.Rmd`: Diagnostics of the marginal distributions.

- **non_stationary_example.Rmd**: Non-stationary bivariate POT analysis using `fbvpot_ns()` (kept at the repo root).

- **functions/**: Utility functions used across the analysis. All `*.R` files are sourced by `scripts/preprocessing.R`. Examples include:
    - `checkDependency.R`: Ensures required R packages are installed and loaded.
    - `auxfun.R`: `run_single_bpot`/`run_single_cbpot` wrappers, variadic `create_result_bpot`/`create_result_cbpot`, and list-based model summaries `summarise_bpot`/`summarise_cbpot`.
    - `fbvpot_ns.R`: Bivariate POT with non-stationary scale (`fbvpot_ns()`) and its S3 methods (`fitted`, `confint`, `plot(num = 1:4)`, `print`).
    - `trunc_bpot_log.R`: Truncated logistic bivariate POT with a single linear-trend scale covariate (`fit_trunc_log()`).
    - `evd_wrappers.R`: Wrappers for `evd`/`ExtremalDep` functions (`mtransform_gp_mk2`, `pb_tvevd`, `db_tvevd`, `pickands_nonpar`, ...) supporting the logistic, asymmetric logistic, negative logistic, asymmetric negative logistic, bilogistic, negative bilogistic, Coles-Tawn and Husler-Reiss models.
    - `copula_wrapper.R`: Copula wrappers and non-parametric copula helpers.
    - `crash_severity.R`: Crash severity densities and integrals (`bpot_dep_args`, `integrate_retry`, `c_bivariate`, `normalize_c_bivariate`, `plot_crash_severity`, ...).
    - `crash_severity_sen.R`: Theoretical density / injury sensitivity curves (`bpot_theoretical_density`, `cbpot_theoretical_density`, `bpot_theoretical_injury`, `cbpot_theoretical_injury`, `plot_theoretical_density`, `plot_theoretical_injury`).
    - `injury_prob.R`: Injury probability computation.
    - `diagnostics_dep.R`: Dependence-structure diagnostics (Pickands function `tbevd_plot`, goodness-of-fit `gof_bpot`).
    - `fevd_wrapper.R`: Non-stationary marginal helpers for `extRemes::fevd` objects (`pit_linear_fevd`, `build_ns_param`).
    - `ciGPprob.R`, `meanExcessFunMk2.R`: Bayesian GP tail probability and mean-excess utilities.

- **scripts/**: Modular scripts for setup and specific statistical methods. Examples include:
    - `preprocessing.R`: Loads datasets and sources all `functions/*.R`.
    - `global_vars.R`: Given a single dataset `Dat` (and crash boundary `x0`) in the environment, defines thresholds, marginal fits, copula data and global variables.
    - `fitModels.R`: Fits all BPOT and CBPOT dependence models for a single dataset.
    - `BPOT_nonpar.R`: Non-parametric BPOT dependence estimation.
    - `RSS2026_result.R`, `synthetic_dependence_test.R`: Sensitivity analyses; depend on objects created by running the notebooks first.

- **data/**: Contains input datasets from different study sites (Sweden and China).

- **tests/**: Testthat infrastructure (`tests/testthat.R` sources all `functions/*.R`). Test files: `test-bpot.R`, `test-cbpot.R`, `test-fbvpot_ns.R`. Run with `testthat::test_dir("tests/testthat")`.

- **plots/**: Output directory for generated figures and visualizations.

---

## Usage

1. Clone the repository
2. Open the R project file (`processing.Rproj`) in RStudio
3. Run `uvr::sync()` to install dependencies
4. Follow the workflow in `results/execution.Rmd` to reproduce the main analysis results

---

## Contact

For questions or issues related to this code, please open an issue in this repository, or contact me: zhankun.chen@tft.lth.se. 
