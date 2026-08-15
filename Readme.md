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

- **Bivariate Peak-Over-Threshold (BPOT) methods** for modeling joint extreme values of conflict proximity and conflict severity
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

- **execution.Rmd**: This file serves as the main script for reproducing the results in the manuscript. It includes sections for loading datasets, applying various statistical models (e.g., BPOT and CBPOT methods), and generating plots. Dependencies are sourced from the `scripts` folder.

- **Procedure.Rmd**: This file includes the separate steps for selecting the models and preconditions which are used in fitting the models.

- **non_stationary_example.Rmd**: Non-stationary bivariate POT analysis using `fbvpot_ns()`.

- **functions/**: Contains utility functions used across the analysis. All `*.R` files are sourced by `scripts/preprocessing.R`. Examples include:
    - `checkDependency.R`: Ensures required R packages are installed and loaded.
    - `auxfun.R`: `run_single_bpot`/`run_single_cbpot` wrappers, `create_result_bpot`/`create_result_cbpot`, and model summaries.
    - `fbvpot_ns.R`: Bivariate POT with non-stationary scale (`fbvpot_ns()`) and its S3 methods (`fitted`, `confint`, `plot(num = 1:4)`, `print`).
    - `trunc_bpot_log.R`: Truncated logistic bivariate POT with a single linear-trend scale covariate (`fit_trunc_log()`).
    - `evd_wrappers.R`: Wrappers for `evd`/`ExtremalDep` functions (`mtransform_gp_mk2`, `pb_tvevd`, `db_tvevd`, `pickands_nonpar`, ...).
    - `copula_wrapper.R`: Copula wrappers and non-parametric copula helpers.
    - `crash_severity.R`, `crash_severity_sen.R`, `injury_prob.R`: Injury probability computation and plotting.
    - `diagnostics_dep.R`: Dependence-structure diagnostics (Pickands function, goodness-of-fit).
    - `fevd_wrapper.R`: Non-stationary marginal helpers for `extRemes::fevd` objects (`pit_linear_fevd`, `build_ns_param`).
    - `ciGPprob.R`, `meanExcessFunMk2.R`: Bayesian GP tail probability and mean-excess utilities.

- **scripts/**: Includes modular scripts for specific statistical methods and copula models. Examples include:
    - `preprocessing.R`: Loads datasets and sources all `functions/*.R`.
    - `global_vars.R`: Defines thresholds and global variables used across the analysis.
    - `BPOT_nonpar.R`, `BPOT_Nonpar_1.R`: Non-parametric BPOT dependence estimation.
    - `Copula_select.R`: Copula model selection.
    - `synthetic_dependence_test.R`, `RSS2026_result.R`: Sensitivity analyses; depend on objects created by running `execution.Rmd` first.

- **data/**: Contains input datasets from different study sites (Sweden and China).

- **tests/**: Testthat infrastructure (`tests/testthat.R` sources all `functions/*.R`). Test files: `test-bpot.R`, `test-cbpot.R`, `test-fbvpot_ns.R`. Run with `testthat::test_dir("tests/testthat")`.

- **plots/**: Output directory for generated figures and visualizations.
> 

---

## Usage

1. Clone the repository
2. Open the R project file (`processing.Rproj`) in RStudio
3. Run `uvr::sync()` to install dependencies
4. Follow the workflow in `execution.Rmd` to reproduce the main analysis results

---

## Contact

For questions or issues related to this code, please open an issue in this repository, or contact me: zhankun.chen@tft.lth.se. 
