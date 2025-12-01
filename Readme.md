# Supplmentary material for "Computing Injury Probability from Microscopic Traffic Data: An Energy-Oriented Approach"

## ⚠️ Disclaimer

**This code is provided for research purposes only.** It is intended solely for academic research, experimentation, and educational use. The code is provided "as is" without warranty of any kind, express or implied. The authors and contributors assume no liability for any consequences arising from the use of this code.

**Please note:**
- This code is not intended for production or commercial use.
- Results should be validated independently before being used in any safety-critical applications.
- Users are responsible for ensuring the code is appropriate for their specific use case.

---

## Project Description

This repository contains the R code for reproducing the results in the STINT (Surrogate Traffic Interaction Model) paper. The project implements bivariate extreme value theory (EVT) methods for estimating crash and injury probabilities from surrogate safety measures (SSM), such as Time-to-Collision (TTC), Post-Encroachment Time (PET), and vehicle speed.

### Key Features

- **Bivariate Peak-Over-Threshold (BPOT) methods** for modeling joint extreme values of conflict proximity and conflict severity
- **Conditional Bivariate POT (CBPOT) methods** using various copula families (Gumbel, Clayton, Gaussian, BB1, Tawn T1)
- **Injury probability estimation** integrating extreme value models with injury severity functions
- **Model diagnostics** including goodness-of-fit tests and tail dependence analysis

---

## Requirements

### R Packages

The following R packages are required. They will be automatically installed if missing when running `functions/checkDependency.R`:

- `readxl`, `dplyr`, `ggplot2`, `ggpubr` - Data manipulation and visualization
- `extRemes`, `evmix`, `evd` - Extreme value analysis
- `copula`, `VC2copula`, `kdecopula`, `VineCopula` - Copula modeling
- `fitdistrplus` - Distribution fitting
- `purrr`, `tidyr` - Data wrangling
- `RColorBrewer`, `ExtremalDep`, `DescTools` - Additional utilities

---

## Workspace Structure

*How to use the repo*

- **execution.Rmd**: This file serves as the main script for reproducing the results in the manuscript. It includes sections for loading datasets, applying various statistical models (e.g., BPOT and CBPOT methods), and generating plots. Dependencies are sourced from the `scripts` folder.

- **Procedure.Rmd**: This file includes the separate steps for selecting the models and preconditions which are used in fitting the models.

- **functions/**: Contains utility functions used across the analysis. Examples include:
    - `checkDependency.R`: Ensures required R packages are installed and loaded. If you encounter any error w.r.t to packages dependencies, you can add the missing packages to the list.
    - `CrashCopula.R`: Implements copula-based crash probability models. Includes also a wrapper of `gofEVCopula` and `gofCopula` from `evd` packages for `VC2Copula` class. 
    - `Truncatedpbeved.R`: A wrapper of `pbvevd` and `dbvevd` from `evd` package for unconditional GP margins. Also includes the computation of injury probability in BPOT approach.
- **scripts/**: Includes modular scripts for specific statistical methods and copula models. Examples include:
    - `BPOT_CT.R`, `BPOT_HR.R`, `BPOT_logistics.R`: Implement bivariate peak-over-threshold methods for different dependence structures, these scripts were generalized as a function thus considered discarded.
    - `CBPOT_clayton.R`, `CBPOT_gumbel.R`, `CBPOT_T1.R`: Implement conditional bivariate peak-over-threshold methods for various copula models, these scripts were generalized as a function thus considered discarded. **Note:** The scripts `CBPOT_*.R` are generally very slow, especially the ones that belong to the VC2copula class.
    - `global_vars.R`: Defines global variables used across the analysis.
    - `theortical_copula.....R`: Conduct the sensitivity analysis in the discussion section 

- **data/**: Contains input datasets from different study sites (Sweden and China).

- **plots/**: Output directory for generated figures and visualizations.
> 

---

## Usage

1. Clone the repository
2. Open the R project file (`processing.Rproj`) in RStudio
3. Follow the workflow in `execution.Rmd` to reproduce the main analysis results

---

## Contact

For questions or issues related to this code, please open an issue in this repository, or contact me: zhankun.chen@tft.lth.se. 
