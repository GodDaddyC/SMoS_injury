# Workspace Structure

*How to use the repo*

- **execution.Rmd**: This file serves as the main script for reproducing the results in the manuscript. It includes sections for loading datasets, applying various statistical models (e.g., BPOT and CBPOT methods), and generating plots. Dependencies are sourced from the `scripts` folder.

- **Procedure.md**: This file includes the separate steps for selecting the models and preconditions which are used in fiting the models

- **functions/**: Contains utility functions used across the analysis. Examples include:
    - `checkDependency.R`: Ensures required R packages are installed and loaded. If you encounter any error w.r.t to packages dependecies, you can add the missing packages to the list.
    - `ciGPprob.R`: Computes confidence intervals for generalized Pareto probabilities.
    - `CrashCopula.R`: Implements copula-based crash probability models. Includes also a wrapper of `gofEVCopula` and `gofCopula` from `evd` packages for `VC2Copula` class. 
    - `meanExcessFunMk2.R`: Calculates mean excess functions for extreme value analysis.
    - `Truncatedpbeved.R`: A wrapper of `pbvevd` and `dbvevd` from `evd` package for unconditional GP margins. Also includes the computation of injury probability

- **scripts/**: Includes modular scripts for specific statistical methods and copula models. Examples include:
    - `BPOT_CT.R`, `BPOT_HR.R`, `BPOT_logistics.R`: Implement bivariate peak-over-threshold methods for different dependence structures.
    - `CBPOT_clayton.R`, `CBPOT_gumbel.R`, `CBPOT_T1.R`: Implement conditional bivariate peak-over-threshold methods for various copula models.
    - `global_vars.R`: Defines global variables used across the analysis.

- **Note:**: The scripts `CBPOT_*.R` are generally very slow. Espeically the ones that belong to VC2copula class. 
