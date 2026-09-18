# Copilot Instructions for SMoS Injury Probability Project

## Build and Environment Setup
This project uses R as the primary language and manages its environment using `uvr` and `renv`.

**Setup Commands:**
To initialize the R environment, run these commands in your R session:
```R
renv::init(bare=TRUE)
uvr sync
```

## High-Level Architecture
The project is structured to reproduce research results related to computing injury probabilities from microscopic traffic data using bivariate extreme value theory (EVT).

- **Data Layer**: Input datasets are stored in `data/`.
- **Functional Core**: Most core logic, including dependency checks (`checkDependency.R`), copula-based crash probability models (`CrashCopula.R`), and injury probability computations (`Truncatedpbeved.R`), is located in the `functions/` directory.
- **Modular Scripts**: The `scripts/` folder contains modular components for various statistical methods (e.g., BPOT, CBPOT). Note that many scripts here are legacy as they have been generalized into functions within `functions/`.
- **Reproducibility Workflow**:
    - `execution.Rmd`: The primary entry point for reproducing the results in the manuscript. It coordinates data loading, model fitting (BPOT/CBPOT), and visualization.
    - `Procedure.Rmd`: Outlines the prerequisite steps for model selection and pre-processing conditions.
- **Output**: Generated figures and visualizations are saved to the `plots/` directory.

## Key Conventions
- **Project Management**: Always open and work within the project context provided by `processing.Rproj`.
- **Dependencies**: Use `functions/checkDependency.R` to ensure all required R packages are available before running analysis scripts.
- **Scripting Style**: Prefer using functions from the `functions/` directory over direct implementation in `.Rmd` files or standalone scripts, unless creating a new modular component.
- **Environment Consistency**: Ensure `uvr sync` has been run to maintain consistency across different machines.
