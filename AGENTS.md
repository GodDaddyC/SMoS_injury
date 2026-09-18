# AGENTS.md

R >= 4.5.0 required (see `uvr.toml`).

## Setup (after clone)

Packages are managed with [uvr](https://github.com/nbafrank/uvr), **not** renv. `.Rprofile` auto-links `.uvr/library/` into `.libPaths()`, but packages are not installed until you run:

```r
uvr::sync()
```

## Architecture

- **functions/** — R functions for EVT-based crash severity models
- **scripts/** — automation for loading variables and setup dependencies
- **results/** — manuscript notebooks (Rmd + knitted `.html`/`.nb.html`)
- **data/** — data and cached output
- **tests/** — unit tests of the key functions

`scripts/preprocessing.R` loads the datasets and sources all `functions/*.R`. `scripts/global_vars.R` is **generic**: given a single dataset `Dat` (and crash boundary `x0`) already in the environment, it defines thresholds, marginal fits, copula data and global variables (`u`, `v`, `v_cbpot`, `pot`, `conseq`, `sq`, `s_un`, `qcrash`, `cop_dat`, ...). To work with several sites, source it once per `Dat` and capture the outputs into named lists (see `scripts/fitModels.R`).

## Entrypoints & execution order

Notebooks live in `results/` (knitted to HTML):
- **`results/execution.Rmd`** — demonstration of the main BPOT/CBPOT functions.
- **`results/estimation.Rmd`** — comparison of estimation methods (censored / Poisson / marginal) for BPOT.
- **`results/dependenceModels.Rmd`** — comparison of BPOT dependence models and CBPOT copulas for the same sampling approach.
- **`results/sensitivity.Rmd`** — sensitivity of crash severity / injury probability to dependence strength.
- **`results/supplimentary.Rmd`** — diagnostics of the marginal distributions.
- **`non_stationary_example.Rmd`** (repo root) — non-stationary bivariate POT via `fbvpot_ns()`.

Scripts that need fitted objects (e.g. `scripts/RSS2026_result.R`, `scripts/synthetic_dependence_test.R`) depend on running the notebooks first.

## BPOT / CBPOT model families

`evd_wrappers.R` provides `pb_tvevd()` and `db_tvevd()` for the bivariate EVT dependence models: `log`, `alog`, `hr`, `neglog`, `aneglog`, `bilog`, `negbilog`, `ct`, plus nonparametric (`nonpar`).

- `run_single_bpot(dat, model, thres, xcrash, speed_ub = 55, estim = "censored")` — BPOT fit + crash probability + plot data.
- `run_single_cbpot(cop_dat, copula, p2, qcrash, pot, pu = pu_default, ...)` — copula-based CBPOT fit (returns `Cop`, `CM`, `CM0`, `x0_un`, `pcrash`, `pu`, `p2`, `s_un`, `qcrash`, `plot_df_q`).

## Result aggregation & summaries

- `create_result_bpot(...)` and `create_result_cbpot(...)` are **variadic**: they simply return `list(...)`. Pass named elements for meaningful labels, e.g. `create_result_bpot(CN = r1, SE = r2)`.
- `summarise_bpot(results, x0 = NULL, severity = NULL, severity_age = NULL, age_mean = 60, labels = NULL, model_names = NULL)` and `summarise_cbpot(results, x0 = NULL, ..., pot = NULL, labels = NULL, model_names = NULL)` take a **list** of single-site results and print each block consecutively. `labels` defaults to `names(results)` (else `"Model i"`). If `model_names` is a character vector of matching length it is pasted into the printout; otherwise a warning is emitted.
- `plot_crash_severity(plot_dfs, injury_df, v = NULL, labels = NULL, legend_title = NULL, ...)` takes a **list** of plot data frames and draws all curves (with an injury-probability secondary axis).

## Dependence-parameter helper

`crash_severity.R` defines `bpot_dep_args(ev_model)` (maps a fitted model to its `dep`/`asy`/`alpha`/`beta` arguments) and `integrate_retry(integrand, x, mar1, thres)` (robust integration with an adaptive upper bound). These are reused across the BPOT probability/density paths (`call_c_bivariate`, `normalize_c_bivariate`, `c_bivariate`, ...).

## Sensitivity analysis (`functions/crash_severity_sen.R`)

- `bpot_dep_change(param, ev_model)` — alter the dependence strength of a fitted BPOT model.
- `bpot_theoretical_density` / `cbpot_theoretical_density` — theoretical crash-severity densities across dependence strengths.
- `bpot_theoretical_injury` / `cbpot_theoretical_injury` — theoretical injury probabilities vs age across dependence strengths (with `filename`/`restart` caching under `data/theoretical_injury/`).
- `plot_theoretical_density(df, ...)` / `plot_theoretical_injury(df, save_file = NULL)` — plot the above; if `save_file` is a string, save the plot under `plots/`.

## fbvpot_ns signature

`fbvpot_ns()` in `functions/fbvpot_ns.R` takes a **data.frame**, not a matrix; the two margins and the scale covariates are specified by column names. Identity/additive link for scale; treatment contrasts for factor covariates (reference level absorbed into the intercept); sequential coefficient names `scale1_0, scale1_1, ...`, `shape1`, `scale2_0, ...`, `shape2`, `dep`.

```r
ns <- fbvpot_ns(data, resp = c("prox3", "maxDV"), threshold = thres.0,
                nsscale1 = c("v_conflict", "movement"), nsscale2 = NULL)
```

S3 methods on the returned `"fbvpot_ns"` object:
- `fitted(ns)` — per-observation fitted scale for each margin (list of two vectors)
- `confint(ns, parm, level)` — normal-approximation intervals from `estimate`/`std.err`
- `plot(ns, num = 1:4)` — four diagnostics: `1` margin QQ, `2` fitted density vs histogram, `3` bivariate quantile curves (sample transformed by the non-stationary scale), `4` Pickands dependence function with CFG estimator and parametric-model CI. `plot(ns)` (no `num`) shows all four with a "Press <Enter>" prompt between plots
- `print(ns)` — concise summary (thresholds, scale models, convergence, deviance)

## Tests

Testthat infrastructure lives in `tests/` (`tests/testthat.R` sources all `functions/*.R`). Test files:
- `tests/testthat/test-bpot.R` — BPOT wrappers and helpers
- `tests/testthat/test-cbpot.R` — CBPOT wrappers and helpers
- `tests/testthat/test-fbvpot_ns.R` — `fbvpot_ns` fitting, factor/numeric covariates, methods

Run with:

```r
testthat::test_dir("tests/testthat")
```

## README file

The `Readme.md` is kept in sync with the code. Trust the code over the README if they ever diverge.
