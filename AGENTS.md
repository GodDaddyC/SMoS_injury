# AGENTS.md

R >= 4.5.0 required (see `uvr.toml`).

## Setup (after clone)

Packages are managed with [uvr](https://github.com/nbafrank/uvr), **not** renv. `.Rprofile` auto-links `.uvr/library/` into `.libPaths()`, but packages are not installed until you run:

```r
uvr::sync()
```

## Architecture

Scripts are organized into three main folders:
- **functions** - R functions for EVT-based crash severity models
- **scripts** — automation for loading variables and setup dependencies
- **results** — contains the main result of the manuscript, organized in different notebooks
- **data** — contains data and cached output
- **tests** - uni-test of the key functions

## Entrypoints & execution order

**`results/execution.Rmd`** —  Examples of BPOT and CBPOT function usage


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

The `Readme.md` may be outdated. Trust the code, not the README.
