# AGENTS.md

R >= 4.5.0 required (see `uvr.toml`).

## Setup (after clone)

Packages are managed with [uvr](https://github.com/nbafrank/uvr), **not** renv. `.Rprofile` auto-links `.uvr/library/` into `.libPaths()`, but packages are not installed until you run:

```r
uvr::sync()
```

## Architecture

Research reproducibility repo (not an R package). Three model families:

- **BPOT** — Bivariate POT via `evd::fbvpot` (logistic, Husler-Reiss, Coles-Tawn)
- **CBPOT** — Copula-based conditional POT via `copula::fitCopula` (gumbel, clayton, normal, BB1, Tawn T1)
- **BPOT-NS** — Bivariate POT with non-stationary scale via `fbvpot_ns()` (`functions/fbvpot_ns.R`): censored likelihood ported from evd's C routine, with per-observation scale driven by covariates (numeric and factor). Also `functions/trunc_bpot_log.R` (`fit_trunc_log()`) for the truncated logistic model with a single linear-trend covariate.

## Entrypoints & execution order

1. **`execution.Rmd`** — primary driver for manuscript results. Chunks must run sequentially:
   - First chunk sources `scripts/preprocessing.R` (loads data + sources all `functions/*.R`)
   - Second chunk sources `scripts/global_vars.R` (defines thresholds, fits GP margins, creates copula data)
   - Later chunks fit BPOT/CBPOT models and produce plots
2. **`Procedure.Rmd`** — model selection diagnostics
3. **`non_stationary_example.Rmd`** — non-stationary bivariate POT analysis
4. **`scripts/synthetic_dependence_test.R`** and **`scripts/RSS2026_result.R`** — sensitivity analyses; depend on objects created by running `execution.Rmd` first

## CBPOT function signatures

`run_single_cbpot()` in `functions/auxfun.R` requires explicit `pot` and `conseq` arguments (no longer reads `POT.1`/`POT.2`/`Conseq.1`/`Conseq.2` from the global env). Pattern:

```r
r1 <- run_single_cbpot(cop_dat_1, copula = gumbelCopula(),
                       p2 = conseq_1, qcrash = qcrash_1, pot = pot_1,
                       method = "mpl")
r2 <- run_single_cbpot(cop_dat_2, copula = rotCopula(gumbelCopula(), flip = c(FALSE, TRUE)),
                       p2 = conseq_2, qcrash = qcrash_2, pot = pot_2,
                       method = "mpl")
result <- create_result_cbpot(CN = r1, SE = r2)
summarise_cbpot(result, pot = list(pot_1, pot_2))
plot_crash_severity(lapply(result, `[[`, "plot_df_q"), injury_df)
```

`summarise_cbpot()` and `plot_crash_severity()` take a **list** of single-site results / plot data frames (labels default to `names(list)`), not a fixed two-input object. `create_result_cbpot()` is variadic: `create_result_cbpot(...)` returns `list(...)`.

Note: the copula data variable is `cop_dat_1` / `cop_dat_2` (lowercase `d`), defined in `scripts/global_vars.R:76-83`.

## BPOT function signatures

```r
r1 <- run_single_bpot(Dat.CN, model = "log", thres = c(u1, v1), xcrash = x0.1)
r2 <- run_single_bpot(Dat.SE, model = "log", thres = c(u2, v2), xcrash = x0.2)
result <- create_result_bpot(CN = r1, SE = r2)
summarise_bpot(result, severity = pis0, severity_age = pis1)
plot_crash_severity(lapply(result, `[[`, "plot_df"), injury_df)
```

`summarise_bpot()` and `plot_crash_severity()` take a **list** of single-site results / plot data frames (labels default to `names(list)`). `create_result_bpot()` is variadic: `create_result_bpot(...)` returns `list(...)`.

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

## README is stale

The `Readme.md` references files that no longer exist (`CrashCopula.R`, `Truncatedpbeved.R`, `BPOT_CT.R`, etc.). Trust the code, not the README.
