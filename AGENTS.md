# AGENTS.md

R >= 4.5.0 required (see `uvr.toml`).

## Setup (after clone)

Packages are managed with [uvr](https://github.com/nbafrank/uvr), **not** renv. `.Rprofile` auto-links `.uvr/library/` into `.libPaths()`, but packages are not installed until you run:

```r
uvr::sync()
```

## Architecture

Research reproducibility repo (not an R package). Two main model families:

- **BPOT** — Bivariate POT via `evd::fbvpot` (logistic, Husler-Reiss, Coles-Tawn)
- **CBPOT** — Copula-based conditional POT via `copula::fitCopula` (gumbel, clayton, normal, BB1, Tawn T1)

## Entrypoints & execution order

1. **`execution.Rmd`** — primary driver for manuscript results. Chunks must run sequentially:
   - First chunk sources `scripts/preprocessing.R` (loads data + sources all `functions/*.R`)
   - Second chunk sources `scripts/global_vars.R` (defines thresholds, fits GP margins, creates copula data)
   - Later chunks fit BPOT/CBPOT models and produce plots
2. **`Procedure.Rmd`** — model selection diagnostics
3. **`scripts/theortical prob *.R`** and **`RSS2026_result.R`** — sensitivity analyses; depend on objects created by running `execution.Rmd` first

## CBPOT function signatures

`runSingleCBPOT()` in `functions/auxfun.R` requires explicit `POT` and `Conseq` arguments (no longer reads `POT.1`/`POT.2`/`Conseq.1`/`Conseq.2` from the global env). Pattern:

```r
r1 <- runSingleCBPOT(Cop.dat.1, copula = gumbelCopula(),
                     P2 = Conseq.1, Qcrash = Qcrash.1, POT = POT.1, Conseq = Conseq.1)
r2 <- runSingleCBPOT(Cop.dat.2, copula = rotCopula(gumbelCopula(), flip = c(FALSE, TRUE)),
                     P2 = Conseq.2, Qcrash = Qcrash.2, POT = POT.2, Conseq = Conseq.2)
result <- create_result_CBPOT(r1, r2)
```

Note: the copula data variable is `Cop.dat.1` / `Cop.dat.2` (lowercase `d`), defined in `scripts/global_vars.R:69-81`.

## BPOT function signatures

```r
r1 <- runSingleBPOT(Dat.CN, model = "log", thres = c(u1, v1), xcrash = x0.1)
r2 <- runSingleBPOT(Dat.SE, model = "log", thres = c(u2, v2), xcrash = x0.2)
result <- create_result_BPOT(r1, r2)
```

## No tests

No test infrastructure exists; functions are validated by running `execution.Rmd` interactively.

## README is stale

The `Readme.md` references files that no longer exist (`CrashCopula.R`, `Truncatedpbeved.R`, `BPOT_CT.R`, etc.). Trust the code, not the README.
