# synthetic_dependence_test.R
# Test whether non-stationarity in copula dependence can be
# explained by non-stationary marginal models.
#
# Two scenarios (select with SCENARIO below):
#   1. Stationary copula: constant tau across seasons.
#      Does fitting stationary margins create a spurious seasonal
#      pattern in the estimated dependence?
#   2. Non-stationary copula: different tau per season.
#      Does the non-stationary marginal model distort or preserve
#      the true seasonal variation in dependence?
#
# Approach:
#   1. Generate bivariate data from a Gumbel copula with either
#      constant or season-varying tau, and non-stationary GEV margins.
#   2. Fit both stationary and non-stationary GEV to each margin.
#   3. PIT-transform to uniform using both sets of models.
#   4. Split pseudo-observations by season; fit Gumbel copula
#      to each season separately.
#   5. Compare estimated tau across seasons.

source("functions/pit_utils.R")

library(extRemes)
library(copula)
library(evd)
library(ggplot2)

# ---- 0. Configuration ----
SEED    <- 2025
N       <- 600
SIGMA   <- 2
XI      <- 0.1

# SCENARIO: 1 = stationary copula, 2 = non-stationary copula (per season)
SCENARIO <- 2

# Base tau and seasonal tau values (used for scenario 1 or 2)
TRUE_TAU_BASE <- 0.5
TRUE_TAU_SEASON <- c(Winter = 0.65, Spring = 0.45, Summer = 0.30, Fall = 0.55)

set.seed(SEED)

# ---- 1. Generate synthetic data ----
t <- rep(1:12, length.out = N)
trend <- seq_len(N) / N
sin_t <- sin(2 * pi * t / 12)
cos_t <- cos(2 * pi * t / 12)

mu1_t <- 10 + 2.5 * sin_t - 1.0 * cos_t + 3.0 * trend
mu2_t <- 12 + 1.5 * sin_t + 2.0 * cos_t + 2.0 * trend

season <- factor(
  dplyr::case_when(
    t %in% c(12, 1, 2) ~ "Winter",
    t %in% 3:5         ~ "Spring",
    t %in% 6:8         ~ "Summer",
    t %in% 9:11        ~ "Fall"
  ),
  levels = c("Winter", "Spring", "Summer", "Fall")
)
seasons <- levels(season)

if (SCENARIO == 1) {
  true_tau_i <- rep(TRUE_TAU_BASE, N)
} else {
  true_tau_i <- TRUE_TAU_SEASON[as.character(season)]
}

true_param_i <- iTau(gumbelCopula(), true_tau_i)

U <- t(sapply(seq_len(N), function(i) {
  rCopula(1, gumbelCopula(param = true_param_i[i], dim = 2))
}))

Y1 <- qgev(U[, 1], loc = mu1_t, scale = SIGMA, shape = XI)
Y2 <- qgev(U[, 2], loc = mu2_t, scale = SIGMA, shape = XI)

df <- data.frame(t, trend, sin_t, cos_t, y1 = Y1, y2 = Y2,
                 season = season, true_tau = true_tau_i)

cat(sprintf("Synthetic data: N = %d, scenario = %d\n", N, SCENARIO))
if (SCENARIO == 1) {
  cat(sprintf("True tau = %.2f (stationary Gumbel copula)\n", TRUE_TAU_BASE))
} else {
  cat("True tau by season (non-stationary Gumbel copula):\n")
  for (s in seasons) {
    cat(sprintf("  %-7s: %.2f\n", s, TRUE_TAU_SEASON[s]))
  }
}

# ---- 2. Fit marginal GEV models ----

cat("\nFitting stationary GEV margins ...\n")
fit_stat_1 <- fevd(y1, data = df, type = "GEV")
fit_stat_2 <- fevd(y2, data = df, type = "GEV")

cat("Fitting non-stationary GEV margins (location ~ seasonal + trend) ...\n")
fit_ns_1 <- fevd(y1, data = df, location.fun = ~ sin_t + cos_t + trend,
                 type = "GEV")
fit_ns_2 <- fevd(y2, data = df, location.fun = ~ sin_t + cos_t + trend,
                 type = "GEV")

# ---- 3. PIT transform to uniform ----

u1_stat <- pit_fevd(fit_stat_1)
u2_stat <- pit_fevd(fit_stat_2)
U_stat <- cbind(u1 = u1_stat, u2 = u2_stat)

cat("Computing per-observation PIT for non-stationary GEV ...\n")
u1_ns <- pit_fevd(fit_ns_1)
u2_ns <- pit_fevd(fit_ns_2)
U_ns <- cbind(u1 = u1_ns, u2 = u2_ns)

# ---- 4. Seasonal copula fitting ----
# Fit a Gumbel copula to the pseudo-observations within each season.

fit_cop_season <- function(U, season_vec) {
  lapply(seasons, function(s) {
    idx <- which(season_vec == s)
    if (length(idx) < 10) {
      return(list(param = NA, tau = NA, n = length(idx), converged = FALSE))
    }
    fit <- tryCatch(
      fitCopula(gumbelCopula(dim = 2), data = U[idx, ],
                method = "ml"),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(param = NA, tau = NA, n = length(idx), converged = FALSE))
    }
    param <- coef(fit)
    tau_val <- tau(gumbelCopula(param = param))
    list(param = param, tau = tau_val, n = length(idx), converged = TRUE)
  })
}

tau_stat_list <- fit_cop_season(U_stat, df$season)
tau_ns_list   <- fit_cop_season(U_ns,   df$season)

# ---- 5. Results table ----

true_tau_season_vec <- if (SCENARIO == 1) {
  rep(TRUE_TAU_BASE, 4)
} else {
  TRUE_TAU_SEASON[seasons]
}

result_table <- data.frame(
  Season            = seasons,
  true_tau          = true_tau_season_vec,
  tau_stationary    = sapply(tau_stat_list, `[[`, "tau"),
  tau_nonstationary = sapply(tau_ns_list,   `[[`, "tau"),
  n                 = sapply(tau_stat_list, `[[`, "n"),
  converged_stat    = sapply(tau_stat_list, `[[`, "converged"),
  converged_ns      = sapply(tau_ns_list,   `[[`, "converged"),
  row.names = NULL
)

cat("\n", paste(rep("=", 75), collapse = ""), "\n")
cat("  Seasonal copula dependence (tau) comparison\n")
cat(paste(rep("=", 75), collapse = ""), "\n\n")

cat(sprintf("%-8s %10s %12s %18s %6s\n",
            "Season", "true_tau", "tau_stat", "tau_nonstat", "n"))
cat(sprintf("%-8s %10s %12s %18s %6s\n",
            "------", "--------", "--------", "-----------", "-----"))
for (i in seq_len(nrow(result_table))) {
  cat(sprintf("%-8s %10.3f %12.3f %18.3f %6d\n",
              result_table$Season[i],
              result_table$true_tau[i],
              result_table$tau_stationary[i],
              result_table$tau_nonstationary[i],
              result_table$n[i]))
}

# ---- 6. Plot: tau by season ----

plot_df <- rbind(
  data.frame(Season = result_table$Season, tau = result_table$true_tau,
             Margins = "True"),
  data.frame(Season = result_table$Season, tau = result_table$tau_stationary,
             Margins = "Stationary GEV"),
  data.frame(Season = result_table$Season, tau = result_table$tau_nonstationary,
             Margins = "Non-stationary GEV")
)
plot_df$Margins <- factor(plot_df$Margins,
                          levels = c("True", "Stationary GEV", "Non-stationary GEV"))

subtitle_str <- if (SCENARIO == 1) {
  paste("True copula: Gumbel, tau =", TRUE_TAU_BASE, "(stationary)")
} else {
  "True copula: Gumbel, varying tau by season (non-stationary)"
}

p_tau <- ggplot(plot_df, aes(x = Season, y = tau, fill = Margins)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           width = 0.75) +
  scale_fill_manual(values = c("True"                  = "#444444",
                                "Stationary GEV"       = "#E69F00",
                                "Non-stationary GEV"   = "#56B4E9")) +
  labs(title    = "Estimated copula dependence (tau) by season",
       subtitle = subtitle_str,
       y        = expression(tau)) +
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_tau)

# ---- 7. PIT histogram diagnostics by season ----

pit_by_season <- function(U, season_vec, margin_name) {
  do.call(rbind, lapply(seasons, function(s) {
    idx <- which(season_vec == s)
    data.frame(Season = s, u = U[idx, 1], Margin = margin_name)
  }))
}

pit_plot_df <- rbind(
  pit_by_season(U_stat, df$season, "Stationary GEV"),
  pit_by_season(U_ns,   df$season, "Non-stationary GEV")
)

expected_bin <- N / 12 / 10  # n per season / 10 bins

p_pit <- ggplot(pit_plot_df, aes(x = u)) +
  facet_grid(rows = vars(Season), cols = vars(Margin)) +
  geom_histogram(boundary = 0, binwidth = 0.1,
                 fill = "steelblue", alpha = 0.7, color = "white") +
  geom_hline(yintercept = expected_bin, linetype = "dashed", linewidth = 0.3) +
  labs(title    = "PIT histograms by season (Margin 1)",
       subtitle = "Uniform distribution expected under correct marginal specification",
       x        = "u", y = "Count") +
  theme_minimal()

print(p_pit)

# ---- 8. Model summary ----

cat("\n", paste(rep("=", 55), collapse = ""), "\n")
cat("  Marginal model summaries\n")
cat(paste(rep("=", 55), collapse = ""), "\n\n")

cat("Stationary GEV (margin 1):\n")
print(summary(fit_stat_1))

cat("\nNon-stationary GEV (margin 1):\n")
print(summary(fit_ns_1))

cat("\nAIC comparison:\n")
cat(sprintf("  Stationary GEV (m1):    %.2f\n", fit_stat_1$results$AIC))
cat(sprintf("  Non-stationary GEV (m1): %.2f\n", fit_ns_1$results$AIC))
cat(sprintf("  Stationary GEV (m2):    %.2f\n", fit_stat_2$results$AIC))
cat(sprintf("  Non-stationary GEV (m2): %.2f\n", fit_ns_2$results$AIC))

cat("\nDone.\n")
