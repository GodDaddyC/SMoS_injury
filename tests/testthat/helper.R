set.seed(42)

# ---- synthetic BPOT data ---------------------------------------------------

n_bpot <- 500
synth_bpot <- data.frame(
  prox  = -rexp(n_bpot, rate = 0.3) - 1,
  Speed = rgamma(n_bpot, shape = 3, rate = 0.08)
)

thres_bpot <- bvtcplot(synth_bpot)
u_synth <- sort(synth_bpot$prox, decreasing = TRUE)[thres_bpot$k]
v_synth <- sort(synth_bpot$Speed, decreasing = TRUE)[thres_bpot$k]

# crash boundary above the prox threshold (avoids density singularity)
x0_synth <- u_synth + 1

# GP margin fits for BPOT tests (both prox and Speed)
mar1_synth <- fevd(x = prox, data = synth_bpot, threshold = u_synth,
                   period.basis = "month", time.units = "months",
                   type = "GP")$results$par
mar2_synth <- fevd(x = Speed, data = synth_bpot, threshold = v_synth,
                   period.basis = "month", time.units = "months",
                   type = "GP")$results$par

eta_synth <- thres_bpot$k / nrow(synth_bpot)

# ---- synthetic CBPOT data --------------------------------------------------

n_cbpot <- 300
cop_gumbel <- gumbelCopula(param = 1.5, dim = 2)
synth_cop_raw <- rCopula(n_cbpot, cop_gumbel)
colnames(synth_cop_raw) <- c("prox", "Speed")

pu_synth <- 0.15
v_cbpot <- quantile(synth_bpot$prox, 1 - pu_synth, na.rm = TRUE)

pot_synth <- fevd(x = prox, data = synth_bpot, threshold = v_cbpot,
                  period.basis = "month", time.units = "months", type = "GP")

conseq0_synth <- synth_bpot %>%
  subset(prox > v_cbpot) %>% {.$Speed} %>%
  fitdistrplus::fitdist(distr = "gamma", method = "mme")

conseq_synth <- synth_bpot %>%
  subset(prox > v_cbpot) %>% {.$Speed} %>%
  fitdistrplus::fitdist(distr = "gamma", method = "mle",
                        start = as.list(conseq0_synth$estimate))

qcrash_synth <- pevd(x0_synth, threshold = v_cbpot,
                     scale = pot_synth$results$par[1],
                     shape = pot_synth$results$par[2],
                     type = "GP", lower.tail = FALSE)

cop_dat_synth <- synth_bpot %>%
  subset(prox > v_cbpot) %>%
  dplyr::select(prox, Speed) %>%
  mutate(
    prox  = pevd(prox, threshold = v_cbpot,
                 scale = pot_synth$results$par[1],
                 shape = pot_synth$results$par[2], type = "GP"),
    Speed = pgamma(Speed, shape = conseq_synth$estimate[1],
                   rate  = conseq_synth$estimate[2])
  ) %>% as.matrix()

sq_synth <- seq(0.5, 60, 0.05)
s_un_synth <- pgamma(sq_synth,
                     shape = conseq_synth$estimate[1],
                     rate  = conseq_synth$estimate[2])

# ---- injury severity functions ---------------------------------------------

pis0_test <- function(speed) {
  1 / (1 + exp(5.261 - 0.104 * speed))
}
pis1_test <- function(speed, age) {
  1 / (1 + exp(5.15 - 0.101 * speed - 0.042 * age))
}

injury_df_test <- data.frame(speed = seq(0, 80, 0.5)) %>%
  mutate(InjuryP = sapply(speed, pis0_test))

# ensure theoretical density output directory exists
dir.create("data/theoretical_density", showWarnings = FALSE, recursive = TRUE)

# temp file helper for theoretical density tests
td_tempfile <- function(prefix) {
  file.path("data", "theoretical_density", paste0(prefix, "_test_", 
             format(Sys.time(), "%H%M%S"), ".csv"))
}
