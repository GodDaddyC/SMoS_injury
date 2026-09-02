# injury severity model
pis0 <- function(speed) {
  return(1 / (1 + exp(5.261 - 0.104 * speed)))
}
pis1 <- function(speed, age) {
  return(1 / (1 + exp(5.15 - 0.101 * speed - 0.042 * age)))
}

# MAIS 3+ for VRU Lubbe et al (2022)
pis2 <- function(speed, age) {
  return(1 / (1 + exp(6.19 - 0.078 * speed - 0.038 * age)))
}
pis3 <- function(speed, age) {
  return(1 / (1 + exp(7.47 - 0.079 * speed - 0.047 * age)))
}

injury_df1 <- data.frame(speed = seq(0, 80, 0.5)) %>%
  mutate(InjuryP = sapply(speed, pis0))

injury_df2 <- data.frame(speed = seq(0, 80, 0.5)) %>%
  mutate(InjuryP = sapply(speed, pis1,age=30))

# shared variables for BPOT
thres_order <- bvtcplot(Dat)$k

u <- sort(Dat$prox, decreasing = TRUE)[thres_order]
v <- sort(Dat$Speed, decreasing = TRUE)[thres_order]

# shared variables for fitting CBPOT
pu_default <- 0.15
v_cbpot <- quantile(Dat$prox, 1 - pu_default, na.rm = TRUE)

pot <- fevd(x = prox, data = Dat, threshold = v_cbpot,
            period.basis = "month", time.units = "months", type = "GP")
pot$results$par
pot$call <- "CBPOT GP margin"
plot(pot)

# finding the parametric distribution for speed | X \leq v
conseq0 <- Dat %>% subset(prox > v_cbpot) %>% {.$Speed} %>%
  fitdistrplus::fitdist(distr = "gamma", method = "mme")
conseq <- Dat %>% subset(prox > v_cbpot) %>% {.$Speed} %>%
  fitdistrplus::fitdist(distr = "gamma", method = "mle",
                        start = as.list(conseq0$estimate))
conseq$estimate
plot(conseq)

sq <- seq(0.5, 60, 0.05)

s_un <- pgamma(sq, shape = conseq$estimate[1],
               rate = conseq$estimate[2])

qcrash <- pevd(x0, threshold = v_cbpot, scale = pot$results$par[1],
               shape = pot$results$par[2], type = "GP", lower.tail = FALSE)

# create data for copula
cop_dat <- Dat %>% subset(prox > v_cbpot) %>%
  dplyr::select(prox, Speed) %>%
  mutate(prox = pevd(prox, threshold = v_cbpot, scale = pot$results$par[1],
                     shape = pot$results$par[2], type = "GP"),
         Speed = pgamma(Speed, shape = conseq$estimate[1],
                        rate = conseq$estimate[2])) %>% as.matrix()
