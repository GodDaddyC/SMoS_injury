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

injury_df <- data.frame(speed = seq(0, 80, 0.5)) %>%
  mutate(InjuryP = sapply(speed, pis0))

# shared variables for BPOT
thres_order1 <- bvtcplot(Dat.CN)$k
thres_order2 <- bvtcplot(Dat.SE)$k

u1 <- sort(Dat.CN$prox, decreasing = TRUE)[thres_order1]
v1 <- sort(Dat.CN$Speed, decreasing = TRUE)[thres_order1]
u2 <- sort(Dat.SE$prox, decreasing = TRUE)[thres_order2]
v2 <- sort(Dat.SE$Speed, decreasing = TRUE)[thres_order2]

# shared variables for fitting CBPOT
pu_default <- 0.15
v_1 <- quantile(Dat.CN$prox, 1 - pu_default, na.rm = TRUE)
v_2 <- quantile(Dat.SE$prox, 1 - pu_default, na.rm = TRUE)

pot_1 <- fevd(x = prox, data = Dat.CN, threshold = v_1,
              period.basis = "month", time.units = "0.5/month", type = "GP")
pot_1$results$par
pot_1$call <- "CBPOT CN GP margin"
plot(pot_1)

pot_2 <- fevd(x = prox, data = Dat.SE, threshold = v_2,
              period.basis = "month", time.units = "months", type = "GP")
pot_2$results$par
pot_2$call <- "CBPOT SE GP margin"
plot(pot_2)

# finding the parametric distribution for speed | X \leq v
conseq0_1 <- Dat.CN %>% subset(prox > v_1) %>% {.$Speed} %>%
  fitdistrplus::fitdist(distr = "gamma", method = "mme")
conseq_1 <- Dat.CN %>% subset(prox > v_1) %>% {.$Speed} %>%
  fitdistrplus::fitdist(distr = "gamma", method = "mle",
                        start = as.list(conseq0_1$estimate))
conseq_1$estimate
plot(conseq_1)

conseq0_2 <- Dat.SE %>% subset(prox > v_2) %>% {.$Speed} %>%
  fitdistrplus::fitdist(distr = "gamma", method = "mme")
conseq_2 <- Dat.SE %>% subset(prox > v_2) %>% {.$Speed} %>%
  fitdistrplus::fitdist(distr = "gamma", method = "mle",
                        start = as.list(conseq0_2$estimate))
conseq_2$estimate
plot(conseq_2)

sq_2 <- seq(0.5, 60, 0.05)

s1_un <- pgamma(sq_2, shape = conseq_1$estimate[1],
                rate = conseq_1$estimate[2])
s2_un <- pgamma(sq_2, shape = conseq_2$estimate[1],
                rate = conseq_2$estimate[2])

qcrash_1 <- pevd(x0_1, threshold = v_1, scale = pot_1$results$par[1],
                 shape = pot_1$results$par[2], type = "GP", lower.tail = FALSE)
qcrash_2 <- pevd(x0_2, threshold = v_2, scale = pot_2$results$par[1],
                 shape = pot_2$results$par[2], type = "GP", lower.tail = FALSE)

# create data for copula
cop_dat_1 <- Dat.CN %>% subset(prox > v_1) %>%
  dplyr::select(prox, Speed) %>%
  mutate(prox = pevd(prox, threshold = v_1, scale = pot_1$results$par[1],
                     shape = pot_1$results$par[2], type = "GP"),
         Speed = pgamma(Speed, shape = conseq_1$estimate[1],
                        rate = conseq_1$estimate[2])) %>% as.matrix()

cop_dat_2 <- Dat.SE %>% subset(prox > v_2) %>%
  dplyr::select(prox, Speed) %>%
  mutate(prox = pevd(prox, threshold = v_2, scale = pot_2$results$par[1],
                     shape = pot_2$results$par[2], type = "GP"),
         Speed = pgamma(Speed, shape = conseq_2$estimate[1],
                        rate = conseq_2$estimate[2])) %>% as.matrix()
