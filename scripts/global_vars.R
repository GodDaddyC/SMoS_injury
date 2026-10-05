# crash and exposure
poisson_ci <- c(1/14 * qchisq(0.025,df=2),1/14 * qchisq(0.975,df=2))
exposure_day <- dim(se_pet)[1]/33
thres_order <- bvtcplot(se_pet)$k

# comment out the section below if no transformation is needed
Dat.PET <- se_pet

# screening for different parameters
# PET.summary <- trans_quickTest(-Dat.PET$prox,type = 'SRP',obsCI = poisson_ci,
#                                expo = exposure_day,prange1=c(0.5,0.8,1,1.2,1.5),
#                                prange2 =c(0.1,0.2,0.5,0.75,1),u.prob = 1-thres_order/length(Dat.PET$prox))
# PET.summary$logical_matrix
# PET.summary$freq_matrix

PET.SRP.optimal <- transParam_select(-Dat.PET$prox,type='SRP',prange1 = c(1.2,1.6),
                                     prange2 = c(0.2,0.7),obsCI = poisson_ci,
                                     u.prob = 1-thres_order/length(Dat.PET$prox),
                                     expo = exposure_day,
                                     user.control = list(ndeps=c(0.05,0.05)),use_phi = FALSE)
PET_SRP <- freq_loss_fun_trans(c(PET.SRP.optimal$par[1],PET.SRP.optimal$par[2]),
                               dat = -Dat.PET$prox,type='SRP',obs = c(0.003,0.52),
                               u.prob =1-thres_order/length(Dat.PET$prox),
                               expo = exposure_day,if.return = TRUE)

alpha_opt <- PET.SRP.optimal$par[1]
#delta_opt <- PET.SRP.optimal$par[2] # can be both 0.2 and 0.7, the study uses 0.2
delta_opt <- 0.2
# change below if you want Dat to be different indicators
# Dat <- se_pet
# x0 <- 0
Dat <- data.frame(prox=SRP_transform(-se_pet$prox, delta = delta_opt, alpha = alpha_opt),Speed=se_pet$Speed)
x0 <- SRP_transform(0, delta = delta_opt, alpha = alpha_opt)


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
u <- sort(Dat$prox, decreasing = TRUE)[thres_order]
v <- sort(Dat$Speed, decreasing = TRUE)[thres_order]

# shared variables for fitting CBPOT
pu_default <- min(0.15,thres_order / length(Dat$prox))
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
