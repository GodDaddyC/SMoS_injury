source('scripts/preprocessing.R')

# CN setup
Dat <- cn_ttc_a
x0 <- 0
source("scripts/global_vars.R")
CN <- list(Dat = Dat, u = u, v = v, thres_order = thres_order, x0 = x0,
           n = nrow(Dat))

# SE setup
Dat <- se_ttc
x0 <- 0
source("scripts/global_vars.R")
SE <- list(Dat = Dat, u = u, v = v, thres_order = thres_order, x0 = x0,
           n = nrow(Dat))

prox_1 <- fevd(x = prox, data = CN$Dat, threshold = CN$u,
               period.basis = "two week", time.units = "0.5/month",
               type = "GP")
conseq_1_np <- fevd(x = Speed, data = CN$Dat, threshold = CN$v,
                    period.basis = "two week", time.units = "0.5/month",
                    type = "GP")
prox_2 <- fevd(x = prox, data = SE$Dat, threshold = SE$u,
               period.basis = "month", time.units = "months", type = "GP")
conseq_2_np <- fevd(x = Speed, data = SE$Dat, threshold = SE$v,
                    period.basis = "month", time.units = "months",
                    type = "GP")

Dat_CN_above <- CN$Dat[CN$Dat$prox > CN$u & CN$Dat$Speed > CN$v, ]
Dat_SE_above <- SE$Dat[SE$Dat$prox > SE$u & SE$Dat$Speed > SE$v, ]

Dat_CN.UF <- mapply(function(col, mar, u, eta, margin) {
  mtransform_gp_mk2(Dat_CN_above[, col], p = mar, thres = u,
                    eta = eta, margin = margin)},
  col = c(1, 2), mar = list(prox_1$results$par, conseq_1_np$results$par),
  u = c(CN$u, CN$v), eta = rep(CN$thres_order / CN$n, 2),
  margin = rep("frechet", 2))
Dat_CN.Exp <- mapply(function(col, mar, u, eta, margin) {
  mtransform_gp_mk2(Dat_CN_above[, col], p = mar, thres = u,
                    eta = eta, margin = margin)},
  col = c(1, 2), mar = list(prox_1$results$par, conseq_1_np$results$par),
  u = c(CN$u, CN$v), eta = rep(CN$thres_order / CN$n, 2),
  margin = rep("exp", 2))
Dat_CN.Un <- mapply(function(col, mar, u, eta, margin) {
  mtransform_gp_mk2(Dat_CN_above[, col], p = mar, thres = u,
                    eta = eta, margin = margin)},
  col = c(1, 2), mar = list(prox_1$results$par, conseq_1_np$results$par),
  u = c(CN$u, CN$v), eta = rep(CN$thres_order / CN$n, 2),
  margin = rep("uniform", 2))


AhatBP.CN <- pickands_nonpar(dat = CN$Dat, thres = c(CN$u, CN$v),
  eta = CN$thres_order / CN$n,
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  k = 10, N = 500, CI = FALSE, ifplot = TRUE)

AhatBP.SE <- pickands_nonpar(dat = SE$Dat, thres = c(SE$u, SE$v),
  eta = SE$thres_order / SE$n,
  mar1 = prox_2$results$par, mar2 = conseq_2_np$results$par,
  k = 10, CI = FALSE, N = 600, ifplot = TRUE)

pcrash_1_np <- pevd(CN$x0, threshold = CN$u, scale = prox_1$results$par[1],
  shape = prox_1$results$par[2], lower.tail = FALSE, type = "GP") *
  CN$thres_order / CN$n

pcrash_2_np <- pevd(SE$x0, threshold = SE$u, scale = prox_2$results$par[1],
  shape = prox_2$results$par[2], lower.tail = FALSE, type = "GP") *
  SE$thres_order / SE$n

pb_tvevd(q1 = CN$x0, q2 = 40, mar1 = prox_1$results$par,
         mar2 = conseq_1_np$results$par, model = "nonpar",
         tail_type = 2, thres = c(CN$u, CN$v),
         eta = CN$thres_order / CN$n, Ahat = AhatBP.CN$beta)

ss_1T <- seq(CN$v, 55, (55 - CN$v) / 150)
ss_2T <- seq(SE$v, 60, (60 - SE$v) / 150)

plot_df_1_np <- create_plot_df_np(ss_1T, x = CN$x0, px = pcrash_1_np,
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  thres = c(CN$u, CN$v), eta = CN$thres_order / CN$n,
  Ahat = AhatBP.CN$beta)

plot_df_2_np <- create_plot_df_np(ss_2T, x = SE$x0, px = pcrash_2_np,
  mar1 = prox_2$results$par, mar2 = conseq_2_np$results$par,
  thres = c(SE$u, SE$v), eta = SE$thres_order / SE$n,
  Ahat = AhatBP.SE$beta)

injury_from_c_bivariate_np(plot_df_1_np$speed, Ahat = AhatBP.CN$beta,
  px = pcrash_1_np, severity = pis0, x0 = CN$x0,
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  thres = c(CN$u, CN$v), eta = CN$thres_order / CN$n)

injury_from_c_bivariate_np(plot_df_2_np$speed, Ahat = AhatBP.SE$beta,
  px = pcrash_2_np, severity = pis0, x0 = SE$x0,
  mar1 = prox_2$results$par, mar2 = conseq_2_np$results$par,
  thres = c(SE$u, SE$v), eta = SE$thres_order / SE$n)

t <- sapply(seq(CN$v, 55, (55 - CN$v) / 150), c_bivariate_np, x = 0,
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  thres = c(CN$u, CN$v), eta = CN$thres_order / CN$n,
  Ahat = AhatBP.CN$beta)
tt <- normalize_c_bivariate_np(CN$x0, Dat = CN$Dat$Speed,
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  thres = c(CN$u, CN$v), eta = CN$thres_order / CN$n,
  Ahat = AhatBP.CN$beta)

t <- sapply(seq(CN$v, 55, (55 - CN$v) / 150), c_bivariate_np, x = 0,
  mar1 = prox_2$results$par, mar2 = conseq_2_np$results$par,
  thres = c(SE$u, SE$v), eta = SE$thres_order / SE$n,
  Ahat = AhatBP.SE$beta)
tt <- normalize_c_bivariate_np(CN$x0, Dat = SE$Dat$Speed,
  mar1 = prox_2$results$par, mar2 = conseq_2_np$results$par,
  thres = c(SE$u, SE$v), eta = SE$thres_order / SE$n,
  Ahat = AhatBP.SE$beta)

db_tvevd(q1 = Dat_CN_above[1, 1], q2 = Dat_CN_above[1, 2], dep = 1.5,
  model = "hr", thres = c(CN$u, CN$v),
  eta = CN$thres_order / CN$n,
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  margin = "frechet")

dExtDep(Dat_CN.UF[1, ], method = "Parametric", model = "HR", par = 2 / 3,
        angular = FALSE, log = FALSE)
dbvevd(x = Dat_CN.UF[1, ], dep = 1.5, model = "hr",
       mar1 = c(1, 1, 1), mar2 = c(1, 1, 1))

dExtDep(Dat_CN.Exp[1, ], method = "Parametric", model = "HR", par = 1.2,
        angular = FALSE, log = FALSE)
dbvevd(x = Dat_CN.Exp[1, ], dep = 0.83, model = "hr",
       mar1 = c(0, 1, 0), mar2 = c(0, 1, 0))

test_diff <- function(tx) {
  rpar <- db_tvevd(q1 = tx[1], q2 = tx[2], dep = 0.83, model = "hr",
    thres = c(CN$u, CN$v), eta = CN$thres_order / CN$n,
    mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par)
  rnonpar <- db_tv_nonpar(q1 = tx[1], q2 = tx[2], thres = c(CN$u, CN$v),
    eta = CN$thres_order / CN$n,
    mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
    Ahat = AhatBP.CN$beta)
  c(rpar, rnonpar)
}
