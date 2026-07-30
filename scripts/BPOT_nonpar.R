
prox_1 <- fevd(x = prox, data = Dat.CN, threshold = u1,
               period.basis = "two week", time.units = "0.5/month",
               type = "GP")
conseq_1_np <- fevd(x = Speed, data = Dat.CN, threshold = v1,
                    period.basis = "two week", time.units = "0.5/month",
                    type = "GP")
prox_2 <- fevd(x = prox, data = Dat.SE, threshold = u2,
               period.basis = "month", time.units = "months", type = "GP")
conseq_2_np <- fevd(x = Speed, data = Dat.SE, threshold = v2,
                    period.basis = "month", time.units = "months",
                    type = "GP")

Dat_CN_above <- Dat.CN[Dat.CN$prox > u1 & Dat.CN$Speed > v1, ]
Dat_SE_above <- Dat.SE[Dat.SE$prox > u2 & Dat.SE$Speed > v2, ]

Dat_CN.UF <- mapply(function(col, mar, u, eta, margin) {
  mtransform_gp_mk2(Dat_CN_above[, col], p = mar, thres = u,
                    eta = eta, margin = margin)},
  col = c(1, 2), mar = list(prox_1$results$par, conseq_1_np$results$par),
  u = c(u1, v1), eta = rep(thres_order1 / dim(Dat.CN)[1], 2),
  margin = rep("frechet", 2))
Dat_CN.Exp <- mapply(function(col, mar, u, eta, margin) {
  mtransform_gp_mk2(Dat_CN_above[, col], p = mar, thres = u,
                    eta = eta, margin = margin)},
  col = c(1, 2), mar = list(prox_1$results$par, conseq_1_np$results$par),
  u = c(u1, v1), eta = rep(thres_order1 / dim(Dat.CN)[1], 2),
  margin = rep("exp", 2))
Dat_CN.Un <- mapply(function(col, mar, u, eta, margin) {
  mtransform_gp_mk2(Dat_CN_above[, col], p = mar, thres = u,
                    eta = eta, margin = margin)},
  col = c(1, 2), mar = list(prox_1$results$par, conseq_1_np$results$par),
  u = c(u1, v1), eta = rep(thres_order1 / dim(Dat.CN)[1], 2),
  margin = rep("uniform", 2))


AhatBP.CN <- pickands_nonpar(dat = Dat.CN, thres = c(u1, v1),
  eta = thres_order1 / dim(Dat.CN)[1],
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  k = 10, N = 500, CI = FALSE, ifplot = TRUE)

AhatBP.SE <- pickands_nonpar(dat = Dat.SE, thres = c(u2, v2),
  eta = thres_order2 / dim(Dat.SE)[1],
  mar1 = prox_2$results$par, mar2 = conseq_2_np$results$par,
  k = 10, CI = FALSE, N = 600, ifplot = TRUE)

pcrash_1_np <- pevd(x0_1, threshold = u1, scale = prox_1$results$par[1],
  shape = prox_1$results$par[2], lower.tail = FALSE, type = "GP") *
  thres_order1 / dim(Dat.CN)[1]

pcrash_2_np <- pevd(x0_2, threshold = u2, scale = prox_2$results$par[1],
  shape = prox_2$results$par[2], lower.tail = FALSE, type = "GP") *
  thres_order2 / dim(Dat.SE)[1]

pb_tvevd(q1 = x0_1, q2 = 40, mar1 = prox_1$results$par,
         mar2 = conseq_1_np$results$par, model = "nonpar",
         tail_type = 2, thres = c(u1, v1),
         eta = thres_order1 / dim(Dat.CN)[1], Ahat = AhatBP.CN$beta)

ss_1T <- seq(v1, 55, (55 - v1) / 150)
ss_2T <- seq(v2, 60, (60 - v2) / 150)

plot_df_1_np <- create_plot_df_np(ss_1T, x = x0_1, px = pcrash_1_np,
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  thres = c(u1, v1), eta = thres_order1 / dim(Dat.CN)[1],
  Ahat = AhatBP.CN$beta)

plot_df_2_np <- create_plot_df_np(ss_2T, x = x0_2, px = pcrash_2_np,
  mar1 = prox_2$results$par, mar2 = conseq_2_np$results$par,
  thres = c(u2, v2), eta = thres_order2 / dim(Dat.SE)[1],
  Ahat = AhatBP.SE$beta)

injury_from_c_bivariate_np(plot_df_1_np$speed, Ahat = AhatBP.CN$beta,
  px = pcrash_1_np, severity = pis0, x0 = x0_1,
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  thres = c(u1, v1), eta = thres_order1 / dim(Dat.CN)[1])

injury_from_c_bivariate_np(plot_df_2_np$speed, Ahat = AhatBP.SE$beta,
  px = pcrash_2_np, severity = pis0, x0 = x0_2,
  mar1 = prox_2$results$par, mar2 = conseq_2_np$results$par,
  thres = c(u2, v2), eta = thres_order2 / dim(Dat.SE)[1])

t <- sapply(seq(v1, 55, (55 - v1) / 150), c_bivariate_np, x = 0,
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  thres = c(u1, v1), eta = thres_order1 / dim(Dat.CN)[1],
  Ahat = AhatBP.CN$beta)
tt <- normalize_c_bivariate_np(x0_1, Dat = Dat.CN$Speed,
  mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
  thres = c(u1, v1), eta = thres_order1 / dim(Dat.CN)[1],
  Ahat = AhatBP.CN$beta)

t <- sapply(seq(v1, 55, (55 - v1) / 150), c_bivariate_np, x = 0,
  mar1 = prox_2$results$par, mar2 = conseq_2_np$results$par,
  thres = c(u2, v2), eta = thres_order2 / dim(Dat.SE)[1],
  Ahat = AhatBP.SE$beta)
tt <- normalize_c_bivariate_np(x0_1, Dat = Dat.SE$Speed,
  mar1 = prox_2$results$par, mar2 = conseq_2_np$results$par,
  thres = c(u2, v2), eta = thres_order2 / dim(Dat.SE)[1],
  Ahat = AhatBP.SE$beta)

db_tvevd(q1 = Dat_CN_above[1, 1], q2 = Dat_CN_above[1, 2], dep = 1.5,
  model = "hr", thres = c(u1, v1),
  eta = thres_order1 / dim(Dat.CN)[1],
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
    thres = c(u1, v1), eta = thres_order1 / dim(Dat.CN)[1],
    mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par)
  rnonpar <- db_tv_nonpar(q1 = tx[1], q2 = tx[2], thres = c(u1, v1),
    eta = thres_order1 / dim(Dat.CN)[1],
    mar1 = prox_1$results$par, mar2 = conseq_1_np$results$par,
    Ahat = AhatBP.CN$beta)
  c(rpar, rnonpar)
}
