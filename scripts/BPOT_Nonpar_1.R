Dat_CN.UF <- mapply(function(col, mar, u, eta, margin) {
    mtransform_gp_mk2(col, p = mar, thres = u, eta = eta, margin = margin)},
    Dat.CN, mar = list(prox_1$results$par, conseq_1_np$results$par),
    u = c(u1, v1), eta = rep(thres_order1 / dim(Dat.CN)[1], 2),
    margin = rep("frechet", 2))
