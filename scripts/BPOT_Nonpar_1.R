Dat_CN.UF <- mapply(function(col,mar, u,eta,margin) {
    mtransform.GPMk2(col, p = mar, thres = u,eta = eta,margin = margin)},
    Dat.CN, mar=list(prox.1$results$par,conseq.1$results$par),
    u = c(u1,v1),eta = rep(thres.order1/dim(Dat.CN)[1], 2),
    margin = rep("frechet", 2))
