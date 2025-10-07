
prox.1 <- fevd(x = prox,data=Dat.CN,threshold = u1,period.basis = "two week",
                time.units = "0.5/month", type = "GP")
conseq.1 <- fevd(x = Speed,data=Dat.CN,threshold = v1,period.basis = "two week",
                  time.units = "0.5/month", type = "GP")
prox.2 <- fevd(x = prox,data=Dat.SE,threshold = u2,period.basis = "month",
                time.units = "months", type = "GP")
conseq.2 <- fevd(x = Speed,data=Dat.SE,threshold = v2,period.basis = "month",
                  time.units = "months", type = "GP")

Dat_CN_above <- Dat.CN[Dat.CN$prox > u1 & Dat.CN$Speed > v1, ]
Dat_CN.UF <- mapply(function(col,mar, u,eta,margin) {
  mtransform.GPMk2(Dat_CN_above[,col], p = mar, thres = u,eta = eta,margin = margin)},
  col=c(1,2), mar=list(prox.1$results$par,conseq.1$results$par),
  u = c(u1,v1),eta = rep(thres.order1/dim(Dat.CN)[1], 2),margin = rep("frechet", 2))

Dat_CN.UF_joint <- mapply(function(col,mar, u,eta,margin) {
  mtransform.GPMk2(Dat_CN_above[,col], p = mar, thres = u,eta = eta,margin = margin)},
  col=c(1,2), mar=list(M.1$estimate[1:2],M.1$estimate[3:4]),
  u = c(u1,v1),eta = rep(thres.order1/dim(Dat.CN)[1], 2),margin = rep("frechet", 2))


tt <- fExtDep(x=Dat_CN.UF,method = "PPP",model="HR")

AhatBP.CN <- pickands.Nonpar(dat=Dat.CN,thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],
                mar1 = prox.1$results$par,mar2 = conseq.1$results$par,k=18,N=500,ifplot=TRUE)

AhatBP.SE <- pickands.Nonpar(dat=Dat.SE,thres=c(u2,v2),eta=thres.order2/dim(Dat.SE)[1],
                             mar1 = prox.2$results$par,mar2 = conseq.2$results$par,k=20,bp=TRUE,N=600,ifplot=TRUE)

Pcrash.1 <- pevd(x0.1,threshold = u1, scale = prox.1$results$par[1],shape = prox.1$results$par[2],
                  lower.tail = FALSE,type = "GP") * thres.order1/dim(Dat.CN)[1]

Pcrash.2 <- pevd(x0.2,threshold = u2, scale = prox.2$results$par[1],shape = prox.2$results$par[2],
                  lower.tail = FALSE,type = "GP") *thres.order2/dim(Dat.SE)[1]

pbTvevd(q1=x0.1,q2=40,mar1 = prox.1$results$par,mar2 = conseq.1$results$par,model = "nonpar",
           tail.type=2,thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],Ahat = AhatBP.CN)

ss.1T <- seq(v1,55,(55- v1)/150)
ss.2T <- seq(v2,60,(60- v2)/150)

plot.df.1 <- create_plot.df.np(ss.1T,x=x0.1,PX=Pcrash.1,mar1 = prox.1$results$par,mar2 = conseq.1$results$par,
                              thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],Ahat = AhatBP.CN)

plot.df.2 <- create_plot.df.np(ss.2T,x=x0.2,PX=Pcrash.2,mar1 = prox.2$results$par,mar2 = conseq.2$results$par,
                               thres=c(u2,v2),eta=thres.order2/dim(Dat.SE)[1],Ahat = AhatBP.SE)
Injury.from_c_bivariate_np(plot.df.1$speed,Ahat = AhatBP.CN,PX = Pcrash.1,severity = PIS0,x0=x0.1,
                           mar1 = prox.1$results$par,mar2 = conseq.1$results$par,
                           thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1])

Injury.from_c_bivariate_np(plot.df.2$speed,Ahat = AhatBP.SE,PX = Pcrash.2,severity = PIS0,x0=x0.2,
                           mar1 = prox.2$results$par,mar2 = conseq.2$results$par,
                           thres=c(u2,v2),eta=thres.order2/dim(Dat.SE)[1])
# debug use
t <- sapply(seq(v1,55,(55- v1)/150), c.bivariate_np, x=x0.1,
            mar1 = prox.1$results$par,mar2 = conseq.1$results$par,
            thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],Ahat = AhatBP.CN)
tt <- normalize_c.bivariate_np(x0.1,Dat=Dat.CN$Speed,mar1 = prox.1$results$par,mar2 = conseq.1$results$par,
                               thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],Ahat = AhatBP.CN)

plot(t/Pcrash.1/tt)

pp1 <- mtransform.GPMk2(x0.2,prox.2$results$par,thres=u2,eta = thres.order2/dim(Dat.SE)[1],margin = "frechet")
pp2 <- mtransform.GPMk2(ss.2T,conseq.2$results$par,thres=v2,eta = thres.order2/dim(Dat.SE)[1],margin = "frechet")
t <- sapply(pp2, c.bivariate1, x=pp1,Ahat = AhatBP.SE)
tt <- normalize_c.bivariate1(x=pp1,dat=pp2,Ahat = AhatBP.SE)
plot(t/Pcrash.2/tt)
c.bivariate1(y=pp2[2], x=pp1,Ahat = AhatBP.SE)
