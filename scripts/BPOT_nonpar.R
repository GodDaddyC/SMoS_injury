
prox.1 <- fevd(x = prox,data=Dat.CN,threshold = u1,period.basis = "two week",
                time.units = "0.5/month", type = "GP")
conseq.1 <- fevd(x = Speed,data=Dat.CN,threshold = v1,period.basis = "two week",
                  time.units = "0.5/month", type = "GP")
prox.2 <- fevd(x = prox,data=Dat.SE,threshold = u2,period.basis = "month",
                time.units = "months", type = "GP")
conseq.2 <- fevd(x = Speed,data=Dat.SE,threshold = v2,period.basis = "month",
                  time.units = "months", type = "GP")

Dat_CN_above <- Dat.CN[Dat.CN$prox > u1 & Dat.CN$Speed > v1, ]
Dat_SE_above <- Dat.SE[Dat.SE$prox > u2 & Dat.SE$Speed > v2, ]

Dat_CN.UF <- mapply(function(col,mar, u,eta,margin) {
  mtransform.GPMk2(Dat_CN_above[,col], p = mar, thres = u,eta = eta,margin = margin)},
  col=c(1,2), mar=list(prox.1$results$par,conseq.1$results$par),
  u = c(u1,v1),eta = rep(thres.order1/dim(Dat.CN)[1], 2),margin = rep("frechet", 2))
Dat_CN.Exp <- mapply(function(col,mar, u,eta,margin) {
  mtransform.GPMk2(Dat_CN_above[,col], p = mar, thres = u,eta = eta,margin = margin)},
  col=c(1,2), mar=list(prox.1$results$par,conseq.1$results$par),
  u = c(u1,v1),eta = rep(thres.order1/dim(Dat.CN)[1], 2),margin = rep("exp", 2))


AhatBP.CN <- pickands.Nonpar(dat=Dat.CN,thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],
                mar1 = prox.1$results$par,mar2 = conseq.1$results$par,k=10,N=500,CI=FALSE,ifplot=TRUE)

AhatBP.SE <- pickands.Nonpar(dat=Dat.SE,thres=c(u2,v2),eta=thres.order2/dim(Dat.SE)[1],
                             mar1 = prox.2$results$par,mar2 = conseq.2$results$par,k=10,CI=FALSE,N=600,ifplot=TRUE)

Pcrash.1 <- pevd(x0.1,threshold = u1, scale = prox.1$results$par[1],shape = prox.1$results$par[2],
                  lower.tail = FALSE,type = "GP") * thres.order1/dim(Dat.CN)[1]

Pcrash.2 <- pevd(x0.2,threshold = u2, scale = prox.2$results$par[1],shape = prox.2$results$par[2],
                  lower.tail = FALSE,type = "GP") *thres.order2/dim(Dat.SE)[1]

pbTvevd(q1=x0.1,q2=40,mar1 = prox.1$results$par,mar2 = conseq.1$results$par,model = "nonpar",
           tail.type=2,thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],Ahat = AhatBP.CN$beta)

ss.1T <- seq(v1,55,(55- v1)/150)
ss.2T <- seq(v2,60,(60- v2)/150)

plot.df.1 <- create_plot.df.np(ss.1T,x=x0.1,PX=Pcrash.1,mar1 = prox.1$results$par,mar2 = conseq.1$results$par,
                              thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],Ahat = AhatBP.CN$beta)

plot.df.2 <- create_plot.df.np(ss.2T,x=x0.2,PX=Pcrash.2,mar1 = prox.2$results$par,mar2 = conseq.2$results$par,
                               thres=c(u2,v2),eta=thres.order2/dim(Dat.SE)[1],Ahat = AhatBP.SE$beta)
Injury.from_c_bivariate_np(plot.df.1$speed,Ahat = AhatBP.CN$beta,PX = Pcrash.1,severity = PIS0,x0=x0.1,
                           mar1 = prox.1$results$par,mar2 = conseq.1$results$par,
                           thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1])

Injury.from_c_bivariate_np(plot.df.2$speed,Ahat = AhatBP.SE$beta,PX = Pcrash.2,severity = PIS0,x0=x0.2,
                           mar1 = prox.2$results$par,mar2 = conseq.2$results$par,
                           thres=c(u2,v2),eta=thres.order2/dim(Dat.SE)[1])
# debug use
t <- sapply(seq(v1,55,(55- v1)/150), c.bivariate_np, x=0,
            mar1 = prox.1$results$par,mar2 = conseq.1$results$par,
            thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],Ahat = AhatBP.CN$beta)
tt <- normalize_c.bivariate_np(x0.1,Dat=Dat.CN$Speed,mar1 = prox.1$results$par,mar2 = conseq.1$results$par,
                               thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],Ahat = AhatBP.CN$beta)

t <- sapply(seq(v1,55,(55- v1)/150), c.bivariate_np, x=0,
            mar1 = prox.2$results$par,mar2 = conseq.2$results$par,
            thres=c(u2,v2),eta=thres.order2/dim(Dat.SE)[1],Ahat = AhatBP.SE$beta)
tt <- normalize_c.bivariate_np(x0.1,Dat=Dat.SE$Speed,mar1 = prox.2$results$par,mar2 = conseq.2$results$par,
                               thres=c(u2,v2),eta=thres.order2/dim(Dat.SE)[1],Ahat = AhatBP.SE$beta)
# 
# plot(t/Pcrash.1/tt)
# 
# pp1 <- mtransform.GPMk2(x0.2,prox.2$results$par,thres=u2,eta = thres.order2/dim(Dat.SE)[1],margin = "frechet")
# pp2 <- mtransform.GPMk2(ss.2T,conseq.2$results$par,thres=v2,eta = thres.order2/dim(Dat.SE)[1],margin = "frechet")
# t <- sapply(pp2, c.bivariate1, x=pp1,Ahat = AhatBP.SE)
# tt <- normalize_c.bivariate1(x=pp1,dat=pp2,Ahat = AhatBP.SE)
# plot(t/Pcrash.2/tt)
# c.bivariate1(y=pp2[2], x=pp1,Ahat = AhatBP.SE)

dbTvevd(q1=Dat_CN_above[1,1],q2=Dat_CN_above[1,2],dep=1.5,model = 'hr',thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],
        mar1 = prox.1$results$par,mar2 = conseq.1$results$par,margin='frechet')

dExtDep(Dat_CN.UF[1,], method="Parametric", model="HR", par=2/3, angular=FALSE, log=FALSE)
dbvevd(x=Dat_CN.UF[1,],dep=1.5,model = 'hr',mar1=c(1,1,1),mar2=c(1,1,1))

dExtDep(Dat_CN.Exp[1,], method="Parametric", model="HR", par=1.2, angular=FALSE, log=FALSE)
dbvevd(x=Dat_CN.Exp[1,],dep=0.83,model = 'hr',mar1=c(0,1,0),mar2=c(0,1,0))

test_diff <- function(tx){
  rpar <- dbTvevd(q1=tx[1],q2=tx[2],dep=0.83,model = 'hr',thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],
          mar1 = prox.1$results$par,mar2 = conseq.1$results$par)
  rnonpar <- dbTvNonpar(q1=tx[1],q2=tx[2],thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1],
             mar1 = prox.1$results$par,mar2 = conseq.1$results$par,Ahat = AhatBP.CN$beta)
  c(rpar,rnonpar)
}

