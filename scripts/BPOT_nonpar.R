
prox.1 <- fevd(x = prox,data=Dat.CN,threshold = u1,period.basis = "two week",
                time.units = "0.5/month", type = "GP")
conseq.1 <- fevd(x = Speed,data=Dat.CN,threshold = v1,period.basis = "two week",
                  time.units = "0.5/month", type = "GP")
prox.2 <- fevd(x = prox,data=Dat.SE,threshold = u2,period.basis = "month",
                time.units = "months", type = "GP")
conseq.2 <- fevd(x = Speed,data=Dat.SE,threshold = v2,period.basis = "month",
                  time.units = "months", type = "GP")

Dat_CN_ex <- Dat.CN %>% filter(prox > u1 & Speed > v1)
Dat_SE_ex <- Dat.SE %>% filter(prox > u2 & Speed > v2)

Acfg <- beed(Dat_CN_ex, x=simplex(2,n=100), 2, est="cfg", margin="emp", k=20,plot = TRUE)
Acfg <- beed(Dat_SE_ex, x=simplex(2,n=100), 2, est="cfg", margin="emp", k=20,plot = TRUE)


M.1 <- abvnonpar(data = Dat.CN,empar=TRUE,method='cfg',convex = TRUE,plot = TRUE)

Pcrash.1 <- pevd(x0.1,threshold = u1, scale = prox.1$results$par[1],shape = prox.1$results$par[2],
                  lower.tail = FALSE,type = "GP") * thres.order1/dim(Dat.CN[1])

Pcrash.2 <- pevd(x0.2,threshold = u2, scale = prox.2$results$par[1],shape = prox.2$results$par[2],
                  lower.tail = FALSE,type = "GP") thres.order2/dim(Dat.SE[1])

pbTvNonpar(q1=x0.1,q2=40,dat=Dat.CN,mar1 = prox.1$results$par,mar2 = conseq.1$results$par,
           tail.type=2,thres=c(u1,v1),eta=thres.order1/dim(Dat.CN)[1])
pbTvNonpar(q1=x0.2,q2=40,dat=Dat.SE,mar1 = prox.2$results$par,mar2 = conseq.2$results$par,
           tail.type=2,thres=c(u2,v2),eta=thres.order2/dim(Dat.SE)[1])

ss.1T <- seq(v1,60,(60- v1)/150)
ss.2T <- seq(v2,60,(60- v2)/150)

plot.df.1 <- create_plot.df(ss.1T,x=x0.1,model=M.1,PX=Pcrash.1)
plot.df.2 <- create_plot.df(ss.2T,x=x0.2,model=M.2,PX=Pcrash.2)




