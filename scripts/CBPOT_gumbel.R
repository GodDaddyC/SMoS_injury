v.1 <- quantile(Dat.CN$prox, 0.7, na.rm =TRUE)
v.2 <- quantile(Dat.SE$prox, 0.8, na.rm =TRUE)


POT.1 <- fevd(x = prox,data=Dat.CN,threshold = v.1,period.basis = "month",
               time.units = "0.5/month", type = "GP")
POT.1$results$par


POT.2 <- fevd(x = prox,data=Dat.SE,threshold = v.2,period.basis = "month",
               time.units = "months", type = "GP")
POT.2$results$par


# finding the parametric distribution for speed | X \leq v
Conseq0.1 <- Dat.CN %>% subset(prox>v.1) %>% {.$Speed}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mme")
Conseq.1<- Dat.CN %>% subset(prox>v.1) %>% {.$Speed}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mle",start=as.list(Conseq0.1$estimate)) # gamma distribution
Conseq.1$estimate
plot(Conseq.1)

Conseq0.2 <- Dat.CN %>% subset(prox>v.1) %>% {.$Speed}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mme")
Conseq.2<- Dat.SE %>% subset(prox>v.2) %>% {.$Speed}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mle",start=as.list(Conseq0.2$estimate)) # gamma distribution
Conseq.2$estimate
plot(Conseq.2)


# create data for copula
Cop.dat.1 <- Dat.CN %>% subset(prox>v.1) %>% 
  dplyr::select(prox,Speed) %>% 
  mutate(prox = pevd(prox,threshold = v.1,scale = POT.1$results$par[1],
                      shape = POT.1$results$par[2], type = "GP"),
         Speed = pgamma(Speed,shape = Conseq.1$estimate[1],
                            rate = Conseq.1$estimate[2])) %>% as.matrix()

Cop.dat.2 <- Dat.SE %>% subset(prox>v.2) %>% 
  dplyr::select(prox,Speed) %>% 
  mutate(prox = pevd(prox,threshold = v.2,scale = POT.2$results$par[1],
                      shape = POT.2$results$par[2], type = "GP"),
         Speed = pgamma(Speed,shape = Conseq.2$estimate[1],
                            rate = Conseq.2$estimate[2])) %>% as.matrix()

# fit copula models
Cop.1 <- fitCopula(gumbelCopula(dim = 2,param = 2), data = Cop.dat.1, method = "mpl")
Cop.1@estimate
Cop.2 <- fitCopula(rotCopula(gumbelCopula(dim = 2,param = 2),flip=c(FALSE,TRUE)), data = Cop.dat.2, method = "mpl")
Cop.2@estimate


# fit multivariate distribution function
CM.1 <- Q.distr.param(Cop.1@copula,mar1=POT.1,mar2 = Conseq.1,type = 4)
CM.2 <- Q.distr.param(Cop.2@copula,mar1=POT.2,mar2 = Conseq.2,type = 4)
CM0.1 <- Q.distr.param(Cop.1@copula,mar1=POT.1,mar2 = Conseq.1,type = 2)
CM0.2 <- Q.distr.param(Cop.2@copula,mar1=POT.2,mar2 = Conseq.2,type = 2)


Qcrash.1 <- pevd(x0.1,threshold = v.1,scale = POT.1$results$par[1],
                  shape = POT.1$results$par[2], type = "GP",lower.tail = FALSE)
Qcrash.2 <- pevd(x0.2,threshold = v.2,scale = POT.2$results$par[1],
                  shape = POT.2$results$par[2], type = "GP",lower.tail = FALSE)

sq.2 <- seq(0,60,0.05)

plot.dfQ.1 <- create_plot.dfQ(sq.2,x=x0.1,model = CM.1,PX=Qcrash.1)
plot.dfQ.2 <- create_plot.dfQ(sq.2,x=x0.2,model = CM.2,PX=Qcrash.2)



