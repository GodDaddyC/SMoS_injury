v.1T <- quantile(CN_TTC_A$N_TTC, 0.7, na.rm =TRUE)
v.2T <- quantile(SE_TTC$N_TTC, 0.8, na.rm =TRUE)
v.1P <- quantile(CN_PET_A$N_PET, 0.7, na.rm =TRUE)
v.2P <- quantile(SE_PET$N_PET, 0.8, na.rm =TRUE)

POT.1T <- fevd(x = N_TTC,data=CN_TTC_A,threshold = v.1T,period.basis = "month",
               time.units = "0.5/month", type = "GP")
POT.1T$results$par
POT.1P <- fevd(x = N_PET,data=CN_PET_A,threshold = v.1P,period.basis = "month",
               time.units = "0.5/month", type = "GP")
POT.1P$results$par

POT.2T <- fevd(x = N_TTC,data=SE_TTC,threshold = v.2T,period.basis = "month",
               time.units = "months", type = "GP")
POT.2T$results$par

POT.2P <- fevd(x = N_PET,data=SE_PET,threshold = v.2P,period.basis = "month",
               time.units = "months", type = "GP")
POT.2P$results$par # no crash

# finding the parametric distribution for speed | X \leq v
Conseq.1T<- CN_TTC_A %>% subset(N_TTC>v.1T) %>% {.$Speed_TTC}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mle",start=list(shape=9.5,rate=0.5)) # gamma distribution
Conseq.1T$estimate
plot(Conseq.1T)

Conseq.1P<- CN_PET_A %>% subset(N_PET>v.1P) %>% {.$Speed_PET}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mle",start=list(shape=8,rate=0.5)) # gamma distribution
Conseq.1P$estimate
plot(Conseq.1P)

Conseq.2T<- SE_TTC %>% subset(N_TTC>v.2T) %>% {.$Speed_TTC}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mle",start=list(shape=4.8,rate=0.4)) # gamma distribution
Conseq.2T$estimate
plot(Conseq.2T)


# create data for copula
Cop.dat.1T <- CN_TTC_A %>% subset(N_TTC>v.1T) %>% 
  dplyr::select(N_TTC,Speed_TTC) %>% 
  mutate(N_TTC = pevd(N_TTC,threshold = v.1T,scale = POT.1T$results$par[1],
                      shape = POT.1T$results$par[2], type = "GP"),
         Speed_TTC = pgamma(Speed_TTC,shape = Conseq.1T$estimate[1],
                            rate = Conseq.1T$estimate[2])) %>% as.matrix()
Cop.dat.1P <- CN_PET_A %>% subset(N_PET>v.1P) %>% 
  dplyr::select(N_PET,Speed_PET) %>% 
  mutate(N_PET = pevd(N_PET,threshold = v.1P,scale = POT.1P$results$par[1],
                      shape = POT.1P$results$par[2], type = "GP"),
         Speed_PET = pgamma(Speed_PET,shape = Conseq.1P$estimate[1],
                            rate = Conseq.1P$estimate[2])) %>% as.matrix()

Cop.dat.2T <- SE_TTC %>% subset(N_TTC>v.2T) %>% 
  dplyr::select(N_TTC,Speed_TTC) %>% 
  mutate(N_TTC = pevd(N_TTC,threshold = v.2T,scale = POT.2T$results$par[1],
                      shape = POT.2T$results$par[2], type = "GP"),
         Speed_TTC = pgamma(Speed_TTC,shape = Conseq.2T$estimate[1],
                            rate = Conseq.2T$estimate[2])) %>% as.matrix()

# fit copula models
Cop.1T.gum <- fitCopula(gumbelCopula(dim = 2), data = Cop.dat.1T, method = "mpl")
Cop.1T.gum@estimate
gofCopula(Cop.1T.gum@copula, Cop.dat.1T, method = "Sn",estim.method='itau',
          N = 1000)

Cop.1P.gum <- fitCopula(normalCopula(dim = 2), data = Cop.dat.1P, method = "mpl")
Cop.1P.gum@estimate
gofCopula(Cop.1P.gum@copula, Cop.dat.1P, method = "Sn", simulation = "pb", N = 500)

Cop.2T.gum <- fitCopula(gumbelCopula(dim = 2), data = Cop.dat.2T, method = "mpl")
Cop.2T.gum@estimate

Cop.2T.clay <- fitCopula(claytonCopula(dim = 2), data = Cop.dat.2T, method = "mpl")
Cop.2T.clay@estimate

# use fitted copula to compute things
CM.1T <- Q.distr.param(Cop.1T.gum@copula,mar1=POT.1T,mar2 = Conseq.1T,type = 4)
pMvdc(c(0,50),CM.1T) # P(X>0,Y<50)

CM.2T <- Q.distr.param(Cop.2T.gum@copula,mar1=POT.2T,mar2 = Conseq.2T,type = 4)
pMvdc(c(0,50),CM.2T) # P(X>0,Y<50)

Qcrash.1T <- pevd(0,threshold = v.1T,scale = POT.1T$results$par[1],
                  shape = POT.1T$results$par[2], type = "GP",lower.tail = FALSE)
Qcrash.2T <- pevd(0,threshold = v.2T,scale = POT.2T$results$par[1],
                  shape = POT.2T$results$par[2], type = "GP",lower.tail = FALSE)

sq.2T <- seq(0,60,0.05)

plot.dfQ.1T <- create_plot.dfQ(sq.2T,x=0,model = CM.1T,PX=Qcrash.1T)
plot.dfQ.2T <- create_plot.dfQ(sq.2T,x=0,model = CM.2T,PX=Qcrash.2T)



