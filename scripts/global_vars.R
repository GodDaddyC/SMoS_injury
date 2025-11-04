# shared variables for BPOT
thres.order1 <- bvtcplot(Dat.CN)$k
thres.order2 <- bvtcplot(Dat.SE)$k


u1 <- sort(Dat.CN$prox,decreasing = TRUE)[thres.order1]
v1 <- sort(Dat.CN$Speed,decreasing = TRUE)[thres.order1]
u2 <- sort(Dat.SE$prox,decreasing = TRUE)[thres.order2]
v2 <- sort(Dat.SE$Speed,decreasing = TRUE)[thres.order2]

# shared variables for fitting CBPOT 
v.1 <- quantile(Dat.CN$prox, 0.85, na.rm =TRUE)
v.2 <- quantile(Dat.SE$prox, 0.85, na.rm =TRUE)


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

Conseq0.2 <- Dat.SE %>% subset(prox>v.2) %>% {.$Speed}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mme")
Conseq.2<- Dat.SE %>% subset(prox>v.2) %>% {.$Speed}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mle",start=as.list(Conseq0.2$estimate)) # gamma distribution
Conseq.2$estimate
plot(Conseq.2)

sq.2 <- seq(0.5,50,0.05) # start from 0.5 for better numerical stability

s1.un <- pgamma(sq.2,shape = Conseq.1$estimate[1],rate = Conseq.1$estimate[2])
s2.un <- pgamma(sq.2,shape = Conseq.2$estimate[1],rate = Conseq.2$estimate[2])

Qcrash.1 <- pevd(x0.1,threshold = v.1,scale = POT.1$results$par[1],
                 shape = POT.1$results$par[2], type = "GP",lower.tail = FALSE) 
Qcrash.2 <- pevd(x0.2,threshold = v.2,scale = POT.2$results$par[1],
                 shape = POT.2$results$par[2], type = "GP",lower.tail = FALSE)

x0.1.un <- 1 - Qcrash.1
x0.2.un <- 1 - Qcrash.2

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
