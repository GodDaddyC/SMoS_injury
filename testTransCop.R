# test the effect of transformation on injury prob from BPOT

tt <- CN_TTC_A
tt.T <- CN_iTTC_A


x0 <- 0
x0.T <- 1

v1 <- quantile(tt$prox, 0.7, na.rm =TRUE)
v1.T <- quantile(tt.T$prox, 0.7, na.rm =TRUE)

POT1 <- fevd(x = prox,data=tt,threshold = v1,period.basis = "month",
              time.units = "0.5/month", type = "GP")

POT1.T <- fevd(x = prox,data=tt.T,threshold = v1.T,period.basis = "month",
              time.units = "months", type = "GP")

Q1 <- pevd(x0,threshold = v1,scale = POT1$results$par[1],
                 shape = POT1$results$par[2], type = "GP",lower.tail = FALSE) 
Q1.T <- pevd(x0.T,threshold = v1.T,scale = POT1.T$results$par[1],
                 shape = POT1.T$results$par[2], type = "GP",lower.tail = FALSE) 

Conseq0 <- tt %>% subset(prox>v1) %>% {.$Speed}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mme")
Conseq<- tt %>% subset(prox>v1) %>% {.$Speed}%>% 
  fitdistrplus::fitdist(distr="gamma",method = "mle",start=as.list(Conseq0$estimate))

Cop.dat <- tt %>% subset(prox>v1) %>% 
  dplyr::select(prox,Speed) %>% 
  mutate(prox = pevd(prox,threshold = v1,scale = POT1$results$par[1],
                     shape = POT1$results$par[2], type = "GP"),
         Speed = pgamma(Speed,shape = Conseq$estimate[1],
                        rate = Conseq$estimate[2])) %>% as.matrix()
Cop.dat.T <- tt.T %>% subset(prox>v1.T) %>% 
  dplyr::select(prox,Speed) %>% 
  mutate(prox = pevd(prox,threshold = v1.T,scale = POT1.T$results$par[1],
                     shape = POT1.T$results$par[2], type = "GP"),
         Speed = pgamma(Speed,shape = Conseq$estimate[1],
                        rate = Conseq$estimate[2])) %>% as.matrix()
Cop <- fitCopula(rotCopula(claytonCopula(0.5)), data = Cop.dat, method = "mpl")
Cop.T <- fitCopula(rotCopula(claytonCopula(0.5)), data = Cop.dat.T, method = "mpl")

CM <- Q.distr.param(Cop@copula,mar1=POT1,mar2 = Conseq,type = 4)
CM.T <- Q.distr.param(Cop.T@copula,mar1=POT1.T,mar2 = Conseq,type = 4)


ss.1T <- seq(0,60,0.05)

df1 <- create_plot.dfQ(ss.1T,x=x0,model = CM,PX=Q1)

df1.T <- create_plot.dfQ(ss.1T,x=x0.T,model=CM.T,PX=Q1.T)

ggplot(df1,aes(x=speed,y=ConditionalD)) + 
  geom_line(aes(colour = "tt")) + 
  geom_line(data = df1.T,aes(x=speed,y=ConditionalD,colour = "tt.T")) +
  scale_colour_manual(name = "Site", values = c("tt" = "red", "tt.T" = "blue")) +
  labs(x = "Speed (km/h)", y = "f(y|TTC<0)") +
  theme(panel.grid.major = element_line(colour = "gray91"),
        panel.grid.minor = element_line(colour = "gray88"),
        panel.background = element_rect(fill = "white",
                                        colour = "white", linetype = "solid"),
        plot.background = element_rect(linetype = "solid"))
Injury.from_cQ_bivariate(df1,CM,PIS0,x0)
Injury.from_cQ_bivariate(df1.T,CM.T,PIS0,x0.T)
Injury.from_cQ_bivariate1(df1,CM,PIS0,x0,Q1)
Injury.from_cQ_bivariate1(df1.T,CM.T,PIS0,x0.T,Q1.T)
normalize_cQ.bivariate(x0,CM)
normalize_cQ.bivariate(x0.T,CM.T)


contourplot2(Cop@copula, dCopula, nlevels = 20, main = "Copula")
contourplot2(Cop.T@copula, dCopula, nlevels = 20, main = "Copula .T")
persp(Cop@copula, dCopula, zlim = c(0, 5), main = "Empirical copula density")
persp(Cop.T@copula, dCopula, zlim = c(0, 5), main = "Empirical copula density .T")

