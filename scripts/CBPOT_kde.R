# semi-parametric approach, copula model is non-parametric while margins are



CopNonpar.1 <- kdecop(Cop.dat.1,mult = 0.2)
CopNonpar.2 <- kdecop(Cop.dat.2,mult = 0.3)



Qcrash.1 <- pevd(x0.1,threshold = v.1,scale = POT.1$results$par[1],
                 shape = POT.1$results$par[2], type = "GP",lower.tail = FALSE)
Qcrash.2 <- pevd(x0.2,threshold = v.2,scale = POT.2$results$par[1],
                 shape = POT.2$results$par[2], type = "GP",lower.tail = FALSE)

x0.1.un <- 1 - Qcrash.1
x0.2.un <- 1 - Qcrash.2



normalize_cQ.bivariate.Nonpar<- function(x,Pcrash, model,lb=0){
  # computes the infinite integral of the conditional density f(y|X >x)
  # use as a nomralization factor for the conditional density
  # lb is the lower bound of the integral, by defalut lb = 0
  c.y <- function(k) {
    temp <- cQ.bivariate.nonpar(v = k, u = x, PX = Pcrash, model = model)
    return(temp)
  }
  
  C <- tryCatch({
    integrate(Vectorize(c.y), lower = lb, upper = 1)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- 0.999
        return(integrate(Vectorize(c.y), lower = lb, upper = ub.alt)$value)
      } 
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(C) 
}
normalize_cQ.bivariate.Nonpar(x0.1.un, Qcrash.1, CopNonpar.1)
normalize_cQ.bivariate.Nonpar(x0.2.un, Qcrash.2, CopNonpar.2)


sq.2 <- seq(0,80,0.05)

s1.un <- pgamma(sq.2,shape = Conseq.1$estimate[1],rate = Conseq.1$estimate[2])
s2.un <- pgamma(sq.2,shape = Conseq.2$estimate[1],rate = Conseq.2$estimate[2])


plot.dfQ.1 <- create_plot.dfQ.nonpar(s1.un,x=x0.1.un,model = CopNonpar.1,PX=Qcrash.1,P2=Conseq.1)
plot.dfQ.2 <- create_plot.dfQ.nonpar(s2.un,x=x0.2.un,model = CopNonpar.2,PX=Qcrash.2,P2=Conseq.2)

scale_factor <- max(plot.dfQ.2$ConditionalD)/max(injury.df$InjuryP)

ggplot(plot.dfQ.2,aes(x=speed,y=ConditionalD)) + 
  geom_line(aes(colour = "SE")) + 
  geom_line(data = plot.dfQ.1,aes(x=speed,y=ConditionalD,colour = "CN")) +
  #geom_line(data = injury.df,aes(x=speed,y=InjuryP * scale_factor) ) +
  # scale_y_continuous(
  #   name = "denisty",  # primary y-axis label
  #   sec.axis = sec_axis(~ . / scale_factor, name = "Injury prob")  # secondary y-axis
  # ) +
  scale_colour_manual(name = "Site", values = c("CN" = "red", "SE" = "blue")) +
  labs(title = "Probability density of impact speed.u",
       x = "speed.u (km/h)", y = "f(y|TTC<0)") +
  theme(panel.grid.major = element_line(colour = "gray91"),
        panel.grid.minor = element_line(colour = "gray88"),
        panel.background = element_rect(fill = "white",
                                        colour = "white", linetype = "solid"),
        plot.background = element_rect(linetype = "solid"))

Injury.from_cQ_bivariate.nonpar(s1.un, CopNonpar.1, PIS0,x0.1.un,Qcrash.1, Conseq.1)/0.3
Injury.from_cQ_bivariate.nonpar(s2.un, CopNonpar.2, PIS0,x0.2.un,Qcrash.2, Conseq.2)/0.2
