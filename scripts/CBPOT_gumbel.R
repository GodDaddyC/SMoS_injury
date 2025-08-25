
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


sq.2 <- seq(0,60,0.05)

plot.dfQ.1 <- create_plot.dfQ(sq.2,x=x0.1,model = CM.1,PX=Qcrash.1)
plot.dfQ.2 <- create_plot.dfQ(sq.2,x=x0.2,model = CM.2,PX=Qcrash.2)



