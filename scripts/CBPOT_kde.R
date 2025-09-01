# semi-parametric approach, copula model is non-parametric while margins are



CopNonpar.1 <- kdecop(Cop.dat.1,mult = 0.4)
CopNonpar.2 <- kdecop(Cop.dat.2,mult = 0.5)


#normalize_cQ.bivariate.Nonpar(x0.1.un, Qcrash.1, CopNonpar.1)
# normalize_cQ.bivariate.Nonpar(x0.2.un, Qcrash.2, CopNonpar.2)



plot.dfQ.1 <- create_plot.dfQ.nonpar(s1.un,x=x0.1.un,model = CopNonpar.1,PX=Qcrash.1,P2=Conseq.1)
plot.dfQ.2 <- create_plot.dfQ.nonpar(s2.un,x=x0.2.un,model = CopNonpar.2,PX=Qcrash.2,P2=Conseq.2)

