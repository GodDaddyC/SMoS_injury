
Cn.CN <- empCopula(X = Cop.dat.1, ties.method="random")
Cn.SE <- empCopula(X = Cop.dat.2, ties.method="random")

Sim1 <- rCopula(n=1e5, copula = Cn.CN)
Sim2 <- rCopula(n=1e5, copula = Cn.SE)

cor(Sim1, method="spearman")
cor(Sim2, method="spearman")
cor(Sim1,method = "kendall")
