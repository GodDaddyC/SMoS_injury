thres.order1 <- bvtcplot(Dat.CN)$k
thres.order2 <- bvtcplot(Dat.SE)$k


u1 <- sort(Dat.CN$prox,decreasing = TRUE)[thres.order1]
v1 <- sort(Dat.CN$Speed,decreasing = TRUE)[thres.order1]
u2 <- sort(Dat.SE$prox,decreasing = TRUE)[thres.order2]
v2 <- sort(Dat.SE$Speed,decreasing = TRUE)[thres.order2]


M.1 <- fbvpot(x = Dat.CN,model = "ct",threshold = c(u1,v1))
M.2 <- fbvpot(x = Dat.SE,model = "ct",threshold = c(u2,v2))


Pcrash.1 <- pevd(x0.1,threshold = u1, scale = M.1$estimate[1],shape = M.1$estimate[2],
                  lower.tail = FALSE,type = "GP") * M.1$nat[1]/M.1$n

Pcrash.2 <- pevd(x0.2,threshold = u2, scale = M.2$estimate[1],shape = M.2$estimate[2],
                  lower.tail = FALSE,type = "GP") * M.2$nat[1]/M.2$n

ss.1T <- seq(v1,60,(60- v1)/150)
ss.2T <- seq(v2,60,(60- v2)/150)

plot.df.1 <- create_plot.df(ss.1T,x=x0.1,model=M.1,PX=Pcrash.1)
plot.df.2 <- create_plot.df(ss.2T,x=x0.2,model=M.2,PX=Pcrash.2)




