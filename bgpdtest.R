thr <- c(u1,v1)

BP <- mgpd_data(Dat.CN, thr=c(u1,v1))
plot(BP[,1],BP[,2])
abline( h=0, v=0 , lty=2)

init <- mgpd_init(BP)
# default initial value (above) doesn't work so we set them manually:
est.log <- fbgpd( c(init,1.2), dat=BP, model="log", fixed=FALSE )

x <-  seq(-3,7,  0.005)
y <-  seq(-6,20,  0.005)
Q <-  c(0.99,0.9,0.75)
# evaluating the density over the grid
z <- outer(x,  y,  dbgpd,  model="log", mar1 = est.log$par[1:3], mar2 = est.log$par[4:6], dep = 1.2)
# calculating the prediction region
reg <- dbgpd_region( x,  y,  z, quant = Q )
# plotting on the original scale (threshold levels are added back)

contour(reg$x+thr[1],  reg$y+thr[2],  reg$z,  levels=reg$q, drawlabels=FALSE, main="Logistic BGPD",  col=c(1,3,4),lwd=1)
abline( h=thr[2],  v=thr[1],  lty=2 )
legend( "bottomright", c(expression(gamma==0.95),  expression(gamma==0.9), expression(gamma==0.75)), lty=1,  col=c(1,3,4), title="Regions", bty="n")
# and the obs
points(BP[,1]+thr[1],BP[,2]+thr[2],cex=0.7)

pbgpd(-u1,40-v1,model="log",mar1=est.log$par[1:3],mar2=est.log$par[4:6],dep=est.log$par[7])
p1 <- pevd(-u1,threshold =est.log$par[1] ,scale=est.log$par[2],shape=est.log$par[3],type="GP")
p2 <- pevd(40-v2,threshold =est.log$par[4],scale=est.log$par[5],shape=est.log$par[6],type="GP")

1- p1 - p2 + pbgpd(-thr[1],40-thr[2],model="log",mar1=est.log$par[1:3],mar2=est.log$par[4:6],dep=est.log$par[7])

mu <- est.log$par[1]
sigma <- est.log$par[2]
xi <- est.log$par[3]