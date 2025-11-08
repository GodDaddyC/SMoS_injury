
x <- x0.1.un
P2 <- Conseq.1
PX <- Qcrash.1
model <- gumbelCopula(param=1.8)

f.y1 <- function(k){
  s <- pgamma(k,shape = P2$estimate[1],rate = P2$estimate[2])
  t <- cQ.bivariate(y=s,x=x,model = model)/PX *0.15 *
  dgamma(k,shape = P2$estimate[1],rate = P2$estimate[2]) * PIS0(k)
  t
}

f.y2 <- function(k){
  s <- qgamma(k,shape = P2$estimate[1],rate = P2$estimate[2])
  t <- cQ.bivariate(y=k,x=x,model = model)/PX * 0.15 *
    dgamma(s,shape = P2$estimate[1],rate = P2$estimate[2]) * PIS0(s)
  t
}

f.y3 <- function(k){
  s <- qgamma(k,shape = P2$estimate[1],rate = P2$estimate[2])
  t <- cQ.bivariate(y=k,x=x,model = model)/PX * PIS0(s) * 0.15
  t
}


tt <- sapply(seq(0,60,0.05), f.y1)
tt2 <- sapply(seq(0,1,0.005), f.y2)
tt3 <- sapply(seq(0,1,0.005), f.y3)

plot(seq(0,60,0.05),tt)
plot(seq(0,1,0.005),tt2)
plot(seq(0,1,0.005),tt3)
integrate(Vectorize(f.y1),lower=0,upper=70)
integrate(Vectorize(f.y2),lower=0,upper=1)
integrate(Vectorize(f.y3),lower=0,upper=1)

