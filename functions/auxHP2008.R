# functions from Heavy tail phenomena Resnick 2008

L1norm<-function (x){
  sum(x)
}

L2norm<-function (x){
  sqrt(sum(x^2))
}

ranktransform<-function (x){
  tr<-(length(x)-rank(x)+1)^(-1)
  #invisible(list(tr=tr))
}

Starica2dPlot<-function (x, y, k, PlotIt = TRUE, norm = L1norm){
  if (length(x) != length(y)) {
    stop("x and y different lengths\n")
  }
  n <- length(x)
  alpha1 <- 1/EstimateEtaHills(k = k, data = x)
  alpha2 <- 1/EstimateEtaHills(k = k, data = y)
  x <- (x/rev(sort(x))[k])^(alpha1)
  y <- (y/rev(sort(y))[k])^(alpha2)
  r <- apply(cbind(x, y), 1, norm)
  u <- rev(sort(r))
  ratio <- (u * (0:(n- 1)))/(length(r[r > 1]))
  if (PlotIt) {
    plot(u[u < 10], ratio[u < 10], xlim = c(0.1, 10),
         type = "l",
         xlab = "scaling constant", ylab = "scaling ratio",
         col = "blue")
    abline(h = 1, col = "red")
    abline(v = 1, lty = 2, lwd = 0.5, col = "red")
    title(paste("k =", k))
  }
  invisible(list(r = u, ratio = ratio,k=k, norm = norm))
}

Starica2dplotrank<-function (x, y, k, PlotIt=TRUE,norm=L1norm)
{
  if (length(x) != length(y)) {
    stop("x and y different lengths\n")
  }
  n <- length(x)
  rx<-ranktransform(x)
  ry<-ranktransform(y)
  R <- apply(cbind(k*rx, k*ry), 1, norm)
  u <- rev(sort(R))
  ratio <- (u * (1:n))/(length(R[R > 1]))
  if (PlotIt) {
    plot(u[u < 5], ratio[u < 5], xlim = c(0.1, 5),
         type = "l",
         xlab = "scaling constant", ylab = "scaling ratio",
         col = "blue")
    abline(h = 1, col = "red")
    abline(v = 1, lty = 2, lwd = 0.5, col = "red")
    title(paste("k =", k))
  }
  invisible(list(r = u, ratio = ratio,k=k, norm = norm))
}

ChooseKrank<-function (x, y, PlotIt = TRUE, norm = L2norm,Lower,Upper){
  if (length(x) != length(y)) {
    stop("x and y different lengths\n")
  }
  n <- length(x)
  rx<-ranktransform(x)
  ry<-ranktransform(y)
  R<-apply(cbind(rx,ry),1,norm) #norms of pairs
  #after rank transform
  R <- rev(sort(R)) #ordered norms; biggest first
  nk <- min(c(500, ceiling(.5*n)))
  Kseq <- (1:nk) #round(exp(seq(log(10), log(n/2), len = nk)))
  dist<-rep(0,nk) #vector of length nk of zeros
  for (i in 1:nk) {
    dist[i]<- L2norm(
      ( i*(1:n)*R/length(i*R[i*R>=1]) )*(i*R>Lower & i*R <=Upper)-
      (i*R>Lower & i*R <=Upper)
    )
  }
  k <- (Kseq[dist == min(dist)])
  if (PlotIt) {
    u<-k*R
    ratio <- k*R*(1:n)/length(k*R[k*R>=1])
    plot(u[u <= Upper & u>Lower], ratio[u <= Upper & u>Lower],
         type = "l", xlab = "scaling constant", ylab = "scaling ratio",
         col="blue")
    abline(h = 1,col="red")
    abline(v = 1,lty=2,lwd=0.5,col="red")
    title(paste("k =",k))
  }
  k
}