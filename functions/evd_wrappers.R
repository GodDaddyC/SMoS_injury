# wrappers for functions from `evd`, `ExtremalDep`

mtransform.GPMk2<- function(x,p,thres,eta, margin="exp"){
  # transform unconditional GP dist for pbvevd of POT1 when inv=F, when inv=T, transform to uniform margin. 
  # p is the vector of scale and shape parameter. this is for atomic x. To use it for a range, use sapply
  if (is.list(p)) {
    if (is.null(dim(x)) && length(x) != length(p)) 
      stop(paste("`p' must have", length(x), "elements"))
    if (!is.null(dim(x)) && ncol(x) != length(p)) 
      stop(paste("`p' must have", ncol(x), "elements"))
    if (is.null(dim(x))) 
      dim(x) <- c(1, length(p))
    for (i in 1:length(p)) 
      x[, i] <- Recall(x[, i], p[[i]],thres[i],eta, inv = inv)
    if (ncol(x) == 1 || (nrow(x) == 1)) 
      x <- drop(x)
    return(x)
  }
  if (is.null(dim(x))) 
    dim(x) <- c(length(x), 1)
  p <- matrix(t(p), nrow = nrow(x), ncol = 2, byrow = TRUE)
  if (min(p[, 1]) <= 0) 
    stop("invalid marginal scale")
  expind <- (p[, 2] == 0)  # check if the margin is exponential 
  nzshapes <- p[!expind, 2]
  
  x <- (x - thres)/p[, 1]
  if (any(x < 0))
    stop("input below thresholds")
  Fx <- ifelse(expind,1-eta * exp(-x),pmax(1-eta * (1 + nzshapes*x)^(-1/nzshapes),0))
  x.t<- switch(margin,
         exp =  -log(Fx),
         frechet = -1/log(Fx),
         uniform = Fx,
         stop("invalid margin type"))
  
  x.t
}


pbTvevd <- function(q1,q2,model = c("log", "alog",
          "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix", "pb","nonpar"),...){
  model <- match.arg(model)
  # m1 <- c("bilog", "negbilog", "ct", "amix")
  # m2 <- c(m1, "log", "hr", "neglog")
  # m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  # if ((model %in% m1) && !missing(dep)) 
  #   warning("ignoring `dep' argument")
  # if ((model %in% m2) && !missing(asy)) 
  #   warning("ignoring `asy' argument")
  # if ((model %in% m3) && !missing(alpha)) 
  #   warning("ignoring `alpha' argument")
  # if ((model %in% m3) && !missing(beta)) 
  #   warning("ignoring `beta' argument")
  switch(model, 
         log = pbTvlog(q1,q2, ...), 
         alog = pbTvalog(q1,q2, ...), 
         hr = pbTvhr(q1,q2, ...),
         #neglog = pbTvneglog(q1,q2, ...), 
         #aneglog = pbTvaneglog(q1,q2, ...),
         #bilog = pbTvbilog(q1,q2, ...),
         #negbilog = pbTvnegbilog(q1,q2,...), 
         ct = pbTvct(q1,q2, ...), 
         #amix = pbTvamix(q1,q2, ...)
         nonpar = pbTvNonpar(q1,q2,...))
}

pbTvlog <- function(q1,q2,dep,mar1,mar2,tail.type,thres,eta){
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 
      0 || dep > 1) 
    stop("invalid argument for `dep'")

  if (length(eta)==1)
    eta <- rep(eta,2)
  x1 <- mtransform.GPMk2(q1, p=mar1,thres=thres[1],eta[1])
  x2 <- mtransform.GPMk2(q2, p=mar2,thres=thres[2],eta[2])
  v <- sum(x1^(1/dep) + x2^(1/dep))^dep
  pp <- exp(-v)
  # P(X>q1,Y>q2)
  if (tail.type==2) {
    pp <- 1- pgev(-log(x1))-pgev(-log(x2)) +pp
  }
  # P(X<=q1,Y>q2)
  else if (tail.type==3) {
    pp <- pgev(-log(x1)) - pp
  }
  # P(X>q1,Y<= q2)
  else if (tail.type==4) {
    pp <- pgev(-log(x2))- pp
  }
  pp
}

pbTvhr <- function(q1,q2,dep,mar1,mar2,tail.type,thres,eta){
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 
      0) 
    stop("invalid argument for `dep'")
  
  if (length(eta)==1)
    eta <- rep(eta,2)
  x1 <- mtransform.GPMk2(q1, p=mar1,thres=thres[1],eta[1])
  x2 <- mtransform.GPMk2(q2, p=mar2,thres=thres[2],eta[2])
  
  fn <- function(x1, x2) {x1 * pnorm(1/dep + dep * log(x1/x2)/2)}
  v <- fn(x1, x2) + fn(x2, x1)
  
  pp <- exp(-v)
  # P(X>q1,Y>q2)
  if (tail.type==2) {
    pp <- 1- pgev(-log(x1))-pgev(-log(x2)) +pp
  }
  # P(X<=q1,Y>q2)
  else if (tail.type==3) {
    pp <- pgev(-log(x1)) - pp
  }
  # P(X>q1,Y<= q2)
  else if (tail.type==4) {
    pp <- pgev(-log(x2)) - pp
  }
  pp
}

pbTvct <- function(q1,q2,alpha, beta,mar1,mar2,tail.type,thres,eta){
  if (length(alpha) != 1 || mode(alpha) != "numeric") 
    stop("invalid argument for `alpha'")
  if (length(beta) != 1 || mode(beta) != "numeric") 
    stop("invalid argument for `beta'")
  if (any(c(alpha, beta) <= 0)) 
    stop("`alpha' and `beta' must be non-negative")
  
  if (length(eta)==1)
    eta <- rep(eta,2)
  x1 <- mtransform.GPMk2(q1, p=mar1,thres=thres[1],eta[1])
  x2 <- mtransform.GPMk2(q2, p=mar2,thres=thres[2],eta[2])
  
  u <- (alpha * x2)/(alpha * x2 + beta * x1)
  v <-  pbeta(u, shape1 = alpha, shape2 = beta + 1) * x2 + 
        pbeta(u, shape1 = alpha + 1, shape2 = beta, 
               lower.tail = FALSE) * x1
  if (x1 + x2 == 0) {
    v <- 0 }
  if (is.infinite(x1) || is.infinite(x2)) {
    v <- Inf }
  
  
  pp <- exp(-v)
  # P(X>q1,Y>q2)
  if (tail.type==2) {
    pp <- 1- pgev(-log(x1))-pgev(-log(x2)) +pp
  }
  # P(X<=q1,Y>q2)
  else if (tail.type==3) {
    pp <- pgev(-log(x1)) - pp
  }
  # P(X>q1,Y<= q2)
  else if (tail.type==4) {
    pp <- pgev(-log(x2))- pp
  }
  pp
}

pickands.Nonpar <- function(dat,mar1,mar2,thres,eta,est="cfg",CI=FALSE,d=2,N=100,k=10,ifplot=FALSE,
                            nboot=500,alpha=0.05){
  if (length(eta)==1){
    eta <- rep(eta,2)
  }
  dat <- dat[dat[,1] > thres[1] & dat[,2] > thres[2], ]
  dat[,1] <- mtransform.GPMk2(dat[,1], p=mar1,thres=thres[1],eta[1],margin="frechet")
  dat[,2] <- mtransform.GPMk2(dat[,2], p=mar2,thres=thres[2],eta[2],margin="frechet")
  S <- simplex(2,N) # the simplex on which spetral measure is applied [0,1]^2
  if (!CI){
    
    bp.est <- beed(data=dat,x=S,d=d,est=est,margin = "frechet", k = k,plot = ifplot)
    return(bp.est)
  }
  else{
    bp.est <- beed.confband(data=dat,x=S,d=d,est=est,margin = "frechet", conf=1-alpha,k = k,plot = ifplot)
    return(bp.est)
  }
}

pbTvNonpar <- function(q1,q2,mar1,mar2,tail.type,thres,eta,Ahat){
  # Ahat is a vector of non-par estimates obatined from pickands.Nonpar
  if (length(eta)==1){
    eta <- rep(eta,2)
  }
  
  x1 <- mtransform.GPMk2(q1, p=mar1,thres=thres[1],eta[1],margin="exp")
  x2 <- mtransform.GPMk2(q2, p=mar2,thres=thres[2],eta[2],margin = "exp")
  w <- x1/(x1+x2)
  v <- A_bp_approx(Ahat,t=w,ord="0")
  
  pp <- exp(-(1/x1+1/x2)*v)
  # P(X>q1,Y>q2)
  if (tail.type==2) {
    pp <- 1- pgev(-log(x1))-pgev(-log(x2)) +pp
  }
  # P(X<=q1,Y>q2)
  else if (tail.type==3) {
    pp <- pgev(-log(x1)) - pp
  }
  # P(X>q1,Y<= q2)
  else if (tail.type==4) {
    pp <- pgev(-log(x2))- pp
  }
  pp
}
## do the same for other evd models.

A_bp_approx <- function(A_bp,t,ord){
  # approximate the derivative of A function estimated by the Bernstein polynomial method.
  # A_bp is the coefficient estimated by `beed()` from `ExtremeDep`
  # ord = 1 or 2 means A' or A''
  beta_hat <- A_bp
  k <- length(beta_hat)
  beta_hat.k <- beta_hat[2:k]
  b_poly <- function(x, j, k) {
    choose(k, j) * x^j * (1 - x)^(k - j)
  }
  A_prox <- switch(ord,
             "0"= sapply(1:(k-1), function(j) b_poly(t,j,k-1)) %*% beta_hat.k + b_poly(t,0,k-1),
             "1"= ExtremalDep:::ph(t,beta_hat) * 2 - 1,
               #sum(sapply(1:(k-2), function(j) {(beta_hat.k[j+1]-beta_hat.k[j]) * dbeta(t,j+1,k-1-j)} ))+
               #(beta_hat.k[1]-beta_hat[1]) * dbeta(t,1,k-1),
             "2"= ExtremalDep:::dh(t,beta_hat) * 2
               #k*(sum(
               #sapply(1:(k-3), function(j) {(beta_hat.k[j+2] - 2*beta_hat.k[j+1] + beta_hat.k[j]) * 
                # dbeta(t,j+1,k-j-2)} ) ) + (beta_hat.k[2]-2*beta_hat.k[1]+1) * dbeta(t,1,k-2))
                )
  return(A_prox)
}

# the density function
dbTvevd <- function(q1,q2,model = c("log", "alog","hr", "neglog", "aneglog", 
                                "bilog", "negbilog", "ct", "amix","nonpar"),...){
  model <- match.arg(model)
  # m1 <- c("bilog", "negbilog", "ct", "amix")
  # m2 <- c(m1, "log", "hr", "neglog")
  # m3 <- c("log", "alog", "hr", "neglog", "aneglog")
  # if ((model %in% m1) && !missing(dep)) 
  #   warning("ignoring `dep' argument")
  # if ((model %in% m2) && !missing(asy)) 
  #   warning("ignoring `asy' argument")
  # if ((model %in% m3) && !missing(alpha)) 
  #   warning("ignoring `alpha' argument")
  # if ((model %in% m3) && !missing(beta)) 
  #   warning("ignoring `beta' argument")
  switch(model, 
         log = dbTvlog(q1, q2, ...), 
         alog = dbTvalog(q1, q2, ...), 
         hr = dbTvhr(q1, q2, ...),
         neglog = dbTvneglog(q1, q2, ...), 
         aneglog = dbTvaneglog(q1, q2, ...),
         bilog = dbTvbilog(q1, q2, ...),
         negbilog = dbTvnegbilog(q1, q2,...), 
         ct = dbTvct(q1, q2, ...), 
         amix = dbTvamix(q1, q2, ...),
         nonpar = dbTvNonpar(q1, q2,...) )
}
dbTvlog <- function(q1,q2,dep,mar1,mar2,thres,eta,log = FALSE){
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 
      0 || dep > 1) 
    stop("invalid argument for `dep'")
  if (length(eta)==1)
    eta <- rep(eta,2)
  x1 <- mtransform.GPMk2(q1, mar1,thres[1],eta[1])
  x2 <- mtransform.GPMk2(q2, mar2,thres[2],eta[2])
  ext <- c((x1 %in% c(0, Inf)),(x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    idep <- 1/dep
    z <- sum(x1^idep + x2^idep)^dep
    lx <- log(c(x1,x2))
    .expr1 <- (idep + mar1[2]) * lx[1] + (idep + mar2[2]) * lx[2] - log(mar1[1] * mar2[1])
    d <- .expr1 + (1 - 2 * idep) * log(z) + log(idep -  1 + z) - z
  }
  if (!log) 
    d <- exp(d)
  d
}

dbTvhr <- function(q1,q2,dep,mar1,mar2,thres,eta,margin="exp",log = FALSE){
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 
      0) 
    stop("invalid argument for `dep'")
  if (length(eta)==1)
    eta <- rep(eta,2)
  x1 <- mtransform.GPMk2(q1, mar1,thres[1],eta[1],margin = margin)
  x2 <- mtransform.GPMk2(q2, mar2,thres[2],eta[2],margin = margin)
  ext <- c((x1 %in% c(0, Inf)),(x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    fn <- function(x1, x2, nm = pnorm) {x1 * nm(1/dep + dep * log(x1/x2)/2)}
    v <- fn(x1, x2) + fn(x2, x1)
    lx <- log(c(x1,x2))
    .expr1 <- fn(x1, x2) * fn(x2, x1) + dep * fn(x1, x2, nm = dnorm)/2
    jac <- mar1[2] * lx[1] + mar2[2] * lx[2] - log(mar1[1] * mar2[1])
    d <- log(.expr1) + jac - v
  }
  if (!log) 
    d <- exp(d)
  d
}

dbTvct <- function(q1,q2,alpha,beta,mar1,mar2,thres,eta,log = FALSE){
  if (length(alpha) != 1 || mode(alpha) != "numeric") 
    stop("invalid argument for `alpha'")
  if (length(beta) != 1 || mode(beta) != "numeric") 
    stop("invalid argument for `beta'")
  if (any(c(alpha, beta) <= 0)) 
    stop("`alpha' and `beta' must be non-negative")
  if (length(eta)==1)
    eta <- rep(eta,2)
  x1 <- mtransform.GPMk2(q1, mar1,thres[1],eta[1],margin = "exp")
  x2 <- mtransform.GPMk2(q2, mar2,thres[2],eta[2],margin = "exp")
  ext <- c((x1 %in% c(0, Inf)),(x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    u <- (alpha * x2)/(alpha * x2 + beta * x1)
    v <- x2 * pbeta(u, shape1 = alpha, shape2 = beta + 
                      1) + x1 * pbeta(u, shape1 = alpha + 1, shape2 = beta, 
                                      lower.tail = FALSE)
    lx <- log(c(x1,x2))
    jac <- (1 + mar1[2]) * lx[1] + (1 + mar2[2]) * 
      lx[2] - log(mar1[1] * mar2[1])
    .c1 <- alpha * beta/(alpha + beta + 1)
    .expr1 <- pbeta(u, shape1 = alpha, shape2 = beta + 1) * 
      pbeta(u, shape1 = alpha + 1, shape2 = beta, lower.tail = FALSE)
    .expr2 <- dbeta(u, shape1 = alpha + 1, shape2 = beta + 
                      1)/(alpha * x2 + beta * x1)
    d <- log(.expr1 + .c1 * .expr2) - v + jac
  }
  if (!log) 
    d <- exp(d)
  d
}


dbTvNonpar <- function(q1,q2,mar1,mar2,thres,eta,Ahat,log = FALSE){
  if (length(eta)==1)
    eta <- rep(eta,2)
  u1 <- mtransform.GPMk2(q1, mar1,thres[1],eta[1],margin = "uniform")
  u2 <- mtransform.GPMk2(q2, mar2,thres[2],eta[2],margin = "uniform")
  
  ext <- c((u1 %in% c(0, Inf)),(u2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    # Z <- x1 + x2
    # X <- x1*x2
    # G <- pbTvevd(q1, q2, model = "nonpar", mar1=mar1, mar2=mar2,
    #                  thres=thres, eta=eta, Ahat=Ahat,tail.type=2)
    # w <- x1/Z
    w <- log(u1) /log(u1 * u2)
    alpha <- A_bp_approx(Ahat,t = w,ord="0")
    print(alpha)
    # A1 <- A_bp_approx(Ahat,t = w,ord="1")
    # A2 <- A_bp_approx(Ahat,t = w,ord="2")
    # d <- G* (A2/Z^3+(A0^2 + X*A0*A1*(x2-x1)/Z^2 - X*A1^2/Z^2)/X^2 )
    # d <- G* (A0^2 + A0*A1*(x2 - x1)/Z - A1^2*X/Z^2 + X*A2/Z^3)
  }
  d
}









