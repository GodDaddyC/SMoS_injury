mtransform.GPMk2<- function(x,p,thres,eta,inv = FALSE, margin="exp"){
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
  if (!inv) {
    x <- (x - thres)/p[, 1]
    if (any(x < 0))
      stop("input below thresholds")
    Fx <- ifelse(expind,1-eta * exp(-x),pmax(1-eta * (1 + nzshapes*x)^(-1/nzshapes),0))
    x.t<- switch(margin,
           exp =  -log(Fx),
           frechet = -1/log(Fx),
           uniform = Fx,
           stop("invalid margin type"))
  }
  else {
    x <- ifelse(expind,thres-log( (1-x[expind,])/eta) * p[expind,1],thres + ((1-x[!expind, ])/eta)^(-nzshapes)*p[!expind,1])
    x[expind, ] <- thres-log( (1-x[expind,])/eta) * p[expind,1]
    x[!expind, ] <- thres + ((1-x[!expind, ])/eta)^(-nzshapes)*p[!expind,1]
  }
  x.t
}


pbTvevd <- function(q1,q2,model = c("log", "alog",
          "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),...){
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
         neglog = pbTvneglog(q1,q2, ...), 
         aneglog = pbTvaneglog(q1,q2, ...),
         bilog = pbTvbilog(q1,q2, ...),
         negbilog = pbTvnegbilog(q1,q2,...), 
         ct = pbTvct(q1,q2, ...), 
         amix = pbTvamix(q1,q2, ...))
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
  if (tail.type==2) {
    pp <- 1- pgev(-log(x1))-pgev(-log(x2)) +pp
  }
  else if (tail.type==3) {
    pp <- pgev(-log(x1)) - pp
  }
  else if (tail.type==4) {
    pp <- pgev(-log(x2))- pp
  }
  pp
}

## do the same for other evd models.


# the density function
dbTvevd <- function(q1,q2,model = c("log", "alog","hr", "neglog", "aneglog", 
                                "bilog", "negbilog", "ct", "amix"),...){
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
         amix = dbTvamix(q1, q2, ...))
}
dbTvlog <- function(q1,q2,dep,mar1,mar2,thres,eta,log = FALSE){
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 
      0 || dep > 1) 
    stop("invalid argument for `dep'")
  if (length(eta)==1)
    eta <- rep(eta,2)
  x1 <- mtransform.GPMk2(q1, mar1,thres[1],eta[1],margin = "exp")
  x2 <- mtransform.GPMk2(q2, mar2,thres[2],eta[2],margin = "exp")
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


c.bivariate <-function(y,x,PX, model, dep, thres, eta, mar1, mar2){
  # approximates the conditional density f(y|X >x) by integrating fxy = f(x=x, y = y) over x
  integrand <- function(q1,...) {
    result <- numeric(length(q1))
    for(i in seq_along(q1)) {
      result[i] <- dbTvevd(q1 = q1[i], q2 = y, model = model, dep = dep, 
                           thres = thres, eta = eta, mar1 = mar1, mar2 = mar2)
    }
    return(result)
  }
  
  # Integrate from x to Inf
  integral_result <- integrate(integrand, lower = x, upper = Inf)$value
  
  # Return conditional density
  return(integral_result / PX)
}

