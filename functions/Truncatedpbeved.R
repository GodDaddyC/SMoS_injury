mtransform.GP <- function(x,p,thres,eta,inv = FALSE, drp = TRUE){
  # transform unconditional GP dist for pbvevd of POT1 when inv=F, when inv=T, transform to uniform margin. 
  # p is the vector of scale and shape parameter. this is for atomic x. To use it for a range, use sappaly
  if (is.list(p)) {
    if (is.null(dim(x)) && length(x) != length(p)) 
      stop(paste("`p' must have", length(x), "elements"))
    if (!is.null(dim(x)) && ncol(x) != length(p)) 
      stop(paste("`p' must have", ncol(x), "elements"))
    if (is.null(dim(x))) 
      dim(x) <- c(1, length(p))
    for (i in 1:length(p)) 
      x[, i] <- Recall(x[, i], p[[i]],thres[i],eta, inv = inv)
    if (ncol(x) == 1 || (nrow(x) == 1 && drp)) 
      x <- drop(x)
    return(x)
  }
  if (is.null(dim(x))) 
    dim(x) <- c(length(x), 1)
  p <- matrix(t(p), nrow = nrow(x), ncol = 2, byrow = TRUE)
  if (min(p[, 1]) <= 0) 
    stop("invalid marginal scale")
  expind <- (p[, 2] == 0)  # check if the margin is expontial 
  nzshapes <- p[!expind, 2]
  if (!inv) {
    x <- (x - thres)/p[, 1]
    if (any(x < 0))
      stop("input below thresholds")
    x[expind, ] <- -1/log(1-eta * exp(-x))
    if (any(!expind)) 
      x[!expind, ] <- pmax(-1/log(1-eta * (1 + nzshapes*x)^(-1/nzshapes)),0)
  }
  else {
    x <- ifelse(expind,thres-log( (1-x[expind,])/eta) * p[expind,1],thres + ((1-x[!expind, ])/eta)^(-nzshapes)*p[!expind,1])
    x[expind, ] <- thres-log( (1-x[expind,])/eta) * p[expind,1]
    x[!expind, ] <- thres + ((1-x[!expind, ])/eta)^(-nzshapes)*p[!expind,1]
  }
  if (ncol(x) == 1 || (nrow(x) == 1 && drp)) 
    x <- drop(x)
  x
}


pbTvevd <- function(q,model = c("log", "alog",
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
         log = pbTvlog(q = q, ...), 
         alog = pbTvalog(q = q, ...), 
         hr = pbTvhr(q = q, ...),
         neglog = pbTvneglog(q = q, ...), 
         aneglog = pbTvaneglog(q = q, ...),
         bilog = pbTvbilog(q = q, ...),
         negbilog = pbTvnegbilog(q = q,...), 
         ct = pbTvct(q = q, ...), 
         amix = pbTvamix(q = q, ...))
}

pbTvlog <- function(q,dep,mar1,mar2,tail.type,thres,eta){
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 
      0 || dep > 1) 
    stop("invalid argument for `dep'")
  if (is.null(dim(q))) 
    dim(q) <- c(1, 2)
  q <- mtransform.GP(q, p=list(mar1, mar2),thres,eta)
  v <- sum(q^(-1/dep))^dep
  pp <- exp(-v)
  if (tail.type==2) {
    pp <- 1 - pgev(q[1],loc=1,scale = 1,shape = 1) - pgev(q[2],loc=1,scale = 1,shape = 1) + pp
  }
  else if (tail.type==3) {
    pp <- pgev(q[1],loc=1,scale = 1,shape = 1) - pp
  }
  else if (tail.type==4) {
    pp <- pgev(q[2],loc=1,scale = 1,shape = 1) - pp
  }
  pp
}

# do the same for other evd models.

