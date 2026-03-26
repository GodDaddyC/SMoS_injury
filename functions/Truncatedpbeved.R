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
    pp <- pgev(-log(x2))- pp
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
  v <- x2 * pbeta(u, shape1 = alpha, shape2 = beta + 1) + 
    x1 * pbeta(u, shape1 = alpha + 1, shape2 = beta, 
               lower.tail = FALSE)
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
  
  pp <- exp(-(x1+x2)*v)
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

c.bivariate <-function(y,x, model, dep,alpha,beta, thres, eta, mar1, mar2){
  # approximates the conditional density f(y|X >x) by integrating fxy = f(x=x, y = y) over x
  # if integral is non-finite, change ulim to a smaller value
  integrand <- function(q1,...) {
    result <- numeric(length(q1))
    for(i in seq_along(q1)) {
      result[i] <- switch(model,
              log=dbTvevd(q1 = q1[i], q2 = y, model = model, dep = dep, 
                          thres = thres, eta = eta, mar1 = mar1, mar2 = mar2),
              hr=dbTvevd(q1 = q1[i], q2 = y, model = model, dep = dep, 
                         thres = thres, eta = eta, mar1 = mar1, mar2 = mar2),
              ct=dbTvevd(q1 = q1[i], q2 = y, model = model, alpha=alpha,beta=beta, 
                         thres = thres, eta = eta, mar1 = mar1, mar2 = mar2)
              )
    }
    return(result)
  }

  # Integrate from x to Inf
  R <- tryCatch({
    integrate(integrand, lower = x, upper = Inf,rel.tol = 1e-3)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered in conditional density. Retrying with finite upper bound for x.")
        if (mar1[2] < 0) {
          ub.alt <- thres[1] - mar1[1] / mar1[2]
          return(integrate(integrand, lower = x, upper = ub.alt, rel.tol = 1e-3)$value)
        } 
        else {
          # It attempts the integration with decreasing multipliers from 5 down to 1.
          for (i in rev(seq(0.05,5,0.05)) ) {
            ub.alt <- x + i * mar1[1]
            message(paste("Attempting integration with upper bound: ", ub.alt, " (multiplier =", i, ")"))
            result <- tryCatch({
              integrate(integrand, lower = x, upper = ub.alt, rel.tol = 1e-3)$value
              }, 
              error = function(e_retry) {
              return(NA)
            })
            # Check if the integration was successful (i.e., didn't return NA).
            if (!is.na(result)) {
              message("Integration successful with a finite upper bound.")
              return(result)
            }
          }
          # If the loop completes without a successful return, it means all attempts failed.
          stop("Failed to find a finite upper bound after multiple retries.", call. = FALSE)
        }
      } 
      else { stop(e)}
    })
        
  return(R)
}

c.bivariate_np <- function(y,x,thres, eta, mar1, mar2,Ahat){
  integrand <- function(q1,...) {
    result <- numeric(length(q1))
    for(i in seq_along(q1)) {
      result[i] <- dbTvNonpar(q1=q1[i],q2=y,Ahat=Ahat,mar1=mar1,mar2=mar2,thres=thres,eta=eta)
    }
    return(result)
  }
  
  # Integrate from x to Inf
  R <- tryCatch({
    integrate(integrand, lower = x, upper = Inf,rel.tol = 1e-3)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered in conditional density. Retrying with finite upper bound for x.")
        if (mar1[2] < 0) {
          # ub.alt <- thres[1] - mar1[1] / mar1[2]
          # return(integrate(integrand, lower = x, upper = ub.alt, rel.tol = 1e-3)$value)
    
          for (j in rev(seq(x,thres[1] - mar1[1] / mar1[2],length.out=10))){
            ub.alt <- j
            message(paste("Attempting integration with upper bound: ", ub.alt))
            result <- tryCatch({
              integrate(integrand, lower = x, upper = ub.alt, rel.tol = 1e-3)$value
            },
            error = function(e_retry) {
              return(NA)
            })
            # Check if the integration was successful (i.e., didn't return NA).
            if (!is.na(result)) {
              message("Integration successful with a finite upper bound.")
              return(result)
            }
           }
        }
        
        else {
          for (i in rev(seq(0.05,5,0.05)) ) {
            ub.alt <- x + i * mar1[1]
            message(paste("Attempting integration with upper bound: ", ub.alt, " (multiplier =", i, ")"))
            result <- tryCatch({
              integrate(integrand, lower = x, upper = ub.alt, rel.tol = 1e-3)$value
            }, 
            error = function(e_retry) {
              return(NA)
            })
            # Check if the integration was successful (i.e., didn't return NA).
            if (!is.na(result)) {
              message("Integration successful with a finite upper bound.")
              return(result)
            }
          }
          # If the loop completes without a successful return, it means all attempts failed.
          stop("Failed to find a finite upper bound after multiple retries.", call. = FALSE)
        }
      }
      else { stop(e)}
    })
  
  return(R)
}

normalize_c.bivariate<- function(x, EVmodel){
  # computes the infinite integral of the conditional density f(y|X >x)
  # use as a nomralization factor for the conditional density
  dat <- EVmodel$data[EVmodel$data[,2]>EVmodel$threshold[2],2]
  c.y <- function(k) {
    temp <- switch(EVmodel$model,
                   log=c.bivariate(y = k, x = x,  model = EVmodel$model, 
                        dep = EVmodel$estimate[5], thres = EVmodel$threshold, 
                        eta = EVmodel$nat[1:2]/EVmodel$n,
                        mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                        mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])),
                   hr=c.bivariate(y = k, x = x, model = EVmodel$model, 
                                  dep = EVmodel$estimate[5], thres = EVmodel$threshold, 
                                  eta = EVmodel$nat[1:2]/EVmodel$n,
                                  mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                  mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])),
                   ct=c.bivariate(y = k, x = x, model = EVmodel$model, 
                                  alpha = EVmodel$estimate[5],beta=EVmodel$estimate[6], 
                                  thres = EVmodel$threshold, 
                                  eta = EVmodel$nat[1:2]/EVmodel$n,
                                  mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                  mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])))
    return(temp)
  }
  
  C <- tryCatch({
    integrate(Vectorize(c.y), lower = min(dat), upper = Inf)$value},
  error = function(e) {
    if (grepl("non-finite function value", e$message) || grepl("Failed to find a finite upper", e$message)) {
        message("Non-finite function value encountered in norming. Retrying with finite upper bound for y.")
        ub.alt <- ifelse(EVmodel$estimate[4]<0,EVmodel$threshold[2] - EVmodel$estimate[3]
                         /EVmodel$estimate[4],55)
        return(integrate(Vectorize(c.y), lower = min(dat), upper = ub.alt)$value)
      } 
      else {
        stop(e)  # rethrow other errors
      }
    })
  nc <- switch(EVmodel$model,
               log = pbTvevd(q1=x,q2=EVmodel$threshold[2],dep=EVmodel$estimate[5],thres=EVmodel$threshold,model=EVmodel$model,
                eta=EVmodel$nat[1:2]/EVmodel$n,mar1=EVmodel$estimate[1:2],mar2=EVmodel$estimate[3:4],tail.type=2),
               hr = pbTvevd(q1=x,q2=EVmodel$threshold[2],dep=EVmodel$estimate[5],thres=EVmodel$threshold,model=EVmodel$model,
                            eta=EVmodel$nat[1:2]/EVmodel$n,mar1=EVmodel$estimate[1:2],mar2=EVmodel$estimate[3:4],tail.type=2),
               ct = pbTvevd(q1=x,q2=EVmodel$threshold[2],alpha=EVmodel$estimate[5],beta=EVmodel$estimate[6],
                            thres=EVmodel$threshold,model=EVmodel$model,eta=EVmodel$nat[1:2]/EVmodel$n,
                            mar1=EVmodel$estimate[1:2],mar2=EVmodel$estimate[3:4],tail.type=2) )
  return(C/nc)
}

Injury.from_c_bivariate <-function(dat,EVmodel,severity,x0,PX,UB=NULL,LB=NULL){
  # computes the injury probability using c.bivariate, severity is a logistic regression model of injury given speed
  # x0 is the crash boundary
  # upper and lower are the integration limits for impact speed, if NULL, use Inf and min(dat)
  f.y <- function(k) {
    temp <- switch(EVmodel$model,
                   log=c.bivariate(y = k, x = x0, model = EVmodel$model, 
                                   dep = EVmodel$estimate[5], thres = EVmodel$threshold, 
                                   eta = EVmodel$nat[1:2]/EVmodel$n,
                                   mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                   mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])),
                   hr=c.bivariate(y = k, x = x0, model = EVmodel$model, 
                                  dep = EVmodel$estimate[5], thres = EVmodel$threshold, 
                                  eta = EVmodel$nat[1:2]/EVmodel$n,
                                  mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                  mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])),
                   ct=c.bivariate(y = k, x = x0, model = EVmodel$model, 
                                  alpha = EVmodel$estimate[5],beta=EVmodel$estimate[6], 
                                  thres = EVmodel$threshold, 
                                  eta = EVmodel$nat[1:2]/EVmodel$n,
                                  mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                  mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])))
    return(temp*severity(k)/PX * EVmodel$nat[3]/EVmodel$n)
  }
  
  
  R <- tryCatch({
    integrate(Vectorize(f.y), lower = ifelse(is.null(LB),min(dat),LB), 
              upper = ifelse(is.null(UB),Inf,UB))$value},
    error = function(e) {
    if (grepl("non-finite function value", e$message)|| grepl("Failed to find a finite upper", e$message)) {
      message("Non-finite function value encountered. Retrying with finite upper bound.")
      ub.alt <- max(dat, na.rm = TRUE)
      return(integrate(Vectorize(f.y), lower = ifelse(is.null(LB),min(dat),LB), upper = ub.alt)$value)
    } 
    else {
      stop(e)  # rethrow other errors
    }
  })
  
  return(R) # injury probability given a crash
}

Injury.from_c_bivariate_E <-function(dat,EVmodel,severity,x0,PX,UB=NULL,LB=NULL,N=500,age.mean=40,age.sd=15){
  # computes the injury probability using c.bivariate, severity is a logistic regression model of injury given speed
  # x0 is the crash boundary
  # upper and lower are the integration limits for impact speed, if NULL, use Inf and min(dat)
  
  f.y <- function(k) {
    temp <- switch(EVmodel$model,
                   log=c.bivariate(y = k, x = x0, model = EVmodel$model, 
                                   dep = EVmodel$estimate[5], thres = EVmodel$threshold, 
                                   eta = EVmodel$nat[1:2]/EVmodel$n,
                                   mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                   mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])),
                   hr=c.bivariate(y = k, x = x0, model = EVmodel$model, 
                                  dep = EVmodel$estimate[5], thres = EVmodel$threshold, 
                                  eta = EVmodel$nat[1:2]/EVmodel$n,
                                  mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                  mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])),
                   ct=c.bivariate(y = k, x = x0, model = EVmodel$model, 
                                  alpha = EVmodel$estimate[5],beta=EVmodel$estimate[6], 
                                  thres = EVmodel$threshold, 
                                  eta = EVmodel$nat[1:2]/EVmodel$n,
                                  mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                  mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])))
    age <- rnorm(N,mean=age.mean,sd=age.sd)
    age <- age[age>0]
    return(mean(temp*severity(k,age)))
  }
  
  
  R <- tryCatch({
    integrate(Vectorize(f.y), lower = ifelse(is.null(LB),min(dat),LB), 
              upper = ifelse(is.null(UB),60,UB))$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)|| grepl("Failed to find a finite upper", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- max(dat, na.rm = TRUE)
        return(integrate(Vectorize(f.y), lower = ifelse(is.null(LB),min(dat),LB), upper = ub.alt)$value)
      } 
      else {
        stop(e)  # rethrow other errors
      }
    })
  
  return(R/PX/normalize_c.bivariate(x0,EVmodel)) # injury probability given a crash
}

create_plot.df <- function(dat,x,model,PX){
  # create a data frame for plotting the conditional density of speed at X
  # model is a fbvpot object
  df <-switch(model$model,
              log = data.frame(speed= dat,
                   JointP = sapply(dat,FUN = pbTvevd,q1=x,model=model$model,dep=model$estimate[5],
                       thres=model$threshold,
                       eta=model$nat[1:2]/model$n,
                       mar1=c(model$estimate[1],model$estimate[2]),
                       mar2=c(model$estimate[3],model$estimate[4]),tail.type=4)) %>%
                    mutate(ConditionP = JointP/PX) %>%na.omit() %>%
                    mutate(ConditionalD = sapply(speed, 
                       function(k) c.bivariate(y = k, x = x, model=model$model,dep=model$estimate[5],
                           thres=model$threshold,
                           eta=model$nat[1:2]/model$n,
                           mar1=c(model$estimate[1],model$estimate[2]),
                           mar2=c(model$estimate[3],model$estimate[4]))
                       )/PX /normalize_c.bivariate(x,model) * model$nat[3]/model$n),
              hr = data.frame(speed= dat,
                    JointP = sapply(dat,FUN = pbTvevd,q1=x,model=model$model,dep=model$estimate[5],
                       thres=model$threshold,
                       eta=model$nat[1:2]/model$n,
                       mar1=c(model$estimate[1],model$estimate[2]),
                       mar2=c(model$estimate[3],model$estimate[4]),tail.type=4)) %>%
                    mutate(ConditionP = JointP/PX) %>% na.omit() %>%
                    mutate(ConditionalD = sapply(speed, 
                       function(k) c.bivariate(y = k, x = x, model=model$model,dep=model$estimate[5],
                           thres=model$threshold,
                           eta=model$nat[1:2]/model$n,
                           mar1=c(model$estimate[1],model$estimate[2]),
                           mar2=c(model$estimate[3],model$estimate[4]))
                       )/PX/normalize_c.bivariate(x,model) * model$nat[3]/model$n),
              ct = data.frame(speed= dat,
                   JointP = sapply(dat,FUN = pbTvevd,q1=x,model=model$model,alpha=model$estimate[5],
                     beta=model$estimate[6],thres=model$threshold,
                     eta=model$nat[1:2]/model$n,
                     mar1=c(model$estimate[1],model$estimate[2]),
                     mar2=c(model$estimate[3],model$estimate[4]),tail.type=4)) %>%
                   mutate(ConditionP = JointP) %>% na.omit() %>%
                   mutate(ConditionalD = sapply(speed, 
                       function(k) c.bivariate(y = k, x = x, model=model$model,alpha=model$estimate[5],
                           beta=model$estimate[6],thres=model$threshold,
                           eta=model$nat[1:2]/model$n,
                           mar1=c(model$estimate[1],model$estimate[2]),
                           mar2=c(model$estimate[3],model$estimate[4]))
                       )/ PX/normalize_c.bivariate(x,model) * model$nat[3]/model$n)
    )
  return(df)
}

create_plot.df.np <- function(dat,x,PX,mar1,mar2,thres,eta,Ahat){
  df <- data.frame(speed= dat,
             JointP = sapply(dat,FUN = pbTvNonpar,q1=x,Ahat=Ahat,
                             thres=thres,eta=eta,
                             mar1=mar1, mar2=mar2,tail.type=4)) %>%
    mutate(ConditionP = JointP/PX) %>%
    na.omit() %>%
    mutate(ConditionalD = sapply(speed, function(k){
               c.bivariate_np(y = k, x = x,Ahat=Ahat,thres=thres,eta=eta,mar1=mar1,mar2=mar2)})/PX) %>%
        mutate(ConditionalD = ifelse(ConditionalD < 0, 0, ConditionalD))
}

TbevdPlot <- function(EVmodel,dat,k=10){
  # produce four plots: the 2*2 plot
  # k used for bernstein polynomial approximation
  # ... pass to pickands.Nonpar
  A.np.conf <- pickands.Nonpar(dat = dat, mar1 = EVmodel$estimate[1:2], mar2 = EVmodel$estimate[3:4],
                          thres = EVmodel$threshold, eta = EVmodel$nat[1:2]/EVmodel$n,
                          est = "cfg", CI = TRUE, d = 2, k = k, ifplot = TRUE)
  A.np<- pickands.Nonpar(dat = dat, mar1 = EVmodel$estimate[1:2], mar2 = EVmodel$estimate[3:4],
                               thres = EVmodel$threshold, eta = EVmodel$nat[1:2]/EVmodel$n,
                               est = "cfg", CI = FALSE, d = 2, k = k, ifplot = FALSE)
  
  # par(mfrow=c(2,2))
  if (EVmodel$model %in% c("log","hr")){
    CB <- confint(EVmodel,parm="dep")
    abvevd(dep=EVmodel$estimate[5],model=EVmodel$model,plot = TRUE,add=TRUE,col="red")
    abvevd(dep=CB[1],model=EVmodel$model,plot = TRUE,add=TRUE,col="red",lty=5)
    abvevd(dep=CB[2],model=EVmodel$model,plot = TRUE,add=TRUE,col="red",lty=5)
    legend("bottomright", 
           legend = c("Nonparametric estimates", "Parametric estimates"), 
           col = c("black", "red"), lty = c(1, 2), 
           bty = "n")
    title(main = sprintf("Dependence diagonstics %s",EVmodel$model))
    
    spec.dens.plot <- data.frame(t = seq(0,1,0.005)) %>% 
      mutate(h.Nonpar = sapply(t,A_bp_approx,A_bp=A.np$beta,ord="2")/2) %>%
      mutate(h.Nonpar.UP = sapply(t,A_bp_approx,A_bp=A.np.conf$up.beta,ord="2")/2) %>%
      mutate(h.Nonpar.LW = sapply(t,A_bp_approx,A_bp=A.np.conf$low.beta,ord="2")/2) %>%
      mutate(h.Par = sapply(t,hbvevd,dep=EVmodel$estimate[5],model=EVmodel$model,half=TRUE)) %>%
      mutate(h.Par.UP = sapply(t,hbvevd,dep=CB[2],model=EVmodel$model,half=TRUE)) %>%
      mutate(h.Par.LW = sapply(t,hbvevd,dep=CB[1],model=EVmodel$model,half=TRUE))
  } 
  if (EVmodel$model %in% c("ct","bilog","negbilog")){
    CB_alpha <- confint(EVmodel,parm="beta")
    CB_beta <- confint(EVmodel,parm="alpha")
    abvevd(alpha=EVmodel$estimate[6],beta=EVmodel$estimate[5],model=EVmodel$model,
           plot = TRUE,add=TRUE,col="red")
    abvevd(alpha=CB_alpha[1],beta=CB_beta[1],model=EVmodel$model,plot = TRUE,add=TRUE,col="red",lty=5)
    abvevd(alpha=CB_alpha[2],beta=CB_beta[2],model=EVmodel$model,plot = TRUE,add=TRUE,col="red",lty=5)
    legend("bottomright", 
           legend = c("Nonparametric estimates", "Parametric estimates"), 
           col = c("black", "red"), lty = c(1, 2), 
           bty = "n")
    title(main = sprintf("Dependence diagonstics %s",EVmodel$model))
    
    spec.dens.plot <- data.frame(t = seq(0,1,0.005)) %>% 
      mutate(h.Nonpar = sapply(t,A_bp_approx,A_bp=A.np$beta,ord="2")/2) %>%
      mutate(h.Nonpar.UP = sapply(t,A_bp_approx,A_bp=A.np.conf$up.beta,ord="2")/2) %>%
      mutate(h.Nonpar.LW = sapply(t,A_bp_approx,A_bp=A.np.conf$low.beta,ord="2")/2) %>%
      mutate(h.Par = sapply(t,hbvevd,alpha=EVmodel$estimate[5],beta=EVmodel$estimate[6],
                            model=EVmodel$model,half=TRUE)) %>%
      mutate(h.Par.UP = sapply(t,hbvevd,alpha=CB_alpha[2],beta=CB_beta[2],model=EVmodel$model,half=TRUE)) %>%
      mutate(h.Par.LW = sapply(t,hbvevd,alpha=CB_alpha[1],beta=CB_beta[1],model=EVmodel$model,half=TRUE))
  }
  
  ggplot(spec.dens.plot, aes(x=t,y=h.Par,color="Model fitted")) + geom_line() +
    geom_line(aes(y=h.Par.LW,color="Model fitted"),linetype="dashed") +
    geom_line(aes(y=h.Par.UP,color="Model fitted"),linetype="dashed") +
    geom_line(aes(y=h.Nonpar.LW,color="Nonparametric"),linetype="dashed") +
    geom_line(aes(y=h.Nonpar.UP,color="Nonparametric"),linetype="dashed") +
    geom_line(aes(y=h.Nonpar,color="Nonparametric")) +
    labs(x="t",y="Spectral density",title = sprintf("Spectral density %s",EVmodel$model)) +
    scale_color_manual(name="",values=c("Model fitted"="red","Nonparametric"="black")) +
    theme_minimal()
  
}

GoF_BPOT <- function(EVmodel,dat){
  # goodness-of-fit test for bivariate POT models using Cramer-von Mises statistic
  n <- nrow(dat)
  dat <- dat %>% filter(prox > EVmodel$threshold[1], Speed >EVmodel$threshold[2])
  U1 <- mtransform.GPMk2(dat[,1],p=EVmodel$estimate[1:2],thres=EVmodel$threshold[1],
                         eta=EVmodel$nat[1]/EVmodel$n,margin="uniform")
  U2 <- mtransform.GPMk2(dat[,2],p=EVmodel$estimate[3:4],thres=EVmodel$threshold[2],
                         eta=EVmodel$nat[2]/EVmodel$n,margin="uniform")
  model <- switch(EVmodel$model,
                  log = gumbelCopula(dim=2,param = 1/EVmodel$estimate[5]),
                  hr = huslerReissCopula(param = EVmodel$estimate[5]))
  Cop <- fitCopula(copula = model,data = cbind(U1, U2),method = "mpl")
  Result <- gofEVCopula(copula = Cop@copula, method = "itau",x = cbind(U1,U2), N = 1000, estimator = "CFG")
  return(Result)
}  
