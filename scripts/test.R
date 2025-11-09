pExtDep1 <- function (q, type, method = "Parametric", model, par, plot = TRUE, 
          main, xlab, cex.lab, cex.axis, lwd, ...) 
{
  if (is.vector(q)) {
    dim <- length(q)
    q0 <- q
    if (is.matrix(par)) {
      npar <- nrow(par)
      p <- vector(length = npar)
      if (any(q == 0)) {
        return(rep(NA, npar))
      }
    }
    else if (is.vector(par)) {
      if (any(q == 0)) {
        return(NA)
      }
    }
  }
  else if (is.matrix(q)) {
    dim <- ncol(q)
    qNA <- apply(q, 1, function(x) any(x == 0))
    q0 <- q
    nq0 <- nrow(q0)
    q <- q[!qNA, ]
    if (is.matrix(q)) {
      nq <- nrow(q)
    }
    else if (is.vector(q)) {
      nq <- 1
    }
    if (is.vector(par)) {
      p <- vector(length = nq)
    }
    else if (is.matrix(par)) {
      p <- matrix(nrow = nrow(par), ncol = nq)
    }
  }
  else {
    stop("q must be a vector or a matrix")
  }
  methods <- c("Parametric", "NonParametric")
  if (!any(method == methods)) {
    stop("Wrong method specified")
  }
  types <- c("lower", "inv.lower", "upper")
  if (method == "Parametric" && !any(type == types)) {
    stop("Wrong probability type specified")
  }
  if (method == "Parametric") {
    if (all(dim != c(2, 3))) {
      stop("Dimension 2 or 3 only")
    }
    models <- c("HR", "ET", "EST")
    if (!any(model == models)) {
      stop("model wrongly specified")
    }
    if (is.vector(par)) {
      if (!dim_ExtDep(model = model, par = par, dim = dim)) {
        stop("Length of 'par' is incorrect")
      }
    }
    else if (is.matrix(par)) {
      if (!dim_ExtDep(model = model, par = par[1, ], dim = dim)) {
        stop("Length of 'par' is incorrect")
      }
    }
    else {
      stop("'par' should be a vector or a matrix")
    }
    if (model == "HR") {
      if (is.vector(q)) {
        if (is.vector(par)) {
          p <- p.hr(q = q, par = par, type = type)
        }
        else if (is.matrix(par)) {
          p <- apply(par, 1, p.hr, q = q, type = type)
        }
      }
      else if (is.matrix(q)) {
        if (is.vector(par)) {
          p <- apply(q, 1, p.hr, par = par, type = type)
        }
        else if (is.matrix(par)) {
          for (i in 1:nq) {
            p[, i] <- apply(par, 1, p.hr, q = q[i, ], 
                            type = type)
          }
        }
      }
    }
    else if (model == "ET") {
      if (is.vector(q)) {
        if (is.vector(par)) {
          p <- p.et(q = q, par = par, type = type)
        }
        else if (is.matrix(par)) {
          p <- apply(par, 1, p.et, q = q, type = type)
        }
      }
      else if (is.matrix(q)) {
        if (is.vector(par)) {
          p <- apply(q, 1, p.et, par = par, type = type)
        }
        else if (is.matrix(par)) {
          for (i in 1:nq) {
            p[, i] <- apply(par, 1, p.et, q = q[i, ], 
                            type = type)
          }
        }
      }
    }
    else if (model == "EST") {
      if (is.vector(q)) {
        if (is.vector(par)) {
          p <- p.est(q = q, par = par, type = type)
        }
        else if (is.matrix(par)) {
          p <- apply(par, 1, p.est, q = q, type = type)
        }
      }
      else if (is.matrix(q)) {
        if (is.vector(par)) {
          p <- apply(q, 1, p.est, par = par, type = type)
        }
        else if (is.matrix(par)) {
          for (i in 1:nq) {
            p[, i] <- apply(par, 1, p.est, q = q[i, ], 
                            type = type)
          }
        }
      }
    }
  }
  else if (method == "NonParametric") {
    if (is.vector(q)) {
      if (is.vector(par)) {
        p <- ph(w = q, beta = par)
      }
      else if (is.matrix(par)) {
        for (i in 1:npar) {
          p[i] <- ph(w = q, beta = par[i, ])
        }
      }
    }
    else if (is.matrix(q)) {
      if (is.vector(par)) {
        for (i in 1:nq) {
          p[i] <- ph(w = q[i, ], beta = par)
        }
      }
      else if (is.matrix(par)) {
        for (i in 1:nq) {
          for (j in 1:npar) {
            p[j, i] <- ph(w = q[i, ], beta = par[j, ])
          }
        }
      }
    }
  }
  if (is.matrix(par)) {
    if (plot) {
      if (missing(main)) {
        main <- ""
      }
      if (missing(cex.lab)) {
        cex.lab <- 1.4
      }
      if (missing(cex.axis)) {
        cex.axis <- 1.4
      }
      if (missing(lwd)) {
        lwd <- 2
      }
      if (is.vector(p)) {
        if (missing(xlab)) {
          if (dim == 2) {
            if (type == "lower") {
              xlab <- paste("P(X<", q[1], ",Y<", q[2], 
                            ")", sep = "")
            }
            else if (type == "inv.lower") {
              xlab <- paste("1-P(X<", q[1], ",Y<", q[2], 
                            ")", sep = "")
            }
            else if (type == "upper") {
              xlab <- paste("P(X>", q[1], ",Y>", q[2], 
                            ")", sep = "")
            }
          }
          else if (dim == 3) {
            if (type == "lower") {
              xlab <- paste("P(X<", q[1], ",Y<", q[2], 
                            ",Z<", q[3], ")", sep = "")
            }
            else if (type == "inv.lower") {
              xlab <- paste("1-P(X<", q[1], ",Y<", q[2], 
                            ",Z<", q[3], ")", sep = "")
            }
            else if (type == "upper") {
              xlab <- paste("P(X>", q[1], ",Y>", q[2], 
                            ",Z>", q[3], ")", sep = "")
            }
          }
        }
        Ke <- density(p)
        Hi <- hist(p, prob = TRUE, col = "lightgrey", 
                   ylim = range(Ke$y), main = main, xlab = xlab, 
                   cex.lab = cex.lab, cex.axis = cex.axis, lwd = lwd, 
                   ...)
        p_ic <- quantile(p, probs = c(0.025, 0.5, 0.975))
        points(x = p_ic, y = c(0, 0, 0), pch = 4, lwd = 4)
        points(x = mean(p), y = 0, pch = 16, lwd = 4)
        lines(Ke, lwd = 2, col = "dimgrey")
      }
      else if (is.matrix(p)) {
        for (i in 1:nq) {
          if (missing(xlab)) {
            if (dim == 2) {
              if (type == "lower") {
                xlab <- paste("P(X<", q[i, 1], ",Y<", 
                              q[i, 2], ")", sep = "")
              }
              else if (type == "inv.lower") {
                xlab <- paste("1-P(X<", q[i, 1], ",Y<", 
                              q[i, 2], ")", sep = "")
              }
              else if (type == "upper") {
                xlab <- paste("P(X>", q[i, 1], ",Y>", 
                              q[i, 2], ")", sep = "")
              }
            }
            else if (dim == 3) {
              if (type == "lower") {
                xlab <- paste("P(X<", q[i, 1], ",Y<", 
                              q[i, 2], ",Z<", q[i, 3], ")", sep = "")
              }
              else if (type == "inv.lower") {
                xlab <- paste("1-P(X<", q[i, 1], ",Y<", 
                              q[i, 2], ",Z<", q[i, 3], ")", sep = "")
              }
              else if (type == "upper") {
                xlab <- paste("P(X>", q[i, 1], ",Y>", 
                              q[i, 2], ",Z>", q[i, 3], ")", sep = "")
              }
            }
          }
          Ke <- density(p[, i])
          Hi <- hist(p[, i], prob = TRUE, col = "lightgrey", 
                     ylim = range(Ke$y), main = main, xlab = xlab, 
                     cex.lab = cex.lab, cex.axis = cex.axis, lwd = lwd, 
                     ...)
          p_ic <- quantile(p[, i], probs = c(0.025, 0.5, 
                                             0.975))
          points(x = p_ic, y = c(0, 0, 0), pch = 4, lwd = 4)
          points(x = mean(p[, i]), y = 0, pch = 16, lwd = 4)
          lines(Ke, lwd = 2, col = "dimgrey")
        }
      }
    }
  }
  if (is.matrix(q0)) {
    if (sum(qNA) != 0) {
      p.temp <- p
      if (is.vector(par)) {
        p <- vector(length = nq0)
        p[!qNA] <- p.temp
        p[qNA] <- rep(NA, sum(qNA))
      }
      else if (is.matrix(par)) {
        p <- matrix(nrow = nrow(par), ncol = nq0)
        p[, !qNA] <- p.temp
      }
    }
  }
  return(p)
}

ph1 <- function (w, beta) 
{
  k <- length(beta) - 1
  j <- 1:k
  res <- 0.5 * (diff(beta) + 1/k) * dbeta(w, j, k - j + 1)
  return(sum(res))
}
dh1 <- function (w, beta, mixture = FALSE) 
{
  k <- length(beta) - 1
  j <- 1:(k - 1)
  const <- 2/(k * (2 - beta[2] - beta[k]))
  res <- diff(diff(beta)) * dbeta(w, j, k - j)
  if (mixture) {
    return(k/2 * const * sum(res))
  }
  else {
    return(k/2 * sum(res))
  }
}

function (x, method = "Parametric", model, par, angular = TRUE, 
          log = FALSE, c = NULL, vectorial = TRUE, mixture = FALSE) 
{
  methods <- c("Parametric", "NonParametric")
  if (!any(method == methods)) {
    stop("Wrong method specified")
  }
  if (method == "Parametric") {
    if (angular) {
      models <- c("PB", "HR", "ET", "EST", "TD", "AL")
      if (!any(model == models)) {
        stop("model wrongly specified")
      }
      if (is.vector(x)) {
        d <- as.integer(length(x))
        if (round(sum(x), 7) != 1) {
          stop("'x' should be a vector in the unit simplex")
        }
      }
      else if (is.matrix(x)) {
        d <- as.integer(ncol(x))
        if (!all(round(rowSums(x), 7) == 1)) {
          stop("'x' should be a matrix with row vectors belonging to the unit simplex")
        }
      }
      else {
        stop(" 'x' should be a vector or a matrix")
      }
      if (model %in% models[-c(1, 2, 5)]) {
        if (d > 3) {
          stop("The angular density is only available for dimension 2 or 3.")
        }
      }
      dim.check <- dim_ExtDep(model = model, par = par, 
                              dim = d)
      if (!dim.check) {
        stop("Wrong length of parameters")
      }
      if (model == "PB") {
        return(dens_pb(x = x, b = par[1:choose(d, 2)], 
                       alpha = par[choose(d, 2) + 1], log = log, vectorial = vectorial))
      }
      if (model == "HR") {
        return(dens_hr(x = x, lambda = par, log = log, 
                       vectorial = vectorial))
      }
      if (model == "TD") {
        return(dens_di(x = x, para = par, log = log, 
                       vectorial = vectorial))
      }
      if (model == "ET") {
        if (is.null(c)) {
          stop("c needs to be specified")
        }
        return(dens_et(x = x, rho = par[1:choose(d, 2)], 
                       mu = par[choose(d, 2) + 1], c = c, log = log, 
                       vectorial = vectorial))
      }
      if (model == "EST") {
        if (is.null(c)) {
          stop("c needs to be specified")
        }
        return(dens_est(x = x, rho = par[1:choose(d, 
                                                  2)], alpha = par[choose(d, 2) + 1:d], mu = par[choose(d, 
                                                                                                        2) + d + 1], c = c, log = log, vectorial = vectorial))
      }
      if (model == "AL") {
        if (is.null(c)) {
          stop("c needs to be specified")
        }
        if (d == 2) {
          return(dens_al(x = x, alpha = par[1], beta = par[2:3], 
                         c = c, log = log, vectorial = vectorial))
        }
        if (d == 3) {
          return(dens_al(x = x, alpha = par[1:4], beta = par[5:13], 
                         c = c, log = log, vectorial = vectorial))
        }
      }
    }
    else {
      models <- c("HR", "ET", "EST")
      if (!any(model == models)) {
        stop("model wrongly specified")
      }
      if (is.vector(x)) {
        d <- as.integer(length(x))
      }
      else if (is.matrix(x)) {
        d <- as.integer(ncol(x))
      }
      else {
        stop(" 'x' should be a vector or a matrix")
      }
      if (d != 2) {
        stop("Density of parametric models only available for d=2")
      }
      dim.check <- dim_ExtDep(model = model, par = par, 
                              dim = d)
      if (!dim.check) {
        stop("Wrong length of parameters")
      }
      if (model == "HR") {
        if (is.vector(x)) {
          out <- dHuslerReiss(x = x, lambda = par)
        }
        else if (is.matrix(x)) {
          out <- apply(x, 1, function(z) dHuslerReiss(x = z, 
                                                      lambda = par))
        }
      }
      if (model == "ET") {
        if (is.vector(x)) {
          out <- dmextst(x = x, scale = par[1], df = par[2])
        }
        else if (is.matrix(x)) {
          out <- apply(x, 1, function(z) dmextst(x = z, 
                                                 scale = par[1], df = par[2]))
        }
      }
      if (model == "EST") {
        if (is.vector(x)) {
          out <- dmextst(x = x, scale = par[1], shape = par[2:3], 
                         df = par[4])
        }
        else if (is.matrix(x)) {
          out <- apply(x, 1, function(z) dmextst(x = z, 
                                                 scale = par[1], shape = par[2:3], df = par[4]))
        }
      }
      if (log) {
        out <- log(out)
      }
      if (!vectorial) {
        if (log) {
          out <- sum(out)
        }
        else {
          out <- prod(out)
        }
      }
      return(out)
    }
  }
  else if (method == "NonParametric") {
    if (angular) {
      if (is.vector(x)) {
        if (length(x) != 2) {
          stop("'x' should of length 2")
        }
        if (sum(x) != 1) {
          stop("'x' should be a vector in the unit simplex")
        }
        dens <- dh(w = x[1], beta = par, mixture = mixture)
      }
      else if (is.matrix(x)) {
        if (ncol(x) != 2) {
          stop("'x' should have 2 columns")
        }
        if (!all(rowSums(x) == 1)) {
          stop("'x' should be a matrix with row vectors belonging to the unit simplex")
        }
        dens <- apply(x, 1, function(y) {
          dh(w = y[1], beta = par, mixture = mixture)
        })
      }
      else {
        stop(" 'x' should be a vector or a matrix")
      }
      if (log) {
        dens <- log(dens)
      }
      if (!vectorial) {
        if (log) {
          dens <- sum(dens)
        }
        else {
          dens <- prod(dens)
        }
      }
      return(dens)
    }
    else {
      stop("Only the angular density available in non-parametric form")
    }
  }
}


dbTvNonpar1 <- function(x1,x2,Ahat){
  
  ext <- c((x1 %in% c(0, Inf)),(x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    Z <- x1 + x2
    X <- x1*x2
    G <- 0.999
    w <- x1/Z
    A0 <- A_bp_approx(Ahat,t = w,ord="0")
    A1 <- A_bp_approx(Ahat,t = w,ord="1")
    A2 <- A_bp_approx(Ahat,t = w,ord="2")
    #du1 <- A0 + (1-w)*A1
    #du2 <- A0 - w*A1
    #du12 <- -A2*(1-w) * w/Z
    # d <- G/X * (du1*du2 - du12) 
    d <- G* (A2/Z^3+(A0^2 + X*A0*A1*(x2-x1)/Z^2 - X*A1^2/Z^2)/X^2 )
  }
  d
}
c.bivariate1 <-function(y,x,Ahat){
  # approximates the conditional density f(y|X >x) by integrating fxy = f(x=x, y = y) over x
  # if integral is non-finite, change ulim to a smaller value
  integrand <- function(q1,...) {
    result <- numeric(length(q1))
    for(i in seq_along(q1)) {
      result[i] <- dbTvNonpar1(x1=q1[i],x2=y,Ahat=Ahat)
    }
    return(result)
  }
  
  # Integrate from x to Inf
  R <- tryCatch({
    integrate(integrand, lower = x, upper = Inf,rel.tol = 1e-3)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered in conditional density. Retrying with finite upper bound for x.")
        ub.alt <- x + 100
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
          
          # If the loop completes without a successful return, it means all attempts failed.
          stop("Failed to find a finite upper bound after multiple retries.", call. = FALSE)
        }
      else { stop(e)}
    })
  
  return(R)
}

normalize_c.bivariate1<- function(x,dat,Ahat){
  # computes the infinite integral of the conditional density f(y|X >x)
  # use as a nomralization factor for the conditional density

  c.y <- function(k) {
    temp <- c.bivariate1(y = k, x = x,Ahat=Ahat)
    return(temp)
  }
  
  nc <- pbTvNonpar1(x1=x,x2=min(dat),Ahat=Ahat,tail.type=2)
  C <- tryCatch({
    integrate(Vectorize(c.y), lower = min(dat), upper = Inf)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message) || grepl("Failed to find a finite upper", e$message)) {
        message("Non-finite function value encountered in norming. Retrying with finite upper bound for y.")
        ub.alt <- max(dat)
        return(integrate(Vectorize(c.y), lower = min(dat), upper = ub.alt)$value)
      } 
      else {
        stop(e)  # rethrow other errors
      }
    })
  
  return(C/nc)
}

pbTvNonpar1 <- function(x1,x2,Ahat,tail.type){
  # Ahat is a vector of non-par estimates obatined from pickands.Nonpar
  w <- x1/(x1+x2)
  v <- A_bp_approx(Ahat,t=w,ord="0")
  
  pp <- exp(-(1/x1+1/x2)*v)
  # P(X>q1,Y>q2)
  if (tail.type==2) {
    pp <- 1- pgev(x1,loc = 1,shape=1)-pgev(x2,loc = 1,shape=1) +pp
  }
  # P(X<=q1,Y>q2)
  else if (tail.type==3) {
    pp <- pgev(x1,loc = 1,shape=1) - pp
  }
  # P(X>q1,Y<= q2)
  else if (tail.type==4) {
    pp <- pgev(x2,loc = 1,shape=1)- pp
  }
  pp
}
