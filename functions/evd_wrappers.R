# wrappers for functions from `evd`, `ExtremalDep`

normalize_eta <- function(eta) {
  if (length(eta) == 1) rep(eta, 2) else eta
}

tail_adjust <- function(v, x1, x2, tail_type) {
  pp <- exp(-v)
  switch(as.character(tail_type),
    "2" = 1 - pgev(-log(x1)) - pgev(-log(x2)) + pp,
    "3" = pgev(-log(x1)) - pp,
    "4" = pgev(-log(x2)) - pp)
}

mtransform_gp_mk2 <- function(x, p, thres, eta, margin = "exp") {
  if (is.list(p)) {
    if (is.null(dim(x)) && length(x) != length(p))
      stop(paste("`p' must have", length(x), "elements"))
    if (!is.null(dim(x)) && ncol(x) != length(p))
      stop(paste("`p' must have", ncol(x), "elements"))
    if (is.null(dim(x)))
      dim(x) <- c(1, length(p))
    for (i in 1:length(p))
      x[, i] <- mtransform_gp_mk2(x[, i], p[[i]], thres[i], eta)
    if (ncol(x) == 1 || (nrow(x) == 1))
      x <- drop(x)
    return(x)
  }
  if (is.null(dim(x)))
    dim(x) <- c(length(x), 1)
  p <- matrix(t(p), nrow = nrow(x), ncol = 2, byrow = TRUE)
  if (min(p[, 1]) <= 0)
    stop("invalid marginal scale")
  expind <- (p[, 2] == 0)
  nzshapes <- p[!expind, 2]

  x <- (x - thres) / p[, 1]
  if (any(x < 0))
    stop("input below thresholds")
  Fx <- ifelse(expind, 1 - eta * exp(-x),
               pmax(1 - eta * (1 + nzshapes * x)^(-1 / nzshapes), 0))
  x_t <- switch(margin,
    exp     = -log(Fx),
    frechet = -1 / log(Fx),
    uniform = Fx,
    stop("invalid margin type"))

  x_t
}

mtransform_gp_mk3 <- function(x, p, thres, eta, margin = "exp") {
  if (is.list(p)) {
    if (is.null(dim(x)) && length(x) != length(p))
      stop(paste("`p' must have", length(x), "elements"))
    if (!is.null(dim(x)) && ncol(x) != length(p))
      stop(paste("`p' must have", ncol(x), "elements"))
    if (is.null(dim(x)))
      dim(x) <- c(1, length(p))
    for (i in 1:length(p))
      x[, i] <- mtransform_gp_mk3(x[, i], p[[i]], thres[i], eta)
    if (ncol(x) == 1 || (nrow(x) == 1))
      x <- drop(x)
    return(x)
  }
  if (is.null(dim(x)))
    dim(x) <- c(length(x), 1)
  p <- matrix(t(p), nrow = nrow(x), ncol = 2, byrow = TRUE)
  if (min(p[, 1]) <= 0)
    stop("invalid marginal scale")
  expind <- (p[, 2] == 0)
  nzshapes <- p[!expind, 2]

  x <- (x - thres) / p[, 1]
  x[x < 0] <- 0
  Fx <- ifelse(expind, 1 - eta * exp(-x),
               pmax(1 - eta * (1 + nzshapes * x)^(-1 / nzshapes), 0))
  x_t <- switch(margin,
    exp     = -log(Fx),
    frechet = -1 / log(Fx),
    uniform = Fx,
    stop("invalid margin type"))

  x_t
}


pb_tvevd <- function(q1, q2, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix", "pb", "nonpar"), ...) {
  model <- match.arg(model)
  switch(model,
    log      = pb_tvlog(q1, q2, ...),
    alog     = pb_tvalog(q1, q2, ...),
    hr       = pb_tvhr(q1, q2, ...),
    neglog   = pb_tvneglog(q1, q2, ...),
    aneglog  = pb_tvaneglog(q1, q2, ...),
    bilog    = pb_tvbilog(q1, q2, ...),
    negbilog = pb_tvnegbilog(q1, q2, ...),
    ct       = pb_tvct(q1, q2, ...),
    nonpar   = pb_tv_nonpar(q1, q2, ...))
}

pb_tvlog <- function(q1, q2, dep, mar1, mar2, tail_type, thres, eta) {
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 || dep > 1)
    stop("invalid argument for `dep'")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, p = mar1, thres = thres[1], eta[1])
  x2 <- mtransform_gp_mk2(q2, p = mar2, thres = thres[2], eta[2])
  v <- (x1^(1 / dep) + x2^(1 / dep))^dep
  tail_adjust(v, x1, x2, tail_type)
}

pb_tvhr <- function(q1, q2, dep, mar1, mar2, tail_type, thres, eta) {
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
    stop("invalid argument for `dep'")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, p = mar1, thres = thres[1], eta[1])
  x2 <- mtransform_gp_mk2(q2, p = mar2, thres = thres[2], eta[2])
  fn <- function(x1, x2) { x1 * pnorm(1 / dep + dep * log(x1 / x2) / 2) }
  v <- fn(x1, x2) + fn(x2, x1)
  tail_adjust(v, x1, x2, tail_type)
}

pb_tvct <- function(q1, q2, alpha, beta, mar1, mar2, tail_type, thres, eta) {
  if (length(alpha) != 1 || mode(alpha) != "numeric")
    stop("invalid argument for `alpha'")
  if (length(beta) != 1 || mode(beta) != "numeric")
    stop("invalid argument for `beta'")
  if (any(c(alpha, beta) <= 0))
    stop("`alpha' and `beta' must be non-negative")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, p = mar1, thres = thres[1], eta[1])
  x2 <- mtransform_gp_mk2(q2, p = mar2, thres = thres[2], eta[2])
  u <- (alpha * x2) / (alpha * x2 + beta * x1)
  v <- pbeta(u, shape1 = alpha, shape2 = beta + 1) * x2 +
        pbeta(u, shape1 = alpha + 1, shape2 = beta,
              lower.tail = FALSE) * x1
  if (x1 + x2 == 0) { v <- 0 }
  if (is.infinite(x1) || is.infinite(x2)) { v <- Inf }
  tail_adjust(v, x1, x2, tail_type)
}

pb_tvneglog <- function(q1, q2, dep, mar1, mar2, tail_type, thres, eta) {
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
    stop("invalid argument for `dep'")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, p = mar1, thres = thres[1], eta[1])
  x2 <- mtransform_gp_mk2(q2, p = mar2, thres = thres[2], eta[2])
  v <- x1 + x2 - (x1^(-dep) + x2^(-dep))^(-1 / dep)
  tail_adjust(v, x1, x2, tail_type)
}

pb_tvalog <- function(q1, q2, dep, asy = c(1, 1), mar1, mar2, tail_type,
                      thres, eta) {
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 || dep > 1)
    stop("invalid argument for `dep'")
  if (length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
      max(asy) > 1)
    stop("invalid argument for `asy'")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, p = mar1, thres = thres[1], eta[1])
  x2 <- mtransform_gp_mk2(q2, p = mar2, thres = thres[2], eta[2])
  v <- ((asy[1] * x1)^(1 / dep) + (asy[2] * x2)^(1 / dep))^dep +
       (1 - asy[1]) * x1 + (1 - asy[2]) * x2
  tail_adjust(v, x1, x2, tail_type)
}

pb_tvaneglog <- function(q1, q2, dep, asy = c(1, 1), mar1, mar2, tail_type,
                         thres, eta) {
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
    stop("invalid argument for `dep'")
  if (length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
      max(asy) > 1)
    stop("invalid argument for `asy'")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, p = mar1, thres = thres[1], eta[1])
  x2 <- mtransform_gp_mk2(q2, p = mar2, thres = thres[2], eta[2])
  v <- x1 + x2 - ((asy[1] * x1)^(-dep) + (asy[2] * x2)^(-dep))^(-1 / dep)
  tail_adjust(v, x1, x2, tail_type)
}

pb_tvbilog <- function(q1, q2, alpha, beta, mar1, mar2, tail_type,
                       thres, eta) {
  if (length(alpha) != 1 || mode(alpha) != "numeric")
    stop("invalid argument for `alpha'")
  if (length(beta) != 1 || mode(beta) != "numeric")
    stop("invalid argument for `beta'")
  if (any(c(alpha, beta) <= 0) || any(c(alpha, beta) >= 1))
    stop("`alpha' and `beta' must be in the open interval (0,1)")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, p = mar1, thres = thres[1], eta[1])
  x2 <- mtransform_gp_mk2(q2, p = mar2, thres = thres[2], eta[2])
  gmafn <- function(z) (1 - alpha) * x1 * (1 - z)^beta -
                        (1 - beta) * x2 * z^alpha
  gma <- uniroot(gmafn, lower = 0, upper = 1,
                 tol = .Machine$double.eps^0.5)$root
  v <- x1 * gma^(1 - alpha) + x2 * (1 - gma)^(1 - beta)
  tail_adjust(v, x1, x2, tail_type)
}

pb_tvnegbilog <- function(q1, q2, alpha, beta, mar1, mar2, tail_type,
                          thres, eta) {
  if (length(alpha) != 1 || mode(alpha) != "numeric")
    stop("invalid argument for `alpha'")
  if (length(beta) != 1 || mode(beta) != "numeric")
    stop("invalid argument for `beta'")
  if (any(c(alpha, beta) <= 0))
    stop("`alpha' and `beta' must be non-negative")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, p = mar1, thres = thres[1], eta[1])
  x2 <- mtransform_gp_mk2(q2, p = mar2, thres = thres[2], eta[2])
  gmafn <- function(z) (1 + alpha) * x1 * z^alpha -
                        (1 + beta) * x2 * (1 - z)^beta
  gma <- uniroot(gmafn, lower = 0, upper = 1,
                 tol = .Machine$double.eps^0.5)$root
  v <- x1 + x2 - x1 * gma^(1 + alpha) - x2 * (1 - gma)^(1 + beta)
  tail_adjust(v, x1, x2, tail_type)
}

pickands_nonpar <- function(dat, mar1, mar2, thres, eta, est = "cfg", CI = FALSE,
                            d = 2, N = 100, k = 10, ifplot = FALSE,
                            nboot = 500, alpha = 0.05) {
  eta <- normalize_eta(eta)
  dat <- dat[dat[, 1] > thres[1] & dat[, 2] > thres[2], ]
  dat[, 1] <- mtransform_gp_mk2(dat[, 1], p = mar1, thres = thres[1], eta[1],
                                margin = "frechet")
  dat[, 2] <- mtransform_gp_mk2(dat[, 2], p = mar2, thres = thres[2], eta[2],
                                margin = "frechet")
  S <- simplex(2, N)
  if (!CI) {
    bp_est <- beed(data = dat, x = S, d = d, est = est, margin = "frechet",
                   k = k, plot = ifplot)
    return(bp_est)
  } else {
    bp_est <- beed.confband(data = dat, x = S, d = d, est = est,
                            margin = "frechet", conf = 1 - alpha,
                            k = k, plot = ifplot)
    return(bp_est)
  }
}

pb_tv_nonpar <- function(q1, q2, mar1, mar2, tail_type, thres, eta, Ahat) {
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, p = mar1, thres = thres[1], eta[1], margin = "exp")
  x2 <- mtransform_gp_mk2(q2, p = mar2, thres = thres[2], eta[2], margin = "exp")
  w <- x1 / (x1 + x2)
  v <- a_bp_approx(Ahat, t = w, ord = "0")
  tail_adjust(v, x1, x2, tail_type)
}

a_bp_approx <- function(A_bp, t, ord) {
  beta_hat <- A_bp
  k <- length(beta_hat)
  beta_hat_k <- beta_hat[2:k]
  b_poly <- function(x, j, k) {
    choose(k, j) * x^j * (1 - x)^(k - j)
  }
  A_prox <- switch(ord,
    "0" = sapply(1:(k - 1), function(j) b_poly(t, j, k - 1)) %*% beta_hat_k +
           b_poly(t, 0, k - 1),
    "1" = ExtremalDep:::ph(t, beta_hat) * 2 - 1,
    "2" = ExtremalDep:::dh(t, beta_hat) * 2
  )
  return(A_prox)
}

# density functions

db_tvevd <- function(q1, q2, model = c("log", "alog", "hr", "neglog", "aneglog",
    "bilog", "negbilog", "ct", "amix", "nonpar"), ...) {
  model <- match.arg(model)
  switch(model,
    log    = db_tvlog(q1, q2, ...),
    alog   = db_tvalog(q1, q2, ...),
    hr     = db_tvhr(q1, q2, ...),
    neglog = db_tvneglog(q1, q2, ...),
    aneglog = db_tvaneglog(q1, q2, ...),
    bilog  = db_tvbilog(q1, q2, ...),
    negbilog = db_tvnegbilog(q1, q2, ...),
    ct     = db_tvct(q1, q2, ...),
    amix   = db_tvamix(q1, q2, ...),
    nonpar = db_tv_nonpar(q1, q2, ...))
}

db_tvlog <- function(q1, q2, dep, mar1, mar2, thres, eta, log = FALSE) {
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 || dep > 1)
    stop("invalid argument for `dep'")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, mar1, thres[1], eta[1])
  x2 <- mtransform_gp_mk2(q2, mar2, thres[2], eta[2])
  ext <- c((x1 %in% c(0, Inf)), (x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    idep <- 1 / dep
    z <- (x1^idep + x2^idep)^dep
    lx <- log(c(x1, x2))
    .expr1 <- (idep + mar1[2]) * lx[1] + (idep + mar2[2]) * lx[2] -
              log(mar1[1] * mar2[1])
    d <- .expr1 + (1 - 2 * idep) * log(z) + log(idep - 1 + z) - z
  }
  if (!log) d <- exp(d)
  d
}

db_tvhr <- function(q1, q2, dep, mar1, mar2, thres, eta, margin = "exp",
                    log = FALSE) {
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
    stop("invalid argument for `dep'")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, mar1, thres[1], eta[1], margin = margin)
  x2 <- mtransform_gp_mk2(q2, mar2, thres[2], eta[2], margin = margin)
  ext <- c((x1 %in% c(0, Inf)), (x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    fn <- function(x1, x2, nm = pnorm) {
      x1 * nm(1 / dep + dep * log(x1 / x2) / 2)
    }
    v <- fn(x1, x2) + fn(x2, x1)
    lx <- log(c(x1, x2))
    .expr1 <- fn(x1, x2) * fn(x2, x1) + dep * fn(x1, x2, nm = dnorm) / 2
    jac <- mar1[2] * lx[1] + mar2[2] * lx[2] - log(mar1[1] * mar2[1])
    d <- log(.expr1) + jac - v
  }
  if (!log) d <- exp(d)
  d
}

db_tvct <- function(q1, q2, alpha, beta, mar1, mar2, thres, eta, log = FALSE) {
  if (length(alpha) != 1 || mode(alpha) != "numeric")
    stop("invalid argument for `alpha'")
  if (length(beta) != 1 || mode(beta) != "numeric")
    stop("invalid argument for `beta'")
  if (any(c(alpha, beta) <= 0))
    stop("`alpha' and `beta' must be non-negative")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, mar1, thres[1], eta[1], margin = "exp")
  x2 <- mtransform_gp_mk2(q2, mar2, thres[2], eta[2], margin = "exp")
  ext <- c((x1 %in% c(0, Inf)), (x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    u <- (alpha * x2) / (alpha * x2 + beta * x1)
    v <- x2 * pbeta(u, shape1 = alpha, shape2 = beta + 1) +
         x1 * pbeta(u, shape1 = alpha + 1, shape2 = beta,
                    lower.tail = FALSE)
    lx <- log(c(x1, x2))
    jac <- (1 + mar1[2]) * lx[1] + (1 + mar2[2]) * lx[2] -
           log(mar1[1] * mar2[1])
    .c1 <- alpha * beta / (alpha + beta + 1)
    .expr1 <- pbeta(u, shape1 = alpha, shape2 = beta + 1) *
              pbeta(u, shape1 = alpha + 1, shape2 = beta,
                    lower.tail = FALSE)
    .expr2 <- dbeta(u, shape1 = alpha + 1, shape2 = beta + 1) /
              (alpha * x2 + beta * x1)
    d <- log(.expr1 + .c1 * .expr2) - v + jac
  }
  if (!log) d <- exp(d)
  d
}

db_tvneglog <- function(q1, q2, dep, mar1, mar2, thres, eta, log = FALSE) {
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
    stop("invalid argument for `dep'")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, mar1, thres[1], eta[1], margin = "exp")
  x2 <- mtransform_gp_mk2(q2, mar2, thres[2], eta[2], margin = "exp")
  ext <- c((x1 %in% c(0, Inf)), (x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    idep <- 1 / dep
    z <- (x1^(-dep) + x2^(-dep))^(-idep)
    v <- x1 + x2 - z
    lx <- log(c(x1, x2))
    fx <- (-dep - 1) * lx
    jac <- (1 + mar1[2]) * lx[1] + (1 + mar2[2]) * lx[2] -
           log(mar1[1] * mar2[1])
    .expr1 <- (1 + dep) * log(z) + log(exp(fx[1]) + exp(fx[2]))
    .expr2 <- fx[1] + fx[2] + (1 + 2 * dep) * log(z) + log(1 + dep + z)
    d <- log(1 - exp(.expr1) + exp(.expr2)) - v + jac
  }
  if (!log) d <- exp(d)
  d
}

db_tvalog <- function(q1, q2, dep, asy = c(1, 1), mar1, mar2, thres, eta,
                      log = FALSE) {
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 0 || dep > 1)
    stop("invalid argument for `dep'")
  if (length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
      max(asy) > 1)
    stop("invalid argument for `asy'")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, mar1, thres[1], eta[1], margin = "exp")
  x2 <- mtransform_gp_mk2(q2, mar2, thres[2], eta[2], margin = "exp")
  ext <- c((x1 %in% c(0, Inf)), (x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    idep <- 1 / dep
    z <- ((asy[1] * x1)^idep + (asy[2] * x2)^idep)^dep
    v <- z + (1 - asy[1]) * x1 + (1 - asy[2]) * x2
    f1asy <- idep * log(asy)
    f2asy <- log(1 - asy)
    lx <- log(c(x1, x2))
    fx <- (idep - 1) * lx
    jac <- (1 + mar1[2]) * lx[1] + (1 + mar2[2]) * lx[2] -
           log(mar1[1] * mar2[1])
    .expr1 <- f2asy[1] + f2asy[2]
    .expr2 <- f2asy[1] + f1asy[2] + fx[2]
    .expr3 <- f2asy[2] + f1asy[1] + fx[1]
    .expr4 <- (1 - idep) * log(z) + log(exp(.expr2) + exp(.expr3))
    .expr5 <- f1asy[1] + f1asy[2] + fx[1] + fx[2] +
              (1 - 2 * idep) * log(z) + log(idep - 1 + z)
    d <- log(exp(.expr1) + exp(.expr4) + exp(.expr5)) - v + jac
  }
  if (!log) d <- exp(d)
  d
}

db_tvaneglog <- function(q1, q2, dep, asy = c(1, 1), mar1, mar2, thres, eta,
                         log = FALSE) {
  if (length(dep) != 1 || mode(dep) != "numeric" || dep <= 0)
    stop("invalid argument for `dep'")
  if (length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
      max(asy) > 1)
    stop("invalid argument for `asy'")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, mar1, thres[1], eta[1], margin = "exp")
  x2 <- mtransform_gp_mk2(q2, mar2, thres[2], eta[2], margin = "exp")
  ext <- c((x1 %in% c(0, Inf)), (x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    idep <- 1 / dep
    z <- ((asy[1] * x1)^(-dep) + (asy[2] * x2)^(-dep))^(-idep)
    v <- x1 + x2 - z
    fasy <- (-dep) * log(asy)
    lx <- log(c(x1, x2))
    fx <- (-dep - 1) * lx
    jac <- (1 + mar1[2]) * lx[1] + (1 + mar2[2]) * lx[2] -
           log(mar1[1] * mar2[1])
    .expr1 <- fasy[1] + fx[1]
    .expr2 <- fasy[2] + fx[2]
    .expr3 <- (1 + dep) * log(z) + log(exp(.expr1) + exp(.expr2))
    .expr4 <- fasy[1] + fasy[2] + fx[1] + fx[2] +
              (1 + 2 * dep) * log(z) + log(1 + dep + z)
    d <- log(1 - exp(.expr3) + exp(.expr4)) - v + jac
  }
  if (!log) d <- exp(d)
  d
}

db_tvbilog <- function(q1, q2, alpha, beta, mar1, mar2, thres, eta,
                       log = FALSE) {
  if (length(alpha) != 1 || mode(alpha) != "numeric")
    stop("invalid argument for `alpha'")
  if (length(beta) != 1 || mode(beta) != "numeric")
    stop("invalid argument for `beta'")
  if (any(c(alpha, beta) <= 0) || any(c(alpha, beta) >= 1))
    stop("`alpha' and `beta' must be in the open interval (0,1)")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, mar1, thres[1], eta[1], margin = "exp")
  x2 <- mtransform_gp_mk2(q2, mar2, thres[2], eta[2], margin = "exp")
  ext <- c((x1 %in% c(0, Inf)), (x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    gmafn <- function(z) (1 - alpha) * x1 * (1 - z)^beta -
                          (1 - beta) * x2 * z^alpha
    gma <- uniroot(gmafn, lower = 0, upper = 1,
                   tol = .Machine$double.eps^0.5)$root
    v <- x1 * gma^(1 - alpha) + x2 * (1 - gma)^(1 - beta)
    lx <- log(c(x1, x2))
    jac <- (1 + mar1[2]) * lx[1] + (1 + mar2[2]) * lx[2] -
           log(mar1[1] * mar2[1])
    .expr1 <- exp((1 - alpha) * log(gma) + (1 - beta) * log(1 - gma))
    .expr2 <- exp(log(1 - alpha) + log(beta) + (beta - 1) * log(1 - gma) +
                  lx[1]) +
              exp(log(1 - beta) + log(alpha) + (alpha - 1) * log(gma) + lx[2])
    d <- log(.expr1 + (1 - alpha) * (1 - beta) / .expr2) - v + jac
  }
  if (!log) d <- exp(d)
  d
}

db_tvnegbilog <- function(q1, q2, alpha, beta, mar1, mar2, thres, eta,
                          log = FALSE) {
  if (length(alpha) != 1 || mode(alpha) != "numeric")
    stop("invalid argument for `alpha'")
  if (length(beta) != 1 || mode(beta) != "numeric")
    stop("invalid argument for `beta'")
  if (any(c(alpha, beta) <= 0))
    stop("`alpha' and `beta' must be non-negative")
  eta <- normalize_eta(eta)
  x1 <- mtransform_gp_mk2(q1, mar1, thres[1], eta[1], margin = "exp")
  x2 <- mtransform_gp_mk2(q2, mar2, thres[2], eta[2], margin = "exp")
  ext <- c((x1 %in% c(0, Inf)), (x2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    gmafn <- function(z) (1 + alpha) * x1 * z^alpha -
                          (1 + beta) * x2 * (1 - z)^beta
    gma <- uniroot(gmafn, lower = 0, upper = 1,
                   tol = .Machine$double.eps^0.5)$root
    v <- x1 + x2 - x1 * gma^(1 + alpha) - x2 * (1 - gma)^(1 + beta)
    lx <- log(c(x1, x2))
    jac <- (1 + mar1[2]) * lx[1] + (1 + mar2[2]) * lx[2] -
           log(mar1[1] * mar2[1])
    .expr1 <- (1 - gma^(1 + alpha)) * (1 - (1 - gma)^(1 + beta))
    .expr2 <- exp(log(1 + alpha) + log(1 + beta) + alpha * log(gma) +
                  beta * log(1 - gma))
    .expr3 <- exp(log(1 + alpha) + log(alpha) + (alpha - 1) * log(gma) +
                  lx[1]) +
              exp(log(1 + beta) + log(beta) + (beta - 1) * log(1 - gma) + lx[2])
    d <- log(.expr1 + .expr2 / .expr3) - v + jac
  }
  if (!log) d <- exp(d)
  d
}

db_tv_nonpar <- function(q1, q2, mar1, mar2, thres, eta, Ahat, log = FALSE) {
  # nonparametric bivariate density via Pickands dependence function
  eta <- normalize_eta(eta)
  u1 <- mtransform_gp_mk2(q1, mar1, thres[1], eta[1], margin = "uniform")
  u2 <- mtransform_gp_mk2(q2, mar2, thres[2], eta[2], margin = "uniform")
  ext <- c((u1 %in% c(0, Inf)), (u2 %in% c(0, Inf)))
  d <- -Inf
  if (all(!ext)) {
    w <- log(u1) / log(u1 * u2)
    alpha <- a_bp_approx(Ahat, t = w, ord = "0")
    print(alpha)
  }
  d
}
