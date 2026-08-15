# truncated bivariate POT with logistic model and non-stationary scale
#   extends evd::fbvpot — censored likelihood + upper truncation
#   only the scale parameters are non-stationary (linear trend)
#   scale on the transformed scale: log-scale = scale_0 + scale_1 * covariate

gp_trunc_to_exp <- function(x, u, sigma, xi, eta) {
  if (any(x < u)) stop("input below threshold")
  sc <- (x - u) / sigma
  zind <- abs(xi) < 1e-7
  Fx <- numeric(length(x))
  if (any(zind)) {
    Fx[zind] <- 1 - eta * exp(-sc[zind])
  }
  if (any(!zind)) {
    t <- 1 + xi[!zind] * sc[!zind]
    ok <- t > 0
    if (!all(ok)) return(rep(NaN, length(x)))
    Fx[!zind] <- pmax(1 - eta * t^(-1 / xi[!zind]), 0)
  }
  -log(Fx)
}

gp_log_jacobian <- function(x_orig, y_exp, u, sigma, xi, eta) {
  sc <- (x_orig - u) / sigma
  zind <- abs(xi) < 1e-7
  lj <- numeric(length(x_orig))
  if (any(zind)) {
    lj[zind] <- log(eta) - log(sigma[zind]) - sc[zind] + y_exp[zind]
  }
  if (any(!zind)) {
    t <- 1 + xi[!zind] * sc[!zind]
    if (any(t <= 0)) return(rep(NaN, length(x_orig)))
    lj[!zind] <- log(eta) - log(sigma[!zind]) +
      (-1 / xi[!zind] - 1) * log(t) + y_exp[!zind]
  }
  lj
}

log_log_cens_contrib <- function(par, data, u1, u2, eta1, eta2,
                                 U1, U2, covar) {
  scale1_0 <- par[1]; scale1_1 <- par[2]; xi1 <- par[3]
  scale2_0 <- par[4]; scale2_1 <- par[5]; xi2 <- par[6]
  alpha  <- par[7]

  if (alpha <= 0 || alpha > 1) return(1e10)
  if (length(par) != 7) return(1e10)

  n <- nrow(data)
  x <- data[, 1]
  y <- data[, 2]

  sigma1 <- exp(scale1_0 + scale1_1 * covar)
  sigma2 <- exp(scale2_0 + scale2_1 * covar)

  xa <- x > u1
  ya <- y > u2

  case <- rep(0L, n)
  case[xa & ya]  <- 1L
  case[xa & !ya] <- 2L
  case[!xa & ya] <- 3L
  case[!xa & !ya] <- 4L

  ll <- numeric(n)

  # --- case 1: both above thresholds (and below upper bounds) ---
  idx <- which(case == 1L)
  if (length(idx) > 0) {
    y1 <- gp_trunc_to_exp(x[idx], u1, sigma1[idx], xi1, eta1)
    if (any(is.nan(y1))) return(1e10)
    y2 <- gp_trunc_to_exp(y[idx], u2, sigma2[idx], xi2, eta2)
    if (any(is.nan(y2))) return(1e10)
    V  <- (y1^(1/alpha) + y2^(1/alpha))^alpha
    if (any(is.nan(V))) return(1e10)

    ldens_exp <- -V + (1/alpha - 1) * (log(y1) + log(y2)) +
                 (1 - 2/alpha) * log(V) + log(V + 1/alpha - 1)

    lj1 <- gp_log_jacobian(x[idx], y1, u1, sigma1[idx], xi1, eta1)
    if (any(is.nan(lj1))) return(1e10)
    lj2 <- gp_log_jacobian(y[idx], y2, u2, sigma2[idx], xi2, eta2)
    if (any(is.nan(lj2))) return(1e10)

    ll[idx] <- ldens_exp + lj1 + lj2
  }

  # --- case 2: x above threshold, y below threshold ---
  idx <- which(case == 2L)
  if (length(idx) > 0) {
    y1 <- gp_trunc_to_exp(x[idx], u1, sigma1[idx], xi1, eta1)
    if (any(is.nan(y1))) return(1e10)
    y2v <- gp_trunc_to_exp(rep(u2, length(idx)), u2, sigma2[idx], xi2, eta2)
    if (any(is.nan(y2v))) return(1e10)
    V  <- (y1^(1/alpha) + y2v^(1/alpha))^alpha

    log_dG_dy1 <- -V + (1 - 1/alpha) * log(V) + (1/alpha - 1) * log(y1)
    lj1 <- gp_log_jacobian(x[idx], y1, u1, sigma1[idx], xi1, eta1)
    if (any(is.nan(lj1))) return(1e10)
    ll[idx] <- log_dG_dy1 + lj1
  }

  # --- case 3: x below threshold, y above threshold ---
  idx <- which(case == 3L)
  if (length(idx) > 0) {
    y1v <- gp_trunc_to_exp(rep(u1, length(idx)), u1, sigma1[idx], xi1, eta1)
    if (any(is.nan(y1v))) return(1e10)
    y2 <- gp_trunc_to_exp(y[idx], u2, sigma2[idx], xi2, eta2)
    if (any(is.nan(y2))) return(1e10)
    V  <- (y1v^(1/alpha) + y2^(1/alpha))^alpha

    log_dG_dy2 <- -V + (1 - 1/alpha) * log(V) + (1/alpha - 1) * log(y2)
    lj2 <- gp_log_jacobian(y[idx], y2, u2, sigma2[idx], xi2, eta2)
    if (any(is.nan(lj2))) return(1e10)
    ll[idx] <- log_dG_dy2 + lj2
  }

  # --- case 4: both below thresholds ---
  idx <- which(case == 4L)
  if (length(idx) > 0) {
    y1v <- gp_trunc_to_exp(rep(u1, length(idx)), u1, sigma1[idx], xi1, eta1)
    if (any(is.nan(y1v))) return(1e10)
    y2v <- gp_trunc_to_exp(rep(u2, length(idx)), u2, sigma2[idx], xi2, eta2)
    if (any(is.nan(y2v))) return(1e10)
    V  <- (y1v^(1/alpha) + y2v^(1/alpha))^alpha
    ll[idx] <- -V
  }

  # --- truncation normalization ---
  if (is.finite(U1) || is.finite(U2)) {
    U1e <- if (is.finite(U1)) U1 else max(x) * 10
    U2e <- if (is.finite(U2)) U2 else max(y) * 10
    y1U <- gp_trunc_to_exp(rep(U1e, n), u1, sigma1, xi1, eta1)
    if (any(is.nan(y1U))) return(1e10)
    y2U <- gp_trunc_to_exp(rep(U2e, n), u2, sigma2, xi2, eta2)
    if (any(is.nan(y2U))) return(1e10)
    VU  <- (y1U^(1/alpha) + y2U^(1/alpha))^alpha
    ll <- ll + VU
  }

  nll <- -sum(ll)
  if (is.na(nll) || is.nan(nll) || is.infinite(nll)) return(1e10)
  nll
}

fit_trunc_log <- function(data, threshold, upper = c(Inf, Inf),
                          covariate = NULL, start = NULL,
                          method = "Nelder-Mead", hessian = TRUE, ...) {
  if (!is.matrix(data)) data <- as.matrix(data)
  if (ncol(data) != 2) stop("data must be a 2-column matrix")
  n <- nrow(data)
  u1 <- threshold[1]
  u2 <- threshold[2]
  U1 <- upper[1]
  U2 <- upper[2]

  if (is.null(covariate)) covariate <- rep(0, n)
  if (length(covariate) != n)
    stop("covariate length must match number of observations")
  covariate <- as.numeric(covariate)

  x <- data[, 1]; y <- data[, 2]
  eta1 <- mean(x > u1)
  eta2 <- mean(y > u2)

  if (is.null(start)) {
    fit1 <- tryCatch(
      evd::fpot(x, threshold = u1, model = "gpd", std.err = FALSE),
      error = function(e) NULL)
    s1 <- if (!is.null(fit1)) fit1$estimate["scale"] else sd(x[x > u1])
    xi1_s <- if (!is.null(fit1)) fit1$estimate["shape"] else -0.1

    fit2 <- tryCatch(
      evd::fpot(y, threshold = u2, model = "gpd", std.err = FALSE),
      error = function(e) NULL)
    s2 <- if (!is.null(fit2)) fit2$estimate["scale"] else sd(y[y > u2])
    xi2_s <- if (!is.null(fit2)) fit2$estimate["shape"] else -0.1

    start <- c(log(s1), 0, xi1_s, log(s2), 0, xi2_s, 0.5)
  }

  names(start) <- c("scale1_0", "scale1_1", "shape1",
                    "scale2_0", "scale2_1", "shape2", "dep")

  opt <- optim(par = start, fn = log_log_cens_contrib,
               data = data, u1 = u1, u2 = u2,
               eta1 = eta1, eta2 = eta2,
               U1 = U1, U2 = U2, covar = covariate,
               method = method, hessian = hessian, ...)

  est <- opt$par
  names(est) <- names(start)
  est_exp <- c(exp(est[1]), est[2], est[3],
               exp(est[4]), est[5], est[6], est[7])
  names(est_exp) <- c("scale1", "scale1_trend", "shape1",
                      "scale2", "scale2_trend", "shape2", "dep")

  vcov <- NULL
  if (hessian) {
    vcov <- tryCatch(solve(opt$hessian), error = function(e) NULL)
  }

  structure(list(
    estimate    = est_exp,
    estimate_raw = est,
    threshold   = threshold,
    upper       = upper,
    model       = "log",
    dep         = est[7],
    n           = n,
    nat         = c(sum(x > u1 & y > u2),
                    sum(x > u1), sum(y > u2)),
    eta         = c(eta1, eta2),
    covariate   = covariate,
    negloglik   = opt$value,
    convergence = opt$convergence,
    counts      = opt$counts,
    vcov        = vcov,
    hessian     = if (hessian) opt$hessian else NULL,
    data        = data,
    call        = match.call()
  ), class = "trunc_bpot_log")
}

# --- methods ---

print.trunc_bpot_log <- function(x, ...) {
  cat("Truncated bivariate POT — logistic model\n")
  cat(sprintf("  Observations: %d\n", x$n))
  cat(sprintf("  Thresholds:  %.4f  %.4f\n", x$threshold[1], x$threshold[2]))
  cat(sprintf("  Upper bounds: %.4f  %.4f\n",
              ifelse(is.finite(x$upper[1]), x$upper[1], Inf),
              ifelse(is.finite(x$upper[2]), x$upper[2], Inf)))
  cat(sprintf("  Convergence:  %d\n", x$convergence))
  cat("\nParameter estimates (on exp-scale for intercept):\n")
  print(x$estimate)
  cat(sprintf("\nNegative log-likelihood: %.4f\n", x$negloglik))
  invisible(x)
}

logLik.trunc_bpot_log <- function(object, ...) {
  val <- -object$negloglik
  attr(val, "df") <- length(object$estimate)
  class(val) <- "logLik"
  val
}

coef.trunc_bpot_log <- function(object, ...) {
  object$estimate
}

vcov.trunc_bpot_log <- function(object, ...) {
  if (is.null(object$vcov)) stop("hessian not available; re-fit with hessian = TRUE")
  object$vcov
}
