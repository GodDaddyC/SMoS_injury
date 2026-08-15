# bivariate POT (logistic) with non-stationary scale
#   extends evd::fbvpot — censored likelihood, but scale is
#   per-observation: scale_i = drop(design %*% coef)
#   (identity/additive link). Shape, dependence remain scalar.
#   The censored likelihood is an R port of evd's C routine
#   nllbvclog (src/bvpot.c), re-implemented so that scale can
#   vary per observation. The GEV-level machinery (param naming,
#   starting values, post-processing) follows evd::fbvlog with
#   nsloc replaced by nsscale.

# R port of evd C function nllbvclog (bvpot.c:5-87) with vectorized scale.
# scale1/scale2 may be vectors of length nn (per observation).
nllbvclog_r <- function(x1, x2, thdi, lambda, n, nn, dep,
                        scale1, shape1, scale2, shape2) {
  if (any(scale1 < 0.01) || any(scale2 < 0.01) || dep < 0.1 || dep > 1)
    return(1e6)

  lambda2 <- (-1 / log(1 - lambda))^(-1 / dep)
  zdn <- -(lambda2[1] + lambda2[2])^dep

  z1 <- x1 / scale1
  z2 <- x2 / scale2

  if (shape1 == 0) {
    t1 <- exp(-z1)
  } else {
    t1 <- 1 + shape1 * z1
    if (any(t1 <= 0)) return(1e6)
    t1 <- t1^(-1 / shape1)
  }
  if (shape2 == 0) {
    t2 <- exp(-z2)
  } else {
    t2 <- 1 + shape2 * z2
    if (any(t2 <= 0)) return(1e6)
    t2 <- t2^(-1 / shape2)
  }

  y1 <- -1 / log(1 - lambda[1] * t1)
  y2 <- -1 / log(1 - lambda[2] * t2)

  tt1 <- lambda[1] * y1^2 * t1^(1 + shape1) / ((1 - lambda[1] * t1) * scale1)
  tt2 <- lambda[2] * y2^2 * t2^(1 + shape2) / ((1 - lambda[2] * t2) * scale2)

  v1 <- y1^(-1 / dep)
  v2 <- y2^(-1 / dep)
  v12 <- (v1 + v2)^(dep - 1)
  v <- v12 * (v1 + v2)
  v1 <- -(v1 / y1) * v12
  v2 <- -(v2 / y2) * v12
  v12 <- (1 - 1 / dep) * v1 * v2 / v

  dvec <- numeric(nn)
  lo1 <- thdi < 1.5
  lo2 <- thdi >= 1.5 & thdi < 2.5
  lo3 <- thdi >= 2.5
  dvec[lo1] <- log(-v1[lo1]) + log(tt1[lo1]) - v[lo1]
  dvec[lo2] <- log(-v2[lo2]) + log(tt2[lo2]) - v[lo2]
  dvec[lo3] <- log(v1[lo3] * v2[lo3] - v12[lo3]) +
               log(tt1[lo3]) + log(tt2[lo3]) - v[lo3]

  -sum(dvec) - (n - nn) * zdn
}

# replicate evd:::sep.bvdata(x, method = "cpot", u = u)
sep_bvdata_cpot <- function(x, u) {
  x1 <- x[, 1]
  x2 <- x[, 2]
  n <- length(x1)
  iau1 <- x1 > u[1] & !is.na(x1)
  iau2 <- x2 > u[2] & !is.na(x2)
  nat <- c(sum(iau1), sum(iau2), sum(iau1 & iau2))
  lambda <- c(sum(iau1) / (n + 1), sum(iau2) / (n + 1))
  x1 <- x1 - u[1]
  x2 <- x2 - u[2]
  x1[!iau1] <- 0
  x2[!iau2] <- 0
  i0 <- iau1 | iau2
  x1 <- x1[i0]
  x2 <- x2[i0]
  nn <- length(x1)
  thdi <- as.logical(x1) + 2 * as.logical(x2)
  list(x1 = x1, x2 = x2, nn = nn, n = n, thdi = thdi,
       lambda = lambda, nat = nat, i0 = i0)
}

# build a non-stationary scale design matrix from covariate columns.
# Intercept column of 1s (coefficient prefix_0); numeric covariate -> one
# column = the value; factor covariate -> one 0/1 indicator per non-reference
# level (reference absorbed into the intercept). Returns the design matrix,
# sequential parameter names and a term-metadata data.frame.
build_ns_design <- function(data, cols, prefix) {
  n <- nrow(data)
  design <- matrix(1, nrow = n, ncol = 1)
  colnames(design) <- "(Intercept)"
  terms <- data.frame(param = paste0(prefix, "_0"),
                      covariate = "(Intercept)",
                      level = NA_character_, stringsAsFactors = FALSE)

  if (!is.null(cols) && length(cols) > 0) {
    for (col in cols) {
      if (!col %in% names(data))
        stop(sprintf("covariate '%s' not found in data", col))
      x <- data[[col]]
      if (is.factor(x)) {
        levs <- levels(x)
        if (length(levs) < 2)
          stop(sprintf("factor '%s' must have at least two levels", col))
        ref <- levs[1]
        for (l in levs[-1]) {
          design <- cbind(design, as.numeric(as.character(x) == l))
          terms <- rbind(terms, data.frame(
            param = paste0(prefix, "_", ncol(design) - 1),
            covariate = col, level = l, stringsAsFactors = FALSE))
        }
      } else if (is.numeric(x)) {
        design <- cbind(design, as.numeric(x))
        terms <- rbind(terms, data.frame(
          param = paste0(prefix, "_", ncol(design) - 1),
          covariate = col, level = NA_character_, stringsAsFactors = FALSE))
      } else {
        stop(sprintf("covariate '%s' must be numeric or a factor", col))
      }
    }
    colnames(design) <- c("(Intercept)", terms$covariate[-1])
    rownames(design) <- NULL
  }

  list(design = design,
       par_names = terms$param,
       terms = terms)
}

# modified evd:::bvstart.vals (method = "pot" branch) with nsscale names
bvstart.vals.ns <- function(x, start, nmdots, param, u, design1, design2,
                            model) {
  scale.param1 <- paste0("scale1_", 0:(ncol(design1) - 1))
  scale.param2 <- paste0("scale2_", 0:(ncol(design2) - 1))

  if (missing(start)) {
    start <- as.list(numeric(length(param)))
    names(start) <- param
    st1 <- fitted(fpot(x[, 1], threshold = u[1], std.err = FALSE))
    st2 <- fitted(fpot(x[, 2], threshold = u[2], std.err = FALSE))
    start[scale.param1] <- c(st1["scale"], rep(0, length(scale.param1) - 1))
    start["shape1"] <- st1["shape"]
    start[scale.param2] <- c(st2["scale"], rep(0, length(scale.param2) - 1))
    start["shape2"] <- st2["shape"]
    if (model == "log") start[["dep"]] <- 0.75
    start <- start[!(param %in% nmdots)]
  }
  if (!is.list(start)) stop("`start' must be a named list")
  start
}

# modified evd:::bvpost.optim (method = "pot" branch) carrying nsscale through
bvpost.optim.ns <- function(x, opt, nm, std.err, corr, sym, cmar, u, nat,
                            data, resp, design1, design2, terms1, terms2,
                            nsscale1 = NULL, nsscale2 = NULL,
                            likelihood = "censored", model, call = NULL) {
  if (opt$convergence != 0) {
    warning(paste("optimization for", model, "may not have succeeded"),
            call. = FALSE)
    if (opt$convergence == 1) opt$convergence <- "iteration limit reached"
  } else {
    opt$convergence <- "successful"
  }

  if (std.err) {
    tol <- .Machine$double.eps^0.5
    var.cov <- qr(opt$hessian, tol = tol)
    if (var.cov$rank != ncol(var.cov$qr))
      stop(paste("observed information matrix for", model,
                 "is singular; use std.err = FALSE"))
    var.cov <- solve(var.cov, tol = tol)
    dimnames(var.cov) <- list(nm, nm)
    std.err <- diag(var.cov)
    if (any(std.err <= 0))
      stop(paste("observed information matrix for", model,
                 "is singular; use std.err = FALSE"))
    std.err <- sqrt(std.err)
    names(std.err) <- nm
    if (corr) {
      .mat <- diag(1 / std.err, nrow = length(std.err))
      corr <- structure(.mat %*% var.cov %*% .mat,
                        dimnames = list(nm, nm))
      diag(corr) <- rep(1, length(std.err))
    } else {
      corr <- NULL
    }
  } else {
    std.err <- var.cov <- corr <- NULL
  }

  dep <- opt$par["dep"]
  dep.sum <- 2 * (1 - abvevd(dep = dep, model = model))

  structure(
    list(estimate    = opt$par,
         std.err     = std.err,
         param       = opt$par,
         deviance    = 2 * opt$value,
         dep.summary = dep.sum,
         corr        = corr,
         var.cov     = var.cov,
         convergence = opt$convergence,
         counts      = opt$counts,
         message     = opt$message,
         data        = data,
         resp        = resp,
         threshold   = u,
         nat         = nat,
         likelihood  = likelihood,
         design1     = design1,
         design2     = design2,
         scale1_terms = terms1,
         scale2_terms = terms2,
         nsscale1    = nsscale1,
         nsscale2    = nsscale2,
         n           = nrow(data),
         sym         = sym,
         model       = model,
         call        = call),
    class = "fbvpot_ns")
}

# fbvpot_ns -------------------------------------------------------------------

fbvpot_ns <- function(data, resp = NULL, threshold = NULL,
                      nsscale1 = NULL, nsscale2 = NULL, model = "log",
                      start, ..., sym = FALSE, std.err = TRUE, corr = FALSE,
                      method = "BFGS", warn.inf = TRUE) {
  model <- match.arg(model, "log")
  call <- match.call()
  if (!is.data.frame(data)) stop("`data' must be a data.frame")
  if (is.null(resp) || length(resp) != 2)
    stop("`resp' must be a character vector of length two naming the two margins")
  if (any(!resp %in% names(data)))
    stop("`resp' columns not found in data")
  if (is.null(threshold) || length(threshold) != 2)
    stop("`threshold' must be a numeric vector of length two")

  x <- as.matrix(data[, resp, drop = FALSE])

  d1 <- build_ns_design(data, nsscale1, "scale1")
  d2 <- build_ns_design(data, nsscale2, "scale2")
  design1 <- d1$design
  design2 <- d2$design
  scale.param1 <- d1$par_names
  scale.param2 <- d2$par_names

  param <- c(scale.param1, "shape1", scale.param2, "shape2", "dep")
  nmdots <- names(list(...))

  start <- bvstart.vals.ns(x, start = start, nmdots = nmdots, param = param,
                           u = threshold, design1 = design1, design2 = design2,
                           model = model)

  fixed.param <- list(...)[nmdots %in% param]
  if (any(!(param %in% c(names(start), names(fixed.param)))))
    stop("unspecified parameters")

  spx <- sep_bvdata_cpot(x, threshold)

  design1_exc <- design1[spx$i0, , drop = FALSE]
  design2_exc <- design2[spx$i0, , drop = FALSE]

  nll <- function(par) {
    scale1 <- drop(design1_exc %*% par[scale.param1])
    scale2 <- drop(design2_exc %*% par[scale.param2])
    nllbvclog_r(spx$x1, spx$x2, spx$thdi, spx$lambda, spx$n, spx$nn,
                par[["dep"]], scale1, par[["shape1"]],
                scale2, par[["shape2"]])
  }

  nm <- names(start)
  nll.opt <- function(p, ...) nll(c(p, unlist(fixed.param)))

  start.arg <- list(p = unlist(start))
  if (warn.inf && do.call(nll.opt, start.arg) == 1e6)
    warning("negative log-likelihood is infinite at starting values")

  dots <- list(...)[!(nmdots %in% param)]
  opt <- do.call(optim, c(list(par = unlist(start), fn = nll.opt,
                               hessian = std.err, method = method), dots))
  if (is.null(names(opt$par))) names(opt$par) <- nm

  bvpost.optim.ns(x, opt, nm, std.err, corr, sym, cmar = c(FALSE, FALSE),
                  u = threshold, nat = spx$nat,
                  data = data, resp = resp,
                  design1 = design1, design2 = design2,
                  terms1 = d1$terms, terms2 = d2$terms,
                  nsscale1 = nsscale1, nsscale2 = nsscale2,
                  likelihood = "censored", model = model,
                  call = call)
}

# fitted.fbvpot_ns -------------------------------------------------------------

fitted.fbvpot_ns <- function(object, ...) {
  margin <- list(
    `1` = list(terms = object$scale1_terms, design = object$design1),
    `2` = list(terms = object$scale2_terms, design = object$design2))
  lapply(margin, function(m) {
    coef <- object$estimate[m$terms$param]
    drop(m$design %*% coef)
  })
}

# print.fbvpot_ns ---------------------------------------------------------------

print.fbvpot_ns <- function(x, ...) {
  cat("Call:\n"); print(x$call)
  cat("\nBivariate POT (logistic) with non-stationary scale\n")
  cat(sprintf("  Observations: %d   Exceedances: %d (margin1), %d (margin2)\n",
              x$n, x$nat[1], x$nat[2]))
  cat(sprintf("  Thresholds:  %.4f  %.4f\n", x$threshold[1], x$threshold[2]))
  cat(sprintf("  Convergence: %s\n", x$convergence))
  cat(sprintf("  Deviance:    %.4f\n", x$deviance))
  cat("\nMarginal scale models\n")
  print_model <- function(terms, est) {
    nm <- terms$param
    tbl <- data.frame(Param = nm, Est = round(est[nm], 4))
    print(tbl, row.names = FALSE)
  }
  if (!is.null(x$nsscale1)) {
    cat("\n  Margin 1 (", x$resp[1], "):\n", sep = "")
    print_model(x$scale1_terms, x$estimate)
  } else {
    cat(sprintf("\n  Margin 1 (%s): scale = %.4f\n", x$resp[1], x$estimate["scale1_0"]))
  }
  if (!is.null(x$nsscale2)) {
    cat("\n  Margin 2 (", x$resp[2], "):\n", sep = "")
    print_model(x$scale2_terms, x$estimate)
  } else {
    cat(sprintf("  Margin 2 (%s): scale = %.4f\n", x$resp[2], x$estimate["scale2_0"]))
  }
  cat(sprintf("\n  Shapes: %.4f, %.4f   Dependence (dep): %.4f\n",
              x$estimate["shape1"], x$estimate["shape2"], x$estimate["dep"]))
  invisible(x)
}

# confint.fbvpot_ns ------------------------------------------------------------

confint.fbvpot_ns <- function(object, parm, level = 0.95, ...) {
  cf <- object$estimate
  pnames <- names(cf)
  if (missing(parm)) {
    parm <- seq_along(pnames)
  } else if (is.character(parm)) {
    parm <- match(parm, pnames, nomatch = 0)
  }
  if (any(!parm)) stop("`parm' contains unknown parameters")

  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  pct <- paste(round(100 * a, 1), "%")
  ci <- array(NA, dim = c(length(parm), 2),
              dimnames = list(pnames[parm], pct))

  if (is.null(object$std.err))
    stop("standard errors not available; refit with std.err = TRUE")
  ses <- object$std.err[parm]
  ci[] <- cf[parm] + ses %o% qnorm(a)
  ci
}

# plot.fbvpot_ns ---------------------------------------------------------------

# per-observation scale for margin m
fbvpot_ns_scale <- function(x, m) {
  prefix <- paste0("scale", m)
  terms  <- x[[paste0(prefix, "_terms")]]
  design <- x[[paste0("design", m)]]
  coef   <- x$estimate[terms$param]
  drop(design %*% coef)
}

# standardized excesses z = (x - u) / sigma_i for exceedances of margin m
fbvpot_ns_excess <- function(x, m) {
  val <- x$data[[x$resp[m]]]
  u <- x$threshold[m]
  scale <- fbvpot_ns_scale(x, m)
  exc <- val > u
  (val[exc] - u) / scale[exc]
}

# transform margin m exceedances to a standard margin using per-observation scale
fbvpot_ns_transform <- function(x, m, margin = "exp") {
  val <- x$data[[x$resp[m]]]
  u <- x$threshold[m]
  scale <- fbvpot_ns_scale(x, m)
  shape <- x$estimate[[paste0("shape", m)]]
  eta <- x$nat[m] / x$n
  exc <- val > u
  p <- cbind(scale = scale[exc], shape = rep(shape, sum(exc)))
  out <- mtransform_gp_mk2(val[exc], p = p, thres = u, eta = eta, margin = margin)
  full <- rep(NA_real_, length(val))
  full[exc] <- out
  list(exceed = exc, value = full)
}

plot.fbvpot_ns <- function(x, num = NULL, ...) {
  if (!inherits(x, "fbvpot_ns"))
    stop("`x' must be an fbvpot_ns object")
  if (!is.null(num) && any(!num %in% 1:4))
    stop("`num' must be a subset of 1:4")

  plots <- list(
    `1` = function() fbvpot_ns_plot_qq(x, ...),
    `2` = function() fbvpot_ns_plot_density(x, ...),
    `3` = function() fbvpot_ns_plot_quantile(x, ...),
    `4` = function() fbvpot_ns_plot_dependence(x, ...))

  if (is.null(num)) {
    for (i in 1:4) {
      print(plots[[as.character(i)]]())
      if (i < 4) readline("Press <Enter> to continue")
    }
    return(invisible(x))
  }

  p <- plots[[as.character(num)]]()
  print(p)
  invisible(p)
}

# num == 1: QQ plot of the margins ------------------------------------------

fbvpot_ns_plot_qq <- function(x, ...) {
  qq_df <- do.call(rbind, lapply(1:2, function(m) {
    z <- fbvpot_ns_excess(x, m)
    k <- length(z)
    pp <- (seq_len(k) - 0.5) / k
    shape <- x$estimate[[paste0("shape", m)]]
    theo <- if (shape == 0) -log(1 - pp) else ((1 - pp)^(-shape) - 1) / shape
    data.frame(margin = m, emp = sort(z), theo = theo)
  }))

  ggplot(qq_df, aes(x = theo, y = emp)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    facet_wrap(~ margin, labeller = label_both) +
    labs(x = "Theoretical GP quantile", y = "Empirical quantile",
         title = "QQ plot of marginal exceedances") +
    theme_minimal()
}

# num == 2: estimated density vs histogram of the margin ---------------------

fbvpot_ns_plot_density <- function(x, nbin = 30, ...) {
  hist_df <- do.call(rbind, lapply(1:2, function(m) {
    data.frame(margin = m, z = fbvpot_ns_excess(x, m))
  }))
  curve_df <- do.call(rbind, lapply(1:2, function(m) {
    shape <- x$estimate[[paste0("shape", m)]]
    z <- fbvpot_ns_excess(x, m)
    ub <- if (shape < 0) min(max(z), -1 / shape) else max(z)
    grid <- seq(0, ub, length.out = 200)
    dens <- if (shape == 0) exp(-grid) else (1 + shape * grid)^(-1 / shape - 1)
    data.frame(margin = m, z = grid, dens = dens)
  }))

  ggplot(hist_df, aes(x = z)) +
    geom_histogram(aes(y = after_stat(density)), bins = nbin,
                   fill = "grey70", colour = "white", alpha = 0.8) +
    geom_line(data = curve_df, aes(x = z, y = dens),
              colour = "red", linewidth = 0.8) +
    facet_wrap(~ margin, labeller = label_both) +
    labs(x = "Standardized excess  z = (x - u)/sigma", y = "Density",
         title = "Fitted GP density vs sample histogram") +
    theme_minimal()
}

# num == 3: quantile plot of the bivariate dependence model ------------------

fbvpot_ns_plot_quantile <- function(x, p = seq(0.75, 0.95, 0.05), ...) {
  dep <- x$estimate[["dep"]]
  tr1 <- fbvpot_ns_transform(x, 1, "exp")
  tr2 <- fbvpot_ns_transform(x, 2, "exp")
  j <- tr1$exceed & tr2$exceed
  pts <- data.frame(z1 = tr1$value[j], z2 = tr2$value[j])
  om <- seq(0.001, 0.999, length.out = 200)
  aom <- (om^(1 / dep) + (1 - om)^(1 / dep))^dep
  curve_df <- do.call(rbind, lapply(p, function(pv) {
    data.frame(z1 = -om / aom * log(pv),
               z2 = -(1 - om) / aom * log(pv),
               p = pv)
  }))

  ggplot() +
    geom_point(data = pts, aes(x = z1, y = z2), alpha = 0.5, size = 1) +
    geom_line(data = curve_df, aes(x = z1, y = z2, colour = factor(p)),
              linewidth = 0.7) +
    scale_colour_brewer(palette = "Set1", name = "p") +
    labs(x = "Margin 1 (exponential scale)", y = "Margin 2 (exponential scale)",
         title = "Quantile curves of the fitted logistic dependence model") +
    theme_minimal()
}

# num == 4: dependence function and CFG estimator ----------------------------

fbvpot_ns_plot_dependence <- function(x, level = 0.95, ...) {
  dep <- x$estimate[["dep"]]
  tr1 <- fbvpot_ns_transform(x, 1, "exp")
  tr2 <- fbvpot_ns_transform(x, 2, "exp")
  dat_exp <- cbind(tr1$value, tr2$value)
  if (sum(is.finite(dat_exp[, 1]) & is.finite(dat_exp[, 2])) < 5)
    stop("not enough joint exceedances for the CFG estimator")

  t <- seq(0, 1, length.out = 200)
  A_cfg <- abvnonpar(t, data = dat_exp, epmar = TRUE, method = "cfg")
  A_par <- abvevd(t, dep = dep, model = "log", plot = FALSE)

  dep_ci <- confint(x, parm = "dep", level = level)[, 1:2]
  dep_lo <- max(min(dep_ci), 0.01)
  dep_hi <- min(max(dep_ci), 1)
  A_lo <- abvevd(t, dep = dep_lo, model = "log", plot = FALSE)
  A_hi <- abvevd(t, dep = dep_hi, model = "log", plot = FALSE)

  dep_df <- data.frame(t, A_cfg, A_par, A_lo, A_hi)
  ggplot(dep_df, aes(x = t)) +
    geom_ribbon(aes(ymin = A_lo, ymax = A_hi), alpha = 0.2) +
    geom_line(aes(y = A_cfg), colour = "black") +
    geom_line(aes(y = A_par), colour = "red") +
    ylim(c(0.5,1)) +
    labs(x = "t", y = "A(t)",
         title = "Dependence function: CFG estimator and fitted model") +
    annotate("text", x = 0.98, y = min(dep_df$A_lo), hjust = 1, vjust = -0.5,
             label = "red: fitted logistic  black: CFG  band: parametric CI",
             size = 3) +
    theme_minimal()
}
