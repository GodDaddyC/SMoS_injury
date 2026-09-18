# wrappers for packages `copula`, `VC2copula`

q_distr_param <- function(Cop, mar1, mar2, type) {
  if (type == 1) {
    Mdistr <- mvdc(Cop, margins = c("evd", mar2$distname),
      paramMargins = list(
        list(scale = mar1$results$par[1], shape = mar1$results$par[2],
             threshold = mar1$threshold, type = "GP"),
        as.list(mar2$estimate)))
  } else if (type == 2) {
    Mdistr <- mvdc(Cop, margins = c("evd", mar2$distname),
      paramMargins = list(
        list(scale = mar1$results$par[1], shape = mar1$results$par[2],
             threshold = mar1$threshold, type = "GP", lower.tail = FALSE),
        c(as.list(mar2$estimate), lower.tail = FALSE)))
  } else if (type == 3) {
    Mdistr <- mvdc(Cop, margins = c("evd", mar2$distname),
      paramMargins = list(
        list(scale = mar1$results$par[1], shape = mar1$results$par[2],
             threshold = mar1$threshold, type = "GP"),
        c(as.list(mar2$estimate), lower.tail = FALSE)))
  } else if (type == 4) {
    Mdistr <- mvdc(Cop, margins = c("evd", mar2$distname),
      paramMargins = list(
        list(scale = mar1$results$par[1], shape = mar1$results$par[2],
             threshold = mar1$threshold, type = "GP", lower.tail = FALSE),
        as.list(mar2$estimate)))
  }
  return(Mdistr)
}

gof_ev_copula_vc2 <- function(copula, x, N = 1000,
    method = c("mpl", "ml", "itau", "irho"),
    estimator = c("CFG", "Pickands"), m = 1000,
    verbose = interactive(),
    ties.method = c("max", "average", "first", "last", "random", "min"),
    fit.ties.meth = eval(formals(rank)$ties.method), ...) {
  stopifnot(is(copula, "copula"), N >= 1L, m >= 100L)
  if (!is.matrix(x)) {
    warning("coercing 'x' to a matrix.")
    stopifnot(is.matrix(x <- as.matrix(x)))
  }
  stopifnot((p <- ncol(x)) > 1, (n <- nrow(x)) > 1, dim(copula) == p)
  method <- match.arg(method)
  estimator <- match.arg(estimator)
  ties.method <- match.arg(ties.method)
  fit.ties.meth <- match.arg(fit.ties.meth)
  if (p != 2)
    stop("The copula and the data should be of dimension two")
  u <- pobs(x, ties.method = ties.method)
  u.fit <- if (ties.method == fit.ties.meth) u
           else pobs(x, ties.method = fit.ties.meth)
  fcop <- fitCopula(copula, u.fit, method, ...)@copula
  g <- seq(0, 1 - 1 / m, by = 1 / m)
  s <- .C("cramer_vonMises_Afun", as.integer(n), as.integer(m),
          as.double(-log(u[, 1])), as.double(-log(u[, 2])),
          as.double(A(fcop, g)), stat = double(2),
          as.integer(estimator == "CFG"))$stat
  s0 <- matrix(NA, N, 2)
  if (verbose) {
    pb <- txtProgressBar(max = N, style = if (isatty(stdout())) 3 else 1)
    on.exit(close(pb))
  }
  for (i in 1:N) {
    u0 <- pobs(rCopula(n, fcop), ties.method = ties.method)
    u0.fit <- if (ties.method == fit.ties.meth) u0
              else pobs(u0, ties.method = fit.ties.meth)
    fcop0 <- fitCopula(copula, u0.fit, method, ...)@copula
    s0[i, ] <- .C("cramer_vonMises_Afun", as.integer(n), as.integer(m),
                  as.double(-log(u0[, 1])), as.double(-log(u0[, 2])),
                  as.double(A(fcop0, g)), stat = double(2),
                  as.integer(estimator == "CFG"))$stat
    if (verbose) setTxtProgressBar(pb, i)
  }
  structure(class = "htest", list(
    method = paste0("Parametric bootstrap based GOF test for EV copulas ",
                    "with argument 'method' set to ", sQuote(method),
                    " and argument 'estimator' set to ", sQuote(estimator)),
    parameter = c(parameter = fcop@parameters),
    statistic = c(statistic = s[1]),
    p.value = (sum(s0[, 1] >= s[1]) + 0.5) / (N + 1),
    data.name = deparse(substitute(x))))
}


gof_copula_vc2 <- function(copula, x, N = 1000,
    method = c("mpl", "ml", "itau", "irho"),
    estimator = c("CFG", "Pickands"), m = 1000,
    verbose = interactive(),
    ties.method = c("max", "average", "first", "last", "random", "min"),
    fit.ties.meth = eval(formals(rank)$ties.method), ...) {
  stopifnot(is(copula, "copula"), N >= 1L, m >= 100L)
  if (!is.matrix(x)) {
    warning("coercing 'x' to a matrix.")
    stopifnot(is.matrix(x <- as.matrix(x)))
  }
  stopifnot((p <- ncol(x)) > 1, (n <- nrow(x)) > 1, dim(copula) == p)
  method <- match.arg(method)
  estimator <- match.arg(estimator)
  ties.method <- match.arg(ties.method)
  fit.ties.meth <- match.arg(fit.ties.meth)
  if (p != 2)
    stop("The copula and the data should be of dimension two")
  u <- pobs(x, ties.method = ties.method)
  u.fit <- if (ties.method == fit.ties.meth) u
           else pobs(x, ties.method = fit.ties.meth)
  fcop <- fitCopula(copula, u.fit, method, ...)@copula
  g <- seq(0, 1 - 1 / m, by = 1 / m)
  s <- .C("cramer_vonMises_Afun", as.integer(n), as.integer(m),
          as.double(-log(u[, 1])), as.double(-log(u[, 2])),
          as.double(A(fcop, g)), stat = double(2),
          as.integer(estimator == "CFG"))$stat
  s0 <- matrix(NA, N, 2)
  if (verbose) {
    pb <- txtProgressBar(max = N, style = if (isatty(stdout())) 3 else 1)
    on.exit(close(pb))
  }
  for (i in 1:N) {
    u0 <- pobs(rCopula(n, fcop), ties.method = ties.method)
    u0.fit <- if (ties.method == fit.ties.meth) u0
              else pobs(u0, ties.method = fit.ties.meth)
    fcop0 <- fitCopula(copula, u0.fit, method, ...)@copula
    s0[i, ] <- .C("cramer_vonMises_Afun", as.integer(n), as.integer(m),
                  as.double(-log(u0[, 1])), as.double(-log(u0[, 2])),
                  as.double(A(fcop0, g)), stat = double(2),
                  as.integer(estimator == "CFG"))$stat
    if (verbose) setTxtProgressBar(pb, i)
  }
  structure(class = "htest", list(
    method = paste0("Parametric bootstrap based GOF test for EV copulas ",
                    "with argument 'method' set to ", sQuote(method),
                    " and argument 'estimator' set to ", sQuote(estimator)),
    parameter = c(parameter = fcop@parameters),
    statistic = c(statistic = s[1]),
    p.value = (sum(s0[, 1] >= s[1]) + 0.5) / (N + 1),
    data.name = deparse(substitute(x))))
}


# Nonparametric approach -------------------------------------------------------

c_q_bivariate_nonpar <- function(v, u, model) {
  integrand <- function(q1, ...) {
    result <- numeric(length(q1))
    for (i in seq_along(q1)) {
      result[i] <- dkdecop(c(q1[i], v), model)
    }
    return(result)
  }

  R <- tryCatch({
    integrate(integrand, lower = u, upper = 1)$value
  }, error = function(e) {
    if (grepl("non-finite function value", e$message)) {
      ub_alt <- 1 - 1e-7
      return(integrate(integrand, lower = u, upper = ub_alt)$value)
    } else { stop(e) }
  })
  return(R)
}

normalize_c_q_bivariate_nonpar <- function(x, pcrash, model, lb = 0) {
  c_y <- function(k) {
    temp <- c_q_bivariate_nonpar(v = k, u = x, px = pcrash, model = model)
    return(temp)
  }

  C <- tryCatch({
    integrate(Vectorize(c_y), lower = lb, upper = 1)$value
  }, error = function(e) {
    if (grepl("non-finite function value", e$message)) {
      ub_alt <- 0.999
      return(integrate(Vectorize(c_y), lower = lb, upper = ub_alt)$value)
    } else { stop(e) }
  })
  return(C)
}

create_plot_df_q_nonpar <- function(dat, x, model, px, p2) {
  df <- data.frame(
    speed_u = dat,
    speed   = qgamma(dat, shape = p2$estimate[1], rate = p2$estimate[2]),
    JointP  = sapply(dat, function(y) { pkdecop(c(x, y), model) })
  ) %>%
    mutate(ConditionP = JointP / px) %>%
    mutate(ConditionalD = sapply(speed_u, c_q_bivariate_nonpar, u = x,
                                 model = model) / px *
             dgamma(speed, shape = p2$estimate[1],
                    rate = p2$estimate[2])) %>% na.omit()
  return(df)
}


injury_from_c_q_bivariate_nonpar <- function(model, severity, x0, px, p2) {
  f_y <- function(k) {
    temp <- c_q_bivariate_nonpar(v = k, u = x0, model = model) *
            severity(qgamma(k, shape = p2$estimate[1],
                            rate = p2$estimate[2]))
    return(temp)
  }

  R <- tryCatch({
    integrate(Vectorize(f_y), lower = 0, upper = 1)$value
  }, error = function(e) {
    if (grepl("non-finite function value", e$message)) {
      ub_alt <- 1 - 1e-5
      return(integrate(Vectorize(f_y), lower = 0, upper = ub_alt)$value)
    } else { stop(e) }
  })
  return(R / px)
}
