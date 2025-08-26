
Q.distr.param <- function(Cop,mar1,mar2, type){
  # creates mvdc model according to the fitted margin. mar1 should be fevd object and mar2 should 
  # be a fitdist object. Cop is a copula object Type refers to which region of events.
  if (type ==1){
    Mdistr<-mvdc(Cop,margins = c("evd",mar2$distname),
         paramMargins = list(list(scale=mar1$results$par[1],
                                  shape=mar1$results$par[2],threshold=mar1$threshold,type="GP"),
                             as.list(mar2$estimate)
         )
    )
  }
  else if (type ==2){
    Mdistr<-mvdc(Cop,margins = c("evd","gamma"),
                 paramMargins = list(list(scale=mar1$results$par[1],
                                          shape=mar1$results$par[2],threshold=mar1$threshold,type="GP",lower.tail=FALSE),
                                     c(as.list(mar2$estimate),lower.tail=FALSE) 
                 )
    )
  }
  else if (type ==3){
    Mdistr<-mvdc(Cop,margins = c("evd","gamma"),
                 paramMargins = list(list(scale=mar1$results$par[1],
                                          shape=mar1$results$par[2],threshold=mar1$threshold,type="GP"),
                                     c(as.list(mar2$estimate),lower.tail=FALSE) 
                 )
    )
  }
  else if (type ==4){
    Mdistr<-mvdc(Cop,margins = c("evd","gamma"),
                 paramMargins = list(list(scale=mar1$results$par[1],
                                          shape=mar1$results$par[2],threshold=mar1$threshold,type="GP",lower.tail=FALSE),
                                     as.list(mar2$estimate)
                 )
    )
  }
  return(Mdistr)
}


cQ.bivariate <-function(y,x,model){
  # approximates the conditional density f(y|X >x) by integrating fxy = f(x=x, y = y) over x
  # model is an mvdc object
  # if integral is non-finite, change ulim to a smaller value
  
  model@paramMargins <- lapply(model@paramMargins, function(param_list) {
    if ("lower.tail" %in% names(param_list)) {
      param_list[setdiff(names(param_list), "lower.tail")]
    } else {
      param_list
    }
  })
  integrand <- function(q1,...) {
    result <- numeric(length(q1))
    for(i in seq_along(q1)) {
      result[i] <- dMvdc(c(q1[i],y),model)*(1-q1[i]^2)
    }
    return(result)
  }
  
  # Integrate from x to Inf
  R <- tryCatch({
    integrate(integrand, lower = x, upper = Inf)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- model@paramMargins[[1]]$threshold - model@paramMargins[[1]]$scale/model@paramMargins[[1]]$shape
        return(integrate(integrand, lower = x, upper = ub.alt)$value)
      }
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(R)
}

normalize_cQ.bivariate<- function(x, model,lb=0){
  # computes the infinite integral of the conditional density f(y|X >x)
  # use as a nomralization factor for the conditional density
  # lb is the lower bound of the integral, by defalut lb = 0
  c.y <- function(k) {
    temp <- cQ.bivariate(y = k, x = x, model = model)
    return(temp)
  }
  
  C <- tryCatch({
    integrate(Vectorize(c.y), lower = lb, upper = Inf)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- 80
        return(integrate(Vectorize(c.y), lower = lb, upper = ub.alt)$value)
      }
      if (grepl("maximum number of subdivisions reached", e$message)) {
        message("Maximum number of subdivisions reached. Retrying with larger tolerence.")
        return(integrate(Vectorize(c.y), lower = lb,upper=Inf,
                         subdivisions = 200,rel.tol = 1e-5)$value)
      }
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(C) 
}

Injury.from_cQ_bivariate <-function(dat,model,severity,x0){
  # computes the injury probability using c.bivariate, severity is a logistic regression model of injury given speed
  f.y <- function(k) {
    temp <- cQ.bivariate(y = k, x = x0, model = model) * severity(k)
    return(temp)
  }
  
  
  R <- tryCatch({
    integrate(Vectorize(f.y), lower = 0, upper = Inf)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- max(dat$speed, na.rm = TRUE)
        return(integrate(Vectorize(f.y), lower = 0, upper = ub.alt)$value)
      } 
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(R/normalize_cQ.bivariate(x = x0, model))
}


create_plot.dfQ <- function(dat,x,model,PX){
  # create a data frame for plotting the conditional density of speed at X
  # model is a mvdc object
  df <- data.frame(speed= dat,
             JointP = sapply(dat,function(y){pMvdc(c(x,y),model)})) %>%
    mutate(ConditionP = JointP/PX) %>% 
    mutate(ConditionalD = sapply(speed,cQ.bivariate,x=x,model = model)/
             normalize_cQ.bivariate(x,model)) %>% na.omit()
  
  return(df)
}

gofEVCopula_VC2 <- function(copula, x, N = 1000,method = c("mpl", "ml", "itau", "irho"),
       estimator = c("CFG", "Pickands"), m = 1000,verbose = interactive(),
       ties.method = c("max", "average", "first", "last", "random", "min"),
       fit.ties.meth = eval(formals(rank)$ties.method), ...){
  # This is a wrapper of copula::gofEVCopula for copula constructed from VC2copula class, the code is basically
  # the same as the original function. 
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
  u.fit <- if (ties.method == fit.ties.meth) 
    u
  else pobs(x, ties.method = fit.ties.meth)
  fcop <- fitCopula(copula, u.fit, method,...)@copula
  g <- seq(0, 1 - 1/m, by = 1/m)
  s <- .C("cramer_vonMises_Afun", as.integer(n), as.integer(m), 
          as.double(-log(u[, 1])), as.double(-log(u[, 2])), 
          as.double(A(fcop, g)), stat = double(2), as.integer(estimator == "CFG"))$stat
  s0 <- matrix(NA, N, 2)
  if (verbose) {
    pb <- txtProgressBar(max = N, style = if (isatty(stdout())) 3 else 1)
    on.exit(close(pb))
  }
  for (i in 1:N) {
    u0 <- pobs(rCopula(n, fcop), ties.method = ties.method)
    u0.fit <- if (ties.method == fit.ties.meth) u0 else pobs(u0, ties.method = fit.ties.meth)
    fcop0 <- fitCopula(copula, u0.fit, method, ...)@copula
    s0[i, ] <- .C("cramer_vonMises_Afun", as.integer(n), as.integer(m), 
                  as.double(-log(u0[, 1])), as.double(-log(u0[, 2])), 
                  as.double(A(fcop0, g)), stat = double(2), 
                  as.integer(estimator == "CFG"))$stat
    if (verbose) 
      setTxtProgressBar(pb, i)
  }
  structure(class = "htest", list(method = paste0("Parametric bootstrap based GOF test for EV copulas with argument 'method' set to ", 
                 sQuote(method), " and argument 'estimator' set to ", 
                 sQuote(estimator)), parameter = c(parameter = fcop@parameters), 
                 statistic = c(statistic = s[1]), 
                 p.value = (sum(s0[,1] >= s[1]) + 0.5)/(N + 1), data.name = deparse(substitute(x))))
}
