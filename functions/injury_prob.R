# functions for computing injury probability using densities


# BPOT --------------------------------------------------------------------

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


# CBPOT -------------------------------------------------------------------

Injury.from_cQ_bivariate <-function(model,severity,x0,PX,P2,LB=NULL,UB=NULL,Pu = 1){
  # computes the injury probability using c.bivariate, severity is a logistic regression model of injury given speed
  f.y <- function(k) {
    S <- pgamma(k,shape = P2$estimate[1],rate = P2$estimate[2])
    temp <- cQ.bivariate(y = S, x = x0, model = model) *severity(k)*
      dgamma(k,shape = P2$estimate[1],rate = P2$estimate[2])
    return(temp)
  }
  
  
  R <- tryCatch({
    integrate(Vectorize(f.y), 
              lower = ifelse(is.null(LB),0, LB),#pgamma(LB,shape = P2$estimate[1],rate = P2$estimate[2])), 
              upper = ifelse(is.null(UB),60 ,UB)#,pgamma(UB,shape = P2$estimate[1],rate = P2$estimate[2]))
    )$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        #ub.alt <- 1 - 1e-3
        ub.alt <- 55
        return(integrate(Vectorize(f.y),
                         lower = ifelse(is.null(LB),0,LB),#pgamma(LB,shape = P2$estimate[1],rate = P2$estimate[2])),
                         upper = ifelse(is.null(UB),min(UB,ub.alt) )#,pgamma(UB,shape = P2$estimate[1],rate = P2$estimate[2]))
        )$value)
      }
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(R/PX * Pu)
}

Injury.from_cQ_bivariate_E <-function(model,severity,x0,PX,P2,LB=0.5,UB=NULL,Pu = 1,N=1000,age.mean=40,
                                      age.sd = 10){
  # computes the injury probability using c.bivariate, severity is contains an assumed distribution
  # for age
  f.y <- function(k) {
    S <- pgamma(k,shape = P2$estimate[1],rate = P2$estimate[2])
    temp <- cQ.bivariate(y = S, x = x0, model = model) *
      dgamma(k,shape = P2$estimate[1],rate = P2$estimate[2])
    age <- rnorm(N,mean=age.mean,sd=age.sd)
    age <- age[age>0]
    return(temp*mean(severity(k,age)))
  }
  
  
  R <- tryCatch({
    integrate(Vectorize(f.y), 
              lower = ifelse(is.null(LB),0, LB),#pgamma(LB,shape = P2$estimate[1],rate = P2$estimate[2])), 
              upper = ifelse(is.null(UB),60 ,UB)#,pgamma(UB,shape = P2$estimate[1],rate = P2$estimate[2]))
    )$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        #ub.alt <- 1 - 1e-3
        ub.alt <- 55
        return(integrate(Vectorize(f.y),
                         lower = ifelse(is.null(LB),0,LB),#pgamma(LB,shape = P2$estimate[1],rate = P2$estimate[2])),
                         upper = ifelse(is.null(UB),min(UB,ub.alt) )#,pgamma(UB,shape = P2$estimate[1],rate = P2$estimate[2]))
        )$value)
      }
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(R /PX * Pu)
}