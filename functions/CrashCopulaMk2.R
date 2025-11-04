cQ.bivariate <-function(y,x,model){
  # approximates the conditional density f(y|X >x) by integrating fxy = f(x=x, y = y) over x
  # model is an mvdc object
  # if integral is non-finite, change ulim to a smaller value
  
  integrand <- function(q1) {
    result <- dCopula(c(q1,y),model) 
    
    return(result)
  }
  
  # Integrate from x to Inf
  R <- tryCatch({
    integrate(Vectorize(integrand), lower = x, upper = 1)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- 1 - 1e-7
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
    integrate(Vectorize(c.y), lower = lb, upper = 1)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- 1-1e-7
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

Injury.from_cQ_bivariate <-function(model,severity,x0,PX,P2,LB=NULL,UB=NULL){
  # computes the injury probability using c.bivariate, severity is a logistic regression model of injury given speed
  f.y <- function(k) {
    S <- qgamma(k,shape = P2$estimate[1],rate = P2$estimate[2])
    temp <- cQ.bivariate(y = k, x = x0, model = model) * 
      dgamma(S,shape = P2$estimate[1],rate = P2$estimate[2])*
      severity(S)
    return(temp)
  }
  
  
  R <- tryCatch({
    integrate(Vectorize(f.y), lower = ifelse(is.null(LB),0,
        pgamma(LB,shape = P2$estimate[1],rate = P2$estimate[2])), 
      upper = ifelse(is.null(UB),1,
        pgamma(UB,shape = P2$estimate[1],rate = P2$estimate[2]))
      )$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- 1 - 1e-5
        return(integrate(Vectorize(f.y), lower = 0, upper = ub.alt)$value)
      } 
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(R/PX)
}

create_plot.dfQ <- function(dat,x,model,PX,P2){
  # create a data frame for plotting the conditional density of speed at X
  # model is a mvdc object
  df <- data.frame(speed.u= dat,speed = qgamma(dat,shape = P2$estimate[1],rate = P2$estimate[2]),
                   JointP = sapply(dat,function(y){pCopula(c(x,y),model)})) %>%
    mutate(ConditionP = JointP/PX) %>% 
    mutate(ConditionalD = sapply(speed.u,cQ.bivariate,x=x,model = model)/PX * 
             dgamma(speed,shape = P2$estimate[1],rate = P2$estimate[2])) %>% na.omit()
  return(df)
}

Kc_test <- function(dat,Copula,n=1e4,m=200){
  Cn <- empCopula(X = dat, ties.method="random")
  Sim <- rCopula(n, copula = Cn)
  Model <- rCopula(n,copula = Copula)
  
  temp1 <- pCopula(Sim,copula = Cn)
  temp2 <- pCopula(Model, copula = Copula)
  
  Kc_emp <- sapply(seq(0,1,1/m),function(k) mean(temp1 < k))
  Kc_model <- sapply(seq(0,1,1/m),function(k) mean(temp2 < k))
  
  result <- CCC(Kc_emp, Kc_model, ci="z-transform")
  return(result$rho.c)
  
}
