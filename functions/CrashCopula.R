
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


cQ.bivariate <-function(y,x,PX, model,ulim.alt=0.45){
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
      result[i] <- dMvdc(c(q1[i],y),model)
    }
    return(result)
  }
  
  # Integrate from x to Inf
  R <- tryCatch({
    integrate(integrand, lower = x, upper = Inf)$value/PX},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- ulim.alt
        return(integrate(integrand, lower = x, upper = ub.alt)$value/PX)
      }
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(R)
}

normalize_cQ.bivariate<- function(x,Pcrash, model,lb=0){
  # computes the infinite integral of the conditional density f(y|X >x)
  # use as a nomralization factor for the conditional density
  # lb is the lower bound of the integral, by defalut lb = 0
  c.y <- function(k) {
    temp <- cQ.bivariate(y = k, x = x, PX = Pcrash, model = model)
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
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(C) 
}

Injury.from_cQ_bivariate <-function(dat,model,Pcrash,severity){
  # computes the injury probability using c.bivariate, severity is a logistic regression model of injury given speed
  f.y <- function(k) {
    temp <- cQ.bivariate(y = k, x = 0, PX = Pcrash, model = model) * severity(k)
    return(temp)
  }
  
  
  R <- tryCatch({
    integrate(Vectorize(f.y), lower = min(dat$speed), upper = Inf)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- max(dat$speed, na.rm = TRUE)
        return(integrate(Vectorize(f.y), lower = min(dat$speed), upper = ub.alt)$value)
      } 
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(R/normalize_cQ.bivariate(x = 0, Pcrash, model))
}


create_plot.df <- function(dat,x,model,PX){
  # create a data frame for plotting the conditional density of speed at X
  # model is a mvdc object
  df <- data.frame(speed= dat,
             JointP = sapply(dat,function(y){pMvdc(c(x,y),model)})) %>%
    mutate(ConditionP = JointP/PX) %>% 
    mutate(ConditionalD = sapply(speed,cQ.bivariate,x=x,PX=PX,model = model)/
             normalize_cQ.bivariate(x,PX,model)) %>% 
    na.omit()
  
  return(df)
}
