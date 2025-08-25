# semi-parametric approach, copula model is non-parametric while margins are

emp.CN <- empCopula(Cop.dat.1,smoothing="beta")
emp.SE <- empCopula(Cop.dat.2,smoothing="beta")
dCopula(c(0.5,0.5),emp.CN)
pCopula(c(0.5,0.5),emp.CN)

tt <- kdecop(Cop.dat.1)
dkdecop(c(0.5,0.5),tt)
pkdecop(c(0.5,0.5),tt)


Qcrash.1 <- pevd(x0.1,threshold = v.1,scale = POT.1$results$par[1],
                 shape = POT.1$results$par[2], type = "GP",lower.tail = FALSE)
Qcrash.2 <- pevd(x0.2,threshold = v.2,scale = POT.2$results$par[1],
                 shape = POT.2$results$par[2], type = "GP",lower.tail = FALSE)

x0.1.un <- 1 - Qcrash.1
x0.2.un <- 1 - Qcrash.2

cQ.bivariate.nonpar <-function(v,u,PX, model,ulim.alt=0.999,mode=0){
  # approximates the conditional density f(y|X >x) by integrating fxy = f(x=x, y = y) over x
  # model is an mvdc object
  # if integral is non-finite, change ulim to a smaller value
  
  integrand <- function(q1,...) {
    result <- numeric(length(q1))
    for(i in seq_along(q1)) {
      #result[i] <- dCopula(c(q1[i],v),model)
      result[i] <- ifelse(mode==0,dCopula(c(q1[i],v),model),dkdecop(c(q1[i],v), model))
    }
    return(result)
  }
  
  # Integrate from x to Inf
  R <- tryCatch({
    integrate(integrand, lower = u, upper = 1)$value/PX},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- ulim.alt
        return(integrate(integrand, lower = u, upper = ub.alt)$value/PX)
      }
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(R)
}

normalize_cQ.bivariate.Nonpar<- function(x,Pcrash, model,lb=0){
  # computes the infinite integral of the conditional density f(y|X >x)
  # use as a nomralization factor for the conditional density
  # lb is the lower bound of the integral, by defalut lb = 0
  c.y <- function(k) {
    temp <- cQ.bivariate.nonpar(v = k, u = x, PX = Pcrash, model = model)
    return(temp)
  }
  
  C <- tryCatch({
    integrate(Vectorize(c.y), lower = lb, upper = 1)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- 0.999
        return(integrate(Vectorize(c.y), lower = lb, upper = ub.alt)$value)
      } 
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(C) 
}
cQ.bivariate.nonpar(0.5,x0.1.un,Qcrash.1,tt,mode=1)
cQ.bivariate.nonpar(0.5,x0.1.un,Qcrash.1,emp.CN,mode=0)/1.72

t1 <- sapply(seq(0.01,0.99,by=0.01),function(k) cQ.bivariate.nonpar(k,x0.1.un,Qcrash.1,tt,mode=1))

create_plot.dfQ.nonpar <- function(dat,x,model,PX){
  # create a data frame for plotting the conditional density of speed at X
  # model is a mvdc object
  df <- data.frame(speed= dat,
                   JointP = sapply(dat,function(y){pMvdc(c(x,y),model)})) %>%
    mutate(ConditionP = JointP/PX) %>% 
    mutate(ConditionalD = sapply(speed,cQ.bivariate,x=x,PX=PX,model = model)/
             normalize_cQ.bivariate(x,PX,model)) %>% na.omit()
  
  return(df)
}
