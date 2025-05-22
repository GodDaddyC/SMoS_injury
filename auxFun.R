# auxFun
library(readxl)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(plotly)
library(extRemes)
library(evd)
library(copula)

#source("C:/Users/zh0815ch/OneDrive - Lund University/Digital assets/R functions and codes/Truncatedpbeved.R")
#source("C:/Users/zh0815ch/OneDrive - Lund University/Digital assets/R functions and codes/copulaForSMoS.R")


#compomentMax <- function(X,Y,block.size){
#   Mn1 <- blockmaxxer(as.data.frame(X),blen = block.size,span = ceiling(length(X) / block.size))
#   Mn2 <- blockmaxxer(as.data.frame(Y),blen = block.size,span = ceiling(length(Y) / block.size))
#   return(cbind(Mn1,Mn2))
# }

empProb <- function(data,x,type,plot=F){
  # emprical probability of the given bivariate data.
  if (type ==1){
    R <- data[which(data[,1] <= x[1] & data[,2] <= x[2]),]
  }
  else if (type ==2){
    R <- data[which(data[,1] >= x[1] & data[,2] >= x[2]),]
  }
  else if (type ==4){
    R <- data[which(data[,1] <= x[1] & data[,2] >= x[2]),]
  }
  else if (type ==3){
    R <- data[which(data[,1] > x[1] & data[,2] <= x[2]),]
  }
  R <- as.data.frame(R)
  if (plot){
    data <- as.data.frame(data)
    colnames(data) = c('V1','V2')
    return(ggplot(data,aes(x = V1, y = V2)) +geom_point(size = 1) + 
             geom_point(data =R ,aes(x= V1,y=V2),size = 1,col ='red')) 
  }
  else{
    return(nrow(R)/nrow(data))
  }
}

#rightEndTest <- function(data,block.size,thres,type,dep,q){
  # check whether the right end points of the model fitted with block size/threshold is greater than q (the first margin). 
#   if (type == 'FML'){
#     if (hasArg(block.size)){
#       block.maxima <- compomentMax(data[,1],data[,2],block.size)
#       fit <- fbvevd(block.maxima,model = dep,std.err = FALSE)
#       param <- fit$estimate
#       return(pgev(0,loc = param[1],scale = param[2],shape = param[3],lower.tail = F) > q)
#     }
#     else if (hasArg(thres)){
#       fit <- fbvpot(data,threshold = thres,model = dep,std.err = FALSE)
#       param <- fit$estimate
#       return(pevd(0,scale = param[1],shape = param[2],lower.tail = F,threshold = thres[1],type = 'GP') > q)
#     }
#   }
#   
#   else if (type == 'IFM'){
#     if (hasArg(block.size)){
#       block.maxima <- compomentMax(data[,1],data[,2],block.size)
#       fit <- fevd(x = Mn1, data = block.maxima,type = 'GEV')
#       param <- fit$results$par
#       return(pevd(0,loc = param[1],scale = param[2],shape = param[3],type = 'GEV',lower.tail = F) > q)
#     }
#     else if (hasArg(thres)){
#       fit <- fevd(x = data[,1],threshold = thres[1],type = 'GP')
#       param <- fit$results$par
#       return(pevd(0,scale = param[1],shape = param[2],lower.tail = F,type = 'GP') > q)
#     }
#   }
# }

#CompareSeparation <- function(Data,u1,u2,x,tail.type){
#   BPD.hr <- B.Extreme.FML(Data,thres = c(u1,u2),dep = 'hr',std.err = F)
#   BPD.log <- B.Extreme.FML(Data,thres = c(u1,u2),dep = 'log',std.err = F)
#   BPD.alog <- B.Extreme.FML(Data,thres = c(u1,u2),dep = 'alog',std.err = F)
#   BPD.neglog <- B.Extreme.FML(Data,thres = c(u1,u2),dep = 'neglog',std.err = F)
#   BPD.aneglog <- B.Extreme.FML(Data,thres = c(u1,u2),dep = 'aneglog',std.err = F)
#   BPD.bilog <- B.Extreme.FML(Data,thres = c(u1,u2),dep = 'bilog',std.err = F)
#   BPD.negbilog <- B.Extreme.FML(Data,thres = c(u1,u2),dep = 'negbilog',std.err = F) 
#   BPD.ct <- B.Extreme.FML(Data,thres = c(u1,u2),dep = 'ct',std.err = F)
#   result.Prob <- c(evalProbPOT1(BPD.hr,x,tail.type)
#                      *empProb(Data,c(u1,u2),2,plot = F),
#                      evalProbPOT1(BPD.log,x,tail.type) *
#                        empProb(Data,c(u1,u2),2,plot = F),
#                      evalProbPOT1(BPD.alog,x,tail.type) *
#                        empProb(Data,c(u1,u2),2,plot = F),
#                      evalProbPOT1(BPD.neglog,x,tail.type) *
#                        empProb(Data,c(u1,u2),2,plot = F),
#                      evalProbPOT1(BPD.aneglog,x,tail.type) *
#                        empProb(Data,c(u1,u2),2,plot = F),
#                      evalProbPOT1(BPD.bilog,x,tail.type) *
#                        empProb(Data,c(u1,u2),2,plot = F),
#                      evalProbPOT1(BPD.negbilog,x,tail.type) *
#                        empProb(Data,c(u1,u2),2,plot = F),
#                      evalProbPOT1(BPD.ct,x,tail.type) *
#                        empProb(Data,c(u1,u2),2,plot = F))
#   return(result.Prob)
# }

#CIprobGEV <- function(x,mu,sig,gamma,V){
  # CI for GEV prob using delta method, theta is a list of mle, V is the covariance matrix
#   d.mu <- -pgev(x,loc=mu,scale=sig,shape=gamma)/sig * (1+gamma*(x-mu)/sig)^(-(1+gamma)/gamma)
#   d.sig <- -pgev(x,loc=mu,scale=sig,shape=gamma)* (x-mu)/sig^2 * (1+gamma*(x-mu)/sig)^(-(1+gamma)/gamma)
#   d.gamma <- -pgev(x,loc=mu,scale=sig,shape=gamma)/gamma^2 * (1+gamma*(x-mu)/sig)^(-(1+gamma)/gamma) * 
#     ((1+gamma*(x-mu)/sig) * log(1+gamma*(x-mu)/sig) - gamma * (x-mu)/sig)
#   grad <- c(d.mu,d.sig,d.gamma)
#   return(c(pgev(x,loc=mu,scale=sig,shape=gamma)-1.96* t(grad) %*% V %*% grad,
#            pgev(x,loc=mu,scale=sig,shape=gamma)+1.96* t(grad) %*% V %*% grad))
# }

CI.prob.GP <- function(x,GP,lower.tail=T,alp=0.05,symmetry="both"){
# CI for the probability of GEV and GPD using delta method, u is the threshold, V is the covariance matrix
    u <- GP$threshold
    sig <- GP$results$par[1]
    if (GP$type == "GP"){
      gamma <- GP$results$par[2]
      d.sig <- (1/gamma) * (x-u)/sig^2 * (1+gamma*(x-u)/sig)^(-(1+gamma)/gamma)
      d.gamma <- (1 + gamma *(x-u)/sig)^(-1/gamma-1) * ((1+gamma*(x-u)/sig) * log(1+gamma*(x-u)/sig)-
        (gamma *(x-u)/sig)) / gamma^2
      grad <- c(d.sig,d.gamma)
      V <- solve(GP$results$hessian,diag(c(1,1)))
      R <- switch(symmetry,
             both = c(pevd(x,scale=sig,shape=gamma,threshold = u,type = 'GP',lower.tail = F)-qnorm(1-alp/2)* sqrt(t(grad) %*% V %*% grad/GP$npy),
               pevd(x,scale=sig,shape=gamma,threshold = u,type = 'GP',lower.tail = F),
               pevd(x,scale=sig,shape=gamma,threshold = u,type = 'GP',lower.tail = F)+qnorm(1-alp/2)* sqrt(t(grad) %*% V %*% grad/GP$npy)),
             lower = c(pevd(x,scale=sig,shape=gamma,threshold = u,type = 'GP',lower.tail = F)-qnorm(1-alp)* sqrt(t(grad) %*% V %*% grad/GP$npy),
                       pevd(x,scale=sig,shape=gamma,threshold = u,type = 'GP',lower.tail = F)),
             upper = c(pevd(x,scale=sig,shape=gamma,threshold = u,type = 'GP',lower.tail = F),
                       pevd(x,scale=sig,shape=gamma,threshold = u,type = 'GP',lower.tail = F)+qnorm(1-alp)* sqrt(t(grad) %*% V %*% grad/GP$npy))
        )
      return(R)
      }
    if (GP$type =="Exponential"){
      d.sig <- 1/sig^2 * (x-u) * exp(-(x-u)/sig)
      V <- invisible(summary(GP)$cov.theta)
      R <- switch(symmetry,
             both = c(pexp(x-u,rate=1/sig,lower.tail = F)-qnorm(1-alp/2)* d.sig^2 * sqrt(V /GP$npy),
               pexp(x-u,rate=1/sig,lower.tail = F),
               pexp(x-u,rate=1/sig,lower.tail = F)+qnorm(1-alp/2)* d.sig^2 * V/sqrt(GP$npy)),
             lower = c(pexp(x-u,rate=1/sig,lower.tail = F)-qnorm(1-alp)* d.sig^2 * sqrt(V /GP$npy),
                       pexp(x-u,rate=1/sig,lower.tail = F)),
             upper = c(pexp(x-u,rate=1/sig,lower.tail = F),
              pexp(x-u,rate=1/sig,lower.tail = F)+qnorm(1-alp)* d.sig^2 * sqrt(V/GP$npy))
             )
    }
    else{
      R <- "Distribution type is not supported"
    }
    return(R)
}

CI.xF <- function(model,alp=0.05,symmetry="both"){
  # sym = "both", "upper" or "lower"
  if (model$type=="Exponential"){
    return("Infinite right end points")
  }
  else{
    if (model$result$par[2]>0){
      return("Infinite right end points")
    }
    x.F <- model$threshold-model$results$par[1]/model$results$par[2]
    Dx.F <- c(-1/model$results$par[2], model$results$par[1]/ model$results$par[2]^2)
    V <- solve(model$results$hessian,diag(c(1,1)))
    switch (symmetry,
             both =  c(x.F-qnorm(1-alp/2) * sqrt(t(Dx.F) %*% V %*% Dx.F/model$npy),x.F,
                        x.F+qnorm(1-alp/2) * sqrt(t(Dx.F) %*% V %*% Dx.F/model$npy)),
             lower = c(x.F-qnorm(1-alp) * sqrt(t(Dx.F) %*% V %*% Dx.F/model$npy),x.F),
             upper = c(x.F,x.F+qnorm(1-alp) * sqrt(t(Dx.F) %*% V %*% Dx.F/model$npy))
                   )
}}


Separation.Check <- function(model1,model2,...){
  # the identifability of null events, ... pass to CI
  x1 <- CI.xF(model1,...)
  x2 <- CI.xF(model1,...)
  if (all(c(is.character(x1),is.character(x2)) ) ){
    return("Fail to separate")
  }
  if (any(c(is.numeric(x1),is.numeric(x2)) ) ){
    C <- if_else(c(is.numeric(x1),is.numeric(x2)),true = "finite",false = "infinite")
    CC<- eval(parse(text=sprintf("max(x%o)",which(C=="finite") )) ) %>%{.<0}
    if (!CC){
      return("Fail to separate")
    }
    else{
      return("succeed")
    }
  }
  else{
    if ((min(x1)>0 & min(x2)>0)){
      return("Fail to separate")
    }
    if (any(max(x1)<0, min(x2)<0)){
      C <- if_else(c(max(x1)<0,max(x2)<0),true = "N",false = "R")
      if (all(C=="R")){
        return("Fail to separate")
      }
      else{
        return("succeed")
      }
    }
  }
}




subsetDataFrame <- function(df, criteria) {
  # Convert the criteria to a logical expression
  expression <- parse(text = paste(criteria, collapse = " & "))
  
  # Subset the data frame based on the logical expression
  subset(df, eval(expression))
}

quickCheck <- function(Dat,criteria){
  Dat.C1 <- subsetDataFrame(Dat,criteria)
  Dat.C2 <- Dat[!rownames(Dat) %in% rownames(Dat.C1),]
  u1.C1 <- quantile(x=Dat1.C1$Mindist,probs =0.1)
  u1.C2 <- quantile(x=Dat1.C2$Mindist,probs =0.1)
  
  M.C1 <- fevd(x=-Dat.C1$Mindist,threshold = -u1.C1,type="GP")
  M.C2 <- fevd(x=-Dat.C2$Mindist,threshold = -u1.C2,type="GP")
  return(Separation.Check(M.C1,M.C2))
}

DetailCheck <- function(Dat,criteria,alp1=0.05,sym1="both",alp2=0.05,sym2="both"){
  #alp1, sym1 for x.f, 2 for prob of zero
  
  Dat1 <- subsetDataFrame(Dat,criteria)
  Dat2 <- Dat[!rownames(Dat) %in% rownames(Dat1),]
  
  u <- min(quantile(x=Dat2$Mindist,probs =0.1),quantile(x=Dat1$Mindist,probs =0.1))
  
  model1 <- fevd(x=-Dat1$Mindist,threshold = -u,type="GP")
  CIpar.1 <- ci(model1,type="parameter",alpha = alp1)
  if (CIpar.1[2,1]<0 & CIpar.1[2,3]>0){
    model1 <- fevd(x=-Dat1$Mindist,threshold = -u,type="Exponential")
  }
  model2 <- fevd(x=-Dat2$Mindist,threshold = -u,type="GP")
  CIpar.2 <- ci(model2,type="parameter",alpha = alp1)
  if (CIpar.2[2,1]<0 & CIpar.2[2,3]>0){
    model2 <- fevd(x=-Dat2$Mindist,threshold = -u,type="Exponential")
  }
  
  CI1.1 <- CI.xF(model1,alp = alp1,symmetry = sym1)
  CI1.2 <- CI.xF(model2,alp = alp1,symmetry = sym1)
  CI2.1 <- CI.prob.GP(0,model1,alp = alp2,symmetry = sym2)
  CI2.2 <- CI.prob.GP(0,model2,alp=alp2,symmetry = sym2)

  cat(sprintf("the %f CI (%s) for the scale parameter of model1 is",alp1,sym1),CIpar.1[1,],"\n")
  cat(sprintf("the %f CI (%s) for the shape parameter of model1 is",alp1,sym1),CIpar.1[2,],"\n")
  cat(sprintf("the %f CI (%s) for the scale parameter of model2 is",alp1,sym1),CIpar.2[1,],"\n")
  cat(sprintf("the %f CI (%s) for the shape parameter of model2 is",alp1,sym1),CIpar.2[2,],"\n")
  cat(sprintf("the %f CI (%s) for the right end point of model1 is",alp1,sym1),CI1.1,"\n")
  cat(sprintf("the %f CI (%s) for the right end point of model2 is",alp1,sym1),CI1.2,"\n")
  cat(sprintf("the %f CI (%s) for the accident prob from model1 is",alp2,sym2),CI2.1,"\n")
  cat(sprintf("the %f CI (%s) for the accident prob from model2 is",alp2,sym2),CI2.2,"\n")
}

Sep.Reward <- function(Dat,criteria,alp1=0.05,sym1="both",alp2=0.05,sym2="both"){
  Dat1 <- subsetDataFrame(Dat,criteria)
  Dat2 <- Dat[!rownames(Dat) %in% rownames(Dat1),]
  
  u <- min(quantile(x=Dat2$Mindist,probs =0.1),quantile(x=Dat1$Mindist,probs =0.1))
  
  model1 <- fevd(x=-Dat1$Mindist,threshold = -u,type="GP")
  CIpar.1 <- ci(model1,type="parameter",alpha = alp1)
  if (CIpar.1[2,1]<0 & CIpar.1[2,3]>0){
    model1 <- fevd(x=-Dat1$Mindist,threshold = -u,type="Exponential")
  }
  model2 <- fevd(x=-Dat2$Mindist,threshold = -u,type="GP")
  CIpar.2 <- ci(model2,type="parameter",alpha = alp1)
  if (CIpar.2[2,1]<0 & CIpar.2[2,3]>0){
    model2 <- fevd(x=-Dat2$Mindist,threshold = -u,type="Exponential")
  }
  CI2.1 <- CI.prob.GP(0,model1,alp = alp2,symmetry = sym2)
  CI2.2 <- CI.prob.GP(0,model2,alp= alp2,symmetry = sym2)
  return(min(CI2.1)-max(CI2.2))
  #return(min(CI2.1)/max(CI2.2))
}
