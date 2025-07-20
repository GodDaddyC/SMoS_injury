# auxFun

lapply(list.files('functions', pattern = "\\.R$", full.names = TRUE), source)


Tail.prob.GP <- function(x,model,lower.tail=FALSE,conditional = FALSE,use.phi=FALSE){
  # model is a fevd object, requires extRemes it assumes use.phi=TRUE, if not, 
  # get rid of the exp() in sigma
  if (conditional){
    if (model$method =='Bayesian'){
      prob.chain <- sapply(1:dim(model$chain.info)[1],function(k){
        pevd(x,threshold = model$threshold, scale = exp(model$result[k,1]), 
             shape = model$result[k,2],lower.tail = lower.tail,type = "GP")
      })
      return(prob.chain)
    }
    if (model$method == "MLE"){
      return(ifelse(use.phi,pevd(x,threshold = model$threshold, scale = exp(model$result$par[1]), 
                                 shape = model$result$par[2],lower.tail = lower.tail,type = "GP"),
                    pevd(x,threshold = model$threshold, scale = model$result$par[1], 
                         shape = model$result$par[2],lower.tail = lower.tail,type = "GP")))
    }
  }
  else{
    if (model$method =='Bayesian'){
      prob.chain <- sum(model$x>model$threshold)/model$n *
        sapply(1:dim(model$chain.info)[1],function(k){
          pevd(x,threshold = model$threshold, scale = exp(model$result[k,1]), 
               shape = model$result[k,2],lower.tail = lower.tail,type = "GP")})
      return(prob.chain)
    }
    if (model$method == "MLE"){
      pu <- sum(model$x>model$threshold)/model$n
      return(ifelse(use.phi,pu*pevd(x,threshold = model$threshold, scale = exp(model$result$par[1]), 
                                    shape = model$result$par[2],lower.tail = lower.tail,type = "GP"),
                    pu*pevd(x,threshold = model$threshold, scale = model$result$par[1], 
                            shape = model$result$par[2],lower.tail = lower.tail,type = "GP")))
    }
  }
}

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
    R <- as.data.frame(R)
    colnames(data) = c('mar1','mar2')
    colnames(R) = c('mar1','mar2')
    return(ggplot(data,aes(x = mar1, y = mar2)) +geom_point(size = 1) + 
             geom_point(data =R ,aes(x= mar1,y=mar2),size = 1,col ='red')+
             geom_vline(xintercept = x[1],col = 'blue') +
             geom_hline(yintercept = x[2],col = 'blue') +
             labs(x = "mar1", y = "mar2",title="empirical prob")) 
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
