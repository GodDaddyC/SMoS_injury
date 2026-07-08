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

empProb <- function(data,x,type,plot=FALSE){
  # emprical probability of the given bivariate data.
  if (type ==1){
    R <- data[which(data[,1] <= x[1] & data[,2] <= x[2]),]
  }
  else if (type ==2){
    R <- data[which(data[,1] >= x[1] & data[,2] >= x[2]),]
  }
  else if (type ==3){
    R <- data[which(data[,1] <= x[1] & data[,2] >= x[2]),]
  }
  else if (type ==4){
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



runExecutionBPOT <- function(Dat.CN,Dat.SE,model,thres1,thres2,speed1.ub =60,speed2.ub=speed1.ub,xcrash=c(x0.1,x0.2)){
  M.1 <- fbvpot(x = Dat.CN,model = model,threshold = thres1)
  M.2 <- fbvpot(x = Dat.SE,model = model,threshold = thres2)
  Pcrash.1 <- pevd(xcrash[1],threshold = thres1[1], scale = M.1$estimate[1],shape = M.1$estimate[2],
                   lower.tail = FALSE,type = "GP") * M.1$nat[1]/M.1$n
  
  Pcrash.2 <- pevd(xcrash[2],threshold = thres2[1], scale = M.2$estimate[1],shape = M.2$estimate[2],
                   lower.tail = FALSE,type = "GP") * M.2$nat[2]/M.2$n
  ss.1T <- seq(thres1[2],speed1.ub,(speed1.ub- thres1[2])/150)
  ss.2T <- seq(thres2[2],speed2.ub,(speed2.ub- thres2[2])/150)
  
  plot.df.1 <- create_plot.df(ss.1T,x=xcrash[1],model=M.1,PX=Pcrash.1)
  plot.df.2 <- create_plot.df(ss.2T,x=xcrash[2],model=M.2,PX=Pcrash.2)
  
  return(list(M.1=M.1,M.2=M.2,Pcrash.1=Pcrash.1,Pcrash.2=Pcrash.2,
              plot.df.1=plot.df.1,plot.df.2=plot.df.2))
  
}

runExecutionCBPOT <- function(Cop.Dat.1,Cop.Dat.2,copula,speed.ub1 =60,speed.ub2=speed.ub1,P2.1,P2.2,
                              Qcrash.1,Qcrash.2,method="mpl",Pu=0.15,...){
  # copula argument should be given as a list
  Cop.1 <- fitCopula(copula[[1]], data = Cop.dat.1,...)
  Cop.2 <- fitCopula(copula[[2]], data = Cop.dat.2,...)
  CM.1 <- Q.distr.param(Cop.1@copula,mar1=POT.1,mar2 = Conseq.1,type = 4)
  CM.2 <- Q.distr.param(Cop.2@copula,mar1=POT.2,mar2 = Conseq.2,type = 4)
  CM0.1 <- Q.distr.param(Cop.1@copula,mar1=POT.1,mar2 = Conseq.1,type = 2)
  CM0.2 <- Q.distr.param(Cop.2@copula,mar1=POT.2,mar2 = Conseq.2,type = 2)
  
  x0.1.un <- 1 - Qcrash.1
  x0.2.un <- 1 - Qcrash.2
  
  sq.1 <- seq (0.5,speed.ub1,0.05)
  sq.2 <- seq(0.5,speed.ub2,0.05) # start from 0.5 for better numerical stability
  s1.un <- pgamma(sq.2,shape = P2.1$estimate[1],rate = P2.1$estimate[2])
  s2.un <- pgamma(sq.2,shape = P2.2$estimate[1],rate = P2.2$estimate[2])
  plot.dfQ.1 <- create_plot.dfQ(s1.un,x=x0.1.un,model = Cop.1@copula,PX=Qcrash.1,P2=P2.1,Pu=Pu)
  plot.dfQ.2 <- create_plot.dfQ(s2.un,x=x0.2.un,model = Cop.2@copula,PX=Qcrash.2,P2=P2.2,Pu=Pu)
  return(list(CM.1 =CM.1,CM.2=CM.2,Cop.1=Cop.1, Cop.2 = Cop.2,Pcrash.1=Qcrash.1 * Pu,Pcrash.2=Qcrash.2 * Pu,
              plot.dfQ.1=plot.dfQ.1,plot.dfQ.2=plot.dfQ.2, CM0.1=CM0.1,CM0.2=CM0.2))
}

GP.param.screening <- function (data, orderlim = NULL, tlim = NULL, alpha = 0.1,min.thres=4,if.CI=FALSE,...) {
  # find the largest threshold that have ci between [-0.5,1], then use this threshold
  # in try.threshold for tcplot
  # if.CI: should the scrrening depending on CI or not
  # 4 is the minimum value for min.thres due to the property of pickand estimator
  if (is.unsorted(data)) {
    data = sort(data, decreasing = TRUE)
  }
  else {
    if (data[1] < data[length(data)]) 
      data = rev(data)
  }
  n = length(data)
  if (!is.null(tlim)) {
    if (length(tlim) != 2) 
      stop("threshold range tlim must be a numeric vector of length 2")
    if (tlim[2] <= tlim[1]) 
      stop("a range of thresholds must be specified by tlim")
    if (is.null(orderlim)) {
      orderlim = c(sum(data >= tlim[2]), max(sum(data >= tlim[1]), 1))
    }
  }
  if (!is.null(orderlim)) {
    if (length(orderlim) != 2 | mode(orderlim) != "numeric") 
      stop("order statistic range orderlim must be an integer vector of length 2")
    if (orderlim[2] <= orderlim[1]) 
      stop("a range of order statistics must be specified by orderlim")
    if (orderlim[2] > floor(n/min.thres)) 
      stop("maximum order statistic in orderlim must be less than floor(n/%d)", min.thres)
  }
  else {
    orderlim = c(3, floor(n/min.thres))
  }
  
  if (max(orderlim) <= 10) 
    stop("must have more than 10 order statistics")
  norder = (diff(orderlim) + 1)
  if (norder < 2) 
    stop("must be more than 2 order statistics considered")
  maxks <- floor(n/min.thres)
  ks <- 1:maxks
  Pick <- log((data[ks] - data[2 * ks])/(data[2 * ks] - data[4 * ks]))/log(2)
  Pickse <- Pick * sqrt((2^(2 * Pick + 1) + 1))/2/(2^Pick - 1)/log(2)/sqrt(ks)
  pickresults <- data.frame(data[ks], ks, Pick, se.H = Pickse)
  if (!is.null(alpha)) {
    Pickci <- cbind(Pick - qnorm(1 - alpha/2) * Pickse, Pick + qnorm(1 - alpha/2) * Pickse)
    pickresults <-  cbind(pickresults, cil.Pick = Pickci[, 1], ciu.Pick = Pickci[, 2])
  }
  if (if.CI){
    thres.range <- pickresults %>% filter(cil.Pick > -0.5 & ciu.Pick < 1)
  }
  else{
    thres.range <- pickresults %>% filter(Pick > -0.5 & Pick < 1)
  }
  
  if (nrow(thres.range) > 0){
    return(thres.range)
  } 
  else{
    return(NA)
  }
}

thres.auto <- function(data,...){
  thres.range <- GP.param.screening(data,...)
  if (!is.data.frame(thres.range)) {
    return(NA)
  } 
  else {
    thres.range <- thres.range$data.ks.
    p_vals <- c()
    for (i in 1:length(thres.range)) {
      p_temp <- tryCatch(
        {gpdAd(data=data[data>thres.range[i]],bootstrap = TRUE,bootnum = 1000)$p.value}, 
        error = function(e) NA
      )
      p_vals <- c(p_vals, p_temp)
    }
    p_vals <- p_vals[!is.na(p_vals)]
    thres.range <- thres.range[!is.na(p_vals)]
    if (sum(p_vals > 0.05) == 0) {
      return(min(thres.range,na.rm = TRUE))
    }
    else{
      return(min(thres.range[p_vals > 0.05], na.rm = TRUE))
    }
  }
}
