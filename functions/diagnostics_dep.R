# diagonstics of the dependence structure

# BPOT --------------------------------------------------------------------

TbevdPlot <- function(EVmodel,dat,k=10){
  # produce four plots: the 2*2 plot
  # k used for bernstein polynomial approximation
  # ... pass to pickands.Nonpar
  A.np.conf <- pickands.Nonpar(dat = dat, mar1 = EVmodel$estimate[1:2], mar2 = EVmodel$estimate[3:4],
                               thres = EVmodel$threshold, eta = EVmodel$nat[1:2]/EVmodel$n,
                               est = "cfg", CI = TRUE, d = 2, k = k, ifplot = TRUE)
  A.np<- pickands.Nonpar(dat = dat, mar1 = EVmodel$estimate[1:2], mar2 = EVmodel$estimate[3:4],
                         thres = EVmodel$threshold, eta = EVmodel$nat[1:2]/EVmodel$n,
                         est = "cfg", CI = FALSE, d = 2, k = k, ifplot = FALSE)
  
  # par(mfrow=c(2,2))
  if (EVmodel$model %in% c("log","hr")){
    CB <- confint(EVmodel,parm="dep")
    abvevd(dep=EVmodel$estimate[5],model=EVmodel$model,plot = TRUE,add=TRUE,col="red")
    abvevd(dep=CB[1],model=EVmodel$model,plot = TRUE,add=TRUE,col="red",lty=5)
    abvevd(dep=CB[2],model=EVmodel$model,plot = TRUE,add=TRUE,col="red",lty=5)
    legend("bottomright", 
           legend = c("Nonparametric estimates", "Parametric estimates"), 
           col = c("black", "red"), lty = c(1, 2), 
           bty = "n")
    title(main = sprintf("Dependence diagonstics %s",EVmodel$model))
    
    spec.dens.plot <- data.frame(t = seq(0,1,0.005)) %>% 
      mutate(h.Nonpar = sapply(t,A_bp_approx,A_bp=A.np$beta,ord="2")/2) %>%
      mutate(h.Nonpar.UP = sapply(t,A_bp_approx,A_bp=A.np.conf$up.beta,ord="2")/2) %>%
      mutate(h.Nonpar.LW = sapply(t,A_bp_approx,A_bp=A.np.conf$low.beta,ord="2")/2) %>%
      mutate(h.Par = sapply(t,hbvevd,dep=EVmodel$estimate[5],model=EVmodel$model,half=TRUE)) %>%
      mutate(h.Par.UP = sapply(t,hbvevd,dep=CB[2],model=EVmodel$model,half=TRUE)) %>%
      mutate(h.Par.LW = sapply(t,hbvevd,dep=CB[1],model=EVmodel$model,half=TRUE))
  } 
  if (EVmodel$model %in% c("ct","bilog","negbilog")){
    CB_alpha <- confint(EVmodel,parm="beta")
    CB_beta <- confint(EVmodel,parm="alpha")
    abvevd(alpha=EVmodel$estimate[6],beta=EVmodel$estimate[5],model=EVmodel$model,
           plot = TRUE,add=TRUE,col="red")
    abvevd(alpha=CB_alpha[1],beta=CB_beta[1],model=EVmodel$model,plot = TRUE,add=TRUE,col="red",lty=5)
    abvevd(alpha=CB_alpha[2],beta=CB_beta[2],model=EVmodel$model,plot = TRUE,add=TRUE,col="red",lty=5)
    legend("bottomright", 
           legend = c("Nonparametric estimates", "Parametric estimates"), 
           col = c("black", "red"), lty = c(1, 2), 
           bty = "n")
    title(main = sprintf("Dependence diagonstics %s",EVmodel$model))
    
    spec.dens.plot <- data.frame(t = seq(0,1,0.005)) %>% 
      mutate(h.Nonpar = sapply(t,A_bp_approx,A_bp=A.np$beta,ord="2")/2) %>%
      mutate(h.Nonpar.UP = sapply(t,A_bp_approx,A_bp=A.np.conf$up.beta,ord="2")/2) %>%
      mutate(h.Nonpar.LW = sapply(t,A_bp_approx,A_bp=A.np.conf$low.beta,ord="2")/2) %>%
      mutate(h.Par = sapply(t,hbvevd,alpha=EVmodel$estimate[5],beta=EVmodel$estimate[6],
                            model=EVmodel$model,half=TRUE)) %>%
      mutate(h.Par.UP = sapply(t,hbvevd,alpha=CB_alpha[2],beta=CB_beta[2],model=EVmodel$model,half=TRUE)) %>%
      mutate(h.Par.LW = sapply(t,hbvevd,alpha=CB_alpha[1],beta=CB_beta[1],model=EVmodel$model,half=TRUE))
  }
  
  ggplot(spec.dens.plot, aes(x=t,y=h.Par,color="Model fitted")) + geom_line() +
    geom_line(aes(y=h.Par.LW,color="Model fitted"),linetype="dashed") +
    geom_line(aes(y=h.Par.UP,color="Model fitted"),linetype="dashed") +
    geom_line(aes(y=h.Nonpar.LW,color="Nonparametric"),linetype="dashed") +
    geom_line(aes(y=h.Nonpar.UP,color="Nonparametric"),linetype="dashed") +
    geom_line(aes(y=h.Nonpar,color="Nonparametric")) +
    labs(x="t",y="Spectral density",title = sprintf("Spectral density %s",EVmodel$model)) +
    scale_color_manual(name="",values=c("Model fitted"="red","Nonparametric"="black")) +
    theme_minimal()
  
}

GoF_BPOT <- function(EVmodel,dat){
  # goodness-of-fit test for bivariate POT models using Cramer-von Mises statistic
  n <- nrow(dat)
  dat <- dat %>% filter(prox > EVmodel$threshold[1], Speed >EVmodel$threshold[2])
  U1 <- mtransform.GPMk2(dat[,1],p=EVmodel$estimate[1:2],thres=EVmodel$threshold[1],
                         eta=EVmodel$nat[1]/EVmodel$n,margin="uniform")
  U2 <- mtransform.GPMk2(dat[,2],p=EVmodel$estimate[3:4],thres=EVmodel$threshold[2],
                         eta=EVmodel$nat[2]/EVmodel$n,margin="uniform")
  model <- switch(EVmodel$model,
                  log = gumbelCopula(dim=2,param = 1/EVmodel$estimate[5]),
                  hr = huslerReissCopula(param = EVmodel$estimate[5]))
  Cop <- fitCopula(copula = model,data = cbind(U1, U2),method = "mpl")
  Result <- gofEVCopula(copula = Cop@copula, method = "itau",x = cbind(U1,U2), N = 1000, estimator = "CFG")
  return(Result)
}  


# CBPOT -------------------------------------------------------------------

Kc_test <- function(dat,Copula,n=1e4,m=200,p.val=0.05){
  Cn <- empCopula(X = dat, ties.method="random")
  Sim <- rCopula(n, copula = Cn)
  Model <- rCopula(n,copula = Copula)
  
  temp1 <- pCopula(Sim,copula = Cn)
  temp2 <- pCopula(Model, copula = Copula)
  
  Kc_emp <- sapply(seq(0,1,1/m),function(k) mean(temp1 < k))
  Kc_model <- sapply(seq(0,1,1/m),function(k) mean(temp2 < k))
  
  result <-Desc::CCC(Kc_emp, Kc_model, ci="z-transform",conf.level = 1-p.val)
  return(result$rho.c)
  
}

Kc_test_two <- function(dat, Copula1,Copula2, n=1e4,m=200,p.val=0.05){
  CI1 <- Kc_test(dat,Copula1,n,m,p.val=p.val)
  CI2 <- Kc_test(dat,Copula2,n,m,p.val=p.val)
  CI1$upr.ci
  CI2$lwr.ci
  ifelse(CI1$upr.ci>CI2$lwr.ci,
         paste("No siginificant difference between Copula1 and Copula2 fits at level",p.val),
         paste("Copula2 is significantly better than Copula1 at level",p.val))
  
}
