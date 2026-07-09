# auxillary functions 


# BPOT --------------------------------------------------------------------

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

create_plot.df <- function(dat,x,model,PX){
  # create a data frame for plotting the conditional density of speed at X
  # model is a fbvpot object
  df <-switch(model$model,
              log = data.frame(speed= dat,
                               JointP = sapply(dat,FUN = pbTvevd,q1=x,model=model$model,dep=model$estimate[5],
                                               thres=model$threshold,
                                               eta=model$nat[1:2]/model$n,
                                               mar1=c(model$estimate[1],model$estimate[2]),
                                               mar2=c(model$estimate[3],model$estimate[4]),tail.type=4)) %>%
                mutate(ConditionP = JointP/PX) %>%na.omit() %>%
                mutate(ConditionalD = sapply(speed, 
                                             function(k) c.bivariate(y = k, x = x, model=model$model,dep=model$estimate[5],
                                                                     thres=model$threshold,
                                                                     eta=model$nat[1:2]/model$n,
                                                                     mar1=c(model$estimate[1],model$estimate[2]),
                                                                     mar2=c(model$estimate[3],model$estimate[4]))
                )/PX /normalize_c.bivariate(x,model) * model$nat[3]/model$n),
              hr = data.frame(speed= dat,
                              JointP = sapply(dat,FUN = pbTvevd,q1=x,model=model$model,dep=model$estimate[5],
                                              thres=model$threshold,
                                              eta=model$nat[1:2]/model$n,
                                              mar1=c(model$estimate[1],model$estimate[2]),
                                              mar2=c(model$estimate[3],model$estimate[4]),tail.type=4)) %>%
                mutate(ConditionP = JointP/PX) %>% na.omit() %>%
                mutate(ConditionalD = sapply(speed, 
                                             function(k) c.bivariate(y = k, x = x, model=model$model,dep=model$estimate[5],
                                                                     thres=model$threshold,
                                                                     eta=model$nat[1:2]/model$n,
                                                                     mar1=c(model$estimate[1],model$estimate[2]),
                                                                     mar2=c(model$estimate[3],model$estimate[4]))
                )/PX/normalize_c.bivariate(x,model) * model$nat[3]/model$n),
              ct = data.frame(speed= dat,
                              JointP = sapply(dat,FUN = pbTvevd,q1=x,model=model$model,alpha=model$estimate[5],
                                              beta=model$estimate[6],thres=model$threshold,
                                              eta=model$nat[1:2]/model$n,
                                              mar1=c(model$estimate[1],model$estimate[2]),
                                              mar2=c(model$estimate[3],model$estimate[4]),tail.type=4)) %>%
                mutate(ConditionP = JointP) %>% na.omit() %>%
                mutate(ConditionalD = sapply(speed, 
                                             function(k) c.bivariate(y = k, x = x, model=model$model,alpha=model$estimate[5],
                                                                     beta=model$estimate[6],thres=model$threshold,
                                                                     eta=model$nat[1:2]/model$n,
                                                                     mar1=c(model$estimate[1],model$estimate[2]),
                                                                     mar2=c(model$estimate[3],model$estimate[4]))
                )/ PX/normalize_c.bivariate(x,model) * model$nat[3]/model$n)
  )
  return(df)
}

create_plot.df.np <- function(dat,x,PX,mar1,mar2,thres,eta,Ahat){
  df <- data.frame(speed= dat,
                   JointP = sapply(dat,FUN = pbTvNonpar,q1=x,Ahat=Ahat,
                                   thres=thres,eta=eta,
                                   mar1=mar1, mar2=mar2,tail.type=4)) %>%
    mutate(ConditionP = JointP/PX) %>%
    na.omit() %>%
    mutate(ConditionalD = sapply(speed, function(k){
      c.bivariate_np(y = k, x = x,Ahat=Ahat,thres=thres,eta=eta,mar1=mar1,mar2=mar2)})/PX) %>%
    mutate(ConditionalD = ifelse(ConditionalD < 0, 0, ConditionalD))
}


summarise_BPOT <- function(result, x0 = c(0,0), severity = NULL,
                           severity.age = NULL, age.mean = 60) {
  # result     : output of runExecutionBPOT
  # x0         : numeric(2), crash TTC thresholds c(x0.CN, x0.SE); required for
  #              severe crash and injury probability sections
  # severity   : severity function for injury prob from impact speed alone (e.g. PIS0)
  # severity.age: severity function including age covariate (e.g. PIS1)
  # age.mean   : mean pedestrian age used in severity.age computation (default 60)

  M1  <- result$M.1
  M2  <- result$M.2
  hdr <- paste0(strrep("=", 55), "\n")
  sec <- paste0(strrep("-", 55), "\n")

  fmt <- function(x, digits = 4) round(x, digits)

  cat(hdr)
  cat("  BPOT Model Summary  (model:", M1$model, ")\n")
  cat(hdr)

  # --- Thresholds ----------------------------------------------------------
  cat("\nThresholds\n", sec)
  cat(sprintf("  %-6s  TTC = %-8.4f  Speed = %.4f\n", "CN:", M1$threshold[1], M1$threshold[2]))
  cat(sprintf("  %-6s  TTC = %-8.4f  Speed = %.4f\n", "SE:", M2$threshold[1], M2$threshold[2]))

  # --- Parameter estimates -------------------------------------------------
  pnames <- c("sigma1", "xi1", "sigma2", "xi2", "dep")
  if (M1$model == "ct") pnames <- c(pnames, "beta")

  cat("\nParameter Estimates\n", sec)
  cat("  CN:", paste(pnames, fmt(M1$estimate), sep = " = ", collapse = "  "), "\n")
  cat("  95% CI (CN):\n"); print(confint(M1))
  cat("  SE:", paste(pnames, fmt(M2$estimate), sep = " = ", collapse = "  "), "\n")
  cat("  95% CI (SE):\n"); print(confint(M2))

  # --- Model fit -----------------------------------------------------------
  cat("\nModel Fit\n", sec)
  cat(sprintf("  Deviance  CN = %.4f\n", deviance(M1)))
  cat(sprintf("  Deviance  SE = %.4f\n", deviance(M2)))

  # --- Crash probability ---------------------------------------------------
  cat("\nCrash Probability  P(TTC < u)\n", sec)
  cat(sprintf("  CN = %.6f\n", result$Pcrash.1))
  cat(sprintf("  SE = %.6f\n", result$Pcrash.2))

  # --- Severe crash and injury (require x0) --------------------------------
  if (!is.null(x0)) {
    cat("\nSevere Crash Probability  P(TTC < u, Speed > 40 km/h)\n", sec)
    p.sev.1 <- pbTvevd(q1 = x0[1], q2 = 40, dep = M1$estimate[5],
                       thres = M1$threshold, model = M1$model,
                       eta   = M1$nat[1] / M1$n,
                       mar1  = M1$estimate[1:2], mar2 = M1$estimate[3:4], tail.type = 2)
    p.sev.2 <- pbTvevd(q1 = x0[2], q2 = 40, dep = M2$estimate[5],
                       thres = M2$threshold, model = M2$model,
                       eta   = M2$nat[1] / M2$n,
                       mar1  = M2$estimate[1:2], mar2 = M2$estimate[3:4], tail.type = 2)
    cat(sprintf("  CN = %.6f\n", p.sev.1))
    cat(sprintf("  SE = %.6f\n", p.sev.2))

    if (!is.null(severity)) {
      cat("\nInjury Probability  (impact speed only)\n", sec)
      ip.1 <- Injury.from_c_bivariate(result$plot.df.1$speed, EVmodel = M1,
                                       severity = severity, x0 = x0[1], PX = result$Pcrash.1)
      ip.2 <- Injury.from_c_bivariate(result$plot.df.2$speed, EVmodel = M2,
                                       severity = severity, x0 = x0[2], PX = result$Pcrash.2)
      cat(sprintf("  CN = %.6f\n", ip.1))
      cat(sprintf("  SE = %.6f\n", ip.2))
    }

    if (!is.null(severity.age)) {
      cat(sprintf("\nInjury Probability  (impact speed + age = %g)\n", age.mean), sec)
      ip.1a <- Injury.from_c_bivariate_E(result$plot.df.1$speed, EVmodel = M1,
                                          severity = severity.age, x0 = x0[1],
                                          PX = result$Pcrash.1, age.mean = age.mean)
      ip.2a <- Injury.from_c_bivariate_E(result$plot.df.2$speed, EVmodel = M2,
                                          severity = severity.age, x0 = x0[2],
                                          PX = result$Pcrash.2, age.mean = age.mean)
      cat(sprintf("  CN = %.6f\n", ip.1a))
      cat(sprintf("  SE = %.6f\n", ip.2a))
    }
  }

  cat(hdr)
  invisible(result)
}


# CBPOT -------------------------------------------------------------------

create_plot.dfQ <- function(dat,x,model,PX,P2,Pu=1){
  # create a data frame for plotting the conditional density of speed at X
  # model is a mvdc object
  df <- data.frame(speed.u= dat,speed = qgamma(dat,shape = P2$estimate[1],rate = P2$estimate[2]),
                   JointP = sapply(dat,function(y){pCopula(c(x,y),model)})) %>%
    mutate(ConditionP = JointP/PX) %>% 
    mutate(ConditionalD = sapply(speed.u,cQ.bivariate,x=x,model = model)/PX * Pu *
             dgamma(speed,shape = P2$estimate[1],rate = P2$estimate[2])) %>% na.omit()
  return(df)
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
