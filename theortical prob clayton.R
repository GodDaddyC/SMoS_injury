copula_params <- c(-0.5,0.5,1,2,3)
if (!all(file.exists("data/theortical_clayton_CN.csv"),file.exists("data/theortical_clayton_SE.csv"))){  
  
  
  # Create named models list
  copula_models.CN <- setNames(
    lapply(copula_params, function(p) {
      Q.distr.param(claytonCopula(param = p), mar1 = POT.1, mar2 = Conseq.1, type = 4)
    }),
    paste0("r=", copula_params)
  )
  
  copula_models.SE <- setNames(
    lapply(copula_params, function(p) {
      Q.distr.param(claytonCopula(param = p), mar1 = POT.2, mar2 = Conseq.2, type = 4)
    }),
    paste0("r=", copula_params)
  )
  
  # Create the base dataframe
  theortical_density.clayton.CN <- data.frame(speed = sq.2) %>%
    mutate(
      original = sapply(speed, dgamma, shape = Conseq.1$estimate[1], rate = Conseq.1$estimate[2]),
      FromData = sapply(speed, cQ.bivariate, x = x0.1, PX = Qcrash.1, model = CM.1) /
        normalize_cQ.bivariate(x0.1, Qcrash.1, CM.1)
    )
  
  theortical_density.clayton.SE <- data.frame(speed = sq.2) %>%
    mutate(
      original = sapply(speed, dgamma, shape = Conseq.2$estimate[1], rate = Conseq.2$estimate[2]),
      FromData = sapply(speed, cQ.bivariate, x = x0.2, PX = Qcrash.2, model = CM.2) /
        normalize_cQ.bivariate(x0.2, Qcrash.2, CM.2)
    )
  
  # Add columns for each copula parameter dynamically
  copula_cols.CN <- map_dfc(names(copula_models.CN), function(name) {
    model <- copula_models.CN[[name]]
    col <- sapply(sq.2, cQ.bivariate, x = x0.1, PX = Qcrash.1, model = model) /
      normalize_cQ.bivariate(x0.1, Qcrash.1, model)
    tibble(!!name := col)
  })
  
  copula_cols.SE <- map_dfc(names(copula_models.SE), function(name) {
    model <- copula_models.SE[[name]]
    col <- sapply(sq.2, cQ.bivariate, x = x0.2, PX = Qcrash.2, model = model) /
      normalize_cQ.bivariate(x0.2, Qcrash.2, model)
    tibble(!!name := col)
  })
  
  
  # Combine the base dataframe and the copula columns
  theortical_density.clayton.CN <- bind_cols(theortical_density.clayton.CN, copula_cols.CN) %>%
    na.omit()
  theortical_density.clayton.SE <- bind_cols(theortical_density.clayton.SE, copula_cols.SE) %>%
    na.omit()
  
  write.table(theortical_density.clayton.CN,file = "data/theortical_clayton_CN.csv",sep=",",row.names = TRUE)
  write.table(theortical_density.clayton.SE,file = "data/theortical_clayton_SE.csv",sep=",",row.names = TRUE)
}

theortical_density.clayton.CN <- read.csv("data/theortical_clayton_CN.csv",sep=",")
theortical_density.clayton.SE <- read.csv("data/theortical_clayton_SE.csv",sep=",")

theortical_plot.clayton.CN <- theortical_density.clayton.CN %>%
  pivot_longer(cols = -speed, names_to = "Dependence", values_to = "Density")
theortical_plot.clayton.SE <- theortical_density.clayton.SE %>%
  pivot_longer(cols = -speed, names_to = "Dependence", values_to = "Density")
for (i in 1:length(copula_params)){
  cat(sprintf("logsitic dependence (%f): the CN injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(plot.dfQ.1,copula_models.CN[[i]],Qcrash.1,PIS0,x0.1),"\n")
  cat(sprintf("logsitic dependence (%f): the SE injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(plot.dfQ.2,copula_models.SE[[i]],Qcrash.2,PIS0,x0.2),"\n")
}
