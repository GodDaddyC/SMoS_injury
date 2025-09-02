copula_params <- c(1.1,1.3,1.5,1.8,2)
if (!all(file.exists("data/theortical_gumbel_CN.csv"),file.exists("data/theortical_gumbel_SE.csv"))){  
  
  
  # Create named models list
  copula_models.CN <- setNames(
    lapply(copula_params, function(p) {
      gumbelCopula(param = p)
    }),
    paste0("r=", copula_params)
  )
  
  copula_models.SE <- setNames(
    lapply(copula_params, function(p) {
      gumbelCopula(param = p)
    }),
    paste0("r=", copula_params)
  )
  
  # Create the base dataframe
  theortical_density.gumbel.CN <- data.frame(speed = sq.2) %>%
    mutate(
      original = sapply(speed, dgamma, shape = Conseq.1$estimate[1], rate = Conseq.1$estimate[2]),
      FromData = plot.dfQ.1$ConditionalD
    )
  
  theortical_density.gumbel.SE <- data.frame(speed = sq.2) %>%
    mutate(
      original = sapply(speed, dgamma, shape = Conseq.2$estimate[1], rate = Conseq.2$estimate[2]),
      FromData = plot.dfQ.2$ConditionalD
    )
  
  # Add columns for each copula parameter dynamically
  copula_cols.CN <- map_dfc(names(copula_models.CN), function(name) {
    model <- copula_models.CN[[name]]
    temp <- create_plot.dfQ(s1.un,x0.1.un,model,Qcrash.1*0.3,Conseq.2)
    col <- temp$ConditionalD
    tibble(!!name := col)
  })
  
  copula_cols.SE <- map_dfc(names(copula_models.SE), function(name) {
    model <- copula_models.SE[[name]]
    temp <- create_plot.dfQ(s2.un,x0.2.un,model,Qcrash.2*0.2,Conseq.2)
    col <- temp$ConditionalD
    tibble(!!name := col)
  })
  
  
  # Combine the base dataframe and the copula columns
  theortical_density.gumbel.CN <- bind_cols(theortical_density.gumbel.CN, copula_cols.CN) %>%
    na.omit()
  theortical_density.gumbel.SE <- bind_cols(theortical_density.gumbel.SE, copula_cols.SE) %>%
    na.omit()
  
  write.table(theortical_density.gumbel.CN,file = "data/theortical_gumbel_CN.csv",sep=",",row.names = TRUE)
  write.table(theortical_density.gumbel.SE,file = "data/theortical_gumbel_SE.csv",sep=",",row.names = TRUE)
}

copula_models.CN <- setNames(
  lapply(copula_params, function(p) {
    gumbelCopula(param = p)
  }),
  paste0("r=", copula_params)
)

copula_models.SE <- setNames(
  lapply(copula_params, function(p) {
    gumbelCopula(param = p)
  }),
  paste0("r=", copula_params)
)
theortical_density.gumbel.CN <- read.csv("data/theortical_gumbel_CN.csv",sep=",")
theortical_density.gumbel.SE <- read.csv("data/theortical_gumbel_SE.csv",sep=",")

theortical_plot.gumbel.CN <- theortical_density.gumbel.CN %>%
  pivot_longer(cols = -speed, names_to = "Dependence", values_to = "Density")
theortical_plot.gumbel.SE <- theortical_density.gumbel.SE %>%
  pivot_longer(cols = -speed, names_to = "Dependence", values_to = "Density")
for (i in 1:length(copula_params)){
  cat(sprintf("logsitic dependence (%f): the CN injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(copula_models.CN[[i]],PIS0,x0.1.un,Qcrash.1*0.3,Conseq.1),"\n")
  cat(sprintf("logsitic dependence (%f): the SE injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(copula_models.SE[[i]],PIS0,x0.2.un,Qcrash.2*0.2,Conseq.2),"\n")
}
