copula_params <- c(1.1,1.3,1.8)

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

if (!file.exists("data/theortical_gumbel.csv")){  
  
  
  # Create the base dataframe
  theortical_density.gumbel<- data.frame(speed =gumbel.result$plot.dfQ.1$speed)%>%
    mutate(
      origin_CN = gumbel.result$plot.dfQ.1$ConditionalD,
      origin_SE = gumbel.result$plot.dfQ.2$ConditionalD
    )
  
  
  # Add columns for each copula parameter
  copula_cols <- mapply(function(nameCN,nameSE) {
    
    col.CN <- create_plot.dfQ(dat=s2.un,x=x0.1.un,model=copula_models.CN[[nameCN]],PX=Qcrash.1,P2=Conseq.1)$ConditionalD
    col.SE <- create_plot.dfQ(dat=s2.un,x=x0.2.un,model=copula_models.SE[[nameSE]],PX=Qcrash.2,P2=Conseq.2)$ConditionalD
    nameCN <- paste0("CN_",nameCN)
    nameSE <- paste0("SE_",nameSE)
    tibble(!!nameCN := col.CN) %>% cbind(tibble(!!nameSE := col.SE) )
  }, names(copula_models.CN), names(copula_models.SE), SIMPLIFY = FALSE)
  
  
  # Combine the base dataframe and the copula columns
  theortical_density.gumbel<- bind_cols(theortical_density.gumbel, copula_cols) %>%
    na.omit()
  
  write.table(theortical_density.gumbel,file = "data/theortical_gumbel.csv",sep=",",row.names = TRUE)
}


theortical_density.gumbel <- read.csv("data/theortical_gumbel.csv",sep=",")
theortical_density.gumbel.CN <- theortical_density.gumbel[,c("speed","origin_CN",
                                                             "CN_r.1.1","CN_r.1.3", "CN_r.1.8")]
colnames(theortical_density.gumbel.CN) <- c("speed","origin_CN",
                                            "r=1.1","r=1.3","r=1.8")
theortical_density.gumbel.SE<- theortical_density.gumbel[,c("speed","origin_SE",
                                                            "SE_r.1.1","SE_r.1.3","SE_r.1.8")]
colnames(theortical_density.gumbel.SE) <- c("speed","origin_SE",
                                            "r=1.1","r=1.3","r=1.8")
theortical_plot.gumbel.CN <- theortical_density.gumbel.CN %>%
  pivot_longer(cols = -speed, names_to = "Dependence strength", values_to = "Density")

theortical_plot.gumbel.SE <- theortical_density.gumbel.SE %>%
  pivot_longer(cols = -speed, names_to = "Dependence strength", values_to = "Density")

ggplot(theortical_plot.gumbel.CN, aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    title = "Theortical density of impact speed:Gumbel copula CN",
    x = "Speed (km/h)",
    y = "f(y|TTC<0)",
    colour = "Dependence strength"
  ) +
  theme_minimal()

ggplot(theortical_plot.gumbel.SE, aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    title = "Theortical density of impact speed: Gumbel copula SE",
    x = "Speed (km/h)",
    y = "f(y|TTC<0)",
    colour = "Dependence strength"
  ) +
  theme_minimal()


for (i in 1:length(copula_params)){
  cat(sprintf("logsitic dependence (%f): the CN injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(copula_models.CN[[i]],PIS0,x0.1.un,Qcrash.1,Conseq.1),"\n")
  cat(sprintf("logsitic dependence (%f): the SE injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(copula_models.SE[[i]],PIS0,x0.2.un,Qcrash.2,Conseq.2),"\n")
  cat(sprintf("logsitic dependence (%f): the CN injury probability (with age) is",copula_params[i]),
      Injury.from_cQ_bivariate_E(copula_models.CN[[i]],PIS1,x0.1.un,Qcrash.1,Conseq.1),"\n")
  cat(sprintf("logsitic dependence (%f): the SE injury probability (with age) is",copula_params[i]),
      Injury.from_cQ_bivariate_E(copula_models.SE[[i]],PIS1,x0.2.un,Qcrash.2,Conseq.2),"\n")
}
