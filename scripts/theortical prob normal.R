copula_params <- c(-0.6,-0.3,0.3,0.6)

copula_models.CN <- setNames(
  lapply(copula_params, function(p) {
    normalCopula(param = p,dim=2,dispstr = "ex")
  }),
  paste0("r=", copula_params)
)

copula_models.SE <- setNames(
  lapply(copula_params, function(p) {
    normalCopula(param = p,dim=2,dispstr = "ex")
  }),
  paste0("r=", copula_params)
)

if (!file.exists("data/theortical_normal.csv")){  
  
  theortical_density.normal <- data.frame(speed =normal.result$plot.dfQ.1$speed)%>%
    mutate(
      origin_CN = normal.result$plot.dfQ.1$ConditionalD,
      origin_SE = normal.result$plot.dfQ.2$ConditionalD
    )
  
  
  # Add columns for each copula parameter
  copula_cols <- mapply(function(nameCN,nameSE) {
    col.CN <- create_plot.dfQ(dat=s2.un,x=x0.1.un,model=copula_models.CN[[nameCN]],PX=Qcrash.1,P2=Conseq.1)$ConditionalD
    col.SE <- create_plot.dfQ(dat=s2.un,x=x0.2.un,model=copula_models.SE[[nameSE]],PX=Qcrash.2,P2=Conseq.2)$ConditionalD
    nameCN <- paste0("CN_",nameCN)
    nameSE <- paste0("SE_",nameSE)
    tibble(!!nameCN := col.CN) %>% cbind(tibble(!!nameSE := col.SE) )
  }, names(copula_models.CN), names(copula_models.SE), SIMPLIFY = FALSE)
  
  theortical_density.normal<- bind_cols(theortical_density.normal, copula_cols) %>%
    na.omit()
  
  write.table(theortical_density.normal,file = "data/theortical_normal.csv",sep=",",row.names = TRUE)
}

theortical_density.normal <- read.csv("data/theortical_normal.csv",sep=",")
theortical_density.normal.CN <- theortical_density.normal[,c("speed","origin_CN",
                                          "CN_r..0.6","CN_r..0.3", "CN_r.0.3","CN_r.0.6")]
colnames(theortical_density.normal.CN) <- c("speed","origin_CN",
                                           "r=-0.6","r=-0.3","r=0.3","r=0.6")
theortical_density.normal.SE<- theortical_density.normal[,c("speed","origin_SE",
                                          "SE_r..0.6","SE_r..0.3", "SE_r.0.3","SE_r.0.6")]
colnames(theortical_density.normal.SE) <- c("speed","origin_SE",
                                            "r=-0.6","r=-0.3","r=0.3","r=0.6")


theortical_plot.normal.CN <- theortical_density.normal.CN %>%
  pivot_longer(cols = -speed, names_to = "Dependence strength", values_to = "Density")
theortical_plot.normal.SE <- theortical_density.normal.SE %>%
  pivot_longer(cols = -speed, names_to = "Dependence strength", values_to = "Density")

ggplot(theortical_plot.normal.CN, aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    title = "Theortical density of impact speed:normal copula CN",
    x = "Speed (km/h)",
    y = "f(y|TTC<0)",
    colour = "Dependence strength"
  ) +
  theme_minimal()

ggplot(theortical_plot.normal.SE, aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    title = "Theortical density of impact speed: normal copula SE",
    x = "Speed (km/h)",
    y = "f(y|TTC<0)",
    colour = "Dependence strength"
  ) +
  theme_minimal()


for (i in 1:length(copula_params)){
  cat(sprintf("normal dependence (%f): the CN injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(copula_models.CN[[i]],PIS0,x0.1.un,Qcrash.1,Conseq.1),"\n")
  cat(sprintf("normal dependence (%f): the SE injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(copula_models.SE[[i]],PIS0,x0.2.un,Qcrash.2,Conseq.2),"\n")
}
