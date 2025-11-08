copula_params <- c(-0.5,-0.2,1,2,4)

copula_models.CN <- setNames(
  lapply(copula_params, function(p) {
    claytonCopula(param = p)
  }),
  paste0("r=", copula_params)
)

copula_models.SE <- setNames(
  lapply(copula_params, function(p) {
    claytonCopula(param = p)
  }),
  paste0("r=", copula_params)
)

if (!file.exists("data/theortical_clayton.csv")){  
  
  theortical_density.clayton <- data.frame(speed =clayton.result$plot.dfQ.1$speed)%>%
    mutate(
      origin_CN = clayton.result$plot.dfQ.1$ConditionalD,
      origin_SE = clayton.result$plot.dfQ.2$ConditionalD
    )
  
  
  # Add columns for each copula parameter
  copula_cols <- mapply(function(nameCN,nameSE) {
    col.CN <- create_plot.dfQ(dat=s2.un,x=x0.1.un,model=copula_models.CN[[nameCN]],PX=Qcrash.1,P2=Conseq.1)$ConditionalD
    col.SE <- create_plot.dfQ(dat=s2.un,x=x0.2.un,model=copula_models.SE[[nameSE]],PX=Qcrash.2,P2=Conseq.2)$ConditionalD
    nameCN <- paste0("CN_",nameCN)
    nameSE <- paste0("SE_",nameSE)
    tibble(!!nameCN := col.CN) %>% cbind(tibble(!!nameSE := col.SE) )
  }, names(copula_models.CN), names(copula_models.SE), SIMPLIFY = FALSE)
  
  theortical_density.clayton<- bind_cols(theortical_density.clayton, copula_cols) %>%
    na.omit()
  
  write.table(theortical_density.clayton,file = "data/theortical_clayton.csv",sep=",",row.names = TRUE)
}

theortical_density.clayton <- read.csv("data/theortical_clayton.csv",sep=",")
theortical_density.clayton.CN <- theortical_density.clayton[,c("speed","origin_CN",
                                          "CN_r..0.5","CN_r..0.2", "CN_r.1","CN_r.2","CN_r.4")]
colnames(theortical_density.clayton.CN) <- c("speed","origin_CN",
                                           "r=-0.5","r=-0.2","r=1","r=2","r=4")
theortical_density.clayton.SE<- theortical_density.clayton[,c("speed","origin_SE",
                                          "SE_r..0.5","SE_r..0.2", "SE_r.1","SE_r.2","SE_r.4")]
colnames(theortical_density.clayton.SE) <- c("speed","origin_SE",
                                             "r=-0.5","r=-0.2","r=1","r=2","r=4")


theortical_plot.clayton.CN <- theortical_density.clayton.CN %>%
  pivot_longer(cols = -speed, names_to = "Dependence strength", values_to = "Density")
theortical_plot.clayton.SE <- theortical_density.clayton.SE %>%
  pivot_longer(cols = -speed, names_to = "Dependence strength", values_to = "Density")

ggplot(theortical_plot.clayton.CN, aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    title = "Theortical density of impact speed:clayton copula CN",
    x = "Speed (km/h)",
    y = "f(y|TTC<0)",
    colour = "Dependence strength"
  ) +
  theme_minimal()

ggplot(theortical_plot.clayton.SE, aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    title = "Theortical density of impact speed: clayton copula SE",
    x = "Speed (km/h)",
    y = "f(y|TTC<0)",
    colour = "Dependence strength"
  ) +
  theme_minimal()


for (i in 1:length(copula_params)){
  cat(sprintf("clayton dependence (%f): the CN injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(copula_models.CN[[i]],PIS0,x0.1.un,Qcrash.1,Conseq.1),"\n")
  cat(sprintf("clayton dependence (%f): the SE injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(copula_models.SE[[i]],PIS0,x0.2.un,Qcrash.2,Conseq.2),"\n")
}
