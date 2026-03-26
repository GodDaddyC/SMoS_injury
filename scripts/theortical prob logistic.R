BPOT_params <- round(1/c(1.1,1.3,1.8),2)

BPOT_models.CN <- setNames(
  lapply(BPOT_params, function(p) {
    BPOT_temp <- logistic.result$M.1
    BPOT_temp$estimate[5] <- p
    return(BPOT_temp)
  }),
  paste0("r=", BPOT_params)
)

BPOT_models.SE <- setNames(
  lapply(BPOT_params, function(p) {
    POT_temp <- logistic.result$M.1
    BPOT_temp$estimate[5] <- p
    return(BPOT_temp)
  }),
  paste0("r=", BPOT_params)
)
  

theortical_density.logistic<- data.frame(speed =logistic.result$plot.df.1$speed)%>%
  mutate(
    origin_CN = logistic.result$plot.df.1$ConditionalD,
    origin_SE = logistic.result$plot.df.2$ConditionalD
  )


# Add columns for each copula parameter
BPOT_cols <- mapply(function(nameCN,nameSE) {
  col.CN <- create_plot.df(dat=logistic.result$plot.df.1$speed,x=x0.1,model=BPOT_models.CN[[nameCN]],
                           PX=logistic.result$Pcrash.1)$ConditionalD
  col.SE <- create_plot.df(dat=logistic.result$plot.df.2$speed,x=x0.1,model=BPOT_models.SE[[nameSE]],
                           PX=logistic.result$Pcrash.2)$ConditionalD
  nameCN <- paste0("CN_",nameCN)
  nameSE <- paste0("SE_",nameSE)
  tibble(!!nameCN := col.CN) %>% cbind(tibble(!!nameSE := col.SE) )
}, names(BPOT_models.CN), names(BPOT_models.SE), SIMPLIFY = FALSE)


# Combine the base dataframe and the copula columns
theortical_density.logistic<- bind_cols(theortical_density.logistic, BPOT_cols) %>%
  na.omit()

theortical_density.logistic.CN <- theortical_density.logistic[,c("speed","origin_CN",
                                                             "CN_r=0.77","CN_r=0.91", "CN_r=0.56")]
colnames(theortical_density.logistic.CN) <- c("speed","origin_CN",
                                            "r=0.77","r=0.91","r=0.56")
theortical_density.logistic.SE<- theortical_density.logistic[,c("speed","origin_SE",
                                                                "SE_r=0.77","SE_r=0.91", "SE_r=0.56")]
colnames(theortical_density.logistic.SE) <- c("speed","origin_SE",
                                              "r=0.77","r=0.91","r=0.56")
theortical_plot.logistic.CN <- theortical_density.logistic.CN %>%
  pivot_longer(cols = -speed, names_to = "Dependence strength", values_to = "Density")

theortical_plot.logistic.SE <- theortical_density.logistic.SE %>%
  pivot_longer(cols = -speed, names_to = "Dependence strength", values_to = "Density")

ggplot(theortical_plot.logistic.CN, aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    title = "Theortical density of impact speed:logistic BPOT CN",
    x = "Speed (km/h)",
    y = "f(y|TTC<0)",
    colour = "Dependence strength"
  ) +
  theme_minimal()

ggplot(theortical_plot.logistic.SE, aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    title = "Crash severity: BPOT",
    x = "Speed (km/h)",
    y = "f(y|TTC<0)",
    colour = "Dependence strength"
  ) +
  theme_minimal()
