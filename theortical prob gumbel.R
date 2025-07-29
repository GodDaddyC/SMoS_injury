copula_params <- c(1.3,1.5, 1.8, 2.5,5)

# Create named models list
copula_models.CN <- setNames(
  lapply(copula_params, function(p) {
    Q.distr.param(gumbelCopula(param = p), mar1 = POT.1T, mar2 = Conseq.1T, type = 4)
  }),
  paste0("r=", copula_params)
)

copula_models.SE <- setNames(
  lapply(copula_params, function(p) {
    Q.distr.param(gumbelCopula(param = p), mar1 = POT.2T, mar2 = Conseq.2T, type = 4)
  }),
  paste0("r=", copula_params)
)

# Create the base dataframe
theortical_density.gumbel.CN <- data.frame(speed = sq.2T) %>%
  mutate(
    original = sapply(speed, dgamma, shape = Conseq.1T$estimate[1], rate = Conseq.1T$estimate[2]),
    FromData = sapply(speed, cQ.bivariate, x = 0, PX = Qcrash.1T, model = CM.1T) /
      normalize_cQ.bivariate(0, Qcrash.1T, CM.1T)
  )

theortical_density.gumbel.SE <- data.frame(speed = sq.2T) %>%
  mutate(
    original = sapply(speed, dgamma, shape = Conseq.2T$estimate[1], rate = Conseq.2T$estimate[2]),
    FromData = sapply(speed, cQ.bivariate, x = 0, PX = Qcrash.2T, model = CM.2T) /
      normalize_cQ.bivariate(0, Qcrash.2T, CM.2T)
  )

# Add columns for each copula parameter dynamically
copula_cols.CN <- map_dfc(names(copula_models.CN), function(name) {
  model <- copula_models.CN[[name]]
  col <- sapply(sq.2T, cQ.bivariate, x = 0, PX = Qcrash.1T, model = model) /
    normalize_cQ.bivariate(0, Qcrash.1T, model)
  tibble(!!name := col)
})

copula_cols.SE <- map_dfc(names(copula_models.SE), function(name) {
  model <- copula_models.SE[[name]]
  col <- sapply(sq.2T, cQ.bivariate, x = 0, PX = Qcrash.2T, model = model) /
    normalize_cQ.bivariate(0, Qcrash.2T, model)
  tibble(!!name := col)
})


# Combine the base dataframe and the copula columns
theortical_density.gumbel.CN <- bind_cols(theortical_density.gumbel.CN, copula_cols.CN) %>%
  na.omit()
theortical_density.gumbel.SE <- bind_cols(theortical_density.gumbel.SE, copula_cols.SE) %>%
  na.omit()

theortical_plot.gumbel.CN <- theortical_density.gumbel.CN %>%
  pivot_longer(cols = -speed, names_to = "Dependence", values_to = "Density")
theortical_plot.gumbel.SE <- theortical_density.gumbel.SE %>%
  pivot_longer(cols = -speed, names_to = "Dependence", values_to = "Density")



for (i in 1:length(copula_params)){
  cat(sprintf("logsitic dependence (%f): the CN injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(plot.dfQ.1T,copula_models.CN[[i]],Qcrash.1T,PIS0),"\n")
  cat(sprintf("logsitic dependence (%f): the SE injury probability is",copula_params[i]),
      Injury.from_cQ_bivariate(plot.dfQ.2T,copula_models.SE[[i]],Qcrash.2T,PIS0),"\n")
}
