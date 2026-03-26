library(reshape2)
age_mean_list <- seq(20,60,by=1)
copula_params <- c(1.1,1.3,1.8)

copula_models.SE <- setNames(
  lapply(copula_params, function(p) {
    gumbelCopula(param = p)
  }),
  paste0("r=", copula_params)
)

CBPOT_result <- data.frame(age_mean = age_mean_list)

for (i in 1:length(copula_params)){
  CBPOT_result <- CBPOT_result %>% cbind(sapply(1:length(age_mean_list),function (k) {
    Injury.from_cQ_bivariate_E(copula_models.SE[[i]],PIS1,x0.2.un,Qcrash.2,Conseq.2,LB=0.5,
                             age.mean = age_mean_list[k])}) )
}
colnames(CBPOT_result) <- c("age_mean", paste0("r=", round(1/copula_params,2)))
CBPOT_result$origin <- sapply(1:length(age_mean_list),function (k) {
      Injury.from_cQ_bivariate_E(gumbel.result$Cop.2@copula,PIS1,x0.2.un,
          Qcrash.2,Conseq.2,LB=0.5,age.mean=age_mean_list[k],)})

CBPOT_result_long <- CBPOT_result %>% 
  pivot_longer(cols = -age_mean, names_to = "Dependence strength", values_to = "Injury")

BPOT_param <- round(1/copula_params,2)
BPOT_result <- data.frame(age_mean = age_mean_list)
for (i in 1:length(BPOT_param)){
  BPOT_temp <- logistic.result$M.2
  BPOT_temp$estimate[5] <- BPOT_param[i]
  BPOT_result <- BPOT_result %>% cbind(sapply(1:length(age_mean_list),function (k) {
    Injury.from_c_bivariate_E(logistic.result$plot.df.2$speed,EVmodel=BPOT_temp,
      severity = PIS1,x0=x0.2,PX=logistic.result$Pcrash.2,age.mean = age_mean_list[k])}) )
}
colnames(BPOT_result) <- c("age_mean", paste0("r=", BPOT_param))
BPOT_result$origin <- sapply(1:length(age_mean_list),function (k){
  Injury.from_c_bivariate_E(logistic.result$plot.df.2$speed,EVmodel=logistic.result$M.2,
                            severity = PIS1,x0=x0.2,PX=logistic.result$Pcrash.2,age.mean = age_mean_list[k])})
BPOT_result_long <- BPOT_result %>% 
  pivot_longer(cols = -age_mean, names_to = "Dependence strength", values_to = "Injury")

ggplot(BPOT_result_long, aes(x = age_mean, y = Injury, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    x = "age",
    y = "AIS3+ probability",
    colour = "Dependence strength"
  ) +
  theme_minimal()

ggplot(CBPOT_result_long, aes(x = age_mean, y = Injury, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    x = "age",
    y = "AIS3+ probability",
    colour = "Dependence strength"
  ) +
  theme_minimal()


injury.df <- data.frame(speed = seq(0,80,0.5)) %>%
  mutate(InjuryP20 = sapply(speed, PIS1,age=20),
         InjuryP60 = sapply(speed, PIS1,age=60))
x_int_20 <- injury.df$speed[which.min(abs(injury.df$InjuryP20 - 0.5))]
x_int_60 <- injury.df$speed[which.min(abs(injury.df$InjuryP60 - 0.5))]

# use the peak of the crash severity distribution at most severe case
injury.df_age <- data.frame(age = age_mean_list) %>%
  mutate(InjuryP20 = sapply(age, PIS1,speed=20),
         InjuryP45 = sapply(age, PIS1,age=60))
x_int_20 <- injury.df$speed[which.min(abs(injury.df$InjuryP20 - 0.5))]
x_int_60 <- injury.df$speed[which.min(abs(injury.df$InjuryP60 - 0.5))]

figure2a <- ggplot(injury.df) + 
  geom_line(aes(x = speed, y = InjuryP20), color = "blue", linewidth = 1) +
  geom_line(aes(x = speed, y = InjuryP60), color = "red", linewidth = 1) +
  # 1. Horizontal and Vertical lines at p = 0.5
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(x_int_20, x_int_60), linetype = "dotted", color = "gray40") +
  annotate("text", x = x_int_20, y = 0.05, label = paste0(round(x_int_20, 1)), color = "blue", angle = 90, vjust = -0.5) +
  annotate("text", x = x_int_60, y = 0.05, label = paste0(round(x_int_60, 1)), color = "red", angle = 90, vjust = -0.5) +
  
  # 2. Bidirectional Arrow (at a representative speed, e.g., 45 km/h)
  geom_segment(aes(x = x_int_20, y = 0.5, xend = x_int_60, yend = 0.5), 
               arrow = arrow(ends = "both", length = unit(0.2, "cm")), 
               color = "black") +
  annotate("text", x = 30, y = 0.6, 
           label = "Sensitivity Range", hjust = 0, size = 2.5) +
  labs(
    x = "Speed (km/h)",
    y = "AIS3+ probability",
    title = "Injury risks for Different Ages"
  ) +
  theme_minimal() +
  annotate("text", x = 60, y = 0.55, label = "Age 20", color = "blue") +
  annotate("text", x = 60, y = 1, label = "Age 60", color = "red")

figure2b <- ggplot(theortical_plot.logistic.SE, aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    title = "Crash severities: BPOT",
    x = "Speed (km/h)",
    y = "f(y|TTC<0)",
    colour = "Dependence strength"
  ) +
  theme_minimal()

figure2c <- ggplot(theortical_plot.gumbel.SE, aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() + 
  labs(
    title = "Cheash severities: CBPOT",
    x = "Speed (km/h)",
    y = "f(y|TTC<0)",
    colour = "Dependence strength"
  ) +
  theme_minimal()
figure2 <- ggarrange(figure2a, figure2b, figure2c, nrow = 1, labels = c("(a)", "(b)", "(c)"))
ggsave(figure2, filename = "plots/RSSfigure2.png", width = 12, height = 4, dpi = 300)


summary_df <- bind_rows(
  BPOT_result_long %>% mutate(Method = "BPOT"),
  CBPOT_result_long %>% mutate(Method = "CBPOT") )

figure3 <- ggplot(summary_df, aes(x = age_mean, y = Injury, color = `Dependence strength`, linetype = Method)) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Sensitivity of injury risk",
    x = "age",
    y = "AIS3+ probability"
  ) +
  scale_color_brewer(palette = "Set1")
ggsave(figure3, filename = "plots/RSSfigure3.png", dpi = 300)
