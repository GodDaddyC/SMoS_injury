age_mean_list <- seq(20, 60, by = 1)
copula_params <- c(1.1, 1.3, 1.8)

copula_models.SE <- setNames(
  lapply(copula_params, function(p) { gumbelCopula(param = p) }),
  paste0("r=", copula_params)
)

cbpot_result <- data.frame(age_mean = age_mean_list)

for (i in seq_along(copula_params)) {
  cbpot_result <- cbpot_result %>%
    cbind(sapply(seq_along(age_mean_list), function(k) {
      injury_from_c_q_bivariate_e(copula_models.SE[[i]], pis1,
        x0_2_un, qcrash_2, conseq_2, lb = 0.5,
        age_mean = age_mean_list[k])
    }))
}
colnames(cbpot_result) <- c("age_mean",
  paste0("r=", round(1 / copula_params, 2)))
cbpot_result$origin <- sapply(seq_along(age_mean_list), function(k) {
  injury_from_c_q_bivariate_e(gumbel_result$Cop_2@copula, pis1,
    x0_2_un, qcrash_2, conseq_2, lb = 0.5,
    age_mean = age_mean_list[k])
})

cbpot_result_long <- cbpot_result %>%
  pivot_longer(cols = -age_mean, names_to = "Dependence strength",
               values_to = "Injury")

bpot_param <- round(1 / copula_params, 2)
bpot_result <- data.frame(age_mean = age_mean_list)
for (i in seq_along(bpot_param)) {
  bpot_temp <- logistic_result$M_2
  bpot_temp$estimate[5] <- bpot_param[i]
  bpot_result <- bpot_result %>%
    cbind(sapply(seq_along(age_mean_list), function(k) {
      injury_from_c_bivariate_e(logistic_result$plot_df_2$speed,
        ev_model = bpot_temp, severity = pis1, x0 = x0_2,
        px = logistic_result$pcrash_2, age_mean = age_mean_list[k])
    }))
}
colnames(bpot_result) <- c("age_mean", paste0("r=", bpot_param))
bpot_result$origin <- sapply(seq_along(age_mean_list), function(k) {
  injury_from_c_bivariate_e(logistic_result$plot_df_2$speed,
    ev_model = logistic_result$M_2, severity = pis1,
    x0 = x0_2, px = logistic_result$pcrash_2,
    age_mean = age_mean_list[k])
})
bpot_result_long <- bpot_result %>%
  pivot_longer(cols = -age_mean, names_to = "Dependence strength",
               values_to = "Injury")

ggplot(bpot_result_long, aes(x = age_mean, y = Injury,
       colour = `Dependence strength`)) +
  geom_line() +
  labs(x = "age", y = "AIS3+ probability",
       colour = "Dependence strength") +
  theme_minimal()

ggplot(cbpot_result_long, aes(x = age_mean, y = Injury,
       colour = `Dependence strength`)) +
  geom_line() +
  labs(x = "age", y = "AIS3+ probability",
       colour = "Dependence strength") +
  theme_minimal()


injury_df2 <- data.frame(speed = seq(0, 80, 0.5)) %>%
  mutate(InjuryP20 = sapply(speed, pis1, age = 20),
         InjuryP60 = sapply(speed, pis1, age = 60))
x_int_20 <- injury_df2$speed[which.min(abs(injury_df2$InjuryP20 - 0.5))]
x_int_60 <- injury_df2$speed[which.min(abs(injury_df2$InjuryP60 - 0.5))]

figure2a <- ggplot(injury_df2) +
  geom_line(aes(x = speed, y = InjuryP20), color = "blue", linewidth = 1) +
  geom_line(aes(x = speed, y = InjuryP60), color = "red", linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(x_int_20, x_int_60), linetype = "dotted",
             color = "gray40") +
  annotate("text", x = x_int_20, y = 0.05,
           label = paste0(round(x_int_20, 1)), color = "blue",
           angle = 90, vjust = -0.5) +
  annotate("text", x = x_int_60, y = 0.05,
           label = paste0(round(x_int_60, 1)), color = "red",
           angle = 90, vjust = -0.5) +
  geom_segment(aes(x = x_int_20, y = 0.5, xend = x_int_60, yend = 0.5),
               arrow = arrow(ends = "both", length = unit(0.2, "cm")),
               color = "black") +
  annotate("text", x = 30, y = 0.6,
           label = "Sensitivity Range", hjust = 0, size = 2.5) +
  labs(x = "Speed (km/h)", y = "AIS3+ probability",
       title = "Injury risks for Different Ages") +
  theme_minimal() +
  annotate("text", x = 60, y = 0.55, label = "Age 20", color = "blue") +
  annotate("text", x = 60, y = 1, label = "Age 60", color = "red")

figure2b <- ggplot(theoretical_plot_logistic_SE,
  aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() +
  labs(title = "Crash severities: BPOT", x = "Speed (km/h)",
       y = "f(y|TTC<0)", colour = "Dependence strength") +
  theme_minimal()

figure2c <- ggplot(theoretical_plot_gumbel_SE,
  aes(x = speed, y = Density, colour = `Dependence strength`)) +
  geom_line() +
  labs(title = "Crash severities: CBPOT", x = "Speed (km/h)",
       y = "f(y|TTC<0)", colour = "Dependence strength") +
  theme_minimal()

figure2 <- ggarrange(figure2a, figure2b, figure2c, nrow = 1,
                     labels = c("(a)", "(b)", "(c)"))
ggsave(figure2, filename = "plots/RSSfigure2.png", width = 12, height = 4,
       dpi = 300)


summary_df <- bind_rows(
  bpot_result_long %>% mutate(Method = "BPOT"),
  cbpot_result_long %>% mutate(Method = "CBPOT"))

figure3 <- ggplot(summary_df,
  aes(x = age_mean, y = Injury, color = `Dependence strength`,
      linetype = Method)) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  labs(title = "Sensitivity of injury risk", x = "age",
       y = "AIS3+ probability") +
  scale_color_brewer(palette = "Set1")
ggsave(figure3, filename = "plots/RSSfigure3.png", dpi = 300)
