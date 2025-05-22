

combined_data <- bind_rows(
  CH_Dat[,c("TTC","PET","MD")] %>% mutate(from = "CH"),
  -SE_Dat[,c("TTC","PET","MD")] %>% mutate(from = "SE")
)

# Histogram + density overlay
g1 <- ggplot(combined_data, aes(x = TTC, fill = from, color = from)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.4, bins = 30) +
  geom_density(linewidth = 1) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  labs(title = "Comparison of TTCs",
       x = "Value",
       y = "Density") +
  theme_minimal()

g2 <- ggplot(combined_data, aes(x = PET, fill = from, color = from)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.4, bins = 30) +
  geom_density(linewidth = 1) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  labs(title = "Comparison of PET",
       x = "Value",
       y = "Density") +
  theme_minimal()

g3 <- ggplot(combined_data, aes(x = MD, fill = from, color = from)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.4, bins = 30) +
  geom_density(linewidth = 1) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  labs(title = "Comparison of MD",
       x = "Value",
       y = "Density") +
  theme_minimal()

ggarrange(g1,g2,g3,ncol=3)