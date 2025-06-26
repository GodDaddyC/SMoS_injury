

combined_data <- bind_rows(
  CN_vru[,c("TTC","PET","MD")] %>% mutate(from = "CN"),
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


# bivariate plotting

# plot CN and SE with also 
ggplot(data = EVT.n,aes(x = prox,y = conseq))+
  geom_point(size = 0.2) + geom_vline(xintercept = u1.n) + geom_vline(xintercept = 0,linetype = "dashed",color="red") + 
  geom_point(data = EVT.n[which(EVT.n$prox>u1.n),],aes(y=conseq),color = 'red',size = 1)  +
  #annotate("text", x = -17, y = 35, label = "threshold for ", color = "red") +
  #annotate("text", x = -30, y = 25, label = "threshold for near interactions")+
  xlab('negated distance (m)') + ylab('Delta-V (m/s)') + theme(axis.line = element_line(linetype = "solid"),
                                                               panel.background = element_rect(fill = NA))