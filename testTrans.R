# test the effect of transformation on injury prob from BPOT

tt <- CN_TTC_A
tt.T <- CN_iTTC_A
t.model <- "ct"

x0 <- 0
x0.T <- 1

k1 <- bvtcplot(tt)$k
k1.T <- bvtcplot(tt.T)$k

u1 <- sort(tt$prox,decreasing = TRUE)[k1]
v1 <- sort(tt$Speed,decreasing = TRUE)[k1]

u1.T <- sort(tt.T$prox,decreasing = TRUE)[k1]
v1.T <- v1

M1 <- fbvpot(x = tt,model = t.model,threshold = c(u1,v1))
M1.T <- fbvpot(x = tt.T,model = t.model,threshold = c(u1.T,v1.T))

P1 <- pevd(x0,threshold = u1, scale = M1$estimate[1],shape = M1$estimate[2],
                 lower.tail = FALSE,type = "GP") * M1$nat[1]/M1$n

P1.T <- pevd(x0.T,threshold = u1.T, scale = M1.T$estimate[1],shape = M1.T$estimate[2],
                 lower.tail = FALSE,type = "GP") * M1.T$nat[1]/M1.T$n

ss.1T <- seq(v1,60,(60- v1)/150)

df1 <- create_plot.df(ss.1T,x=x0,model=M1,PX=P1)
df1.T <- create_plot.df(ss.1T,x=x0.T,model=M1.T,PX=P1.T)

ggplot(df1,aes(x=speed,y=ConditionalD)) + 
  geom_line(aes(colour = "tt")) + 
  geom_line(data = df1.T,aes(x=speed,y=ConditionalD,colour = "tt.T")) +
  geom_vline(xintercept = v1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = v1.T, linetype = "dashed", color = "blue") +
  annotate("text", x = 20, y = 0.02, 
           label = paste("u=", round(v1,3)), color = "red")+
  scale_colour_manual(name = "Site", values = c("tt" = "red", "tt.T" = "blue")) +
  labs(x = "Speed (km/h)", y = "f(y|TTC<0)") +
  theme(panel.grid.major = element_line(colour = "gray91"),
        panel.grid.minor = element_line(colour = "gray88"),
        panel.background = element_rect(fill = "white",
                                        colour = "white", linetype = "solid"),
        plot.background = element_rect(linetype = "solid"))
Injury.from_c_bivariate(df1,EVmodel=M1,severity = PIS0,x0=x0)
Injury.from_c_bivariate(df1.T,EVmodel=M1.T,severity = PIS0,x0=x0.T)
