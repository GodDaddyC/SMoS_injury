thres.order1T <- bvtcplot(CN_TTC_A)$k
thres.order2T <- bvtcplot(SE_TTC)$k
thres.order1P <- bvtcplot(CN_PET_A)$k
thres.order2P <- bvtcplot(SE_PET)$k

u.T1 <- sort(CN_TTC_A$N_TTC,decreasing = TRUE)[thres.order1T]
u.ST1 <- sort(CN_TTC_A$Speed_TTC,decreasing = TRUE)[thres.order1T]
u.T2 <- sort(SE_TTC$N_TTC,decreasing = TRUE)[thres.order2T]
u.ST2 <- sort(SE_TTC$Speed_TTC,decreasing = TRUE)[thres.order2T]

u.P1 <- sort(CN_PET_A$N_PET,decreasing = TRUE)[thres.order1P]
u.SP1 <- sort(CN_PET_A$Speed_PET,decreasing = TRUE)[thres.order1P]
u.P2 <- sort(SE_PET$N_PET,decreasing = TRUE)[thres.order2P]
u.SP2 <- sort(SE_PET$Speed_PET,decreasing = TRUE)[thres.order2P]

M.1T <- fbvpot(x = CN_TTC_A,model = "log",threshold = c(u.T1,u.ST1))
M.2T <- fbvpot(x = SE_TTC,model = "log",threshold = c(u.T2,u.ST2))
M.1P <- fbvpot(x = CN_PET_A,model = "log",threshold = c(u.P1,u.SP1))
M.2P <- fbvpot(x = SE_PET,model = "log",threshold = c(u.P2,u.SP2))

Pcrash.1T <- pevd(0,threshold = u.T1, scale = M.1T$estimate[1],shape = M.1T$estimate[2],
                  lower.tail = FALSE,type = "GP") * thres.order1T/dim(CN_TTC_A)[1]
Pcrash.2T <- pevd(0,threshold = u.T2, scale = M.2T$estimate[1],shape = M.2T$estimate[2],
                  lower.tail = FALSE,type = "GP") * thres.order2T/dim(SE_TTC)[1]

ss.1T <- seq(u.ST1,60,(60- u.ST1)/150)
ss.2T <- seq(u.ST2,60,(60- u.ST2)/150)

plot.df.1T <- data.frame(speed= ss.1T,
                         JointP = sapply(ss.1T,FUN = pbTvevd,q1=0,model="log",dep=M.1T$estimate[5],thres=c(u.T1,u.ST1),
                                         eta=thres.order1T/dim(CN_TTC_A)[1],
                                         mar1=c(M.1T$estimate[1],M.1T$estimate[2]),
                                         mar2=c(M.1T$estimate[3],M.1T$estimate[4]),tail.type=4)) %>%
  mutate(ConditionP = JointP/Pcrash.1T) %>% na.omit()
plot.df.2T <- data.frame(speed= ss.2T,
                         JointP = sapply(ss.2T,FUN = pbTvevd,q1=0,model="log",dep=M.2T$estimate[5],thres=c(u.T2,u.ST2),
                                         eta=thres.order2T/dim(SE_TTC)[1],
                                         mar1=c(M.2T$estimate[1],M.2T$estimate[2]),
                                         mar2=c(M.2T$estimate[3],M.2T$estimate[4]),tail.type=4)) %>%
  mutate(ConditionP = JointP/Pcrash.2T)

plot.df.1T <- plot.df.1T %>% 
  mutate(ConditionalD = sapply(speed, 
                               function(k) c.bivariate(y = k, x = 0, PX = Pcrash.1T, model = "log", 
                                                       dep = M.1T$estimate[5], thres = M.1T$threshold, 
                                                       eta = thres.order1T/dim(CN_TTC_A)[1],
                                                       mar1 = c(M.1T$estimate[1], M.1T$estimate[2]), 
                                                       mar2 = c(M.1T$estimate[3], M.1T$estimate[4])) )/
                                 nomralize_c.bivariate(x = 0, Pcrash = Pcrash.1T, EVmodel = M.1T)
  )

plot.df.2T <- plot.df.2T %>% 
  mutate(ConditionalD = sapply(speed, 
                               function(k) c.bivariate(y = k, x = 0, PX = Pcrash.2T, model = "log", 
                                                       dep = M.2T$estimate[5], thres = M.2T$threshold, 
                                                       eta = thres.order2T/dim(SE_TTC)[1],
                                                       mar1 = c(M.2T$estimate[1], M.2T$estimate[2]), 
                                                       mar2 = c(M.2T$estimate[3], M.2T$estimate[4])) )/
           nomralize_c.bivariate(x = 0, Pcrash = Pcrash.2T, EVmodel = M.2T)
  )


