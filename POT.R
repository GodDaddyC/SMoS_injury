# POT models


# change the mrlplot to the one you used in the other paper, 
evmix::mrlplot(SE_Dat$N_TTC)
evmix::mrlplot(SE_Dat$N_PET)
evmix::mrlplot(SE_Dat$maxDV_PET)
evmix::mrlplot(SE_Dat$maxDV_TTC)

evmix::mrlplot(CN_vru$TTC)
# Coordinates for the circle and the arrow
# x_circle <- -4.7  # Adjust based on your data
# y_circle <- 1.4  # Adjust based on your data
# x_arrow_from <- -4.7  # Arrow start x
# y_arrow_from <- 2.5  # Arrow start y
# x_arrow_to <- x_circle  # Arrow end x
# y_arrow_to <- 1.9  # Arrow end y

# Add a circle at (x_circle, y_circle)
# symbols(x_circle, y_circle, circles = 1, inches = 0.45,, add = TRUE)
# # Add an arrow pointing to the circle
# arrows(x_arrow_from, y_arrow_from, x_arrow_to, y_arrow_to, col = "red", lwd = 2)
# text(x_arrow_from, y_arrow_from + 0.1, labels = "Select threshold from this region", col = "red", pos = 3)

evd::tcplot(CN_vru$TTC,tlim = c(-2.5,-1.5))
evd::tcplot(SE_Dat$N_PET,tlim = c(-3,-2))




# run the threshold selection for CN and SE


# fit the univariate POT model to the CN and SE data

POT.CN.TTC <- fevd(x = TTC,data=CN_vru,threshold = quantile(CN_vru$TTC,0.8), type = "GP",period.basis = "two weeks",
                   time.units = "2/month")
POT.CN.TTC$results
POT.CN.PET <- fevd(x = PET,data=CN_vru,threshold = quantile(CN_vru$TTC,0.8), type = "GP",period.basis = "two weeks",
                   time.units = "2/month")
POT.CN.PET$results
POT.SE.TTC <- fevd(x = N_TTC,data=SE_TTC,threshold = quantile(SE_TTC$N_TTC,0.8),period.basis = "month",
                  time.units = "months", type = "GP")
POT.SE.TTC$results
POT.SE.PET<- fevd(x = N_PET,data=SE_PET,threshold = quantile(SE_PET$N_PET,0.8), period.basis = "month",
                  time.units = "months",type = "GP")
POT.SE.PET$results

POT.SE.DV_PET <- fevd(x = maxDV_PET,data=SE_PET,threshold = quantile(SE_PET$maxDV_PET,0.8), period.basis = "month",
                      time.units = "months", type = "GP")
POT.SE.DV_PET$results
POT.SE.DV_TTC <- fevd(x = maxDV_TTC,data=SE_TTC,threshold = quantile(SE_TTC$maxDV_TTC,0.8), period.basis = "month",
                      time.units = "months", type = "GP")
POT.SE.DV_TTC$results # need to remove DV >15 to prevent very heavy tail in the estimation.

plot(POT.CN.TTC)
plot(POT.CN.PET)
plot(POT.SE.TTC)
plot(POT.SE.PET)
plot(POT.SE.DV_PET)
plot(POT.SE.DV_TTC)

# fitting bivariate model
thres.order <- bvtcplot(SE_TTC)$k
u.TTC <- sort(SE_TTC$N_TTC,decreasing = TRUE)[thres.order]
u.S_TTC <- sort(SE_TTC$Speed_TTC,decreasing = TRUE)[thres.order]

M.test <- fbvpot(x = SE_TTC,model = "log",threshold = c(u.TTC,u.S_TTC))
M.test$estimate

# does the marginal fitting do well in the 
M.prox <- fevd(x = N_TTC,data=SE_TTC,threshold = u.TTC,period.basis = "month",
               time.units = "months", type = "GP")
plot(M.prox)
M.conseq <- fevd(x = Speed_TTC,data=SE_TTC,threshold = u.S_TTC,period.basis = "month",
               time.units = "months", type = "GP")
plot(M.conseq)


# compute the crash proability at different level of DV, this computes P(X>0,Y>6)
pbTvlog(q1=0,q2=16,dep=M.test$estimate[5],thres=c(u.TTC,u.S_TTC),eta=thres.order/dim(SE_TTC)[1],
           mar1=c(M.test$estimate[1],M.test$estimate[2]),
           mar2=c(M.test$estimate[3],M.test$estimate[4]),tail.type=3)

# the conditional distribution P(Y <= y | X>0), y> u.y, cutoff at y = 16
ss <- seq(u.S_TTC,50,(50- u.S_TTC)/100)
JointP <- sapply(ss,FUN = pbTvevd,q1=0,model="log",dep=M.test$estimate[5],thres=c(u.TTC,u.S_TTC),
       eta=thres.order/dim(SE_TTC)[1],
       mar1=c(M.test$estimate[1],M.test$estimate[2]),
       mar2=c(M.test$estimate[3],M.test$estimate[4]),tail.type=4)
names(JointP) <- NULL
Pcrash <- pevd(0,threshold = u.TTC, scale = M.test$estimate[1],shape = M.test$estimate[2],
               lower.tail = FALSE,type = "GP") * thres.order/dim(SE_TTC)[1]
ConditionP <- JointP/Pcrash
plot(ss,ConditionP) # plot the conditional probability P(Y <= y | X>0)

# compute the conditional density f(y|X>0).
# c.bivariate(y = 5, x = 0, PX = Pcrash, model = "log", 
#             dep = M.test$estimate[5], thres = c(u.TTC,u.DV_TTC), eta = thres.order/dim(SE_TTC)[1],
#             mar1 = c(M.test$estimate[1], M.test$estimate[2]), 
#             mar2 = c(M.test$estimate[3], M.test$estimate[4]))
Conditionf<- sapply(ss, function(k) c.bivariate(y = k, x = 0, PX = Pcrash, model = "log", 
               dep = M.test$estimate[5], thres = c(u.TTC,u.S_TTC), eta = thres.order/dim(SE_TTC)[1],
               mar1 = c(M.test$estimate[1], M.test$estimate[2]), 
               mar2 = c(M.test$estimate[3], M.test$estimate[4])))
plot(ss,Conditionf)


injury.intergrand <- function(y){
  f.y <- c.bivariate(y = y, x = 0, PX = Pcrash, model = "log", 
                     dep = M.test$estimate[5], thres = c(u.TTC,u.S_TTC), eta = thres.order/dim(SE_TTC)[1],
                     mar1 = c(M.test$estimate[1], M.test$estimate[2]), 
                     mar2 = c(M.test$estimate[3], M.test$estimate[4]))
  return(PIS0(y) * f.y)
}

#AIS 3+ injury probability | crash
integrate(injury.intergrand,lower = min(ss),upper=max(ss))$value

# prob of injury crash
integrate(injury.intergrand,lower = min(ss),upper=max(ss))$value * Pcrash

# yearly frequency of injury crash
integrate(injury.intergrand,lower = min(ss),upper=max(ss))$value * Pcrash *
  dim(SE_TTC)[1] * thres.order/dim(SE_TTC)[1] /33 * 365 
