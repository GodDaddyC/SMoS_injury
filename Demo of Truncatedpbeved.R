# example of using Truncatedpbeved.R
# Simulated data
# simulate some GP distribution
test_m1 <- revd(500,threshold = 5,scale = 1.5,shape = 0.2,type = "GP")
test_m2 <- revd(500,threshold = 16,scale = 0.8,shape = -0.15,type = "GP")

test_t <- as.matrix(data.frame(test_m1,test_m2))


# How to fit a bivariate extreme value distribution truncated below threshold and 
# calulate the probability?

# directly fit
test_bi <- fbvpot(x = test_t,model = "log",threshold = c(5,16))
test_bi$estimate

# compute the probability of the bivariate extreme value distribution truncated below threshold
pbTvevd(q1=10,q2=18,model = "log",dep = test_bi$estimate[5],thres = c(5,16),eta = 1,tail.type = 2,
        mar1 = c(test_bi$estimate[1],test_bi$estimate[2]),mar2 = c(test_bi$estimate[3],test_bi$estimate[4]))

# what does the tail type mean:
empProb(test_t,x = c(10,17),type =1 ,plot = TRUE)
empProb(test_t,x = c(10,17),type =2 ,plot = TRUE)
empProb(test_t,x = c(10,17),type =3 ,plot = TRUE)
empProb(test_t,x = c(10,17),type =4 ,plot = TRUE)

# get the conditional density

dbTvevd(q1=10,q2=17,model = "log",dep = test_bi$estimate[5],thres = c(5,16),eta = 1,
        mar1 = c(test_bi$estimate[1],test_bi$estimate[2]),mar2 = c(test_bi$estimate[3],test_bi$estimate[4]))
c.bivariate(y = 17, x = 10, PX = 0.077, model = "log", 
                      dep = test_bi$estimate[5], thres = c(5,16), eta = 1,
                      mar1 = c(test_bi$estimate[1], test_bi$estimate[2]), 
                      mar2 = c(test_bi$estimate[3], test_bi$estimate[4]))

sapply(seq(16,18,0.1), function(k) c.bivariate(y = k, x = 10, PX = 0.077, model = "log", 
                      dep = test_bi$estimate[5], thres = c(5,16), eta = 1,
                      mar1 = c(test_bi$estimate[1], test_bi$estimate[2]), 
                      mar2 = c(test_bi$estimate[3], test_bi$estimate[4])))

# real example:
source("auxfun.R")

SE_Dat <- read.table("data/VehicleVRU_v5.csv", sep = ",", header = TRUE) 

SE_Dat$N_MD <- -SE_Dat$MD
#SE_Dat$N_MDc <- -SE_Dat$MDc
SE_Dat$N_PET <- -SE_Dat$PET
SE_Dat$N_TTC <- -SE_Dat$TTC
SE_Dat$maxDV_PET <- apply(SE_Dat[,c("DV1_PET","DV2_PET")],1,max)
SE_Dat$maxDV_TTC <- apply(SE_Dat[,c("DV1_TTC","DV2_TTC")],1,max)
SE_Dat <- SE_Dat[sapply(SE_Dat$maxDV_TTC, function(x) all(is.finite(x)) ), ]
SE_Dat <- SE_Dat %>% mutate(Speed_TTC = case_when(type1=="vru"~ Speed2_TTC,
                                                  type2=="vru"~ Speed1_TTC,
                                                  TRUE ~ NA_real_) * 3.6) %>%
  mutate(Speed_PET = case_when(type1=="vru"~ Speed2_PET,
                               type2=="vru"~ Speed1_PET,
                               TRUE ~ NA_real_) * 3.6) 

#create dataset for bivariate of TTC and PET
SE_TTC <- SE_Dat %>% select(N_TTC,Speed_TTC) %>% subset(Speed_TTC < 50) # otherwise the tail is too heavy
SE_PET <- SE_Dat %>% select(N_PET,Speed_PET) %>% subset(Speed_PET < 50)

# fitting bivariate model
thres.order <- bvtcplot(SE_TTC)$k
u.TTC <- sort(SE_TTC$N_TTC,decreasing = TRUE)[thres.order]
u.S_TTC <- sort(SE_TTC$Speed_TTC,decreasing = TRUE)[thres.order]

M.test <- fbvpot(x = SE_TTC,model = "log",threshold = c(u.TTC,u.S_TTC))
M.test$estimate

# does the marginal fitting do well in the chosen threshold
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
# if an non-finite error exists, change the `upper` in `c.bivariate` to a smaller value (e.g. 2)

Pcrash <- pevd(0,threshold = u.TTC, scale = M.test$estimate[1],shape = M.test$estimate[2],
               lower.tail = FALSE,type = "GP") * thres.order/dim(SE_TTC)[1]
ConditionP <- JointP/Pcrash
plot(ss,ConditionP) # plot the conditional probability P(Y <= y | X>0)

# compute injury
# AIS 3+ for pedetrain Kong and Yang (2010)
PIS0 <- function(speed){return(1/(1+exp(5.261 - 0.104*speed))) } # Eq.6 

injury.intergrand <- function(y){
  f.y <- c.bivariate(y = y, x = 0, PX = Pcrash, model = "log", 
                     dep = M.test$estimate[5], thres = c(u.TTC,u.S_TTC), eta = thres.order/dim(SE_TTC)[1],
                     mar1 = c(M.test$estimate[1], M.test$estimate[2]), 
                     mar2 = c(M.test$estimate[3], M.test$estimate[4]))
  return(PIS0(y) * f.y)
}

#AIS 3+ injury probability | crash
integrate(injury.intergrand,lower = min(ss),upper=Inf)$value

# prob of injury crash
integrate(injury.intergrand,lower = min(ss),upper=Inf)$value * Pcrash

# yearly frequency of injury crash
integrate(injury.intergrand,lower = min(ss),upper=max(ss))$value * Pcrash *
  dim(SE_TTC)[1] * thres.order/dim(SE_TTC)[1] /33 * 365 
