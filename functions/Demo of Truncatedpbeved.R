# example of using Truncatedpbeved.R

# simulate some GP distribution
test_m1 <- revd(500,threshold = 5,scale = 1.5,shape = 0.2,type = "GP")
test_m2 <- revd(500,threshold = 16,scale = 0.8,shape = -0.15,type = "GP")
test_truth <- rgev(500,1,1,0.4)
test_t <- as.matrix(data.frame(test_m1,test_m2))

data(damage)
fit1 <- fevd(Dam, damage, threshold=6, type="GP", time.units="2.05/year")
fit1$results$par
fit11 <- mtransform(subset(damage,Dam>6) %>% select(Dam),c(6,4.589,0.51))
fit12 <- mtransform.GP(subset(damage,Dam>6) %>% select(Dam),c(4.589,0.51),thres = 6,eta = 0.125)

# Transform it to unit Frechet margin
UF <- mtransform.GP(test_t,p=list(c(1.5,0.2),c(0.8,-0.15)),thres=c(5,16))
# when the margins are GEV, mtransform maps the margins to standard Exp
# when the margins are GP, mtransform maps the margins to uniform
# mtransform.GP maps unconditional GP margins to unit Frechet
UF2 <- mtransform(test_t,p=list(c(5,1.5,0.2),c(16,0.8,-0.15)))

# is it unit Frechet? Almost..
UF.t <- fevd(x = UF[,1],type = "GEV")
UF.t$results$par
UF.t2 <- fevd(x = UF[,2],type = "GEV")
UF.t2$results$par


# How to fit a bivariate extreme value distribution truncated below threshold and 
# calulate the probability?

# directly fit
test_bi <- fbvpot(x = test_t,model = "log",threshold = c(5,16))
test_bi$estimate



pbTvevd(q = c(10,18),model = "log",dep = test_bi$estimate[5],thres = c(5,16),eta = 1,tail.type = 1,
        mar1 = c(test_bi$estimate[1],test_bi$estimate[2]),mar2 = c(test_bi$estimate[3],test_bi$estimate[4]))

# what does the tail type mean:
empProb(test_t,x = c(10,18),type =1 ,plot = TRUE)
empProb(test_t,x = c(10,18),type =2 ,plot = TRUE)
empProb(test_t,x = c(10,18),type =3 ,plot = TRUE)
empProb(test_t,x = c(10,18),type =4 ,plot = TRUE)

# get the conditional density

dbTvevd(q1=7,q2=17,model = "log",dep = test_bi$estimate[5],thres = c(5,16),eta = 1,
        mar1 = c(test_bi$estimate[1],test_bi$estimate[2]),mar2 = c(test_bi$estimate[3],test_bi$estimate[4]))
c.bivariate(y = 18, x = 10, PX = 0.077, q2 = 18, model = "log", 
                      dep = test_bi$estimate[5], thres = c(5,16), eta = 1,
                      mar1 = c(test_bi$estimate[1], test_bi$estimate[2]), 
                      mar2 = c(test_bi$estimate[3], test_bi$estimate[4]))

c.bivariate(y = 18,x = 10,PX = 0.077,q2=18,model = "log",dep = test_bi$estimate[5],thres = c(5,16),eta = 1,
            mar1 = c(test_bi$estimate[1],test_bi$estimate[2]),mar2 = c(test_bi$estimate[3],test_bi$estimate[4]))
c.bivariate(y = 18, x = 10, PX = 0.077, model = "log", 
                      dep = test_bi$estimate[5], thres = c(5,16), eta = 1,
                      mar1 = c(test_bi$estimate[1], test_bi$estimate[2]), 
                      mar2 = c(test_bi$estimate[3], test_bi$estimate[4]))
sapply(seq(16,20,0.1), function(k) c.bivariate(y = k, x = 10, PX = 0.077, q2 = k, model = "log", 
                      dep = test_bi$estimate[5], thres = c(5,16), eta = 1,
                      mar1 = c(test_bi$estimate[1], test_bi$estimate[2]), 
                      mar2 = c(test_bi$estimate[3], test_bi$estimate[4])))
