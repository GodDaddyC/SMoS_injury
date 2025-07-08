# example of using Truncatedpbeved.R

# simulate some GP distribution
test_m1 <- revd(500,threshold = 5,scale = 1.5,shape = 0.2,type = "GP")
test_m2 <- revd(500,threshold = 16,scale = 0.8,shape = -0.15,type = "GP")
test_t <- as.matrix(data.frame(test_m1,test_m2))

# Transform it to unit Frechet margin
UF <- mtransform.GP(test_t,p=list(c(1.5,0.2),c(0.8,-0.15)),thres=c(5,16))

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

# note that the pbvevd use Exponential margin for parameterization, however the parameterization 
# for bivariate POT in Coles 2001 is unit Frechet. Note that the result will be the same if eta = 1,
# meaning the entire sample are beyond threshold. If eta is not the same, you have to use pbTvevd.

# if you use unit Frechet margin:
pbTvevd(q = c(10,18),model = "log",dep = test_bi$estimate[5],thres = c(5,16),eta = 1,tail.type = 2,
        mar1 = c(test_bi$estimate[1],test_bi$estimate[2]),mar2 = c(test_bi$estimate[3],test_bi$estimate[4]))
# if you use exponential margin:
pbvevd(q = c(10,18),model = "log",dep = test_bi$estimate[5],
       mar1 = c(5,test_bi$estimate[1],test_bi$estimate[2]),mar2 = c(16,test_bi$estimate[3],test_bi$estimate[4]))



# what does the tail type mean:
empProb(test_t,x = c(10,18),type =1 ,plot = TRUE)
empProb(test_t,x = c(10,18),type =2 ,plot = TRUE)
empProb(test_t,x = c(10,18),type =3 ,plot = TRUE)
empProb(test_t,x = c(10,18),type =4 ,plot = TRUE)

# get the conditional density

