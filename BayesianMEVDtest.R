# test the baeysian estimate of dependence parameter

library(ExtremalDep)
data(pollution)

f.hr <- fExtDep(x=PNS, method="BayesianPPP", model="HR",Hpar = list(mean.lambda=0.7,sd.lambda=0.05),
                MCpar = list(sdq=0.05,sdw=0.02,sdw=0.02),
                Nbin=2000,Nsim=10000,
                par.start = rep(0.5, 3))

plot(x=f.hr, type="pickands",
     labels=c(expression(PM[10]), expression(NO), expression(SO[2])), 
     cex.lab=2)


data(WindSpeedGust)

years <- format(ParcayMeslay$time, format="%Y")
attach(ParcayMeslay[which(years %in% c(2004:2013)),])

# Marginal quantiles
WS_th <- quantile(WS,.9)
DP_th <- quantile(DP,.9)

# Standardisation to unit Frechet (requires evd package)
pars.WS <- evd::fpot(WS, WS_th, model="pp")$estimate
pars.DP <- evd::fpot(DP, DP_th, model="pp")$estimate

# transform the marginal distribution to common unit Frechet:
data_uf <- trans2UFrechet(cbind(WS,DP), type="Empirical")

# compute exceedances 
rdata <- rowSums(data_uf)
r0 <- quantile(rdata, probs=.90)
extdata_WSDP <- data_uf[rdata>=r0,]

# Fit
SP_mle <- fExtDep.np(x=extdata_WSDP, method="Frequentist", k0=10, type="maxima")

# Plot
plot(x=SP_mle, type="summary")

####################################################
### Example 2 - Pollution levels in Milan, Italy ###
####################################################

## Not run: 

### Here we will only model the dependence structure	
data(MilanPollution)

data <- Milan.winter[,c("NO2","SO2")] 
data <- as.matrix(data[complete.cases(data),])

# Thereshold
u <- apply(data, 2, function(x) quantile(x, prob=0.9, type=3))

# Hyperparameters
hyperparam <- list(mu.nbinom = 6, var.nbinom = 8, a.unif=0, b.unif=0.2)

### Standardise data to univariate Frechet margins

f1 <- fGEV(data=data[,1], method="Bayesian", sig0 = 0.0001, nsim = 5e+4)
diagnostics(f1)
burn1 <- 1:30000
gev.pars1 <- apply(f1$param_post[-burn1,],2,mean)
sdata1 <- trans2UFrechet(data=data[,1], pars=gev.pars1, type="GEV")

f2 <- fGEV(data=data[,2], method="Bayesian", sig0 = 0.0001, nsim = 5e+4)
diagnostics(f2)
burn2 <- 1:30000
gev.pars2 <- apply(f2$param_post[-burn2,],2,mean)
sdata2 <- trans2UFrechet(data=data[,2], pars=gev.pars2, type="GEV")

sdata <- cbind(sdata1,sdata2)

### Bayesian estimation using Bernstein polynomials

pollut1 <- fExtDep.np(x=sdata, method="Bayesian", u=TRUE,
                      mar.fit=FALSE, k0=5, hyperparam = hyperparam, nsim=5e+4)

diagnostics(pollut1)
pollut1_sum <- summary(object=pollut1, burn=3e+4, plot=TRUE)
pl1 <- plot(x=pollut1, type="Qsets", summary.mcmc=pollut1_sum, 
            mar1=gev.pars1, mar2=gev.pars2, P = 1/c(600, 1200, 2400),
            dep=TRUE, data=data, xlim=c(0,400), ylim=c(0,400))

pl1b <- plot(x=pollut1, type="Qsets", summary.mcmc=pollut1_sum, est.out=pl1$est.out, 
             mar1=gev.pars1, mar2=gev.pars2, P = 1/c(1200),
             dep=FALSE, data=data, xlim=c(0,400), ylim=c(0,400))

### Frequentist estimation using Bernstein polynomials

pollut2 <- fExtDep.np(x=sdata, method="Frequentist", mar.fit=FALSE, type="rawdata", k0=8)
plot(x=pollut2, type = "summary", CEX=1.5)

pl2 <- plot(x=pollut2, type="Qsets", mar1=gev.pars1, mar2=gev.pars2, 
            P = 1/c(600, 1200, 2400),
            dep=TRUE, data=data, xlim=c(0,400), ylim=c(0,400),
            xlab=expression(NO[2]), ylab=expression(SO[2]),
            col.Qfull = c("red", "green", "blue"))

### Frequentist estimation using EKdH estimator

pollut3 <- fExtDep.np(x=data, method="Empirical")
plot(x=pollut3, type = "summary", CEX=1.5)

pl3 <- plot(x=pollut3, type="Qsets", mar1=gev.pars1, mar2=gev.pars2,
            P = 1/c(600, 1200, 2400),
            dep=TRUE, data=data, xlim=c(0,400), ylim=c(0,400),
            xlab=expression(NO[2]), ylab=expression(SO[2]),
            col.Qfull = c("red", "green", "blue"))


## End(Not run)