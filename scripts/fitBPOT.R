source('../scripts/preprocessing.R')
Dat <- se_ttc
x0 <- 0
source("../scripts/global_vars.R")

# adjust the speed_ub to be lower if the error non-finite value exists (default = 80)
log_M <- run_single_bpot(Dat, model = "log", thres = c(u, v), xcrash = x0,estim = 'margins')
neglog_M <- run_single_bpot(Dat, model = "neglog", thres = c(u, v), xcrash = x0,estim = 'margins') 
hr_M <- run_single_bpot(Dat, model = "hr", thres = c(u, v), xcrash = x0,estim = 'margins')
ct_M <- run_single_bpot(Dat, model = "ct", thres = c(u, v), xcrash = x0,estim = 'margins')
alog_M <- run_single_bpot(Dat, model = "alog", thres = c(u, v),xcrash = x0,estim = 'margins')
aneglog_M <- run_single_bpot(Dat, model = "aneglog", thres = c(u, v), xcrash = x0,estim = 'margins')
bilog_M <- run_single_bpot(Dat, model = "bilog", thres = c(u, v), xcrash = x0,estim = 'margins')



all_bpot_results <- list("logistic"= log_M,
                         "neglogistic"= neglog_M,
                         "hr"= hr_M,
                         "ct"= ct_M,
                         "alogistic"= alog_M,
                         #"aneglogistic"= aneglog_M,
                         "bilogistic"= bilog_M)
