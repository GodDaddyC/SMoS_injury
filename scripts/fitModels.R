source('../scripts/preprocessing.R')
Dat <- se_ttc
x0 <- 0
source("../scripts/global_vars.R")

# adjust the speed_ub to be lower if the error non-finite value exists (default = 60)
log_M <- run_single_bpot(Dat, model = "log", thres = c(u, v), xcrash = x0)
neglog_M <- run_single_bpot(Dat, model = "neglog", thres = c(u, v), xcrash = x0) 
hr_M <- run_single_bpot(Dat, model = "hr", thres = c(u, v), xcrash = x0)
ct_M <- run_single_bpot(Dat, model = "ct", thres = c(u, v), xcrash = x0)
alog_M <- run_single_bpot(Dat, model = "alog", thres = c(u, v),xcrash = x0)
aneglog_M <- run_single_bpot(Dat, model = "aneglog", thres = c(u, v), xcrash = x0)
bilog_M <- run_single_bpot(Dat, model = "bilog", thres = c(u, v), xcrash = x0)

gumbel_M <- run_single_cbpot(cop_dat, copula = gumbelCopula(),
                             p2 = conseq, qcrash = qcrash, pot = pot,method = "ml")
clayton_M <- run_single_cbpot(cop_dat, copula = claytonCopula(),
                              p2 = conseq, qcrash = qcrash, pot = pot,method = "ml")
gaussian_M <- run_single_cbpot(cop_dat,copula = normalCopula(dispstr = "ex"),
                               p2 = conseq, qcrash = qcrash, pot = pot,
                               method = "ml")