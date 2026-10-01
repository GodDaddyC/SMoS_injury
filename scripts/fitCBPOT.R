source('../scripts/preprocessing.R')
source("../scripts/global_vars.R")

frank_M <- run_single_cbpot(cop_dat, copula = frankCopula(),
                             p2 = conseq, qcrash = qcrash, pot = pot,method = "ml")
clayton_M <- run_single_cbpot(cop_dat, copula = claytonCopula(),
                              p2 = conseq, qcrash = qcrash, pot = pot,method = "ml")
gaussian_M <- run_single_cbpot(cop_dat,copula = normalCopula(dispstr = "ex"),
                               p2 = conseq, qcrash = qcrash, pot = pot,
                               method = "ml")
gumbel_M <- run_single_cbpot(cop_dat, copula = gumbelCopula(),
                            p2 = conseq, qcrash = qcrash, pot = pot,method = "ml")

all_cbpot_results <- list("frank"= frank_M,
                          "clayton"= clayton_M,
                          "gaussian"= gaussian_M,
                          "gumbel"= gumbel_M)
