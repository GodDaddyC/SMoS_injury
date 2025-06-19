# POT models

POT.CN.TTC <- fevd(x = TTC,data=CN_vru,threshold = quantile(CN_vru$TTC,0.8), type = "GP")
POT.CN.TTC$results
plot(POT.CN.TTC)
