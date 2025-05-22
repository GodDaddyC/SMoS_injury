M.CH1 <- fevd(x=MD,data=BM.CH, type = "GEV")
M.CH1.1 <- fevd(x=MD,data=BM.CH,location.fun =~expo, type = "GEV")
M.CH2 <- fevd(x=PET,data=BM.CH, type = "GEV")
M.CH2.1 <- fevd(x=PET,data=BM.CH,location.fun =~expo, type = "GEV")
M.CH3 <- fevd(x=TTC,data=BM.CH, type = "GEV")
M.CH3.1 <- fevd(x=TTC,data=BM.CH,location.fun =~expo, type = "GEV")

M.SE1 <- fevd(x=MD,data=BM.SE, type = "GEV")
M.SE1.1 <- fevd(x=MD,data=BM.SE,location.fun =~expo,type = "GEV")
M.SE2 <- fevd(x=PET,data=BM.SE, type = "GEV")
M.SE2.1 <- fevd(x=PET,data=BM.SE,location.fun =~expo,type = "GEV")
M.SE3 <- fevd(x=TTC,data=BM.SE, type = "GEV")
M.SE3.1 <- fevd(x=TTC,data=BM.SE,location.fun =~expo,type = "GEV")

combined_data <- bind_rows(
  BM.CH[,c("TTC","PET","MD","expo")] %>% mutate(from = "CH"),
  BM.SE[,c("TTC","PET","MD","expo")] %>% mutate(from = "SE")
)

M.join.1 <- fevd(x=MD,data=combined_data,location.fun =~expo,type = "GEV")
M.join.2 <- fevd(x=MD,data=combined_data,location.fun =~expo + from,type = "GEV")
lr.test(M.join.1,M.join.2)
