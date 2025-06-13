
B_CH.PET0 <- fevd(x=PET,data=BM.CH, type = "GEV")
B_CH.PET1 <- fevd(x=PET,data=BM.CH,location.fun =~expo, type = "GEV")
B_CH.TTC0 <- fevd(x=TTC,data=BM.CH, type = "GEV")
B_CH.TTC1 <- fevd(x=TTC,data=BM.CH,location.fun =~expo, type = "GEV")

B_SE.PET0 <- fevd(x=PET,data=BM.SE, type = "GEV")
B_SE.PET1 <- fevd(x=PET,data=BM.SE,location.fun =~expo,type = "GEV")
B_SE.TTC0 <- fevd(x=TTC,data=BM.SE, type = "GEV")
B_SE.TTC1 <- fevd(x=TTC,data=BM.SE,location.fun =~expo,type = "GEV")

combined_data <- bind_rows(
  BM.CH[,c("TTC","PET","MD","expo")] %>% mutate(from = "CH"),
  BM.SE[,c("TTC","PET","MD","expo")] %>% mutate(from = "SE")
)

M.join.1 <- fevd(x=MD,data=combined_data,location.fun =~expo,type = "GEV")
M.join.2 <- fevd(x=MD,data=combined_data,location.fun =~expo + from,type = "GEV")
lr.test(M.join.1,M.join.2)
