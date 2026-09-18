
lapply(list.files("../functions", pattern = "\\.R$", full.names = TRUE), source)

# cn_ttc <- read.table("data/TTC_CH.csv", sep = ",", header = TRUE)
# cn_pet <- read.table("data/PET_CH.csv", sep = ",", header = TRUE)
# cn_md  <- read.table("data/MD_CH.csv", sep = ",", header = TRUE)
# 
# names(cn_ttc)[c(6, 12)] <- c("Speed1_TTC", "Speed2_TTC")
# names(cn_pet)[c(6, 13)] <- c("Speed1_PET", "Speed2_PET")
# 
# cn_dat <- cn_md %>%
#   inner_join(cn_pet %>%
#     dplyr::select(track_id1, track_id2, PET, Speed1_PET, Speed2_PET),
#     by = c("track_id1", "track_id2")) %>%
#   inner_join(cn_ttc %>%
#     dplyr::select(track_id1, track_id2, TTC, Speed1_TTC, Speed2_TTC),
#     by = c("track_id1", "track_id2"))
# 
# cn_vru <- cn_dat %>%
#   subset(., (type1 == "PED" & type2 == "CAR") |
#             (type1 == "BIC" & type2 == "CAR") |
#             (type2 == "PED" & type1 == "CAR") |
#             (type2 == "BIC" & type1 == "CAR"))
# 
# cn_vru$N_PET <- -cn_vru$PET
# cn_vru$MD    <- -cn_vru$MD
# cn_vru$N_TTC <- -cn_vru$TTC
# cn_vru <- cn_vru %>%
#   mutate(Speed_TTC = case_when(type1 == "CAR" ~ Speed1_TTC,
#                                type2 == "CAR" ~ Speed2_TTC,
#                                TRUE ~ NA_real_) * 3.6) %>%
#   mutate(Speed_PET = case_when(type1 == "CAR" ~ Speed1_PET,
#                                type2 == "CAR" ~ Speed2_PET,
#                                TRUE ~ NA_real_) * 3.6)
# 
# cn_ttc_a <- cn_vru %>% mutate(prox = N_TTC, Speed = Speed_TTC) %>%
#   dplyr::select(prox, Speed) %>% subset(prox < 0, Speed < 50)
# 
# cn_pet_a <- cn_vru %>% mutate(prox = N_PET, Speed = Speed_PET) %>%
#   dplyr::select(prox, Speed) %>% subset(Speed < 50)
# 
# cn_ttc_i <- cn_ttc %>%
#   mutate(Speed_TTC = case_when(type1 == "CAR" ~ Speed1_TTC,
#                                type2 == "CAR" ~ Speed2_TTC,
#                                type1 == "TRUCK_BUS" ~ Speed1_TTC,
#                                type2 == "TRUCK_BUS" ~ Speed2_TTC,
#                                TRUE ~ NA_real_) * 3.6) %>%
#   dplyr::select(TTC, Speed_TTC, cycle, X15.min) %>%
#   subset(Speed_TTC < 50)
# 
# cn_pet_i <- cn_pet %>%
#   mutate(Speed_PET = case_when(type1 == "CAR" ~ Speed1_PET,
#                                type2 == "CAR" ~ Speed2_PET,
#                                type1 == "TRUCK_BUS" ~ Speed1_PET,
#                                type2 == "TRUCK_BUS" ~ Speed2_PET,
#                                TRUE ~ NA_real_) * 3.6) %>%
#   dplyr::select(PET, Speed_PET, cycle, X15.min) %>%
#   subset(Speed_PET < 50)

# SE data
se_dat <- read.table("../data/VehicleVRU_v5.csv", sep = ",", header = TRUE)

se_dat$N_PET <- -se_dat$PET
se_dat$N_TTC <- -se_dat$TTC
se_dat$maxDV_PET <- apply(se_dat[, c("DV1_PET", "DV2_PET")], 1, max)
se_dat$maxDV_TTC <- apply(se_dat[, c("DV1_TTC", "DV2_TTC")], 1, max)
se_dat <- se_dat[sapply(se_dat$maxDV_TTC, function(x) all(is.finite(x))), ]
se_dat <- se_dat %>%
  mutate(Speed_TTC = case_when(type1 == "vru" ~ Speed2_TTC,
                               type2 == "vru" ~ Speed1_TTC,
                               TRUE ~ NA_real_) * 3.6) %>%
  mutate(Speed_PET = case_when(type1 == "vru" ~ Speed2_PET,
                               type2 == "vru" ~ Speed1_PET,
                               TRUE ~ NA_real_) * 3.6)

se_ttc <- se_dat %>% mutate(prox = N_TTC, Speed = Speed_TTC) %>%
  dplyr::select(prox, Speed) %>% subset(Speed < 50)

se_pet <- se_dat %>% mutate(prox = N_PET, Speed = Speed_PET) %>%
  dplyr::select(prox, Speed) %>% subset(Speed < 50)
