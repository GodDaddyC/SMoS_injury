
source("auxfun.R")

CN_TTC <- read.table("data/TTC_CH.csv", sep = ",", header = TRUE)
CN_PET <- read.table("data/PET_CH.csv", sep = ",", header = TRUE) 
CN_MD <- read.table("data/MD_CH.csv", sep = ",", header = TRUE)

names(CN_TTC)[c(6,12)] <- c("Speed1_TTC","Speed2_TTC")
names(CN_PET)[c(6,13)] <- c("Speed1_PET","Speed2_PET")


CN_Dat <- CN_MD %>%
  inner_join(CN_PET %>% dplyr::select(track_id1, track_id2, PET,Speed1_PET,Speed2_PET),  
             by = c("track_id1","track_id2")) %>%
  inner_join(CN_TTC %>% dplyr::select(track_id1, track_id2, TTC,Speed1_TTC,Speed2_TTC),
             by = c("track_id1","track_id2")) 

CN_vru <- CN_Dat %>%
  subset(.,(type1=="PED" & type2=="CAR") | (type1=="BIC" & type2=="CAR") | 
           (type2=="PED" & type1=="CAR") | (type2=="BIC" & type1=="CAR") )

# CN_car <- CN_Dat %>%
#   subset(.,(type1=="PED" & type2=="CAR") | (type1=="BIC" & type2=="CAR") | 
#            (type2=="PED" & type1=="CAR") | (type2=="BIC" & type1=="CAR") )
# CN data


CN_vru$N_PET <- -CN_vru$PET
CN_vru$MD <- -CN_vru$MD
CN_vru$N_TTC <- -CN_vru$TTC
CN_vru <- CN_vru %>% 
  mutate(Speed_TTC = case_when(type1=="CAR"~ Speed1_TTC,
                               type2=="CAR"~ Speed2_TTC,
                               TRUE ~ NA_real_)*3.6) %>%
  mutate(Speed_PET = case_when(type1=="CAR"~ Speed1_PET,
                               type2=="CAR"~ Speed2_PET,
                               TRUE ~ NA_real_)*3.6)
CN_TTC_A <- CN_vru %>% mutate(prox = N_TTC,Speed=Speed_TTC) %>%
  dplyr::select(prox,Speed) %>% subset(Speed < 50) # aggregated data
CN_PET_A <- CN_vru %>% mutate(prox = N_PET,Speed=Speed_PET) %>%
  dplyr::select(prox,Speed) %>% subset(Speed < 50)
CN_iTTC_A <- CN_vru %>% mutate(prox = 1/(TTC+1),Speed=Speed_TTC) %>%
  dplyr::select(prox,Speed) %>% subset(Speed < 50)
CN_iPET_A <- CN_vru %>% mutate(prox = 1/(PET+5),Speed=Speed_PET) %>%
  dplyr::select(prox,Speed) %>% subset(Speed < 50)
  


CN_TTC_I <- CN_TTC %>% mutate(Speed_TTC = case_when(type1=="CAR"~ Speed1_TTC,  #separated data.
               type2=="CAR"~ Speed2_TTC,
               type1=="TRUCK_BUS" ~Speed1_TTC,
               type2=="TRUCK_BUS" ~Speed2_TTC,TRUE ~ NA_real_)*3.6) %>%
  dplyr::select(TTC,Speed_TTC,cycle,X15.min) %>% subset(Speed_TTC < 50)


CN_PET_I <- CN_PET %>% mutate(Speed_PET = case_when(type1=="CAR"~ Speed1_PET,
                type2=="CAR"~ Speed2_PET,
                type1=="TRUCK_BUS" ~Speed1_PET,
                type2=="TRUCK_BUS" ~Speed2_PET,TRUE ~ NA_real_)*3.6) %>%
  dplyr::select(PET,Speed_PET,cycle,X15.min) %>% subset(Speed_PET < 50)


# CN_vru$X30.mins <- ceiling(CN_vru$X15.min/2)

## BM for each 30 mins
# BM.CH <- data.frame(MD = blockmaxxer(CN_vru,which = "MD",blocks = CN_vru$X30.mins)$MD,
#                     PET = blockmaxxer(CN_vru,which = "PET",blocks = CN_vru$X30.mins)$PET,
#                     TTC = blockmaxxer(CN_vru,which = "TTC",blocks = CN_vru$X30.mins)$TTC,
#                     expo = sapply(unique(CN_vru$X30.mins),
#                                   function (x) {dim(subset(CN_vru,X30.mins==x))[1]})) %>% 
#   filter(.,PET>-3,TTC>-3,TTC<0)


# SE data
SE_Dat <- read.table("data/VehicleVRU_v5.csv", sep = ",", header = TRUE) 

#SE_Dat$N_MDc <- -SE_Dat$MDc
SE_Dat$N_PET <- -SE_Dat$PET
SE_Dat$N_TTC <- -SE_Dat$TTC
SE_Dat$maxDV_PET <- apply(SE_Dat[,c("DV1_PET","DV2_PET")],1,max)
SE_Dat$maxDV_TTC <- apply(SE_Dat[,c("DV1_TTC","DV2_TTC")],1,max)
SE_Dat <- SE_Dat[sapply(SE_Dat$maxDV_TTC, function(x) all(is.finite(x)) ), ]
SE_Dat <- SE_Dat %>% mutate(Speed_TTC = case_when(type1=="vru"~ Speed2_TTC,
                                                  type2=="vru"~ Speed1_TTC,
                                                  TRUE ~ NA_real_) * 3.6) %>%
  mutate(Speed_PET = case_when(type1=="vru"~ Speed2_PET,
                               type2=="vru"~ Speed1_PET,
                               TRUE ~ NA_real_) * 3.6) 

#create dataset for bivariate of TTC and PET
SE_TTC <- SE_Dat %>% mutate(prox=N_TTC,Speed=Speed_TTC) %>% subset(Speed < 50) # otherwise the tail is too heavy
SE_PET <- SE_Dat %>% mutate(prox=N_PET,Speed=Speed_PET) %>% subset(Speed < 50)
SE_iTTC <- SE_Dat %>% mutate(prox = 1/(TTC+1),Speed=Speed_TTC) %>%
  dplyr::select(prox,Speed) %>% subset(Speed < 50) 
SE_iPET <- SE_Dat %>% mutate(prox = 1/(PET+5),Speed=Speed_PET) %>% 
  dplyr::select(prox,Speed) %>% subset(Speed < 50) 
  

# SE_Dat$X30_mins <- ceiling(SE_Dat$X15_mins/2)


# BM.SE <- data.frame(MD=blockmaxxer(SE_Dat,which = "N_MD",blocks = SE_Dat$X30_mins)$N_MD,
#                    PET=blockmaxxer(SE_Dat,which = "N_PET",blocks = SE_Dat$X30_mins)$N_PET,
#                    TTC=blockmaxxer(SE_Dat,which = "N_TTC",blocks = SE_Dat$X30_mins)$N_TTC,
#                    expo = sapply(unique(SE_Dat$X30_mins),
#                                  function (x) {dim(subset(SE_Dat,X30_mins==x))[1]})) %>%
#   filter(.,PET>-3,TTC>-3,TTC<0)


# clean the variable space a bit
# rm(list = ls(pattern = "Dat."))
