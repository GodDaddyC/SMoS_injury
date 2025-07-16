
source("auxfun.R")

CN_TTC <- read.table("data/TTC_CH.csv", sep = ",", header = TRUE)
CN_PET <- read.table("data/PET_CH.csv", sep = ",", header = TRUE) 
CN_MD <- read.table("data/MD_CH.csv", sep = ",", header = TRUE) 

CN_Dat <- CN_MD %>%
  inner_join(CN_PET %>% select(track_id1, track_id2, PET),  by = c("track_id1","track_id2")) %>%
  inner_join(CN_TTC %>% select(track_id1, track_id2, TTC),  by = c("track_id1","track_id2")) 

CN_vru <- CN_Dat %>%
  subset(.,(type1=="PED" & type2=="CAR") | (type1=="BIC" & type2=="CAR") | 
           (type2=="PED" & type1=="CAR") | (type2=="BIC" & type1=="CAR") )

CN_car <- CN_Dat %>%
  subset(.,(type1=="PED" & type2=="CAR") | (type1=="BIC" & type2=="CAR") | 
           (type2=="PED" & type1=="CAR") | (type2=="BIC" & type1=="CAR") )
# CN data


CN_vru$PET <- -CN_vru$PET
CN_vru$MD <- -CN_vru$MD
CN_vru$TTC <- -CN_vru$TTC
CN_vru$X30.mins <- ceiling(CN_vru$X15.min/2)

## BM for each 30 mins
BM.CH <- data.frame(MD = blockmaxxer(CN_vru,which = "MD",blocks = CN_vru$X30.mins)$MD,
                    PET = blockmaxxer(CN_vru,which = "PET",blocks = CN_vru$X30.mins)$PET,
                    TTC = blockmaxxer(CN_vru,which = "TTC",blocks = CN_vru$X30.mins)$TTC,
                    expo = sapply(unique(CN_vru$X30.mins),
                                  function (x) {dim(subset(CN_vru,X30.mins==x))[1]})) %>% 
  filter(.,PET>-3,TTC>-3,TTC<0)


# SE data
SE_Dat <- read.table("data/data_SWE22.csv", sep = ",", header = TRUE) 

SE_Dat$N_MD <- -SE_Dat$MD
#SE_Dat$N_MDc <- -SE_Dat$MDc
SE_Dat$N_PET <- -SE_Dat$PET
SE_Dat$N_TTC <- -SE_Dat$TTC
SE_Dat$maxDV_PET <- apply(SE_Dat[,c("DV1_PET","DV2_PET")],1,max)
SE_Dat$maxDV_TTC <- apply(SE_Dat[,c("DV1_TTC","DV2_TTC")],1,max)
SE_Dat <- SE_Dat[sapply(SE_Dat$maxDV_TTC, function(x) all(is.finite(x)) ), ]

#create dataset for bivariate of TTC and PET
SE_TTC <- SE_Dat %>% select(N_TTC,maxDV_TTC) %>% subset(maxDV_TTC < 16) # otherwise the tail is too heavy
SE_PET <- SE_Dat %>% select(N_PET,maxDV_PET) %>% subset(maxDV_PET < 16)

# SE_Dat$X30_mins <- ceiling(SE_Dat$X15_mins/2)


# BM.SE <- data.frame(MD=blockmaxxer(SE_Dat,which = "N_MD",blocks = SE_Dat$X30_mins)$N_MD,
#                    PET=blockmaxxer(SE_Dat,which = "N_PET",blocks = SE_Dat$X30_mins)$N_PET,
#                    TTC=blockmaxxer(SE_Dat,which = "N_TTC",blocks = SE_Dat$X30_mins)$N_TTC,
#                    expo = sapply(unique(SE_Dat$X30_mins),
#                                  function (x) {dim(subset(SE_Dat,X30_mins==x))[1]})) %>%
#   filter(.,PET>-3,TTC>-3,TTC<0)


# clean the variable space a bit
rm(list = ls(pattern = "Dat."))
