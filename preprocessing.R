source("auxfun.R")


# CN data
Dat.TTC <- read.table("data/TTC_CH.csv", sep = ",", header = TRUE) %>%
  subset(.,(type1=="PED" & type2=="CAR") | (type1=="BIC" & type2=="CAR") | 
           (type2=="PED" & type1=="CAR") | (type2=="BIC" & type1=="CAR") )
Dat.TTC$X30.mins <- ceiling(Dat.TTC$X15.min/2)

Dat.PET <- read.table("data/PET_CH.csv", sep = ",", header = TRUE) %>%
  subset(.,(type1=="PED" & type2=="CAR") | (type1=="BIC" & type2=="CAR") | 
           (type2=="PED" & type1=="CAR") | (type2=="BIC" & type1=="CAR") )

Dat.MD  <- read.table("data/MD_CH.csv", sep = ",", header = TRUE) %>%
  subset(.,(type1=="PED" & type2=="CAR") | (type1=="BIC" & type2=="CAR") | 
           (type2=="PED" & type1=="CAR") | (type2=="BIC" & type1=="CAR") )

CH_Dat <- Dat.MD %>%
  inner_join(Dat.PET %>% select(track_id1, track_id2, PET),  by = c("track_id1","track_id2")) %>%
  inner_join(Dat.TTC %>% select(track_id1, track_id2, TTC),  by = c("track_id1","track_id2")) %>%
  select(TTC,PET,MD,cycle,X15.min)

CH_Dat$PET <- -CH_Dat$PET
CH_Dat$MD <- -CH_Dat$MD
CH_Dat$TTC <- -CH_Dat$TTC
CH_Dat$X30.mins <- ceiling(CH_Dat$X15.min/2)

## BM for each 30 mins
BM.CH <- data.frame(MD = blockmaxxer(CH_Dat,which = "MD",blocks = CH_Dat$X30.mins)$MD,
                    PET = blockmaxxer(CH_Dat,which = "PET",blocks = CH_Dat$X30.mins)$PET,
                    TTC = blockmaxxer(CH_Dat,which = "TTC",blocks = CH_Dat$X30.mins)$TTC,
                    expo = sapply(unique(CH_Dat$X30.mins),
                                  function (x) {dim(subset(CH_Dat,X30.mins==x))[1]})) %>% 
  filter(.,MD>-4,PET>-4,TTC>-4,TTC<0)


# SE data
SE_Dat <- read.table("data/dat_SWE.csv", sep = ",", header = TRUE) 

SE_Dat$N_MD <- -SE_Dat$MD
SE_Dat$N_MDc <- -SE_Dat$MDc
SE_Dat$N_PET <- -SE_Dat$PET
SE_Dat$N_TTC <- -SE_Dat$TTC

BM.SE <- data.frame(MD=blockmaxxer(SE_Dat,which = "N_MD",blocks = SE_Dat$X30_mins)$N_MD,
                   PET=blockmaxxer(SE_Dat,which = "N_PET",blocks = SE_Dat$X30_mins)$N_PET,
                   TTC=blockmaxxer(SE_Dat,which = "N_TTC",blocks = SE_Dat$X30_mins)$N_TTC,
                   expo = sapply(unique(SE_Dat$X30_mins),
                                 function (x) {dim(subset(SE_Dat,X30_mins==x))[1]})) %>%
  filter(.,MD>-5,PET>-5,TTC>-4,TTC<0)


# clean the variable space a bit
rm(list = ls(pattern = "Dat."))
