pkgs <- c("readxl", "dplyr", "ggplot2", "ggpubr", "extRemes", "evmix",
          "purrr", "tidyr", "evd", "eva", "copula", "RColorBrewer",
          "fitdistrplus", "VC2copula", "kdecopula", "VineCopula",
          "ExtremalDep", "DescTools")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

lapply(list.files("../functions", pattern = "\\.R$", full.names = TRUE), source)

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
