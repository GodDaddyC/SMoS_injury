pkgs <- c("readxl", "dplyr", "ggplot2","ggpubr","extRemes","evmix","purrr","tidyr",
          "evd","eva", "copula", "RColorBrewer","fitdistrplus","VC2copula","kdecopula","VineCopula",
          "ExtremalDep","DescTools")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

