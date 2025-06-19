pkgs <- c("readxl", "dplyr", "ggplot2","ggpubr","extRemes","evmix",
          "evd", "copula", "RColorBrewer")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}


library(readxl)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(extRemes)
library(evd)
library(copula)