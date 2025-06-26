pkgs <- c("readxl", "dplyr", "ggplot2","ggpubr","extRemes","evmix",
          "evd", "copula", "RColorBrewer")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(pkg, character.only = TRUE)
}

