# Select copula models

C.gumbel <- gumbelCopula(1.1)
contourplot2(C.gumbel, dCopula, nlevels = 20,
             main = "Theoretical Gumbel Copula")
contourplot2(rotCopula(C.gumbel, flip = c(TRUE, FALSE)), dCopula,
             nlevels = 20,
             main = "Theoretical Gumbel Copula (Rot 270)")
contourplot2(rotCopula(C.gumbel, flip = c(FALSE, TRUE)), dCopula,
             nlevels = 20,
             main = "Theoretical Gumbel Copula (Rot 90)")
persp(C.gumbel, dCopula, zlim = c(0, 3), main = "Gumbel Copula Density")

C.hr <- huslerReissCopula(1)
contourplot2(C.hr, dCopula, nlevels = 20,
             main = "Theoretical Husler-Reiss Copula")
contourplot2(rotCopula(C.hr, flip = c(TRUE, FALSE)), dCopula,
             nlevels = 20,
             main = "Theoretical Husler-Reiss Copula (Rot 270)")
contourplot2(rotCopula(C.hr, flip = c(FALSE, TRUE)), dCopula,
             nlevels = 20,
             main = "Theoretical Husler-Reiss Copula (Rot 90)")
persp(C.hr, dCopula, zlim = c(0, 5),
      main = "Husler-Reiss Copula Density")

C.clay <- claytonCopula(3)
plot(rCopula(500, C.clay))
contourplot2(C.clay, dCopula, nlevels = 20,
             main = "dCopula(<rotCopula>)")
contourplot2(rotCopula(C.clay, flip = c(TRUE, FALSE)), dCopula,
             nlevels = 20, main = "Clayton Copula (Rot 270)")
contourplot2(rotCopula(C.clay, flip = c(FALSE, TRUE)), dCopula,
             nlevels = 20, main = "Clayton Copula (Rot 90)")
persp(C.clay, dCopula, zlim = c(0, 10),
      main = "Clayton Copula Density")


# Empirical copula
emp.CN <- empCopula(cop_dat_1, smoothing = "beta")
emp.SE <- empCopula(cop_dat_2, smoothing = "beta")
contourplot2(emp.CN, dCopula, nlevels = 20,
             main = "Empirical Copula CN")
contourplot2(emp.SE, dCopula, nlevels = 20,
             main = "Empirical Copula SE")
persp(emp.CN, dCopula, zlim = c(0, 5),
      main = "Empirical copula density CN")
persp(emp.SE, dCopula, zlim = c(0, 5),
      main = "Empirical copula density SE")

plot(kdecop(cop_dat_1), type = "surface", margin = "unif",
     main = "Empirical Copula Density CN")
plot(kdecop(cop_dat_2), type = "surface", margin = "unif",
     main = "Empirical Copula Density SE")
