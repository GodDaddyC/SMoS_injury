# diagnostics of the dependence structure

# BPOT --------------------------------------------------------------------

tbevd_plot <- function(ev_model, dat, k = 10) {
  A_np_conf <- pickands_nonpar(dat = dat,
    mar1 = ev_model$estimate[1:2], mar2 = ev_model$estimate[3:4],
    thres = ev_model$threshold,
    eta = ev_model$nat[1:2] / ev_model$n,
    est = "cfg", CI = TRUE, d = 2, k = k, ifplot = TRUE)
  A_np <- pickands_nonpar(dat = dat,
    mar1 = ev_model$estimate[1:2], mar2 = ev_model$estimate[3:4],
    thres = ev_model$threshold,
    eta = ev_model$nat[1:2] / ev_model$n,
    est = "cfg", CI = FALSE, d = 2, k = k, ifplot = FALSE)

  spec_dens_plot <- data.frame(t = seq(0, 1, 0.005)) %>%
    mutate(h.Nonpar    = sapply(t, a_bp_approx, A_bp = A_np$beta,
                                ord = "2") / 2) %>%
    mutate(h.Nonpar.UP = sapply(t, a_bp_approx,
                                A_bp = A_np_conf$up.beta, ord = "2") / 2) %>%
    mutate(h.Nonpar.LW = sapply(t, a_bp_approx,
                                A_bp = A_np_conf$low.beta, ord = "2") / 2)

  if (ev_model$model %in% c("log", "hr", "neglog")) {
    CB <- confint(ev_model, parm = "dep")
    abvevd(dep = ev_model$estimate[5], model = ev_model$model,
           plot = TRUE, add = TRUE, col = "red")
    abvevd(dep = CB[1], model = ev_model$model,
           plot = TRUE, add = TRUE, col = "red", lty = 5)
    abvevd(dep = CB[2], model = ev_model$model,
           plot = TRUE, add = TRUE, col = "red", lty = 5)

    spec_dens_plot <- spec_dens_plot %>%
      mutate(h.Par    = sapply(t, hbvevd, dep = ev_model$estimate[5],
                               model = ev_model$model, half = TRUE)) %>%
      mutate(h.Par.UP = sapply(t, hbvevd, dep = CB[2],
                               model = ev_model$model, half = TRUE)) %>%
      mutate(h.Par.LW = sapply(t, hbvevd, dep = CB[1],
                               model = ev_model$model, half = TRUE))
  }

  if (ev_model$model %in% c("ct", "bilog", "negbilog")) {
    CB_alpha <- confint(ev_model, parm = "beta")
    CB_beta  <- confint(ev_model, parm = "alpha")
    abvevd(alpha = ev_model$estimate[6], beta = ev_model$estimate[5],
           model = ev_model$model, plot = TRUE, add = TRUE, col = "red")
    abvevd(alpha = CB_alpha[1], beta = CB_beta[1],
           model = ev_model$model, plot = TRUE, add = TRUE,
           col = "red", lty = 5)
    abvevd(alpha = CB_alpha[2], beta = CB_beta[2],
           model = ev_model$model, plot = TRUE, add = TRUE,
           col = "red", lty = 5)

    spec_dens_plot <- spec_dens_plot %>%
      mutate(h.Par    = sapply(t, hbvevd,
                               alpha = ev_model$estimate[5],
                               beta  = ev_model$estimate[6],
                               model = ev_model$model, half = TRUE)) %>%
      mutate(h.Par.UP = sapply(t, hbvevd,
                               alpha = CB_alpha[2], beta = CB_beta[2],
                               model = ev_model$model, half = TRUE)) %>%
      mutate(h.Par.LW = sapply(t, hbvevd,
                               alpha = CB_alpha[1], beta = CB_beta[1],
                               model = ev_model$model, half = TRUE))
  }

  if (ev_model$model %in% c("alog", "aneglog")) {
    CB_asy1 <- confint(ev_model, parm = "asy1")
    CB_asy2 <- confint(ev_model, parm = "asy2")
    CB <- confint(ev_model, parm = "dep")
    abvevd(dep = ev_model$estimate[7], asy = ev_model$estimate[5:6],
           model = ev_model$model, plot = TRUE, add = TRUE, col = "red")
    abvevd(dep = CB[1], asy = c(CB_asy1[1], CB_asy2[1]),
           model = ev_model$model, plot = TRUE, add = TRUE, col = "red",
           lty = 5)
    abvevd(dep = CB[2], asy = c(CB_asy1[1], CB_asy2[1]),
           model = ev_model$model, plot = TRUE, add = TRUE, col = "red",
           lty = 5)

    spec_dens_plot <- spec_dens_plot %>%
      mutate(h.Par    = sapply(t, hbvevd, dep = ev_model$estimate[7],
                               asy = ev_model$estimate[5:6],
                               model = ev_model$model, half = TRUE)) %>%
      mutate(h.Par.UP = sapply(t, hbvevd, dep = CB[2],
                               asy = c(CB_asy1[2], CB_asy2[2]),
                               model = ev_model$model, half = TRUE)) %>%
      mutate(h.Par.LW = sapply(t, hbvevd, dep = CB[1],
                               asy = c(CB_asy1[1], CB_asy2[1]),
                               model = ev_model$model, half = TRUE))
  }

  legend("bottomright",
         legend = c("Nonparametric estimates", "Parametric estimates"),
         col = c("black", "red"), lty = c(1, 2), bty = "n")
  title(main = sprintf("Dependence diagnostics %s", ev_model$model))

  ggplot(spec_dens_plot, aes(x = t, y = h.Par, color = "Model fitted")) +
    geom_line() +
    geom_line(aes(y = h.Par.LW, color = "Model fitted"), linetype = "dashed") +
    geom_line(aes(y = h.Par.UP, color = "Model fitted"), linetype = "dashed") +
    geom_line(aes(y = h.Nonpar.LW, color = "Nonparametric"),
              linetype = "dashed") +
    geom_line(aes(y = h.Nonpar.UP, color = "Nonparametric"),
              linetype = "dashed") +
    geom_line(aes(y = h.Nonpar, color = "Nonparametric")) +
    labs(x = "t", y = "Spectral density",
         title = sprintf("Spectral density %s", ev_model$model)) +
    scale_color_manual(name = "",
      values = c("Model fitted" = "red", "Nonparametric" = "black")) +
    theme_minimal()
}

gof_bpot <- function(ev_model, dat) {
  n <- nrow(dat)
  dat <- dat %>% filter(prox > ev_model$threshold[1],
                        Speed > ev_model$threshold[2])
  U1 <- mtransform_gp_mk2(dat[, 1], p = ev_model$estimate[1:2],
                          thres = ev_model$threshold[1],
                          eta = ev_model$nat[1] / ev_model$n,
                          margin = "uniform")
  U2 <- mtransform_gp_mk2(dat[, 2], p = ev_model$estimate[3:4],
                          thres = ev_model$threshold[2],
                          eta = ev_model$nat[2] / ev_model$n,
                          margin = "uniform")
  model <- switch(ev_model$model,
    log = gumbelCopula(dim = 2, param = 1 / ev_model$estimate[5]),
    hr  = huslerReissCopula(param = ev_model$estimate[5]))
  Cop <- fitCopula(copula = model, data = cbind(U1, U2), method = "mpl")
  Result <- gofEVCopula(copula = Cop@copula, method = "itau",
                        x = cbind(U1, U2), N = 1000, estimator = "CFG")
  return(Result)
}


# CBPOT -------------------------------------------------------------------

kc_test <- function(dat, Copula, n = 1e4, m = 200, p_val = 0.05) {
  Cn <- empCopula(X = dat, ties.method = "random")
  Sim <- rCopula(n, copula = Cn)
  Model <- rCopula(n, copula = Copula)

  temp1 <- pCopula(Sim, copula = Cn)
  temp2 <- pCopula(Model, copula = Copula)

  Kc_emp   <- sapply(seq(0, 1, 1 / m), function(k) mean(temp1 < k))
  Kc_model <- sapply(seq(0, 1, 1 / m), function(k) mean(temp2 < k))

  result <- Desc::CCC(Kc_emp, Kc_model, ci = "z-transform",
                      conf.level = 1 - p_val)
  return(result$rho.c)
}

kc_test_two <- function(dat, Copula1, Copula2, n = 1e4, m = 200,
                        p_val = 0.05) {
  CI1 <- kc_test(dat, Copula1, n, m, p_val = p_val)
  CI2 <- kc_test(dat, Copula2, n, m, p_val = p_val)
  CI1$upr.ci
  CI2$lwr.ci
  ifelse(CI1$upr.ci > CI2$lwr.ci,
    paste("No significant difference between Copula1 and Copula2 fits at level",
          p_val),
    paste("Copula2 is significantly better than Copula1 at level", p_val))
}
