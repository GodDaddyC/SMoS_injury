# diagnostics of the dependence structure

# BPOT --------------------------------------------------------------------

tbevd_plot <- function(ev_model, dat, k = 10, filename = NULL) {
  A_np_conf <- pickands_nonpar(dat = dat,
    mar1 = ev_model$estimate[1:2], mar2 = ev_model$estimate[3:4],
    thres = ev_model$threshold,
    eta = ev_model$nat[1:2] / ev_model$n,
    est = "cfg", CI = FALSE, d = 2, k = k, ifplot = FALSE,alpha = 0.05)
  A_np <- pickands_nonpar(dat = dat,
    mar1 = ev_model$estimate[1:2], mar2 = ev_model$estimate[3:4],
    thres = ev_model$threshold,
    eta = ev_model$nat[1:2] / ev_model$n,
    est = "cfg", CI = FALSE, d = 2, k = k, ifplot = FALSE)

  tgrid <-  seq(0, 1, 0.005)
  spec_dens_plot <- data.frame(t = tgrid) %>%
    mutate(h.Nonpar    = sapply(t, a_bp_approx, A_bp = A_np$beta,
                                ord = "2")) %>%
    mutate(h.Nonpar.UP = sapply(t, a_bp_approx,
                                A_bp = A_np_conf$up.beta, ord = "2")) %>%
    mutate(h.Nonpar.LW = sapply(t, a_bp_approx,
                                A_bp = A_np_conf$low.beta, ord = "2"))
  
  abvnonpar(x = tgrid,data=dat,method='pot',k=ev_model$nat[3],plot=TRUE)
  if (ev_model$model %in% c("log", "hr", "neglog")) {
    abvevd(dep = ev_model$estimate[5], model = ev_model$model,add=TRUE, col = "red")
    spec_dens_plot <- spec_dens_plot %>%
      mutate(h.Par    = sapply(t, hbvevd, dep = ev_model$estimate[5],
                               model = ev_model$model, half = TRUE))
  }

  if (ev_model$model %in% c("ct", "bilog", "negbilog")) {
    abvevd(alpha = ev_model$estimate[5], beta = ev_model$estimate[6],
           model = ev_model$model, add=TRUE,col = "red")
    spec_dens_plot <- spec_dens_plot %>%
      mutate(h.Par    = sapply(t, hbvevd,
                               alpha = ev_model$estimate[6],
                               beta  = ev_model$estimate[5],
                               model = ev_model$model, half = TRUE))
  }

  if (ev_model$model %in% c("alog", "aneglog")) {
    abvevd(dep = ev_model$estimate[7], asy = ev_model$estimate[5:6],
           model = ev_model$model,add=TRUE, col = "red")
    spec_dens_plot <- spec_dens_plot %>%
      mutate(h.Par    = sapply(t, hbvevd, dep = ev_model$estimate[7],
                               asy = ev_model$estimate[c(6,5)],
                               model = ev_model$model, half = TRUE))
  }

  legend("bottomright",
         legend = c("Nonparametric estimates", "Parametric estimates"),
         col = c("black", "red"), lty = c(1, 2), bty = "n")
  title(main = sprintf("Dependence diagnostics %s", ev_model$model))

  p <- ggplot(spec_dens_plot, aes(x = t, y = h.Par, color = "Model fitted")) +
    geom_line() +
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

  if (!is.null(filename)) save_plot(p, filename)
  return(p)
}

bpot_ab_values <- function(ev_model, t) {
  model <- ev_model$model
  if (model %in% c("log", "hr", "neglog")) {
    abvevd(x = t, dep = ev_model$estimate[5], model = model, plot = FALSE)
  } else if (model %in% c("ct", "bilog", "negbilog")) {
    abvevd(x = t, alpha = ev_model$estimate[5], beta = ev_model$estimate[6],
           model = model, plot = FALSE)
  } else if (model %in% c("alog", "aneglog")) {
    abvevd(x = t, dep = ev_model$estimate[7], asy = ev_model$estimate[5:6],
           model = model, plot = FALSE)
  } else {
    stop(sprintf("Model '%s' not supported by bpot_ab_values", model))
  }
}

bpot_hb_values <- function(ev_model, t) {
  model <- ev_model$model
  if (model %in% c("log", "hr", "neglog")) {
    hbvevd(x = t, dep = ev_model$estimate[5], model = model,
           half = TRUE, plot = FALSE)
  } else if (model %in% c("ct", "bilog", "negbilog")) {
    hbvevd(x = t, alpha = ev_model$estimate[6], beta = ev_model$estimate[5],
           model = model, half = TRUE, plot = FALSE)
  } else if (model %in% c("alog", "aneglog")) {
    hbvevd(x = t, dep = ev_model$estimate[7], asy = ev_model$estimate[c(6, 5)],
           model = model, half = TRUE, plot = FALSE)
  } else {
    stop(sprintf("Model '%s' not supported by bpot_hb_values", model))
  }
}

tbevd_plot_multi <- function(ev_models, dat, k = 10,
                             type = c("spectral", "pickands"),
                             labels = NULL, filename = NULL, ...) {
  type <- match.arg(type)

  if (!is.list(ev_models) || length(ev_models) == 0)
    stop("`ev_models` must be a non-empty list of fitted BPOT models")
  if (is.null(names(ev_models)))
    names(ev_models) <- paste0("Model ", seq_along(ev_models))
  if (!is.null(labels)) {
    if (length(labels) != length(ev_models))
      stop("`labels` must have the same length as `ev_models`")
    names(ev_models) <- labels
  }

  ref <- ev_models[[1]]
  tgrid <- seq(0, 1, 0.005)
  A_np_conf <- pickands_nonpar(dat = dat,
                               mar1 = ref$estimate[1:2], mar2 = ref$estimate[3:4],
                               thres = ref$threshold,
                               eta = ref$nat[1:2] / ref$n,
                               est = "cfg", CI = TRUE, d = 2, k = k,
                               ifplot = FALSE)
  if (type == "spectral") {
    np_est   <- sapply(tgrid, a_bp_approx, A_bp = (A_np_conf$up.beta+A_np_conf$low.beta)/2, ord = "2")
    param_fn <- bpot_hb_values
    np_band <- data.frame(t = tgrid, 
                          ymin = sapply(tgrid, a_bp_approx,A_bp = A_np_conf$low.beta, ord = "2"), 
                          ymax = sapply(tgrid, a_bp_approx,A_bp = A_np_conf$up.beta, ord = "2"))
    np_line <- data.frame(t = tgrid, y = np_est, Model = "Nonparametric")
    param_df <- do.call(rbind, lapply(seq_along(ev_models), function(i) {
      data.frame(t = tgrid, y = param_fn(ev_models[[i]], tgrid),
                 Model = names(ev_models)[i])
    }))
  } else {
    np_est <- abvnonpar(x = tgrid,data=dat,method='pot',k=ref$nat[3])
    np_band <- NULL
    param_fn <- bpot_ab_values
    np_line <- data.frame(t = tgrid, y = np_est, Model = "Nonparametric")
    param_df <- do.call(rbind, lapply(seq_along(ev_models), function(i) {
      data.frame(t = tgrid, y = param_fn(ev_models[[i]], tgrid),
                 Model = names(ev_models)[i])
    }))
  }

  lvls <- c("Nonparametric", names(ev_models))
  plot_df <- rbind(np_line, param_df)
  plot_df$Model <- factor(plot_df$Model, levels = lvls)
  ltys <- setNames(c(2, rep(1, length(ev_models))), lvls)

  ylab <- if (type == "spectral") "Spectral density" else "Dependence function A(t)"

  p <- ggplot(plot_df, aes(x = t, y = y, color = Model, linetype = Model)) +
    geom_line(na.rm = TRUE) +
    scale_linetype_manual(values = ltys) +
    labs(x = "t", y = ylab, color = "", linetype = "") +
    theme_minimal()

  if (!is.null(np_band)) {
    p <- p + geom_ribbon(data = np_band,
                         aes(x = t, ymin = ymin, ymax = ymax),
                         inherit.aes = FALSE, alpha = 0.15, fill = "grey50")
  }

  if (!is.null(filename)) save_plot(p, filename)
  return(p)
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
