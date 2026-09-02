# auxiliary functions

# saving plots ------------------------------------------------------------

save_plot <- function(p, filename, ...) {
  if (is.null(filename)) return(invisible(p))
  out_dir <- file.path("../", "plots")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, filename), plot = p, ...)
  invisible(p)
}

# BPOT --------------------------------------------------------------------

run_single_bpot <- function(dat, model, thres, xcrash, speed_ub = 80,
                            estim="joint",likelihood = "censored") {
  if (!estim %in% c("joint", "margins")) stop("estim must be either 'joint','margins'")
  if (!likelihood %in% c("censored", "poisson", "copula")) stop("likelihood must be either 'censored','poisson','copula'")
  
  ss <- seq(thres[2], speed_ub, (speed_ub - thres[2]) / 150)
  m1 <- m2 <- NULL
  if (estim != "margins" && likelihood != "copula"){
    M <- fbvpot(x = dat, model = model, threshold = thres,likelihood = likelihood)
    pcrash <- pevd(xcrash, threshold = thres[1], scale = M$estimate[1],
                   shape = M$estimate[2], lower.tail = FALSE, type = "GP") *M$nat[2] / M$n
  }
  else{
    m1 <- fevd(x = dat[,1], threshold = thres[1],type = "GP")
    pcrash <- pevd(xcrash, threshold = thres[1], scale = m1$results$par[1],
                   shape = m1$results$par[2], lower.tail = FALSE, type = "GP") * m1$rate
    m2 <- fevd(x = dat[,2], threshold = thres[2],type = "GP")
    
    if (likelihood=="copula"){
      dat <- dat[dat[,1] > thres[1] & dat[,2] > thres[2],]
      data_t <-cbind(mtransform_gp_mk2(dat[,1], p = m1$results$par, thres = thres[1],eta=m1$rate,margin = "uniform"),
                     mtransform_gp_mk2(dat[,2], p = m2$results$par, thres = thres[2],eta=m2$rate,margin = "uniform") )
      cop <- switch(model,
                    "log" = fitCopula(gumbelCopula(), data = data_t, method = "ml"),
                    "hr" = fitCopula(huslerReissCopula(), data = data_t, method = "ml"),
                    "neglog" = fitCopula(galambosCopula(), data = data_t, method = "ml"),
                    otherwise = stop("Unsupported model for copula likelihood"))
      dep<- switch(model,
                   "log" = 1/cop@estimate,
                   "hr" = cop@estimate,
                   "neglog" = cop@estimate)
      M <- list(estimate = c(m1$results$par, m2$results$par, dep),
                threshold = thres,
                data = dat,
                model = model,
                nat = c(length(m1$x[m1$x>m1$threshold]), length(m2$x[m2$x>m2$threshold]),nrow(dat)),
                n = length(m1$x))
    }
    
    else{
      # summarise_bpot does not apply to marginal approach, there will be a dimension error
      # CI needs to be computed separately for each margin (use extRemes::ci()), and the dependence parameter (confint(M))
      M <- fbvpot(x = dat, model = model, threshold = thres,scale1 = m1$results$par[1], scale2 = m2$results$par[1],
                   shape1 = m1$results$par[2], shape2 = m2$results$par[2],likelihood = likelihood)
      M$estimate <- c(m1$results$par, m2$results$par, M$estimate)
    }
  }
  plot_df <- create_plot_df(ss, x = xcrash, model = M, px = pcrash)
  list(M = M, pcrash = pcrash, plot_df = plot_df, ss = ss, xcrash = xcrash,
       estim = estim, m1 = m1, m2 = m2)
}

create_plot_df <- function(dat, x, model, px) {
  .args <- c(list(q1 = x, model = model$model, thres = model$threshold,
                  eta = model$nat[1:2] / model$n,
                  mar1 = model$estimate[1:2], mar2 = model$estimate[3:4],
                  tail_type = 4),
             bpot_dep_args(model))
  .args2 <- c(list(q1 = x, model = model$model, thres = model$threshold,
                  eta = model$nat[1:2] / model$n,
                  mar1 = model$estimate[1:2], mar2 = model$estimate[3:4],
                  tail_type = 2),
             bpot_dep_args(model))

  joint_p <- sapply(dat, function(q2) {
    do.call(pb_tvevd, c(.args2, q2 = q2))
  })

  cond_d <- sapply(dat, function(k) {
    do.call(c_bivariate, c(list(y = k, x = x), .args[intersect(
      names(.args), c("model", "dep", "alpha", "beta", "asy", "thres",
                      "eta", "mar1", "mar2"))]))
  }) / px / normalize_c_bivariate(x, model)

  df <- data.frame(speed = dat, JointP = joint_p) %>%
    mutate(ConditionP = JointP / px) %>% na.omit() %>%
    mutate(ConditionalD = cond_d)
  return(df)
}

create_plot_df_np <- function(dat, x, px, mar1, mar2, thres, eta, Ahat) {
  df <- data.frame(speed = dat,
    JointP = sapply(dat, FUN = pb_tv_nonpar, q1 = x, Ahat = Ahat,
                    thres = thres, eta = eta,
                    mar1 = mar1, mar2 = mar2, tail_type = 4)) %>%
    mutate(ConditionP = JointP / px) %>%
    na.omit() %>%
    mutate(ConditionalD = sapply(speed, function(k) {
      c_bivariate_np(y = k, x = x, Ahat = Ahat, thres = thres,
                     eta = eta, mar1 = mar1, mar2 = mar2)}) / px) %>%
    mutate(ConditionalD = ifelse(ConditionalD < 0, 0, ConditionalD))
}


pb_tvevd_sev <- function(M, q1, q2 = 40, tail_type = 2) {
  .args <- c(list(q1 = q1, q2 = q2, model = M$model,
                  thres = M$threshold,
                  eta = M$nat[1] / M$n,
                  mar1 = M$estimate[1:2], mar2 = M$estimate[3:4],
                  tail_type = tail_type),
             bpot_dep_args(M))
  do.call(pb_tvevd, .args)
}

summarise_bpot <- function(results, x0 = NULL, severity_thresholds = 40,
                           severity = NULL,severity_age = NULL, age_mean = 30,
                           labels = NULL) {
  # results are a list of output from run_single_bpot
  n <- length(results)
  if (n == 0) stop("`results' must contain at least one model")

  if (is.null(labels))
    labels <- if (!is.null(names(results))) names(results)
              else paste0("Model ", seq_len(n))
  if (length(labels) != n) stop("`labels' must have length ", n)

  lw <- max(nchar(labels))
  hdr <- paste0(strrep("=", 55), "\n")
  sec <- paste0(strrep("-", 55), "\n")
  
  pnames <- c("sigma1", "xi1", "sigma2", "xi2")
  

  fmt <- function(x, digits = 4) round(x, digits)

  cat(hdr)
  cat("  BPOT Model Summary  (model:", results[[1]]$M$model, ")\n")
  cat(hdr)

  cat("\nThresholds\n", sec)
  for (i in seq_len(n)) {
    M <- results[[i]]$M
    cat(sprintf("  %-*s  TTC = %-8.4f  Speed = %.4f\n",
                lw, labels[i], M$threshold[1], M$threshold[2]))
  }

  cat("\nParameter Estimates\n", sec)
  for (i in seq_len(n)) {
    M <- results[[i]]$M
    if (M$model %in% c("ct", "bilog", "negbilog"))
      pnames <- c(pnames[1:4], "alpha", "beta")
    if (M$model %in% c("alog","aneglog"))
      pnames <- c(pnames[1:4], "asy1", "asy2", "dep")
    else
      pnames <- c(pnames[1:4],"dep")
    cat(sprintf("  %-*s  %s\n", lw, labels[i],
                paste(pnames, fmt(M$estimate), sep = " = ",
                      collapse = "  ")))
    cat(sprintf("  95%% CI (%s):\n", labels[i])); print(confint(M))
  }

  cat("\nModel Fit\n", sec)
  for (i in seq_len(n))
    cat(sprintf("  Deviance  %-*s = %.4f\n", lw, labels[i],
                deviance(results[[i]]$M)))
  for (i in seq_len(n))
    cat(sprintf("  AIC  %-*s = %.4f\n", lw, labels[i],
                AIC(results[[i]]$M)))


  cat("\nCrash Probability  P(TTC < 0)\n", sec)
  for (i in seq_len(n))
    cat(sprintf("  %-*s = %.6f\n", lw, labels[i], results[[i]]$pcrash))

  if (!is.null(x0)) {
    for (i in 1:length(severity_thresholds)){
      cat(sprintf("\nSevere Crash Probability  P(TTC < 0, Speed > %s km/h)\n",severity_thresholds[i]), sec)
      for (j in seq_len(n))
        cat(sprintf("  %-*s = %.6f\n", lw, labels[j],
                    pb_tvevd_sev(results[[j]]$M, q1 = x0,q2=severity_thresholds[i])))
    }
    

    if (!is.null(severity)) {
      cat("\nInjury Probability  (impact speed only)\n", sec)
      for (i in seq_len(n)) {
        ip <- injury_from_c_bivariate(results[[i]]$plot_df$speed,
                                      ev_model = results[[i]]$M,
                                      severity = severity,
                                      x0 = x0, px = results[[i]]$pcrash)
        cat(sprintf("  %-*s = %.6f\n", lw, labels[i], ip))
      }
    }

    if (!is.null(severity_age)) {
      cat(sprintf("\nInjury Probability  (impact speed + age = %g)\n",
                  age_mean), sec)
      for (i in seq_len(n)) {
        ip <- injury_from_c_bivariate_e(results[[i]]$plot_df$speed,
                                        ev_model = results[[i]]$M,
                                        severity = severity_age,
                                        x0 = x0, px = results[[i]]$pcrash,
                                        age_mean = age_mean)
        cat(sprintf("  %-*s = %.6f\n", lw, labels[i], ip))
      }
    }
  }

  cat(hdr)
  invisible(results)
}


# CBPOT -------------------------------------------------------------------

create_plot_df_q <- function(dat, x, model,CM0, px, p2, pu) {
  df <- data.frame(
    speed_u = dat,
    speed   = qgamma(dat, shape = p2$estimate[1], rate = p2$estimate[2]),
    JointP  = sapply(dat, function(y) { pMvdc(c(x, y), CM0) })
  ) %>%
    mutate(ConditionP = JointP / px) %>%
    mutate(ConditionalD = sapply(speed_u, c_q_bivariate, x = x,
                                 model = model) / px * pu *
             dgamma(speed, shape = p2$estimate[1],
                    rate = p2$estimate[2])) %>% na.omit()
  return(df)
}

run_single_cbpot <- function(cop_dat, copula, p2, qcrash, pot,
                             pu = pu_default, ...) {
  Cop <- fitCopula(copula, data = cop_dat, ...)
  CM  <- q_distr_param(Cop@copula, mar1 = pot, mar2 = p2, type = 4)
  CM0 <- q_distr_param(Cop@copula, mar1 = pot, mar2 = p2, type = 2)

  x0_un <- 1 - qcrash

  sq   <- seq(0.5, 60, 0.05)
  s_un <- pgamma(sq, shape = p2$estimate[1], rate = p2$estimate[2])
  plot_df_q <- create_plot_df_q(s_un, x = x0_un, model = Cop@copula,CM0=CM0,
                                px = qcrash, p2 = p2, pu = pu)

  list(Cop = Cop, CM = CM, CM0 = CM0,
       x0_un = x0_un, pcrash = qcrash * pu, pu = pu,
       p2 = p2, pot = pot, s_un = s_un, qcrash = qcrash,
       plot_df_q = plot_df_q)
}

summarise_cbpot <- function(results, x0 = NULL,severity_thresholds=40, 
                            severity = NULL,severity_age = NULL, age_mean = 30,
                            pot = NULL, labels = NULL) {
  n <- length(results)
  if (n == 0) stop("`results' must contain at least one model")

  if (is.null(labels))
    labels <- if (!is.null(names(results))) names(results)
              else paste0("Model ", seq_len(n))
  if (length(labels) != n) stop("`labels' must have length ", n)

  if (!is.null(pot) && length(pot) != n)
    stop("`pot' must have length ", n)

  lw <- max(nchar(labels))
  hdr <- paste0(strrep("=", 55), "\n")
  sec <- paste0(strrep("-", 55), "\n")

  fmt <- function(x, digits = 4) round(x, digits)

  cat(hdr)
  cat(" CBPOT Model Summary  (model:",
      results[[1]]$Cop@copula@fullname, ")\n")
  cat(hdr)

  pnames <- c("sigma", "xi", "gamma", "r", "dep")

  cat("\nThresholds\n", sec)
  for (i in seq_len(n)) {
    th <- results[[i]]$CM@paramMargins %>% unlist() %>% {.["threshold"]}
    cat(sprintf("  %-*s  TTC = %s\n", lw, labels[i], th))
  }

  cat("\nParameter Estimates\n", sec)
  for (i in seq_len(n)) {
    r <- results[[i]]
    param_keep <- r$CM@paramMargins %>% unlist() %>%
      {.[c("scale.scale", "shape.shape", "shape", "rate")]} %>%
      c(., r$CM@copula@parameters) %>% as.numeric()
    cat(sprintf("  %-*s  %s\n", lw, labels[i],
                paste(pnames, fmt(param_keep), sep = " = ",
                      collapse = "  ")))
    cat(sprintf("  95%% CI (%s):\n", labels[i]))
    if (!is.null(pot)) ci(pot[[i]], type = "parameter")
    cat(sprintf("CBPOT %s gamma margin\n", labels[i]))
    confint(r$p2)
  }

  cat("\nDependence parameter\n", sec)
  for (i in seq_len(n)) {
    r <- results[[i]]
    cat(sprintf("  %-*s  lb = %.4f  estim = %.4f  ub = %.4f\n", lw, labels[i],
                r$Cop@estimate - 1.96 * sqrt(r$Cop@var.est),
                r$Cop@estimate,
                r$Cop@estimate + 1.96 * sqrt(r$Cop@var.est)))
  }

  cat("\nModel Fit\n", sec)
  for (i in seq_len(n))
    cat(sprintf("  Deviance  %-*s = %.4f\n", lw, labels[i],
                -2 * results[[i]]$Cop@loglik))

  cat("\nCrash Probability  P(TTC < 0)\n", sec)
  for (i in seq_len(n))
    cat(sprintf("  %-*s = %.6f\n", lw, labels[i], results[[i]]$pcrash))

  if (!is.null(x0)) {
    for (i in 1:length(severity_thresholds)){
      cat(sprintf("\nSevere Crash Probability  Q(TTC < 0, Speed > %s km/h)\n",severity_thresholds[i]), sec)
      for (j in seq_len(n))
        cat(sprintf("  %-*s = %.6f\n", lw, labels[j],
                    pMvdc(c(x0, severity_thresholds[i]),results[[j]]$CM0)))
    }
    
  }

  if (!is.null(severity)) {
    cat("\nInjury Probability  (impact speed only)\n", sec)
    for (i in seq_len(n)) {
      r <- results[[i]]
      ip <- injury_from_c_q_bivariate(r$Cop@copula, severity,
                                      r$x0_un, r$qcrash, r$p2)
      cat(sprintf("  %-*s = %.6f\n", lw, labels[i], ip))
    }
  }

  if (!is.null(severity_age)) {
    cat(sprintf("\nInjury Probability  (impact speed + age = %g)\n",
                age_mean), sec)
    for (i in seq_len(n)) {
      r <- results[[i]]
      ip <- injury_from_c_q_bivariate_e(r$Cop@copula, severity_age,
                                        r$x0_un, r$qcrash, r$p2,
                                        age_mean = age_mean)
      cat(sprintf("  %-*s = %.6f\n", lw, labels[i], ip))
    }
  }

  cat(hdr)
  invisible(results)
}
