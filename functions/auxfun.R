# auxiliary functions

# BPOT --------------------------------------------------------------------

run_single_bpot <- function(dat, model, thres, xcrash, speed_ub = 60) {
  M <- fbvpot(x = dat, model = model, threshold = thres)
  pcrash <- pevd(xcrash, threshold = thres[1], scale = M$estimate[1],
                 shape = M$estimate[2], lower.tail = FALSE, type = "GP") *
            M$nat[2] / M$n
  ss <- seq(thres[2], speed_ub, (speed_ub - thres[2]) / 150)
  plot_df <- create_plot_df(ss, x = xcrash, model = M, px = pcrash)
  list(M = M, pcrash = pcrash, plot_df = plot_df, ss = ss, xcrash = xcrash)
}

create_result_bpot <- function(r1, r2) {
  list(M_1 = r1$M, M_2 = r2$M,
       pcrash_1 = r1$pcrash, pcrash_2 = r2$pcrash,
       plot_df_1 = r1$plot_df, plot_df_2 = r2$plot_df)
}

create_plot_df <- function(dat, x, model, px) {
  mod <- model$model
  .args <- list(q1 = x, model = mod, thres = model$threshold,
                eta = model$nat[1:2] / model$n,
                mar1 = model$estimate[1:2], mar2 = model$estimate[3:4],
                tail_type = 4)
  if (mod == "ct") {
    .args$alpha <- model$estimate[5]
    .args$beta  <- model$estimate[6]
  } else {
    .args$dep <- model$estimate[5]
  }

  joint_p <- sapply(dat, function(q2) {
    do.call(pb_tvevd, c(.args, q2 = q2))
  })

  cond_d <- sapply(dat, function(k) {
    do.call(c_bivariate, c(list(y = k, x = x), .args[intersect(
      names(.args), c("model", "dep", "alpha", "beta", "thres", "eta",
                      "mar1", "mar2"))]))
  }) / px / normalize_c_bivariate(x, model) * model$nat[3] / model$n

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


summarise_bpot <- function(result, x0 = c(0, 0), severity = NULL,
                           severity_age = NULL, age_mean = 60) {
  M1 <- result$M_1
  M2 <- result$M_2
  hdr <- paste0(strrep("=", 55), "\n")
  sec <- paste0(strrep("-", 55), "\n")

  fmt <- function(x, digits = 4) round(x, digits)

  cat(hdr)
  cat("  BPOT Model Summary  (model:", M1$model, ")\n")
  cat(hdr)

  cat("\nThresholds\n", sec)
  cat(sprintf("  %-6s  TTC = %-8.4f  Speed = %.4f\n",
              "CN:", M1$threshold[1], M1$threshold[2]))
  cat(sprintf("  %-6s  TTC = %-8.4f  Speed = %.4f\n",
              "SE:", M2$threshold[1], M2$threshold[2]))

  pnames <- c("sigma1", "xi1", "sigma2", "xi2", "dep")
  if (M1$model == "ct") pnames <- c(pnames[1:4], "alpha", "beta")

  cat("\nParameter Estimates\n", sec)
  cat("  CN:", paste(pnames, fmt(M1$estimate), sep = " = ",
                     collapse = "  "), "\n")
  cat("  95% CI (CN):\n"); print(confint(M1))
  cat("  SE:", paste(pnames, fmt(M2$estimate), sep = " = ",
                     collapse = "  "), "\n")
  cat("  95% CI (SE):\n"); print(confint(M2))

  cat("\nModel Fit\n", sec)
  cat(sprintf("  Deviance  CN = %.4f\n", deviance(M1)))
  cat(sprintf("  Deviance  SE = %.4f\n", deviance(M2)))

  cat("\nCrash Probability  P(TTC < u)\n", sec)
  cat(sprintf("  CN = %.6f\n", result$pcrash_1))
  cat(sprintf("  SE = %.6f\n", result$pcrash_2))

  if (!is.null(x0)) {
    cat("\nSevere Crash Probability  P(TTC < u, Speed > 40 km/h)\n", sec)
    if (M1$model != "ct") {
      p_sev_1 <- pb_tvevd(q1 = x0[1], q2 = 40, dep = M1$estimate["dep"],
                          thres = M1$threshold, model = M1$model,
                          eta = M1$nat[1] / M1$n,
                          mar1 = M1$estimate[1:2],
                          mar2 = M1$estimate[3:4], tail_type = 2)
      p_sev_2 <- pb_tvevd(q1 = x0[2], q2 = 40, dep = M2$estimate["dep"],
                          thres = M2$threshold, model = M2$model,
                          eta = M2$nat[1] / M2$n,
                          mar1 = M2$estimate[1:2],
                          mar2 = M2$estimate[3:4], tail_type = 2)
    } else {
      p_sev_1 <- pb_tvevd(q1 = x0[1], q2 = 40, alpha = M1$estimate["alpha"],
                          beta = M1$estimate["beta"],
                          thres = M1$threshold, model = M1$model,
                          eta = M1$nat[1] / M1$n,
                          mar1 = M1$estimate[1:2],
                          mar2 = M1$estimate[3:4], tail_type = 2)
      p_sev_2 <- pb_tvevd(q1 = x0[2], q2 = 40, alpha = M2$estimate["alpha"],
                          beta = M2$estimate["beta"],
                          thres = M2$threshold, model = M2$model,
                          eta = M2$nat[1] / M2$n,
                          mar1 = M2$estimate[1:2],
                          mar2 = M2$estimate[3:4], tail_type = 2)
    }

    cat(sprintf("  CN = %.6f\n", p_sev_1))
    cat(sprintf("  SE = %.6f\n", p_sev_2))

    if (!is.null(severity)) {
      cat("\nInjury Probability  (impact speed only)\n", sec)
      ip_1 <- injury_from_c_bivariate(result$plot_df_1$speed,
                                      ev_model = M1, severity = severity,
                                      x0 = x0[1], px = result$pcrash_1)
      ip_2 <- injury_from_c_bivariate(result$plot_df_2$speed,
                                      ev_model = M2, severity = severity,
                                      x0 = x0[2], px = result$pcrash_2)
      cat(sprintf("  CN = %.6f\n", ip_1))
      cat(sprintf("  SE = %.6f\n", ip_2))
    }

    if (!is.null(severity_age)) {
      cat(sprintf("\nInjury Probability  (impact speed + age = %g)\n",
                  age_mean), sec)
      ip_1a <- injury_from_c_bivariate_e(result$plot_df_1$speed,
                                         ev_model = M1,
                                         severity = severity_age,
                                         x0 = x0[1], px = result$pcrash_1,
                                         age_mean = age_mean)
      ip_2a <- injury_from_c_bivariate_e(result$plot_df_2$speed,
                                         ev_model = M2,
                                         severity = severity_age,
                                         x0 = x0[2], px = result$pcrash_2,
                                         age_mean = age_mean)
      cat(sprintf("  CN = %.6f\n", ip_1a))
      cat(sprintf("  SE = %.6f\n", ip_2a))
    }
  }

  cat(hdr)
  invisible(result)
}


# CBPOT -------------------------------------------------------------------

create_plot_df_q <- function(dat, x, model, px, p2, pu) {
  df <- data.frame(
    speed_u = dat,
    speed   = qgamma(dat, shape = p2$estimate[1], rate = p2$estimate[2]),
    JointP  = sapply(dat, function(y) { pCopula(c(x, y), model) })
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
  plot_df_q <- create_plot_df_q(s_un, x = x0_un, model = Cop@copula,
                                px = qcrash, p2 = p2, pu = pu)

  list(Cop = Cop, CM = CM, CM0 = CM0,
       x0_un = x0_un, pcrash = qcrash * pu, pu = pu,
       p2 = p2, s_un = s_un,
       plot_df_q = plot_df_q)
}

create_result_cbpot <- function(r1, r2) {
  list(CM_1 = r1$CM, CM_2 = r2$CM,
       Cop_1 = r1$Cop, Cop_2 = r2$Cop,
       x0_un_1 = r1$x0_un, x0_un_2 = r2$x0_un,
       pcrash_1 = r1$pcrash, pcrash_2 = r2$pcrash,
       plot_df_q_1 = r1$plot_df_q, plot_df_q_2 = r2$plot_df_q,
       CM0_1 = r1$CM0, CM0_2 = r2$CM0)
}

summarise_cbpot <- function(result, x0 = c(0, 0), severity = NULL,
                            severity_age = NULL, age_mean = 60,
                            pot_1, pot_2, conseq_1, conseq_2) {
  M1 <- result$CM_1; M2 <- result$CM_2
  C1 <- result$Cop_1; C2 <- result$Cop_2
  hdr <- paste0(strrep("=", 55), "\n")
  sec <- paste0(strrep("-", 55), "\n")

  fmt <- function(x, digits = 4) round(x, digits)

  cat(hdr)
  cat(" CBPOT Model Summary  (model:", C1@copula@fullname, ")\n")
  cat(hdr)

  cat("\nThresholds\n", sec)
  cat(sprintf("  %-6s  TTC = %s", "CN:",
              M1@paramMargins %>% unlist() %>% {.["threshold"]}))
  cat(sprintf("  %-6s  TTC = %s", "SE:",
              M2@paramMargins %>% unlist() %>% {.["threshold"]}))

  pnames <- c("sigma", "xi", "gamma", "r", "dep")
  param_keep <- M1@paramMargins %>% unlist() %>%
    {.[c("scale.scale", "shape.shape", "shape", "rate")]} %>%
    c(., M1@copula@parameters) %>% as.numeric()

  cat("\nParameter Estimates\n", sec)
  cat("  CN:", paste(pnames, fmt(param_keep), sep = " = ",
                     collapse = "  "), "\n")
  cat("  95% CI (CN):\n")
  ci(pot_1, type = "parameter")
  cat("CBPOT CN gamma margin")
  confint(conseq_1)
  cat("  SE:", paste(pnames, fmt(param_keep), sep = " = ",
                     collapse = "  "), "\n")
  cat("  95% CI (SE):\n")
  ci(pot_2, type = "parameter")
  cat("CBPOT SE gamma margin")
  confint(conseq_2)

  cat("\ Dependence parameter \n", sec)
  cat("CN:", "lb", C1@estimate - 1.96 * sqrt(C1@var.est),
      "estim", C1@estimate, "ub", C1@estimate + 1.96 * sqrt(C1@var.est), "\n")
  cat("SE:", "lb", C2@estimate - 1.96 * sqrt(C2@var.est),
      "estim", C2@estimate, "ub", C2@estimate + 1.96 * sqrt(C2@var.est), "\n")

  cat(sprintf("  Deviance  CN = %.4f\n", -2 * C1@loglik))
  cat(sprintf("  Deviance  SE = %.4f\n", -2 * C2@loglik))

  cat("\nCrash Probability  P(TTC < u)\n", sec)
  cat(sprintf("  CN = %.6f\n", result$pcrash_1))
  cat(sprintf("  SE = %.6f\n", result$pcrash_2))

  if (!is.null(x0)) {
    cat("\nSevere Crash Probability  P(TTC < u, Speed > 40 km/h)\n", sec)
    p_sev_1 <- pMvdc(c(x0[1], 40), result$CM0_1)
    p_sev_2 <- pMvdc(c(x0[2], 40), result$CM0_2)
    cat(sprintf("  CN = %.6f\n", p_sev_1))
    cat(sprintf("  SE = %.6f\n", p_sev_2))

    if (!is.null(severity)) {
      cat("\nInjury Probability  (impact speed only)\n", sec)
      ip_1 <- injury_from_c_q_bivariate(C1@copula, severity,
                                        result$x0_un_1, qcrash_1, conseq_1)
      ip_2 <- injury_from_c_q_bivariate(C2@copula, severity,
                                        result$x0_un_2, qcrash_2, conseq_2)
      cat(sprintf("  CN = %.6f\n", ip_1))
      cat(sprintf("  SE = %.6f\n", ip_2))
    }

    if (!is.null(severity_age)) {
      cat(sprintf("\nInjury Probability  (impact speed + age = %g)\n",
                  age_mean), sec)
      ip_1a <- injury_from_c_q_bivariate_e(C1@copula, severity_age,
                                           result$x0_un_1, qcrash_1, conseq_1,
                                           age_mean = age_mean)
      ip_2a <- injury_from_c_q_bivariate_e(C2@copula, severity_age,
                                           result$x0_un_2, qcrash_2, conseq_2,
                                           age_mean = age_mean)
      cat(sprintf("  CN = %.6f\n", ip_1a))
      cat(sprintf("  SE = %.6f\n", ip_2a))
    }
  }

  cat(hdr)
  invisible(result)
}
