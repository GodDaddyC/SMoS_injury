Tail.prob.GP <- function(x,model,lower.tail=FALSE,conditional = FALSE,use.phi=FALSE){
  # model is a fevd object, requires extRemes it assumes use.phi=TRUE, if not, 
  # get rid of the exp() in sigma
  if (conditional){
    if (model$method =='Bayesian'){
      prob.chain <- sapply(1:dim(model$chain.info)[1],function(k){
        pevd(x,threshold = model$threshold, scale = exp(model$result[k,1]), 
             shape = model$result[k,2],lower.tail = lower.tail,type = "GP")
      })
      return(prob.chain)
    }
    if (model$method == "MLE"){
      return(ifelse(use.phi,pevd(x,threshold = model$threshold, scale = exp(model$result$par[1]), 
                                 shape = model$result$par[2],lower.tail = lower.tail,type = "GP"),
                    pevd(x,threshold = model$threshold, scale = model$result$par[1], 
                         shape = model$result$par[2],lower.tail = lower.tail,type = "GP")))
    }
  }
  else{
    if (model$method =='Bayesian'){
      prob.chain <- sum(model$x>model$threshold)/model$n *
        sapply(1:dim(model$chain.info)[1],function(k){
          pevd(x,threshold = model$threshold, scale = exp(model$result[k,1]), 
               shape = model$result[k,2],lower.tail = lower.tail,type = "GP")})
      return(prob.chain)
    }
    if (model$method == "MLE"){
      pu <- sum(model$x>model$threshold)/model$n
      return(ifelse(use.phi,pu*pevd(x,threshold = model$threshold, scale = exp(model$result$par[1]), 
                                    shape = model$result$par[2],lower.tail = lower.tail,type = "GP"),
                    pu*pevd(x,threshold = model$threshold, scale = model$result$par[1], 
                            shape = model$result$par[2],lower.tail = lower.tail,type = "GP")))
    }
  }
}


# CI crash probability delta method
ci_delta <- function(x, GP, alpha=0.05) {
  u <- GP$threshold
  sig <- exp(GP$results$par[1])
  pu <- sum(GP$x > GP$threshold) / GP$n
  if (GP$type == "GP") {
    gamma <- GP$results$par[2]
    aa <- 1 + gamma * (x - u) / sig
    prob <- Tail.prob.GP(x, model = GP)
    d_sig <- pu * (x - u) / sig^2 * aa^(-1 / gamma - 1)
    d_gamma <- prob * (log(aa) / gamma^2 - x / (gamma * sig * aa))
    grad <- c(d_sig, d_gamma)
    V <- solve(GP$results$hessian, diag(c(1, 1)))
    return(c(prob - qnorm(1 - alpha / 2) * sqrt(t(grad) %*% V %*% grad),
             prob,
             prob + qnorm(1 - alpha / 2) * sqrt(t(grad) %*% V %*% grad))) }
  else if (GP$type == "Exponential") {
    d_sig <- -exp(-(x - u) / sig) * (x - u) / sig^2
    grad <- c(d_sig)
    V <- solve(GP$results$hessian, diag(1))
    return(c(
      exp(-(x - u) / sig) - qnorm(1 - alpha / 2) * sqrt(t(grad) %*% V %*% grad),
      exp(-(x - u) / sig),
      exp(-(x - u) / sig) + qnorm(1 - alpha / 2) * sqrt(t(grad) %*% V %*% grad)))
  } else {
    stop("Unsupported GP$type. Please specify either 'GP' or 'Exponential'.")
  }
}

# CI crash probability profile likelihood method
ci_prolik <- function(x0, GP, alpha = 0.05, expand_factor = 2,
                      nint = 1000, if_plot = FALSE, use_phi = FALSE) {
  u <- GP$threshold
  data <- GP$x[GP$x > u]
  p0 <- Tail.prob.GP(x0, model = GP, use.phi = use_phi)
  p_guess <- c(p0 * 1e-2, p0 * 10)
  if (p0 == 0) { stop("target probability is 0") }
  max_iter <- 8
  iter <- 0
  repeat {
    p <- seq(p_guess[1], p_guess[2],
             length = ifelse(p_guess[2] / p_guess[1] < 1e-4, nint,
                             10 * nint))

    theta <- GP$results$par[2]
    pu <- sum(GP$x > GP$threshold) / GP$n

    l_p_gamma <- function(theta, x0, data, pu, p0, u) {
      sig <- theta * (x0 - u) / ((p0 / pu)^(-theta) - 1)
      y <- 1 + theta * (data - u) / sig
      if (any(y <= 0)) { return(1e-6) }
      else {
        return(ifelse(use_phi,
          length(data - u) * (log(pu) - sig) -
            (theta + 1) / theta * sum(log(y)),
          length(data - u) * (log(pu) - log(sig)) -
            (theta + 1) / theta * sum(log(y))))
      }
    }

    l <- sapply(p, function(p0) {
      output <- optim(theta, fn = l_p_gamma, method = "BFGS",
                      p0 = p0, x0 = x0, data = data, u = u, pu = pu,
                      control = list(fnscale = -5))
      output$value
    })
    p_mle <- p[which.max(l)]
    l_max <- max(l)
    q <- l_max - 0.5 * qchisq(1-alpha, 1)
    p_lb <- l[l >= q][1]
    p_ub <- l[l <= q & p > p_mle][1]
    ci_lb <- p[match(p_lb, l)]
    ci_ub <- p[match(p_ub, l)]
    if (!any(is.na(c(ci_lb, ci_ub))) || iter >= max_iter) { break }
    if (is.na(ci_lb) || ci_lb == p_guess[1]) {
      p_guess[1] <- p_guess[1] / expand_factor
    }
    if (is.na(ci_ub) || ci_ub == p_guess[2]) {
      p_guess[2] <- p_guess[2] * expand_factor
    }
    iter <- iter + 1
  }

  if (if_plot) {
    plot(p, l, type = "l", xlab = "p", ylab = "Profile Log-likelihood")
    abline(h = q, col = "red")
    abline(v = p_mle, col = 4)
    abline(v = ci_lb, col = "red")
    abline(v = ci_ub, col = "red")
  }
  return(c(ci_lb, p_mle, ci_ub))
}

# CI crash probability from posterior distribution
ci_bayes <- function(x0, GP, alpha=0/05, if_plot = FALSE) {
  pp <- Tail.prob.GP(x0, GP)
  pp_density <- density(pp)
  MAP <- pp_density$x[which.max(pp_density$y)]
  MP <- quantile(pp, 0.5)
  BCI <- quantile(pp, probs = c(alpha,1 - alpha))
  if (if_plot) {
    plot(pp_density, main = sprintf("Posterior distribution of P(X>%f)", x0),
         xlab = "p")
    abline(v = MAP, col = "red")
    abline(v = MP, col = "green")
    abline(v = BCI[1], col = "blue")
    abline(v = BCI[2], col = "blue")
  }
  return(c(BCI[1], MP, MAP, BCI[2]))
}


# CI crash probability parametric bootstrap
ci_boot <- function(x0,GP,alpha=0.05,B=1000){
  u <- GP$threshold
  
  param_est <- GP$results$par
  V <- solve(GP$results$hessian, diag(c(1, 1)))
  V <- (V + t(V)) / 2
  R <- tryCatch(chol(V), error = function(e) {
    tryCatch(chol(V + diag(1e-8, 2)), error = function(e2) NULL)
  })
  if (is.null(R))
    stop("Covariance matrix is not positive definite")
  
  boots <- numeric(B)
  i <- 0L
  attempts <- 0L
  max_attempts <- 100L * B
  while (i < B && attempts < max_attempts) {
    attempts <- attempts + 1L
    param_b <- as.vector(param_est + t(R) %*% rnorm(2))
    
    GP_temp <- GP
    GP_temp$results$par <- param_b
    pb <- tryCatch(Tail.prob.GP(x0,GP_temp),
                   error = function(e) NA_real_)
    if (!is.finite(pb)) next
    
    i <- i + 1L
    boots[i] <- pb
  }
  
  if (i < B)
    warning(sprintf("only %d of %d bootstrap replicates were valid", i, B))
  if (i < 10)
    stop("too few valid bootstrap replicates to form a confidence interval")
  
  boots <- boots[seq_len(i)]
  probs <- c(alpha / 2, 1 - alpha / 2)
  ci <- quantile(boots, probs = probs, na.rm = TRUE)
  out <- c(lower = as.numeric(ci[1]),
           estimate = Tail.prob.GP(x0,GP),
           upper = as.numeric(ci[2]))
  out
}

ci_boot2 <- function(x0, GP, alpha = 0.05, B = 1000) {
  
  u <- GP$threshold
  
  boots <- numeric(B)
  i <- 0L
  attempts <- 0L
  max_attempts <- 100L * B
  
  while (i < B && attempts < max_attempts) {
    
    attempts <- attempts + 1L
    
    # Parametric bootstrap sample
    boot_sample <- revd(
      floor(GP$rate * dim(GP$cov.data)[1]),
      threshold = u,
      scale = GP$results$par[1],
      shape = GP$results$par[2],
      type = "GP"
    )
    
    # Re-estimate GP parameters
    GP_temp <- tryCatch(
      fevd(
        x = boot_sample,
        threshold = u,
        type = "GP"
      ),
      error = function(e) NULL
    )
    
    if (is.null(GP_temp))
      next
    
    # Recalculate target probability
    pb <- tryCatch(
      Tail.prob.GP(x0, GP_temp),
      error = function(e) NA_real_
    )
    
    if (!is.finite(pb))
      next
    
    i <- i + 1L
    boots[i] <- pb
  }
  
  if (i < B)
    warning(sprintf(
      "Only %d of %d bootstrap replicates were valid",
      i, B
    ))
  
  if (i < 10)
    stop("Too few valid bootstrap replicates.")
  
  boots <- boots[seq_len(i)]
  
  # Correct percentile confidence interval
  probs <- c(alpha / 2, 1 - alpha / 2)
  
  ci <- quantile(
    boots,
    probs = probs,
    na.rm = TRUE
  )
  
  c(
    lower = as.numeric(ci[1]),
    estimate = Tail.prob.GP(x0, GP),
    upper = as.numeric(ci[2])
  )
}

# bivariate POT — parametric bootstrap CI ----------------------------------

valid_bpot_params <- function(theta, model) {
  scale1 <- theta[1]; scale2 <- theta[3]
  if (!is.finite(scale1) || !is.finite(scale2) || scale1 <= 0 || scale2 <= 0)
    return(FALSE)

  if (model %in% c("ct", "bilog", "negbilog")) {
    return(is.finite(theta[5]) && is.finite(theta[6]) &&
           theta[5] > 0 && theta[6] > 0)
  }

  if (model %in% c("alog", "aneglog")) {
    dep <- theta[7]
    asy <- theta[5:6]
    valid_dep <- is.finite(dep) && dep > 0 &&
                 (!model %in% "alog" || dep <= 1)
    return(valid_dep && all(is.finite(asy)) &&
           all(asy >= 0 & asy <= 1))
  }

  dep <- theta[5]
  if (!is.finite(dep) || dep <= 0) return(FALSE)
  if (model %in% c("log") && dep > 1) return(FALSE)
  TRUE
}

bpot_prob <- function(ev_model, q1, q2, tail_type) {
  .args <- c(list(q1 = q1, q2 = q2, model = ev_model$model,
                  thres = ev_model$threshold,
                  eta = ev_model$nat[1:2] / ev_model$n,
                  mar1 = ev_model$estimate[1:2],
                  mar2 = ev_model$estimate[3:4],
                  tail_type = tail_type),
             bpot_dep_args(ev_model))
  do.call(pb_tvevd, .args)
}

# parametric bootstrap for joint approach BPOT
ci_bpot_boot_joint <- function(ev_model, q1, q2, tail_type = 2, B = 1000,
                         alpha = 0.05, seed = 1098) {
  if (!inherits(ev_model, "bvpot"))
    stop("`ev_model` must be a fitted bivariate POT model (from `fbvpot`)")

  supported <- c("log", "alog", "hr", "neglog", "aneglog", "bilog",
                 "negbilog", "ct")
  if (!ev_model$model %in% supported)
    stop(sprintf("Model '%s' not supported; use one of: %s",
                 ev_model$model, paste(supported, collapse = ", ")))

  
  V <- ev_model$var.cov
  param_init <- length(ev_model$estimate)-dim(V)[1]+1
  theta_hat <- ev_model$estimate[param_init:length(ev_model$estimate)]
  if (is.null(V) || any(!is.finite(V)))
    stop("Covariance matrix not available; refit `fbvpot` with `std.err = TRUE`")

  p <- dim(V)[1]
  if (!is.null(seed)) set.seed(seed)

  p_hat <- bpot_prob(ev_model, q1 = q1, q2 = q2, tail_type = tail_type)

  V <- (V + t(V)) / 2
  R <- tryCatch(chol(V), error = function(e) {
    tryCatch(chol(V + diag(1e-8, p)), error = function(e2) NULL)
  })
  if (is.null(R))
    stop("Covariance matrix is not positive definite")

  boots <- numeric(B)
  i <- 0L
  attempts <- 0L
  max_attempts <- 100L * B
  while (i < B && attempts < max_attempts) {
    attempts <- attempts + 1L
    theta_b <- as.vector(theta_hat + t(R) %*% rnorm(p))
    if (!valid_bpot_params(theta_b, ev_model$model)) next

    m_b <- ev_model
    m_b$estimate <- theta_b

    pb <- tryCatch(bpot_prob(m_b, q1 = q1, q2 = q2, tail_type = tail_type),
                   error = function(e) NA_real_)
    if (!is.finite(pb)) next

    i <- i + 1L
    boots[i] <- pb
  }

  if (i < B)
    warning(sprintf("only %d of %d bootstrap replicates were valid", i, B))
  if (i < 10)
    stop("too few valid bootstrap replicates to form a confidence interval")

  boots <- boots[seq_len(i)]
  probs <- c(alpha / 2, 1 -  alpha / 2)
  ci <- quantile(boots, probs = probs, na.rm = TRUE)

  out <- c(lower = as.numeric(ci[1]),
           estimate = as.numeric(p_hat),
           upper = as.numeric(ci[2]))
  attr(out, "bootstrap") <- boots
  attr(out, "alpha") <- alpha
  attr(out, "tail_type") <- tail_type
  attr(out, "n_valid") <- i
  class(out) <- c("ci_bpot_boot", "numeric")
  out
}

# parametric bootstrap for marginal approach BPOT
ci_bpot_boot_marginal <- function(result, q1, q2, tail_type = 2, B = 1000,
                                  alpha = 0.05, seed = 1098) {
  if (!is.list(result) || is.null(result$M) || is.null(result$m1) ||
      is.null(result$m2))
    stop("`result` must be the output of `run_single_bpot(estim = 'margins')`")

  ev_model <- result$M
  mar1 <- result$m1
  mar2 <- result$m2

  supported <- c("log", "alog", "hr", "neglog", "aneglog", "bilog",
                 "negbilog", "ct")
  if (!ev_model$model %in% supported)
    stop(sprintf("Model '%s' not supported; use one of: %s",
                 ev_model$model, paste(supported, collapse = ", ")))

  # --- marginal parameters: separate covariance per margin ---
  theta_mar1 <- mar1$results$par
  theta_mar2 <- mar2$results$par
  V_mar1 <- tryCatch(solve(mar1$results$hessian), error = function(e) NULL)
  V_mar2 <- tryCatch(solve(mar2$results$hessian), error = function(e) NULL)
  if (is.null(V_mar1) || is.null(V_mar2))
    stop("Marginal covariance matrices not available (check `mar1`/`mar2` hessian)")

  # --- dependence parameters: separate covariance ---
  V_dep <- ev_model$var.cov
  if (is.null(V_dep) || any(!is.finite(V_dep)))
    stop("Dependence covariance matrix not available; refit `fbvpot` with `std.err = TRUE`")
  n_dep <- ncol(V_dep)
  start_dep <- length(ev_model$estimate) - n_dep + 1
  theta_dep <- ev_model$estimate[start_dep:length(ev_model$estimate)]

  chol_psd <- function(V) {
    V <- (V + t(V)) / 2
    tryCatch(chol(V), error = function(e) {
      tryCatch(chol(V + diag(1e-8, ncol(V))), error = function(e2) NULL)
    })
  }
  R1   <- chol_psd(V_mar1)
  R2   <- chol_psd(V_mar2)
  Rdep <- chol_psd(V_dep)
  if (is.null(R1) || is.null(R2) || is.null(Rdep))
    stop("One of the covariance matrices is not positive definite")

  if (!is.null(seed)) set.seed(seed)

  p_hat <- bpot_prob(ev_model, q1 = q1, q2 = q2, tail_type = tail_type)

  p1 <- length(theta_mar1)
  p2 <- length(theta_mar2)

  boots <- numeric(B)
  i <- 0L
  attempts <- 0L
  max_attempts <- 100L * B
  while (i < B && attempts < max_attempts) {
    attempts <- attempts + 1L
    theta_b <- c(
      as.vector(theta_mar1 + t(R1) %*% rnorm(p1)),
      as.vector(theta_mar2 + t(R2) %*% rnorm(p2)),
      as.vector(theta_dep  + t(Rdep) %*% rnorm(n_dep))
    )
    if (!valid_bpot_params(theta_b, ev_model$model)) next

    m_b <- ev_model
    m_b$estimate <- theta_b

    pb <- tryCatch(bpot_prob(m_b, q1 = q1, q2 = q2, tail_type = tail_type),
                   error = function(e) NA_real_)
    if (!is.finite(pb)) next

    i <- i + 1L
    boots[i] <- pb
  }

  if (i < B)
    warning(sprintf("only %d of %d bootstrap replicates were valid", i, B))
  if (i < 10)
    stop("too few valid bootstrap replicates to form a confidence interval")

  boots <- boots[seq_len(i)]
  probs <- c(alpha / 2, 1 - alpha / 2)
  ci <- quantile(boots, probs = probs, na.rm = TRUE)

  out <- c(lower = as.numeric(ci[1]),
           estimate = as.numeric(p_hat),
           upper = as.numeric(ci[2]))
  attr(out, "bootstrap") <- boots
  attr(out, "alpha") <- alpha
  attr(out, "tail_type") <- tail_type
  attr(out, "n_valid") <- i
  class(out) <- c("ci_bpot_boot", "numeric")
  out
}

print.ci_bpot_boot <- function(x, ...) {
  cat("Parametric bootstrap CI for bivariate probability\n")
  cat(sprintf("  lower    = %.6g\n", x[["lower"]]))
  cat(sprintf("  estimate = %.6g\n", x[["estimate"]]))
  cat(sprintf("  upper    = %.6g\n", x[["upper"]]))
  cat(sprintf("  (alpha = %g, tail_type = %d, %d valid replicates)\n",
              attr(x, "alpha"), attr(x, "tail_type"), attr(x, "n_valid")))
  invisible(x)
}


# box plot of probability estimates with confidence intervals across models
plot_prob_ci_box <- function(ci_list, labels = NULL,
                             xlab = "Model", ylab = "Probability",
                             title = NULL, filename = NULL, ...) {
  if (!is.list(ci_list) || length(ci_list) == 0)
    stop("`ci_list` must be a non-empty list of CI results")

  if (is.null(names(ci_list)))
    names(ci_list) <- paste0("Model ", seq_along(ci_list))
  if (!is.null(labels)) {
    if (length(labels) != length(ci_list))
      stop("`labels` must have the same length as `ci_list`")
    names(ci_list) <- labels
  }

  get_val <- function(nm) vapply(ci_list, function(x) {
    v <- x[[nm]]
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, numeric(1))

  df <- data.frame(
    model    = factor(names(ci_list), levels = names(ci_list)),
    estimate = get_val("estimate"),
    lower    = get_val("lower"),
    upper    = get_val("upper")
  )

  p <- ggplot(df, aes(x = model, y = estimate)) +
    geom_crossbar(aes(ymin = lower, ymax = upper), width = 0.25,
                  fill = "grey85", colour = "black") +
    geom_point(size = 2.5) +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal()

  if (!is.null(filename)) save_plot(p, filename)
  p
}


# box plot of probability CIs for several estimation methods on one figure
plot_prob_ci_box_grouped <- function(ci_groups, labels = NULL,
                                     xlab = "Model", ylab = "Probability",
                                     title = NULL, legend_title = "Method",
                                     dodge_width = 0.8,
                                     filename = NULL, ...) {
  if (!is.list(ci_groups) || length(ci_groups) == 0)
    stop("`ci_groups` must be a non-empty list of CI-result lists")
  if (any(!vapply(ci_groups, is.list, logical(1))))
    stop("each element of `ci_groups` must be a list of CI results")

  if (is.null(names(ci_groups)))
    names(ci_groups) <- paste0("Method ", seq_along(ci_groups))
  if (!is.null(labels)) {
    if (length(labels) != length(ci_groups))
      stop("`labels` must have the same length as `ci_groups`")
    names(ci_groups) <- labels
  }

  get_val <- function(ci, nm) {
    v <- ci[[nm]]
    if (is.null(v)) NA_real_ else as.numeric(v)
  }

  rows <- do.call(rbind, lapply(seq_along(ci_groups), function(g) {
    grp <- ci_groups[[g]]
    gname <- names(ci_groups)[g]
    if (is.null(names(grp)))
      names(grp) <- paste0("Model ", seq_along(grp))
    do.call(rbind, lapply(seq_along(grp), function(k) {
      ci <- grp[[k]]
      data.frame(model    = names(grp)[k],
                 method   = gname,
                 estimate = get_val(ci, "estimate"),
                 lower    = get_val(ci, "lower"),
                 upper    = get_val(ci, "upper"))
    }))
  }))

  model_lvls <- unique(rows$model)
  rows$model  <- factor(rows$model, levels = model_lvls)
  rows$method <- factor(rows$method, levels = names(ci_groups))

  n_grp <- length(ci_groups)
  dodge <- position_dodge(width = dodge_width)

  p <- ggplot(rows, aes(x = model, y = estimate,
                        colour = method, fill = method)) +
    geom_crossbar(aes(ymin = lower, ymax = upper),
                  width = dodge_width / n_grp * 0.9,
                  position = dodge, alpha = 0.35, linewidth = 0.4) +
    geom_point(position = dodge, size = 2.5) +
    labs(x = xlab, y = ylab, title = title,
         colour = legend_title, fill = legend_title) +
    theme_minimal()

  if (!is.null(filename)) save_plot(p, filename)
  p
}


# CBPOT parametric bootstrap CI -------------------------------------------

cbpot_copula_rebuild <- function(cop, theta) {
  cls <- class(cop)[1]
  switch(cls,
    gumbelCopula      = gumbelCopula(param = theta),
    huslerReissCopula = huslerReissCopula(param = theta),
    galambosCopula    = galambosCopula(param = theta),
    frankCopula       = frankCopula(param=theta),
    claytonCopula     = claytonCopula(param = theta),
    normalCopula      = normalCopula(param = theta),
    stop(sprintf("Unsupported copula family '%s' in CBPOT bootstrap", cls)))
}

cbpot_sev_prob <- function(result, x0, y0,
                           cop_theta = NULL, gp_par = NULL, gam_par = NULL) {
  if (is.null(cop_theta)) cop_theta <- result$Cop@estimate
  if (is.null(gp_par))    gp_par    <- result$pot$results$par
  if (is.null(gam_par))   gam_par   <- result$p2$estimate

  S1 <- pevd(x0, threshold = result$pot$threshold, scale = gp_par[1],
             shape = gp_par[2], type = "GP",)
  S2 <- pgamma(y0, shape = gam_par[1], rate = gam_par[2])
  cop <- cbpot_copula_rebuild(result$Cop@copula, cop_theta)
  (1- S1- S2 + pCopula(c(S1, S2), cop) )* result$pu
}

ci_cbpot_boot <- function(result, x0, y0, B = 1000, alpha = 0.05,
                          seed = NULL) {
  if (!is.list(result) || is.null(result$Cop) || is.null(result$pot) ||
      is.null(result$p2))
    stop("`result` must be the output of `run_single_cbpot()`")

  # --- copula (dependence) block: parameters handled separately ---
  theta_cop <- result$Cop@estimate
  V_cop <- result$Cop@var.est
  if (length(V_cop) != length(theta_cop) || any(!is.finite(V_cop)))
    stop("Copula variance (`Cop@var.est`) not available for bootstrap")

  # --- GP marginal block ---
  gp_par <- result$pot$results$par
  V_gp <- tryCatch(solve(result$pot$results$hessian), error = function(e) NULL)
  if (is.null(V_gp))
    stop("GP marginal covariance (`pot$results$hessian`) not available")

  # --- gamma marginal block ---
  gam_par <- result$p2$estimate
  V_gam <- result$p2$vcov
  if (is.null(V_gam))
    stop("Gamma marginal covariance (`p2$vcov`) not available")

  chol_psd <- function(V) {
    V <- (V + t(V)) / 2
    tryCatch(chol(V), error = function(e) {
      tryCatch(chol(V + diag(1e-8, ncol(V))), error = function(e2) NULL)
    })
  }
  R_cop <- chol_psd(V_cop)
  R_gp  <- chol_psd(V_gp)
  R_gam <- chol_psd(V_gam)
  if (is.null(R_cop) || is.null(R_gp) || is.null(R_gam))
    stop("One of the covariance matrices is not positive definite")

  if (!is.null(seed)) set.seed(seed)

  p_hat <- cbpot_sev_prob(result, x0 = x0, y0 = y0)

  nc <- length(theta_cop)
  np <- length(gp_par)
  ng <- length(gam_par)

  boots <- numeric(B)
  i <- 0L
  attempts <- 0L
  max_attempts <- 100L * B
  while (i < B && attempts < max_attempts) {
    attempts <- attempts + 1L

    theta_cop_b <- as.vector(theta_cop + t(R_cop) %*% rnorm(nc))
    cop_b <- tryCatch(cbpot_copula_rebuild(result$Cop@copula, theta_cop_b),
                      error = function(e) NULL)
    if (is.null(cop_b)) next

    gp_b <- as.vector(gp_par + t(R_gp) %*% rnorm(np))
    if (!is.finite(gp_b[1]) || gp_b[1] <= 0) next

    gam_b <- as.vector(gam_par + t(R_gam) %*% rnorm(ng))
    if (!is.finite(gam_b[1]) || !is.finite(gam_b[2]) ||
        gam_b[1] <= 0 || gam_b[2] <= 0) next

    pb <- tryCatch(cbpot_sev_prob(result, x0 = x0, y0 = y0,
                                  cop_theta = theta_cop_b,
                                  gp_par = gp_b, gam_par = gam_b),
                   error = function(e) NA_real_)
    if (!is.finite(pb)) next

    i <- i + 1L
    boots[i] <- pb
  }

  if (i < B)
    warning(sprintf("only %d of %d bootstrap replicates were valid", i, B))
  if (i < 10)
    stop("too few valid bootstrap replicates to form a confidence interval")

  boots <- boots[seq_len(i)]
  probs <- c(alpha / 2, 1 - alpha / 2)
  ci <- quantile(boots, probs = probs, na.rm = TRUE)

  out <- c(lower = as.numeric(ci[1]),
           estimate = as.numeric(p_hat),
           upper = as.numeric(ci[2]))
  attr(out, "bootstrap") <- boots
  attr(out, "alpha") <- alpha
  attr(out, "n_valid") <- i
  class(out) <- c("ci_bpot_boot", "numeric")
  out
}
