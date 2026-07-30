ci_delta <- function(x, GP, alpha) {
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
      exp(-(x - u) / sig) - 1.96 * sqrt(t(grad) %*% V %*% grad),
      exp(-(x - u) / sig),
      exp(-(x - u) / sig) + 1.96 * sqrt(t(grad) %*% V %*% grad)))
  } else {
    stop("Unsupported GP$type. Please specify either 'GP' or 'Exponential'.")
  }
}

ci_prolik <- function(x0, GP, alpha = 0.95, expand_factor = 2,
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
    q <- l_max - 0.5 * qchisq(alpha, 1)
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

ci_bayes <- function(x0, GP, alpha, if_plot = FALSE) {
  pp <- Tail.prob.GP(x0, GP)
  pp_density <- density(pp)
  MAP <- pp_density$x[which.max(pp_density$y)]
  MP <- quantile(pp, 0.5)
  BCI <- quantile(pp, probs = c(1 - alpha, alpha))
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
