# functions for computing injury probability using densities

# BPOT --------------------------------------------------------------------

injury_from_c_bivariate <- function(dat, ev_model, severity, x0, px,
                                    ub = NULL, lb = NULL) {
  f_y <- function(k) {
    temp <- call_c_bivariate(y = k, x = x0, ev_model = ev_model)
    return(temp * severity(k) / px * ev_model$nat[3] / ev_model$n)
  }

  R <- tryCatch({
    integrate(Vectorize(f_y),
              lower = ifelse(is.null(lb), min(dat), lb),
              upper = ifelse(is.null(ub), Inf, ub))$value
  }, error = function(e) {
    if (grepl("non-finite function value|Failed to find a finite upper",
              e$message)) {
      ub_alt <- max(dat, na.rm = TRUE)
      return(integrate(Vectorize(f_y),
                       lower = ifelse(is.null(lb), min(dat), lb),
                       upper = ub_alt)$value)
    } else { stop(e) }
  })

  return(R)
}

injury_from_c_bivariate_e <- function(dat, ev_model, severity, x0, px,
                                      ub = NULL, lb = NULL, N = 500,
                                      age_mean = 40, age_sd = 15) {
  f_y <- function(k) {
    temp <- call_c_bivariate(y = k, x = x0, ev_model = ev_model)
    age <- rnorm(N, mean = age_mean, sd = age_sd)
    age <- age[age > 0]
    return(mean(temp * severity(k, age)))
  }

  R <- tryCatch({
    integrate(Vectorize(f_y),
              lower = ifelse(is.null(lb), min(dat), lb),
              upper = ifelse(is.null(ub), 60, ub))$value
  }, error = function(e) {
    if (grepl("non-finite function value|Failed to find a finite upper",
              e$message)) {
      ub_alt <- max(dat, na.rm = TRUE)
      return(integrate(Vectorize(f_y),
                       lower = ifelse(is.null(lb), min(dat), lb),
                       upper = ub_alt)$value)
    } else { stop(e) }
  })

  return(R / px / normalize_c_bivariate(x0, ev_model))
}


# CBPOT -------------------------------------------------------------------

injury_from_c_q_bivariate <- function(model, severity, x0, px, p2,
                                      lb = NULL, ub = NULL, pu = 1) {
  f_y <- function(k) {
    S <- pgamma(k, shape = p2$estimate[1], rate = p2$estimate[2])
    temp <- c_q_bivariate(y = S, x = x0, model = model) * severity(k) *
            dgamma(k, shape = p2$estimate[1], rate = p2$estimate[2])
    return(temp)
  }

  R <- tryCatch({
    integrate(Vectorize(f_y),
              lower = ifelse(is.null(lb), 0, lb),
              upper = ifelse(is.null(ub), 60, ub))$value
  }, error = function(e) {
    if (grepl("non-finite function value", e$message)) {
      ub_alt <- 55
      return(integrate(Vectorize(f_y),
                       lower = ifelse(is.null(lb), 0, lb),
                       upper = ifelse(is.null(ub), min(ub, ub_alt),
                                      ub_alt))$value)
    } else { stop(e) }
  })
  return(R / px * pu)
}

injury_from_c_q_bivariate_e <- function(model, severity, x0, px, p2,
                                        lb = 0.5, ub = NULL, pu = 1,
                                        N = 1000, age_mean = 40,
                                        age_sd = 10) {
  f_y <- function(k) {
    S <- pgamma(k, shape = p2$estimate[1], rate = p2$estimate[2])
    temp <- c_q_bivariate(y = S, x = x0, model = model) *
            dgamma(k, shape = p2$estimate[1], rate = p2$estimate[2])
    age <- rnorm(N, mean = age_mean, sd = age_sd)
    age <- age[age > 0]
    return(temp * mean(severity(k, age)))
  }

  R <- tryCatch({
    integrate(Vectorize(f_y),
              lower = ifelse(is.null(lb), 0, lb),
              upper = ifelse(is.null(ub), 60, ub))$value
  }, error = function(e) {
    if (grepl("non-finite function value", e$message)) {
      ub_alt <- 55
      return(integrate(Vectorize(f_y),
                       lower = ifelse(is.null(lb), 0, lb),
                       upper = ifelse(is.null(ub), min(ub, ub_alt),
                                      ub_alt))$value)
    } else { stop(e) }
  })
  return(R / px * pu)
}
