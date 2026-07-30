# crash severity densities derived from different approaches

# shared helpers ---------------------------------------------------------------

integrate_safe <- function(f, lower, upper, ..., rel.tol = 1e-3) {
  tryCatch({
    integrate(f, lower = lower, upper = upper, ..., rel.tol = rel.tol)$value
  }, error = function(e) {
    if (grepl("non-finite|Failed to find a finite upper", e$message)) {
      NA
    } else {
      stop(e)
    }
  })
}

call_c_bivariate <- function(y, x, ev_model) {
  .args <- list(y = y, x = x, model = ev_model$model,
                thres = ev_model$threshold,
                eta   = ev_model$nat[1:2] / ev_model$n,
                mar1  = ev_model$estimate[1:2],
                mar2  = ev_model$estimate[3:4])
  if (ev_model$model == "ct") {
    .args$alpha <- ev_model$estimate[5]
    .args$beta  <- ev_model$estimate[6]
  } else {
    .args$dep <- ev_model$estimate[5]
  }
  do.call(c_bivariate, .args)
}


# BPOT --------------------------------------------------------------------

c_bivariate <- function(y, x, model, dep, alpha, beta, thres, eta, mar1, mar2) {
  integrand <- function(q1, ...) {
    result <- numeric(length(q1))
    for (i in seq_along(q1)) {
      result[i] <- switch(model,
        log = db_tvevd(q1 = q1[i], q2 = y, model = model, dep = dep,
                       thres = thres, eta = eta, mar1 = mar1, mar2 = mar2),
        hr  = db_tvevd(q1 = q1[i], q2 = y, model = model, dep = dep,
                       thres = thres, eta = eta, mar1 = mar1, mar2 = mar2),
        ct  = db_tvevd(q1 = q1[i], q2 = y, model = model, alpha = alpha,
                       beta = beta, thres = thres, eta = eta,
                       mar1 = mar1, mar2 = mar2)
      )
    }
    return(result)
  }

  R <- tryCatch({
    integrate(integrand, lower = x, upper = Inf, rel.tol = 1e-3)$value
  }, error = function(e) {
    if (grepl("non-finite function value", e$message)) {
      if (mar1[2] < 0) {
        ub_alt <- thres[1] - mar1[1] / mar1[2]
        return(integrate(integrand, lower = x, upper = ub_alt,
                         rel.tol = 1e-3)$value)
      } else {
        for (i in rev(seq(0.05, 5, 0.05))) {
          ub_alt <- x + i * mar1[1]
          result <- tryCatch({
            integrate(integrand, lower = x, upper = ub_alt,
                      rel.tol = 1e-3)$value
          }, error = function(e_retry) { return(NA) })
          if (!is.na(result)) return(result)
        }
        stop("Failed to find a finite upper bound after multiple retries.",
             call. = FALSE)
      }
    } else { stop(e) }
  })

  return(R)
}

c_bivariate_np <- function(y, x, thres, eta, mar1, mar2, Ahat) {
  integrand <- function(q1, ...) {
    result <- numeric(length(q1))
    for (i in seq_along(q1)) {
      result[i] <- db_tv_nonpar(q1 = q1[i], q2 = y, Ahat = Ahat,
                                mar1 = mar1, mar2 = mar2,
                                thres = thres, eta = eta)
    }
    return(result)
  }

  R <- tryCatch({
    integrate(integrand, lower = x, upper = Inf, rel.tol = 1e-3)$value
  }, error = function(e) {
    if (grepl("non-finite function value", e$message)) {
      if (mar1[2] < 0) {
        for (j in rev(seq(x, thres[1] - mar1[1] / mar1[2],
                          length.out = 10))) {
          ub_alt <- j
          result <- tryCatch({
            integrate(integrand, lower = x, upper = ub_alt,
                      rel.tol = 1e-3)$value
          }, error = function(e_retry) { return(NA) })
          if (!is.na(result)) return(result)
        }
      } else {
        for (i in rev(seq(0.05, 5, 0.05))) {
          ub_alt <- x + i * mar1[1]
          result <- tryCatch({
            integrate(integrand, lower = x, upper = ub_alt,
                      rel.tol = 1e-3)$value
          }, error = function(e_retry) { return(NA) })
          if (!is.na(result)) return(result)
        }
        stop("Failed to find a finite upper bound after multiple retries.",
             call. = FALSE)
      }
    } else { stop(e) }
  })

  return(R)
}

normalize_c_bivariate <- function(x, ev_model) {
  dat <- ev_model$data[ev_model$data[, 2] > ev_model$threshold[2], 2]
  c_y <- function(k) { call_c_bivariate(y = k, x = x, ev_model = ev_model) }

  C <- tryCatch({
    integrate(Vectorize(c_y), lower = min(dat), upper = Inf)$value
  }, error = function(e) {
    if (grepl("non-finite function value|Failed to find a finite upper",
              e$message)) {
      ub_alt <- ifelse(ev_model$estimate[4] < 0,
                       ev_model$threshold[2] -
                         ev_model$estimate[3] / ev_model$estimate[4],
                       55)
      return(integrate(Vectorize(c_y), lower = min(dat),
                       upper = ub_alt)$value)
    } else { stop(e) }
  })

  nc <- switch(ev_model$model,
    log = pb_tvevd(q1 = x, q2 = ev_model$threshold[2],
                   dep = ev_model$estimate[5], thres = ev_model$threshold,
                   model = ev_model$model,
                   eta = ev_model$nat[1:2] / ev_model$n,
                   mar1 = ev_model$estimate[1:2],
                   mar2 = ev_model$estimate[3:4], tail_type = 2),
    hr  = pb_tvevd(q1 = x, q2 = ev_model$threshold[2],
                   dep = ev_model$estimate[5], thres = ev_model$threshold,
                   model = ev_model$model,
                   eta = ev_model$nat[1:2] / ev_model$n,
                   mar1 = ev_model$estimate[1:2],
                   mar2 = ev_model$estimate[3:4], tail_type = 2),
    ct  = pb_tvevd(q1 = x, q2 = ev_model$threshold[2],
                   alpha = ev_model$estimate[5],
                   beta  = ev_model$estimate[6],
                   thres = ev_model$threshold, model = ev_model$model,
                   eta = ev_model$nat[1:2] / ev_model$n,
                   mar1 = ev_model$estimate[1:2],
                   mar2 = ev_model$estimate[3:4], tail_type = 2))
  return(C / nc)
}


# CBPOT -------------------------------------------------------------------

c_q_bivariate <- function(y, x, model) {
  integrand <- function(q1) { dCopula(c(q1, y), model) }
  R <- tryCatch({
    integrate(Vectorize(integrand), lower = x, upper = 1)$value
  }, error = function(e) {
    if (grepl("non-finite function value", e$message)) {
      ub_alt <- 1 - 1e-7
      return(integrate(integrand, lower = x, upper = ub_alt)$value)
    } else { stop(e) }
  })
  return(R)
}

normalize_c_q_bivariate <- function(x, model, lb = 0) {
  c_y <- function(k) { c_q_bivariate(y = k, x = x, model = model) }

  C <- tryCatch({
    integrate(Vectorize(c_y), lower = lb, upper = 1)$value
  }, error = function(e) {
    if (grepl("non-finite function value", e$message)) {
      ub_alt <- 1 - 1e-7
      return(integrate(Vectorize(c_y), lower = lb, upper = ub_alt)$value)
    }
    if (grepl("maximum number of subdivisions reached", e$message)) {
      return(integrate(Vectorize(c_y), lower = lb, upper = Inf,
                       subdivisions = 200, rel.tol = 1e-5)$value)
    } else { stop(e) }
  })
  return(C)
}


# plotting ----------------------------------------------------------------

plot_crash_severity <- function(plot_df_1, plot_df_2, injury_df,
                                v1 = NULL, v2 = NULL,
                                label1 = "CN", label2 = "SE",
                                legend_title = "Site",
                                sec_axis_label = "Injury prob",
                                save_path = NULL) {
  scale_factor <- max(plot_df_1$ConditionalD) / max(injury_df$InjuryP)
  col1 <- "red"
  col2 <- "blue"

  p <- ggplot(plot_df_1, aes(x = speed, y = ConditionalD)) +
    geom_line(aes(colour = label1)) +
    geom_line(data = plot_df_2, aes(x = speed, y = ConditionalD,
                                    colour = label2)) +
    geom_line(data = injury_df, aes(x = speed, y = InjuryP * scale_factor))

  if (!is.null(v1)) {
    p <- p +
      geom_vline(xintercept = v1, linetype = "dashed", color = col1) +
      annotate("text", x = v1 + 2, y = 0.02,
               label = paste("u=", round(v1, 3)), color = col1)
  }
  if (!is.null(v2)) {
    p <- p +
      geom_vline(xintercept = v2, linetype = "dashed", color = col2) +
      annotate("text", x = v2 + 2, y = 0.015,
               label = paste("u=", round(v2, 3)), color = col2)
  }

  col_values <- setNames(c(col1, col2), c(label1, label2))

  p <- p +
    scale_y_continuous(
      name = "density",
      sec.axis = sec_axis(~ . / scale_factor, name = sec_axis_label)
    ) +
    scale_colour_manual(name = legend_title, values = col_values) +
    labs(x = "y", y = "f(y|crash)") +
    theme(
      panel.grid.major   = element_line(colour = "gray91"),
      panel.grid.minor   = element_line(colour = "gray88"),
      panel.background   = element_rect(fill = "white", colour = "white",
                                        linetype = "solid"),
      plot.background    = element_rect(linetype = "solid")
    )

  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, dpi = 300, device = "png")
  }
  return(p)
}


plot_crash_severity_single <- function(plot_df, injury_df,
                                       v = NULL, label = NULL,
                                       show_legend = FALSE,
                                       sec_axis_label = "Injury prob",
                                       save_path = NULL) {
  scale_factor <- max(plot_df$ConditionalD) / max(injury_df$InjuryP)

  if (!is.null(label)) {
    p <- ggplot(plot_df, aes(x = speed, y = ConditionalD, colour = label))
  } else {
    p <- ggplot(plot_df, aes(x = speed, y = ConditionalD))
  }

  p <- p +
    geom_line(show.legend = show_legend) +
    geom_line(data = injury_df,
              aes(x = speed, y = InjuryP * scale_factor),
              colour = "black", inherit.aes = FALSE)

  if (!is.null(v)) {
    p <- p +
      geom_vline(xintercept = v, linetype = "dashed", color = "red") +
      annotate("text", x = v + 2, y = 0.02,
               label = paste("u=", round(v, 3)), color = "red")
  }

  p <- p +
    scale_y_continuous(
      name = "density",
      sec.axis = sec_axis(~ . / scale_factor, name = sec_axis_label)
    ) +
    labs(x = "Speed (km/h)", y = "f(y|TTC<0)") +
    theme(
      panel.grid.major   = element_line(colour = "gray91"),
      panel.grid.minor   = element_line(colour = "gray88"),
      panel.background   = element_rect(fill = "white", colour = "white",
                                        linetype = "solid"),
      plot.background    = element_rect(linetype = "solid")
    )

  if (!is.null(save_path)) {
    ggsave(save_path, plot = p)
  }
  return(p)
}
