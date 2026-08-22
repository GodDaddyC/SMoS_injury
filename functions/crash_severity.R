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

bpot_dep_args <- function(ev_model) {
  if (ev_model$model %in% c("alog", "aneglog")) {
    list(dep = ev_model$estimate[7], asy = ev_model$estimate[5:6])
  } else if (ev_model$model %in% c("ct", "bilog", "negbilog")) {
    list(alpha = ev_model$estimate[5], beta = ev_model$estimate[6])
  } else {
    list(dep = ev_model$estimate[5])
  }
}

integrate_retry <- function(integrand, x, mar1, thres) {
  ub <- ifelse(mar1[2] < 0, thres[1] - mar1[1] / mar1[2], Inf)
  tryCatch({
    integrate(integrand, lower = x, upper = ub, rel.tol = 1e-3)$value
  }, error = function(e) {
    if (grepl("non-finite function value", e$message)) {
      for (i in rev(seq(0.05, 5, 0.05))) {
        ub_alt <- x + i * mar1[1]
        if (ub_alt > ub) next
        result <- tryCatch({
          integrate(integrand, lower = x, upper = ub_alt,
                    rel.tol = 1e-3)$value
        }, error = function(e_retry) return(NA))
        if (!is.na(result)) return(result)
      }
      stop("Failed to find a finite upper bound after multiple retries.",
           call. = FALSE)
    } else stop(e)
  })
}

call_c_bivariate <- function(y, x, ev_model) {
  .args <- c(list(y = y, x = x, model = ev_model$model,
                  thres = ev_model$threshold,
                  eta   = ev_model$nat[1:2] / ev_model$n,
                  mar1  = ev_model$estimate[1:2],
                  mar2  = ev_model$estimate[3:4]),
             bpot_dep_args(ev_model))
  do.call(c_bivariate, .args)
}


# BPOT --------------------------------------------------------------------

c_bivariate <- function(y, x, model, dep, alpha, beta, asy, thres, eta,
                        mar1, mar2) {
  if (model %in% c("alog", "aneglog")) {
    extra <- list(dep = dep, asy = asy)
  } else if (model %in% c("ct", "bilog", "negbilog")) {
    extra <- list(alpha = alpha, beta = beta)
  } else {
    extra <- list(dep = dep)
  }
  base <- list(q2 = y, model = model, thres = thres, eta = eta,
               mar1 = mar1, mar2 = mar2)

  integrand <- function(q1, ...) {
    vapply(q1, function(q) {
      do.call(db_tvevd, c(list(q1 = q), base, extra))
    }, numeric(1))
  }

  integrate_retry(integrand, x, mar1, thres)
}

c_bivariate_np <- function(y, x, thres, eta, mar1, mar2, Ahat) {
  integrand <- function(q1, ...) {
    vapply(q1, function(q) {
      db_tv_nonpar(q1 = q, q2 = y, Ahat = Ahat, mar1 = mar1, mar2 = mar2,
                   thres = thres, eta = eta)
    }, numeric(1))
  }

  integrate_retry(integrand, x, mar1, thres)
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

  nc <- do.call(pb_tvevd, c(list(q1 = x, q2 = ev_model$threshold[2],
                                  model = ev_model$model,
                                  thres = ev_model$threshold,
                                  eta = ev_model$nat[1:2] / ev_model$n,
                                  mar1 = ev_model$estimate[1:2],
                                  mar2 = ev_model$estimate[3:4],
                                  tail_type = 2),
                             bpot_dep_args(ev_model)))
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

plot_crash_severity <- function(plot_dfs, injury_df, v = NULL, labels = NULL,
                                legend_title = NULL,
                                sec_axis_label = "Injury prob",
                                save_path = NULL) {
  n <- length(plot_dfs)
  if (n == 0) stop("`plot_dfs' must contain at least one data frame")

  if (is.null(labels))
    labels <- if (!is.null(names(plot_dfs))) names(plot_dfs)
              else paste0("Model ", seq_len(n))
  if (length(labels) != n) stop("`labels' must have length ", n)
  if (!is.null(v) && length(v) != n) stop("`v' must have length ", n)

  plot_df <- do.call(rbind, Map(function(d, lbl) {
    data.frame(speed = d$speed, ConditionalD = d$ConditionalD, model = lbl)
  }, plot_dfs, labels))

  scale_factor <- max(vapply(plot_dfs, function(d) max(d$ConditionalD),
                             numeric(1))) / max(injury_df$InjuryP)

  p <- ggplot(plot_df, aes(x = speed, y = ConditionalD, colour = model)) +
    geom_line() +
    geom_line(data = injury_df,
              aes(x = speed, y = InjuryP * scale_factor),
              colour = "black", inherit.aes = FALSE)

  if (!is.null(v)) {
    vline_df <- data.frame(x = v, model = labels)
    p <- p +
      geom_vline(data = vline_df, aes(xintercept = x, colour = model),
                 linetype = "dashed")
  }

  p <- p +
    scale_y_continuous(
      name = "density",
      sec.axis = sec_axis(~ . / scale_factor, name = sec_axis_label)
    ) +
    scale_colour_discrete(name = legend_title) +
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
