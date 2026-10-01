# sensitivity analysis of dependence strength

#-----------------------BPOT-------------------------------------------
bpot_dep_change <- function(param, ev_model) {
  M_name <- ev_model$model
  if (M_name %in% c("log", "hr", "neglog")) ev_model$estimate[5] <- param
  if (M_name %in% c("alog", "aneglog")) ev_model$estimate[5:7] <- param
  if (M_name %in% c("ct", "bilog", "negbilog", "amix")) {
    if (length(param) != 2) return("wrong number of dependence parameters")
    else ev_model$estimate[c(5, 6)] <- param
  }
  return(ev_model)
}

# Invert extremal-dependence to the parameters of pickands dependence function
bpot_dep_from_chi <- function(chi,
                              model = c("log", "alog", "hr", "neglog",
                                        "aneglog", "bilog", "negbilog",
                                        "ct", "amix"),
                              asy = c(1, 1), tol = 1e-8) {
  model <- match.arg(model)
  identifiable_models <- c("log", "alog", "hr", "neglog", "aneglog")

  if (!is.numeric(chi) || any(!is.finite(chi)) || any(chi < 0 | chi > 1)) {
    stop("`chi` must contain finite values in [0, 1]", call. = FALSE)
  }
  if (length(tol) != 1 || !is.numeric(tol) || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive finite number", call. = FALSE)
  }
  if (!model %in% identifiable_models) {
    stop(sprintf("model '%s' has two dependence parameters; `chi` alone is not identifying",
                 model), call. = FALSE)
  }

  asymmetric <- model %in% c("alog", "aneglog")
  if (asymmetric) {
    if (length(asy) != 2 || !is.numeric(asy) || any(!is.finite(asy))) {
      stop("`asy` must be a numeric vector of length 2", call. = FALSE)
    }
    if (model == "alog" && (min(asy) < 0 || max(asy) > 1)) {
      stop("`asy` must lie in [0, 1] for the asymmetric logistic model",
           call. = FALSE)
    }
    if (model == "aneglog" && (min(asy) <= 0 || max(asy) > 1)) {
      stop("`asy` must lie in (0, 1] for the asymmetric negative logistic model",
           call. = FALSE)
    }
  }

  pickands <- function(dep, t) {
    args <- list(x = t, dep = dep, model = model, plot = FALSE)
    if (asymmetric) args$asy <- asy
    do.call(evd::abvevd, args)
  }

  chi_at_dep <- function(dep) {
    t_star <- if (asymmetric) {
      optimize(function(t) pickands(dep, t), interval = c(0, 1),
               tol = tol)$minimum
    } else {
      0.5
    }
    2 * (1 - pickands(dep, t_star))
  }

  invert_one <- function(target) {
    if (target == 0 && model == "log") return(1)
    if (target == 0 && model == "alog") {
      if (any(asy == 0)) {
        stop("`chi` = 0 does not identify `dep` when `asy` contains zero",
             call. = FALSE)
      }
      return(1)
    }
    if (target <= 0 || target >= 1) {
      stop(sprintf("`chi` = %g has no finite parameter value for model '%s'",
                   target, model), call. = FALSE)
    }

    lower <- .Machine$double.eps^0.5
    upper <- 1
    bounded <- model %in% c("log", "alog")
    f <- function(dep) chi_at_dep(dep) - target
    f_lower <- f(lower)
    f_upper <- f(upper)

    if (!is.finite(f_lower) || !is.finite(f_upper)) {
      stop(sprintf("`chi` = %g is outside the attainable range for model '%s'",
                   target, model), call. = FALSE)
    }
    while (!bounded && is.finite(f_upper) &&
           sign(f_lower) == sign(f_upper) && upper < 2^52) {
      upper <- upper * 2
      f_upper <- f(upper)
    }
    if (!is.finite(f_lower) || !is.finite(f_upper) ||
        sign(f_lower) == sign(f_upper)) {
      stop(sprintf("`chi` = %g is outside the attainable range for model '%s'",
                   target, model), call. = FALSE)
    }

    uniroot(f, interval = c(lower, upper), tol = tol)$root
  }

  vapply(chi, invert_one, numeric(1))
}

bpot_theoretical_severe <- function(params, result, severity_boundary) {
  if (length(params) == 0) {
    stop("`params` must contain at least one value", call. = FALSE)
  }
  if (length(severity_boundary) != 1 || !is.numeric(severity_boundary) ||
      !is.finite(severity_boundary) || severity_boundary < 0) {
    stop("`severity_boundary` must be one non-negative finite number",
         call. = FALSE)
  }
  if (!is.list(result) || is.null(result$M) || is.null(result$xcrash)) {
    stop("`result` must be the output of `run_single_bpot()`", call. = FALSE)
  }

  vapply(params, function(param) {
    ev_model <- bpot_dep_change(param, result$M)
    if (is.character(ev_model)) stop(ev_model, call. = FALSE)
    pb_tvevd_sev(ev_model, q1 = result$xcrash,
                 q2 = severity_boundary)
  }, numeric(1))
}

bpot_theoretical_density <- function(params, result) {
  models_list <- setNames(
    lapply(params, function(p) {
      bpot_temp <- bpot_dep_change(p, result$M)
      return(bpot_temp)
    }),
    paste0("dependence param(s):", params)
  )

  theoretical_density <- data.frame(speed = result$plot_df$speed) %>%
    mutate(origin = result$plot_df$ConditionalD)

  bpot_cols <- lapply(seq_along(models_list), function(k) {
    M_temp <- models_list[[k]]
    dd <- create_plot_df(dat = result$plot_df$speed, x = result$xcrash,
                         model = M_temp, px = result$pcrash)$ConditionalD
  })
  bpot_cols <- do.call(bind_cols, bpot_cols)
  names(bpot_cols) <- names(models_list)

  theoretical_density <- bind_cols(theoretical_density, bpot_cols) %>% na.omit()
  return(theoretical_density)
}

bpot_theoretical_injury <- function(params,result, age_mean_list,filename = NULL, restart = FALSE){
  # slow, thus prepared results are loaded
  model_names <- result$M$model
  if (is.null(filename)) {
    saved_file <- sprintf("../data/theoretical_injury/theoretical_%s.csv",
                          model_names)
  } else {
    saved_file <- paste0("../data/theoretical_injury/", filename)
  }
  
  if (restart && !is.null(filename)) file.remove(saved_file)
  
  if (file.exists(saved_file)) {
    theoretical_injury <- read.csv(saved_file, sep = ",")
    return(theoretical_injury)
  }
  
  bpot_result <- data.frame(age_mean = age_mean_list)
  
  models_list <- setNames(
    lapply(params, function(p) {
      bpot_temp <- bpot_dep_change(p, result$M)
      return(bpot_temp)
    }),
    paste0("dependence param(s):", params)
  )
  for (i in 1:length(models_list)) {
    bpot_result <- bpot_result %>%
      cbind(sapply(seq_along(age_mean_list), function(k) {
        injury_from_c_bivariate_e(result$plot_df$speed,
                                  ev_model = models_list[[i]], severity = pis1, x0 = x0,
                                  px = result$pcrash, age_mean = age_mean_list[k])}))
  }
  colnames(bpot_result) <- c("age_mean", names(models_list))
  write.table(bpot_result, file = saved_file, sep = ",",row.names = TRUE)
  return(bpot_result)
}

# -----CBPOT-------------------------------------------------------
cbpot_theoretical_density <- function(params, cop_name, result) {
  models_list <- setNames(
    lapply(params, function(p) {
      do.call(paste0(cop_name, "Copula"), list(param = p))
    }),
    paste0("r=", params)
  )

  theoretical_density <- data.frame(speed = result$plot_df_q$speed) %>%
    mutate(origin = result$plot_df_q$ConditionalD)

  cbpot_cols <- lapply(seq_along(models_list), function(k) {
    M_temp <- models_list[[k]]
    dd <- create_plot_df_q(dat = result$s_un, x = result$x0_un,
                           model = M_temp, px = 1 - result$x0_un,
                           p2 = result$p2, pu = result$pu)$ConditionalD
  })

  cbpot_cols <- do.call(bind_cols, cbpot_cols)
  names(cbpot_cols) <- names(models_list)
  theoretical_density <- bind_cols(theoretical_density, cbpot_cols) %>%
    na.omit()
  return(theoretical_density)
}

cbpot_theoretical_severe <- function(params, cop_name, result,
                                     severity_boundary) {
  if (length(severity_boundary) != 1 || !is.numeric(severity_boundary) ||
      !is.finite(severity_boundary) || severity_boundary < 0) {
    stop("`severity_boundary` must be one non-negative finite number",
         call. = FALSE)
  }
  if (!is.character(cop_name) || length(cop_name) != 1 ||
      !nzchar(cop_name)) {
    stop("`cop_name` must be one non-empty character string", call. = FALSE)
  }
  if (!is.list(result) || is.null(result$qcrash) || is.null(result$p2) ||
      is.null(result$Cop)) {
    stop("`result` must be the output of `run_single_cbpot()`", call. = FALSE)
  }
  if (!inherits(result$Cop@copula, paste0(cop_name, "Copula"))) {
    stop("`cop_name` does not match the copula in `result`", call. = FALSE)
  }

  speed_tail <- pgamma(
    severity_boundary,
    shape = result$p2$estimate[1],
    rate = result$p2$estimate[2],
    lower.tail = FALSE
  )

  vapply(params, function(param) {
    copula <- do.call(paste0(cop_name, "Copula"), list(param = param))
    pCopula(c(result$qcrash, speed_tail), copula)
  }, numeric(1))
}

cbpot_theoretical_injury <- function(params, cop_name, result,age_mean_list,conseq,
                                     filename = NULL, restart = FALSE) {
  if (is.null(filename)) {
    saved_file <- sprintf("../data/theoretical_injury/theoretical_%s.csv",
                          cop_name)
  } else {
    saved_file <- paste0("../data/theoretical_injury/", filename)
  }
  
  if (restart && !is.null(filename)) file.remove(saved_file)
  
  if (file.exists(saved_file)) {
    theoretical_injury <- read.csv(saved_file, sep = ",")
    return(theoretical_injury)
  }
  
  models_list <- setNames(
    lapply(params, function(p) {
      do.call(paste0(cop_name, "Copula"), list(param = p))
    }),
    paste0("r=", params)
  )
  
  cbpot_result <- data.frame(age_mean = age_mean_list)
  
  for (i in 1:length(params)) {
    cbpot_result <- cbpot_result %>%
      cbind(sapply(seq_along(age_mean_list), function(k) {
        injury_from_c_q_bivariate_e(models_list[[i]], pis1,
                                    result$x0_un, result$qcrash, conseq, lb = 0.5,
                                    age_mean = age_mean_list[k])
      }))
  }
  colnames(cbpot_result) <- c("age_mean",paste0("dependence:", params))
  write.table(cbpot_result, file = saved_file, sep = ",",row.names = TRUE)
  return(cbpot_result)
}

plot_theoretical_density <- function(theoretical_df,
    title = "Theoretical density of crash severity",
    xlab = "Speed (km/h)", ylab = "f(y|TTC<0)",
    include_origin = TRUE, filename = NULL) {
  if (!include_origin) {
    theoretical_df <- theoretical_df %>% dplyr::select(-any_of("origin"))
  }
  plot_df <- theoretical_df %>%
    pivot_longer(cols = -speed, names_to = "Dependence strength",
                 values_to = "Density")

  p <- ggplot(plot_df, aes(x = speed, y = Density,
                           colour = `Dependence strength`)) +
    geom_line() +
    labs(title = title, x = xlab, y = ylab,
         colour = "Dependence strength") +
    theme_minimal()

  if (!is.null(filename)) save_plot(p, filename, dpi = 300)
  return(p)
}

plot_theoretical_injury <- function(theoretical_df, filename = NULL) {
  plot_df <- theoretical_df %>%
    pivot_longer(cols = -age_mean, names_to = "Dependence strength",
                 values_to = "Injury")

  p <- ggplot(plot_df, aes(x = age_mean, y = Injury,
                           colour = `Dependence strength`)) +
    geom_line() +
    labs(x = "age", y = "AIS3+ probability",
         colour = "Dependence strength") +
    theme_minimal()

  if (!is.null(filename)) save_plot(p, filename, dpi = 300)

  return(p)
}
