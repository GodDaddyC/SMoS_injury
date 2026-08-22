# sensitivity analysis of dependence strength

bpot_dep_change <- function(param, ev_model) {
  M_name <- ev_model$model
  if (M_name %in% c("log", "hr", "neglog")) ev_model$estimate[5] <- param
  if (M_name == "alog") ev_model$estimate[7] <- param
  if (M_name %in% c("ct", "bilog", "negbilog", "amix")) {
    if (length(param) != 2) return("wrong number of dependence parameters")
    else ev_model$estimate[c(5, 6)] <- param
  }
  return(ev_model)
}

bpot_theoretical_density <- function(params, result, filename = NULL) {
  model_names <- result$M$model
  if (is.null(filename)) {
    saved_file <- sprintf("data/theoretical_density/theoretical_%s.csv",
                          model_names)
  } else {
    saved_file <- paste0("data/theoretical_density/", filename)
  }
  if (file.exists(saved_file)) {
    theoretical_density <- read.csv(saved_file, sep = ",")
    return(theoretical_density)
  }

  models_list <- setNames(
    lapply(params, function(p) {
      bpot_temp <- bpot_dep_change(p, result$M)
      return(bpot_temp)
    }),
    paste0("r=", params)
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
  write.table(theoretical_density, file = saved_file, sep = ",",
              row.names = TRUE)
  return(theoretical_density)
}

cbpot_theoretical_density <- function(params, cop_name, result,
                                      filename = NULL, restart = FALSE) {
  if (is.null(filename)) {
    saved_file <- sprintf("data/theoretical_density/theoretical_%s.csv",
                          cop_name)
  } else {
    saved_file <- paste0("data/theoretical_density/", filename)
  }

  if (restart && !is.null(filename)) file.remove(saved_file)

  if (file.exists(saved_file)) {
    theoretical_density <- read.csv(saved_file, sep = ",")
    return(theoretical_density)
  }

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
  write.table(theoretical_density, file = saved_file, sep = ",",
              row.names = TRUE)
  return(theoretical_density)
}


plot_theoretical_density <- function(theoretical_df,
    title = "Theoretical density of impact speed: BPOT",
    xlab = "Speed (km/h)", ylab = "f(y|TTC<0)",
    include_origin = TRUE) {
  if (!include_origin) {
    theoretical_df <- theoretical_df %>% select(-any_of("origin"))
  }
  plot_df <- theoretical_df %>%
    pivot_longer(cols = -speed, names_to = "Dependence strength",
                 values_to = "Density")

  ggplot(plot_df, aes(x = speed, y = Density, colour = `Dependence strength`)) +
    geom_line() +
    labs(title = title, x = xlab, y = ylab,
         colour = "Dependence strength") +
    theme_minimal()
}
