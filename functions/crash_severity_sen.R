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
    include_origin = TRUE) {
  if (!include_origin) {
    theoretical_df <- theoretical_df %>% dplyr::select(-any_of("origin"))
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

plot_theoretical_injury <- function(theoretical_df, save_file = NULL) {
  plot_df <- theoretical_df %>%
    pivot_longer(cols = -age_mean, names_to = "Dependence strength",
                 values_to = "Injury")

  p <- ggplot(plot_df, aes(x = age_mean, y = Injury,
                           colour = `Dependence strength`)) +
    geom_line() +
    labs(x = "age", y = "AIS3+ probability",
         colour = "Dependence strength") +
    theme_minimal()

  if (is.character(save_file)) {
    dir.create("plots", showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path("plots", save_file), plot = p, dpi = 300)
  }

  return(p)
}
