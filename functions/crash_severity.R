# crash severity densities derived from difference approaches


# BPOT --------------------------------------------------------------------

c.bivariate <-function(y,x, model, dep,alpha,beta, thres, eta, mar1, mar2){
  # approximates the conditional density f(y|X >x) by integrating fxy = f(x=x, y = y) over x
  # if integral is non-finite, change ulim to a smaller value
  integrand <- function(q1,...) {
    result <- numeric(length(q1))
    for(i in seq_along(q1)) {
      result[i] <- switch(model,
                          log=dbTvevd(q1 = q1[i], q2 = y, model = model, dep = dep, 
                                      thres = thres, eta = eta, mar1 = mar1, mar2 = mar2),
                          hr=dbTvevd(q1 = q1[i], q2 = y, model = model, dep = dep, 
                                     thres = thres, eta = eta, mar1 = mar1, mar2 = mar2),
                          ct=dbTvevd(q1 = q1[i], q2 = y, model = model, alpha=alpha,beta=beta, 
                                     thres = thres, eta = eta, mar1 = mar1, mar2 = mar2)
      )
    }
    return(result)
  }
  
  # Integrate from x to Inf
  R <- tryCatch({
    integrate(integrand, lower = x, upper = Inf,rel.tol = 1e-3)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered in conditional density. Retrying with finite upper bound for x.")
        if (mar1[2] < 0) {
          ub.alt <- thres[1] - mar1[1] / mar1[2]
          return(integrate(integrand, lower = x, upper = ub.alt, rel.tol = 1e-3)$value)
        } 
        else {
          # It attempts the integration with decreasing multipliers from 5 down to 1.
          for (i in rev(seq(0.05,5,0.05)) ) {
            ub.alt <- x + i * mar1[1]
            message(paste("Attempting integration with upper bound: ", ub.alt, " (multiplier =", i, ")"))
            result <- tryCatch({
              integrate(integrand, lower = x, upper = ub.alt, rel.tol = 1e-3)$value
            }, 
            error = function(e_retry) {
              return(NA)
            })
            # Check if the integration was successful (i.e., didn't return NA).
            if (!is.na(result)) {
              message("Integration successful with a finite upper bound.")
              return(result)
            }
          }
          # If the loop completes without a successful return, it means all attempts failed.
          stop("Failed to find a finite upper bound after multiple retries.", call. = FALSE)
        }
      } 
      else { stop(e)}
    })
  
  return(R)
}

c.bivariate_np <- function(y,x,thres, eta, mar1, mar2,Ahat){
  integrand <- function(q1,...) {
    result <- numeric(length(q1))
    for(i in seq_along(q1)) {
      result[i] <- dbTvNonpar(q1=q1[i],q2=y,Ahat=Ahat,mar1=mar1,mar2=mar2,thres=thres,eta=eta)
    }
    return(result)
  }
  
  # Integrate from x to Inf
  R <- tryCatch({
    integrate(integrand, lower = x, upper = Inf,rel.tol = 1e-3)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered in conditional density. Retrying with finite upper bound for x.")
        if (mar1[2] < 0) {
          # ub.alt <- thres[1] - mar1[1] / mar1[2]
          # return(integrate(integrand, lower = x, upper = ub.alt, rel.tol = 1e-3)$value)
          
          for (j in rev(seq(x,thres[1] - mar1[1] / mar1[2],length.out=10))){
            ub.alt <- j
            message(paste("Attempting integration with upper bound: ", ub.alt))
            result <- tryCatch({
              integrate(integrand, lower = x, upper = ub.alt, rel.tol = 1e-3)$value
            },
            error = function(e_retry) {
              return(NA)
            })
            # Check if the integration was successful (i.e., didn't return NA).
            if (!is.na(result)) {
              message("Integration successful with a finite upper bound.")
              return(result)
            }
          }
        }
        
        else {
          for (i in rev(seq(0.05,5,0.05)) ) {
            ub.alt <- x + i * mar1[1]
            message(paste("Attempting integration with upper bound: ", ub.alt, " (multiplier =", i, ")"))
            result <- tryCatch({
              integrate(integrand, lower = x, upper = ub.alt, rel.tol = 1e-3)$value
            }, 
            error = function(e_retry) {
              return(NA)
            })
            # Check if the integration was successful (i.e., didn't return NA).
            if (!is.na(result)) {
              message("Integration successful with a finite upper bound.")
              return(result)
            }
          }
          # If the loop completes without a successful return, it means all attempts failed.
          stop("Failed to find a finite upper bound after multiple retries.", call. = FALSE)
        }
      }
      else { stop(e)}
    })
  
  return(R)
}

normalize_c.bivariate<- function(x, EVmodel){
  # computes the infinite integral of the conditional density f(y|X >x)
  # use as a nomralization factor for the conditional density
  dat <- EVmodel$data[EVmodel$data[,2]>EVmodel$threshold[2],2]
  c.y <- function(k) {
    temp <- switch(EVmodel$model,
                   log=c.bivariate(y = k, x = x,  model = EVmodel$model, 
                                   dep = EVmodel$estimate[5], thres = EVmodel$threshold, 
                                   eta = EVmodel$nat[1:2]/EVmodel$n,
                                   mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                   mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])),
                   hr=c.bivariate(y = k, x = x, model = EVmodel$model, 
                                  dep = EVmodel$estimate[5], thres = EVmodel$threshold, 
                                  eta = EVmodel$nat[1:2]/EVmodel$n,
                                  mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                  mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])),
                   ct=c.bivariate(y = k, x = x, model = EVmodel$model, 
                                  alpha = EVmodel$estimate[5],beta=EVmodel$estimate[6], 
                                  thres = EVmodel$threshold, 
                                  eta = EVmodel$nat[1:2]/EVmodel$n,
                                  mar1 = c(EVmodel$estimate[1], EVmodel$estimate[2]), 
                                  mar2 = c(EVmodel$estimate[3], EVmodel$estimate[4])))
    return(temp)
  }
  
  C <- tryCatch({
    integrate(Vectorize(c.y), lower = min(dat), upper = Inf)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message) || grepl("Failed to find a finite upper", e$message)) {
        message("Non-finite function value encountered in norming. Retrying with finite upper bound for y.")
        ub.alt <- ifelse(EVmodel$estimate[4]<0,EVmodel$threshold[2] - EVmodel$estimate[3]
                         /EVmodel$estimate[4],55)
        return(integrate(Vectorize(c.y), lower = min(dat), upper = ub.alt)$value)
      } 
      else {
        stop(e)  # rethrow other errors
      }
    })
  nc <- switch(EVmodel$model,
               log = pbTvevd(q1=x,q2=EVmodel$threshold[2],dep=EVmodel$estimate[5],thres=EVmodel$threshold,model=EVmodel$model,
                             eta=EVmodel$nat[1:2]/EVmodel$n,mar1=EVmodel$estimate[1:2],mar2=EVmodel$estimate[3:4],tail.type=2),
               hr = pbTvevd(q1=x,q2=EVmodel$threshold[2],dep=EVmodel$estimate[5],thres=EVmodel$threshold,model=EVmodel$model,
                            eta=EVmodel$nat[1:2]/EVmodel$n,mar1=EVmodel$estimate[1:2],mar2=EVmodel$estimate[3:4],tail.type=2),
               ct = pbTvevd(q1=x,q2=EVmodel$threshold[2],alpha=EVmodel$estimate[5],beta=EVmodel$estimate[6],
                            thres=EVmodel$threshold,model=EVmodel$model,eta=EVmodel$nat[1:2]/EVmodel$n,
                            mar1=EVmodel$estimate[1:2],mar2=EVmodel$estimate[3:4],tail.type=2) )
  return(C/nc)
}


# CBPOT -------------------------------------------------------------------


cQ.bivariate <-function(y,x,model){
  # approximates the conditional density f(y|X >x) by integrating fxy = f(x=x, y = y) over x
  # model is an mvdc object
  # if integral is non-finite, change ulim to a smaller value
  
  integrand <- function(q1) {
    result <- dCopula(c(q1,y),model) 
    
    return(result)
  }
  
  # Integrate from x to Inf
  R <- tryCatch({
    integrate(Vectorize(integrand), lower = x, upper = 1)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- 1 - 1e-7
        return(integrate(integrand, lower = x, upper = ub.alt)$value)
      }
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(R)
}

normalize_cQ.bivariate<- function(x, model,lb=0){
  # computes the infinite integral of the conditional density f(y|X >x)
  # use as a nomralization factor for the conditional density
  # lb is the lower bound of the integral, by defalut lb = 0
  c.y <- function(k) {
    temp <- cQ.bivariate(y = k, x = x, model = model)
    return(temp)
  }
  
  C <- tryCatch({
    integrate(Vectorize(c.y), lower = lb, upper = 1)$value},
    error = function(e) {
      if (grepl("non-finite function value", e$message)) {
        message("Non-finite function value encountered. Retrying with finite upper bound.")
        ub.alt <- 1-1e-7
        return(integrate(Vectorize(c.y), lower = lb, upper = ub.alt)$value)
      }
      if (grepl("maximum number of subdivisions reached", e$message)) {
        message("Maximum number of subdivisions reached. Retrying with larger tolerence.")
        return(integrate(Vectorize(c.y), lower = lb,upper=Inf,
                         subdivisions = 200,rel.tol = 1e-5)$value)
      }
      else {
        stop(e)  # rethrow other errors
      }
    })
  return(C) 
}


# plotting ----------------------------------------------------------------

plot_crash_severity <- function(plot.df.1, plot.df.2, injury.df,
                                v1 = NULL, v2 = NULL,
                                label1 = "CN", label2 = "SE",
                                legend.title = "Site",
                                sec.axis.label = "Injury prob",
                                save.path = NULL) {
  # plot.df.1, plot.df.2: data frames with columns speed and ConditionalD (one per site)
  # injury.df: data frame with columns speed and InjuryP
  # SAVE: if TRUE, saves the plot to save.path
  # v1, v2: optional thresholds for site 1 and site 2 drawn as vertical reference lines
  # label1, label2: legend labels for the two sites
  # legend.title: title of the colour legend
  # sec.axis.label: label for the secondary y-axis

  scale_factor <- max(plot.df.1$ConditionalD) / max(injury.df$InjuryP)

  col1 <- "red"
  col2 <- "blue"

  p <- ggplot(plot.df.1, aes(x = speed, y = ConditionalD)) +
    geom_line(aes(colour = label1)) +
    geom_line(data = plot.df.2, aes(x = speed, y = ConditionalD, colour = label2)) +
    geom_line(data = injury.df, aes(x = speed, y = InjuryP * scale_factor))

  if (!is.null(v1)) {
    p <- p +
      geom_vline(xintercept = v1, linetype = "dashed", color = col1) +
      annotate("text", x = v1 + 2, y = 0.02, label = paste("u=", round(v1, 3)), color = col1)
  }
  if (!is.null(v2)) {
    p <- p +
      geom_vline(xintercept = v2, linetype = "dashed", color = col2) +
      annotate("text", x = v2 + 2, y = 0.015, label = paste("u=", round(v2, 3)), color = col2)
  }

  col.values <- setNames(c(col1, col2), c(label1, label2))

  p <- p +
    scale_y_continuous(
      name = "density",
      sec.axis = sec_axis(~ . / scale_factor, name = sec.axis.label)
    ) +
    scale_colour_manual(name = legend.title, values = col.values) +
    labs(x = "y", y = "f(y|crash)") +
    theme(
      panel.grid.major = element_line(colour = "gray91"),
      panel.grid.minor = element_line(colour = "gray88"),
      panel.background = element_rect(fill = "white", colour = "white", linetype = "solid"),
      plot.background = element_rect(linetype = "solid")
    )

  if (!is.null(save.path)) {
    ggsave(save.path, plot = p,dpi = 300,device='png')
  }
  return(p)
}


plot_crash_severity_single <- function(plot.df, injury.df,
                                       v = NULL,
                                       label = NULL,
                                       show.legend = FALSE,
                                       sec.axis.label = "Injury prob",
                                       save.path = NULL) {
  # plot.df: data frame with columns speed and ConditionalD (single site)
  # injury.df: data frame with columns speed and InjuryP
  # SAVE: if TRUE, saves the plot to save.path
  # v: optional threshold drawn as a vertical reference line
  # label: optional site label shown in the legend; if NULL no colour legend is drawn
  # show.legend: if FALSE (default) suppresses the colour legend
  # sec.axis.label: label for the secondary y-axis

  scale_factor <- max(plot.df$ConditionalD) / max(injury.df$InjuryP)

  if (!is.null(label)) {
    p <- ggplot(plot.df, aes(x = speed, y = ConditionalD, colour = label))
  } else {
    p <- ggplot(plot.df, aes(x = speed, y = ConditionalD))
  }

  p <- p +
    geom_line(show.legend = show.legend) +
    geom_line(data = injury.df, aes(x = speed, y = InjuryP * scale_factor),
              colour = "black", inherit.aes = FALSE)

  if (!is.null(v)) {
    p <- p +
      geom_vline(xintercept = v, linetype = "dashed", color = "red") +
      annotate("text", x = v + 2, y = 0.02, label = paste("u=", round(v, 3)), color = "red")
  }

  p <- p +
    scale_y_continuous(
      name = "density",
      sec.axis = sec_axis(~ . / scale_factor, name = sec.axis.label)
    ) +
    labs(x = "Speed (km/h)", y = "f(y|TTC<0)") +
    theme(
      panel.grid.major = element_line(colour = "gray91"),
      panel.grid.minor = element_line(colour = "gray88"),
      panel.background = element_rect(fill = "white", colour = "white", linetype = "solid"),
      plot.background = element_rect(linetype = "solid")
    )

  if (!is.null(save.path)) {
    ggsave(save.path, plot = p)
  }

  return(p)
}
