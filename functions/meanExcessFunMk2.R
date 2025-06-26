meanExcessFunMk2 <- function(data,plottitle=NULL,tlim=NULL,u.prob=NULL,
                             nt = min(120,length(data)), alpha=0.05,p.or.n =FALSE,...){
  # modified from the meanExcess function in the evmix package
  if (is.null(plottitle)){
    plottitle <- "Mean Excess Plot with Fitted Line"
  }
  try.threshold <- ifelse(is.null(u.prob), quantile(data, probs = 0.85), quantile(data,probs=u.prob))
  if (is.unsorted(data)) {
    data = sort(data)
  }
  else {
    if (data[1] > data[length(data)]) 
      data = rev(data)
  }
  
  if (is.null(tlim)) {
    thresholds = seq(median(data) - 2 * .Machine$double.eps, 
                     data[length(data) - 6], length.out = nt)
  } else {
    thresholds = seq(tlim[1], tlim[2], length.out = nt)
  }
  
  n <- length(data)
  data <- data[data > min(thresholds)]
  
  me.calc <- function(u, x, alpha) {
    excesses = x[x > u] - u
    nxs = length(excesses)
    meanxs = mean(excesses)
    sdxs = ifelse(nxs <= 5, NA, sd(excesses))
    results = c(u, nxs, meanxs, sdxs)
    if (!is.null(alpha)) {
      results = c(results, meanxs + qnorm(c(alpha/2, 1 - alpha/2)) * sdxs/sqrt(nxs))
    }
    return(results)
  }
  me <- t(sapply(thresholds, FUN = me.calc, x = data, alpha = alpha)) %>% as.data.frame()
  
  # compute the gradient
  fitresults <- fgpd(data, try.threshold, std.err = FALSE)
  mleparams <- fitresults$mle
  mrlint <- (mleparams[1] - mleparams[2] * try.threshold) / (1 - mleparams[2])
  mrlgrad <- mleparams[2] / (1 - mleparams[2])
  
  line_data <- data.frame(
    x = c(try.threshold, max(thresholds)),
    y = mrlint + mrlgrad * c(try.threshold, max(thresholds)))
  
  dashed_line_data <- data.frame(
    x = c(min(thresholds), try.threshold),
    y = mrlint + mrlgrad * c(min(thresholds), try.threshold))
  
  # plottomg
  if (!is.null(alpha)) {
    names(me) = c("u", "nu", "mean.excess", "sd.excess", 
                  "cil.excess", "ciu.excess")
    # ggplot for visualization
    p <- ggplot() +
      # Add mean excess and CI lines
      geom_line(data = me, aes(x = u, y = mean.excess, color = "Mean Excess"), 
                linetype = 'solid', linewidth = 1) +
      geom_line(data = me, aes(x = u, y = cil.excess), linetype = 'dashed', 
                linewidth = 0.6) + 
      geom_line(data = me, aes(x = u, y = ciu.excess), linetype = 'dashed', 
                linewidth = 0.6) + 
      geom_line(data = line_data, aes(x = x, y = y, color = "Gradient Line"), 
                linetype = "solid", linewidth = 0.8) +
      geom_line(data = dashed_line_data, aes(x = x, y = y, color = "Gradient Line"), 
                linetype = "dashed", linewidth = 0.8) +
      geom_vline(aes(xintercept = try.threshold, color = "Threshold Line"), 
                 linetype = "dotted") +
      scale_color_manual(
        values = c(
          "Mean Excess" = "blue",
          "Gradient Line" = "red",
          "Threshold Line" = "red"
        ),
        breaks = c("Mean Excess", "Gradient Line", "Threshold Line"),
        labels = c("Mean Excess",paste0("Gradient Line (Slope: ", round(mrlgrad, 2), ")"),
                   sprintf("Threshold (%.2f percentile)",u.prob))
      ) +
      labs(color = plottitle) +
      theme_minimal()
    p
  }
  else {
    names(me) = c("u", "nu", "mean.excess", "sd.excess")
    p <- ggplot() +
      geom_line(data = me, aes(x = u, y = mean.excess, color = "Mean Excess"), 
                linetype = 'solid', linewidth = 1) +
      geom_line(data = line_data, aes(x = x, y = y, color = "Gradient Line"), 
                linetype = "solid", linewidth = 0.8) +
      geom_line(data = dashed_line_data, aes(x = x, y = y, color = "Gradient Line"), 
                linetype = "dashed", linewidth = 0.8) +
      geom_vline(aes(xintercept = try.threshold, color = "Threshold Line"), 
                 linetype = "dotted") +
      # Custom legend labels
      scale_color_manual(
        values = c(
          "Mean Excess" = "blue",
          "Gradient Line" = "red",
          "Threshold Line" = "red"
        ),
        breaks = c("Mean Excess", "Gradient Line", "Threshold Line"),
        labels = c(
          "Mean Excess (u)",
          paste0("Gradient Line (Slope: ", round(mrlgrad, 2), ")"),
          sprintf("Threshold (%.2f percentile)",u.prob)
        )
      ) +
      labs(title = plottitle) +  # Legend title
      theme_minimal()
    p
  }
  
  
}

