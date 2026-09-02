mean_excess_fun_mk2 <- function(data, plottitle = NULL, tlim = NULL,
                                   u_prob = NULL,
                                   nt = min(120, length(data)),
                                   alpha = 0.05, filename = NULL, ...) {
  if (is.null(plottitle)) {
    plottitle <- "Mean Excess Plot with Fitted Line"
  }
  try_threshold <- ifelse(is.null(u_prob),
    quantile(data, probs = 0.85), quantile(data, probs = u_prob))
  if (is.unsorted(data)) {
    data <- sort(data)
  } else {
    if (data[1] > data[length(data)]) data <- rev(data)
  }

  if (is.null(tlim)) {
    thresholds <- seq(median(data) - 2 * .Machine$double.eps,
                      data[length(data) - 6], length.out = nt)
  } else {
    thresholds <- seq(tlim[1], tlim[2], length.out = nt)
  }

  n <- length(data)
  data <- data[data > min(thresholds)]

  me_calc <- function(u, x, alpha) {
    excesses <- x[x > u] - u
    nxs <- length(excesses)
    meanxs <- mean(excesses)
    sdxs <- ifelse(nxs <= 5, NA, sd(excesses))
    results <- c(u, nxs, meanxs, sdxs)
    if (!is.null(alpha)) {
      results <- c(results,
        meanxs + qnorm(c(alpha / 2, 1 - alpha / 2)) * sdxs / sqrt(nxs))
    }
    return(results)
  }
  me <- t(sapply(thresholds, FUN = me_calc, x = data, alpha = alpha)) %>%
    as.data.frame()

  fitresults <- fgpd(data, try_threshold, std.err = FALSE)
  mleparams <- fitresults$mle
  mrlint <- (mleparams[1] - mleparams[2] * try_threshold) /
            (1 - mleparams[2])
  mrlgrad <- mleparams[2] / (1 - mleparams[2])

  line_data <- data.frame(
    x = c(try_threshold, max(thresholds)),
    y = mrlint + mrlgrad * c(try_threshold, max(thresholds)))

  dashed_line_data <- data.frame(
    x = c(min(thresholds), try_threshold),
    y = mrlint + mrlgrad * c(min(thresholds), try_threshold))

  if (!is.null(alpha)) {
    names(me) <- c("u", "nu", "mean.excess", "sd.excess",
                   "cil.excess", "ciu.excess")
    p <- ggplot() +
      geom_line(data = me, aes(x = u, y = mean.excess, color = "Mean Excess"),
                linetype = "solid", linewidth = 1) +
      geom_line(data = me, aes(x = u, y = cil.excess),
                linetype = "dashed", linewidth = 0.6) +
      geom_line(data = me, aes(x = u, y = ciu.excess),
                linetype = "dashed", linewidth = 0.6) +
      geom_line(data = line_data,
                aes(x = x, y = y, color = "Gradient Line"),
                linetype = "solid", linewidth = 0.8) +
      geom_line(data = dashed_line_data,
                aes(x = x, y = y, color = "Gradient Line"),
                linetype = "dashed", linewidth = 0.8) +
      geom_vline(aes(xintercept = try_threshold, color = "Threshold Line"),
                 linetype = "dotted") +
      scale_color_manual(
        values = c("Mean Excess" = "blue", "Gradient Line" = "red",
                   "Threshold Line" = "red"),
        breaks = c("Mean Excess", "Gradient Line", "Threshold Line"),
        labels = c("Mean Excess",
                   paste0("Gradient Line (Slope: ", round(mrlgrad, 2), ")"),
                   sprintf("Threshold (%.2f percentile)", u_prob))
      ) +
      labs(color = plottitle) +
      theme_minimal()
    if (!is.null(filename)) save_plot(p, filename)
    p
  } else {
    names(me) <- c("u", "nu", "mean.excess", "sd.excess")
    p <- ggplot() +
      geom_line(data = me, aes(x = u, y = mean.excess, color = "Mean Excess"),
                linetype = "solid", linewidth = 1) +
      geom_line(data = line_data,
                aes(x = x, y = y, color = "Gradient Line"),
                linetype = "solid", linewidth = 0.8) +
      geom_line(data = dashed_line_data,
                aes(x = x, y = y, color = "Gradient Line"),
                linetype = "dashed", linewidth = 0.8) +
      geom_vline(aes(xintercept = try_threshold, color = "Threshold Line"),
                 linetype = "dotted") +
      scale_color_manual(
        values = c("Mean Excess" = "blue", "Gradient Line" = "red",
                   "Threshold Line" = "red"),
        breaks = c("Mean Excess", "Gradient Line", "Threshold Line"),
        labels = c("Mean Excess (u)",
                   paste0("Gradient Line (Slope: ", round(mrlgrad, 2), ")"),
                   sprintf("Threshold (%.2f percentile)", u_prob))
      ) +
      labs(title = plottitle) +
      theme_minimal()
    if (!is.null(filename)) save_plot(p, filename)
    p
  }
}
