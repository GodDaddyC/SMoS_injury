# POT models


# change the mrlplot to the one you used in the other paper, 
evmix::mrlplot(SE_Dat$N_TTC)
evmix::mrlplot(SE_Dat$N_PET)
evmix::mrlplot(SE_Dat$maxDV_PET)
evmix::mrlplot(SE_Dat$maxDV_TTC)

evmix::mrlplot(CN_vru$TTC)
# Coordinates for the circle and the arrow
# x_circle <- -4.7  # Adjust based on your data
# y_circle <- 1.4  # Adjust based on your data
# x_arrow_from <- -4.7  # Arrow start x
# y_arrow_from <- 2.5  # Arrow start y
# x_arrow_to <- x_circle  # Arrow end x
# y_arrow_to <- 1.9  # Arrow end y

# Add a circle at (x_circle, y_circle)
# symbols(x_circle, y_circle, circles = 1, inches = 0.45,, add = TRUE)
# # Add an arrow pointing to the circle
# arrows(x_arrow_from, y_arrow_from, x_arrow_to, y_arrow_to, col = "red", lwd = 2)
# text(x_arrow_from, y_arrow_from + 0.1, labels = "Select threshold from this region", col = "red", pos = 3)

evd::tcplot(CN_vru$TTC,tlim = c(-2.5,-1.5))
evd::tcplot(SE_Dat$N_PET,tlim = c(-3,-2))




# run the threshold selection for CN and SE


# fit the univariate POT model to the CN and SE data

POT.CN.TTC <- fevd(x = TTC,data=CN_vru,threshold = quantile(CN_vru$TTC,0.8), type = "GP",
                   time.units = "2/month")
POT.CN.TTC$results
POT.CN.PET <- fevd(x = PET,data=CN_vru,threshold = quantile(CN_vru$TTC,0.8), type = "GP")
POT.CN.PET$results
POT.SE.TTC <- fevd(x = N_TTC,data=SE_Dat,threshold = quantile(SE_Dat$N_TTC,0.8),period.basis = "month",
                  time.units = "months", type = "GP")
POT.SE.TTC$results
POT.SE.PET<- fevd(x = N_PET,data=SE_Dat,threshold = quantile(SE_Dat$N_PET,0.8), type = "GP")
POT.SE.PET$results

POT.SE.DV_PET <- fevd(x = maxDV_PET,data=SE_Dat,threshold = quantile(SE_Dat$maxDV_PET,0.8), type = "GP")
POT.SE.DV_PET$results
POT.SE.DV_TTC <- fevd(x = maxDV_TTC,data=SE_Dat,threshold = quantile(SE_Dat$maxDV_TTC,0.8), type = "GP")
POT.SE.DV_TTC$results

plot(POT.CN.TTC)
plot(POT.CN.PET)
plot(POT.SE.TTC)
plot(POT.SE.PET)

# figure out the correct unit for return level plot
return_periods_to_plot <- c(10, 50, 100, 200, 300, 365)
estimated_return_levels <- return.level(POT.SE.TTC, return.period = return_periods_to_plot)

# Combine into a data frame for easier plotting
return_level_df <- data.frame(
  ReturnPeriod_Days = return_periods_to_plot,
  ReturnLevel = estimated_return_levels
)

message("\nReturn levels for various periods (in days):")
print(return_level_df)
plot(POT.SE.TTC, type = "rl")
points(return_level_df$ReturnPeriod_Days,return_level_df$ReturnLevel)

data(Tphap)
fit <- fevd(-MinT ~1, Tphap, threshold=-73, type="GP", units="deg F",
            time.units="62/year", verbose=TRUE)

fit
plot(fit)
plot(fit, "trace")



tt <- fevd(x = PET,data=CN_PET,threshold = quantile(CN_PET$PET,0.8), type = "GP",
                   threshold.fun ~ v1)
CN_PET$PET <- -CN_PET$PET
