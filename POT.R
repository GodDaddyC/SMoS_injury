# POT models


# change the mrlplot to the one you used in the other paper, 
mrlplot(SE_Dat$N_TTC)
# Coordinates for the circle and the arrow
x_circle <- -4.7  # Adjust based on your data
y_circle <- 1.4  # Adjust based on your data
x_arrow_from <- -4.7  # Arrow start x
y_arrow_from <- 2.5  # Arrow start y
x_arrow_to <- x_circle  # Arrow end x
y_arrow_to <- 1.9  # Arrow end y

# Add a circle at (x_circle, y_circle)
symbols(x_circle, y_circle, circles = 1, inches = 0.45,, add = TRUE)
# Add an arrow pointing to the circle
arrows(x_arrow_from, y_arrow_from, x_arrow_to, y_arrow_to, col = "red", lwd = 2)
text(x_arrow_from, y_arrow_from + 0.1, labels = "Select threshold from this region", col = "red", pos = 3)

meanExcessFunMk2(data = -CN_vru$TTC,u.prob=0.9)
# run the threshold selection for CN and SE


# fit the univariate POT model to the CN and SE data

POT.CN.TTC <- fevd(x = TTC,data=CN_vru,threshold = quantile(CN_vru$TTC,0.8), type = "GP")
POT.CN.TTC$results
plot(POT.CN.TTC)


# figure out the correct unit for return level plot