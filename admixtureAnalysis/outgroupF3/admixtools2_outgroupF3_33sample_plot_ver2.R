# ------------------------------------------------------------------
# admixtools2_outgroupF3_33sample_plot.R
#
# Plot the Outgroup F3 statistics
# Date: 2025/04/09
# Author: Shingo Fujimoto
# ------------------------------------------------------------------

# Import the table of outgroup F3
F3statDF <- read.csv("outgroupF3.csv")

# PDF format output
pdf(file = "OutgroupF3_medakaWGS_HSOK_sak_lat_33samples.pdf", width = 6, height = 4)
par(mar=c(5, 8, 1, 1))
par(mfrow = c(1,2))
plot(F3statDF$est, seq(1, length(F3statDF$est), 1),
     ylab = "", xlab = "Outgroup F3 (X, Y, HSOK)", yaxt = "n",
     xlim = c(0.33, 0.35), ylim = c(0.5, length(F3statDF$est) + 0.5),
     main = "", pch = 16, col = "gray", cex.axis = 0.75)

# add the three population combinations at y-axis
axis(side = 2, seq(1, length(F3statDF$est), 1),
     labels = paste0(F3statDF$pop3, ",", F3statDF$pop2), las = 1, cex.axis = 0.75
     )

# SE error bar of Outgroup F3
segments(y0 = seq(1, length(F3statDF$est), 1),
         y1 = seq(1, length(F3statDF$est), 1),
         x0 = F3statDF$est - F3statDF$se,
         x1 = F3statDF$est + F3statDF$se, lwd = 2, col = "gray")

dev.off()

