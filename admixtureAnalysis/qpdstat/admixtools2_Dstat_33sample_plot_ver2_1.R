# ------------------------------------------------------------------
# admixtools2_interCladeAdmix_dxyDirectionTestt.R
#
# admixture graph estimation using admixtools2
# Date: 2023/01/26
# update: 2024/05/09, Figure remake using the latest dataset
# update: 2024/10/17, Combined the two Dstat tables
# Author: Shingo Fujimoto
# ------------------------------------------------------------------
# Merge the D stat tables 
dstatDF <- rbind(read.csv("Dstat_EJWJWK_Osak_HSOK.csv"),
                 read.csv("Dstat_EJ_WJWK_Osak-sep_HSOK.csv"),
                 read.csv("Dstat_EJsep_WJWK_Osak_HSOK.csv"),
                 read.csv("Dstat_iso_EJsep_WJOsak_HSOK.csv"),
                 read.csv("Dstat_WJsep_WK_Osak_HSOK.csv"),
                 read.csv("Dstat_WJ_WKsep_Osak_HSOK.csv")
                 )
FigOrder <- -1 * as.numeric(rownames(dstatDF)) # plot order in the figure

pdf(file = "Dstat_medakaWGS_HSOK_sak_lat_33samples_eJsep_plot.pdf", width = 6, height = 8)
par(mfrow = c(1,2))
par(mar=c(3, 6, 1, 1))
plot(dstatDF$est, -1 * as.numeric(rownames(dstatDF)),
     ylab = "", xlab = "", yaxt = "n", xlim = c(-0.4, 0.75),
     main = "Patterson's D (W, X; Y, HSOK)",
     pch = 1, type = "n")

# SE error bar of D stat
segments(y0 = FigOrder, y1 = FigOrder,
         x0 = dstatDF$est - dstatDF$se, x1 = dstatDF$est + dstatDF$se, lwd = 3, col = "gray")
# plot the D stat
points(dstatDF$est, FigOrder, pch = 1, col = "gray", cex = 0.75)


# Significant level of positive D value(Z > 3)
segments(y0 = min(c(1:length(dstatDF$est))), y1 = -1 * max(c(1:length(dstatDF$est))) - 4,
         x0 = 0, x1 = 0, lwd = 1, lty = "dotted")
segments(y0 = min(c(1:length(dstatDF$est))), y1 = -1 * max(c(1:length(dstatDF$est))) - 4,
         x0 = 0.033, x1 = 0.033, lwd = 1, lty = "dotted", col = "gray")
segments(y0 = min(c(1:length(dstatDF$est))), y1 = -1 * max(c(1:length(dstatDF$est))) - 4,
         x0 = - 0.033, x1 = - 0.033, lwd = 1, lty = "dotted", col = "gray")

# Label of yaxis 
axis(side = 2, at = FigOrder,
     labels = paste0(dstatDF$pop1, ",", dstatDF$pop2, ";", dstatDF$pop3, ", O" ), las = 1, cex.axis = 0.5)



dev.off()

