# ------------------------------------------------------------
# population_statistics_summary.R
# Tajima's D statistics
# 2021/05/13
# update: 2021/07/26
# update: 2025/01/24, O. sakaizumii (N = 40), Pi, Tajima's D, sexFst, Fd
# update: 2025/12/08, centromere, QTL overlap
# Shingo Fujimoto
# ------------------------------------------------------------
library("tidyverse")
library("ggplot2")
library("mgcv")
source("../scripts/population_statistics_summary.R") # self-made functions to make the 

setwd("../Fig2-3_genomeScan_formatOutput/sak40_pi_tajimaD_sexFst_f_d")
# population statistics merged table
sakStatDF <- read.csv("O_sak40_pi_tajimaD_sexFst_f_d.csv")
# threshold of each pop-stat
thresholdDF <- read.csv("threshold0005.csv")

# ------------------------------------------------------------------------------
# Anal fin length, QTL region on chromosome 1
# import and modify Tajima's D, inter-sexual Fst files
# ------------------------------------------------------------------------------
# 1, import Tajima's D dataset (calculated by vcf-kit)
pdf(file = "Chr1_TajimaD_sexFst.pdf", width = 6.75, height = 6.75)

extractChr <- subset(sakStatDF, chromosome == 1)
extractChr <- dplyr::select(extractChr, chromosome, window_start, window_end, TajimaD) %>% distinct()
tempThreshold <- filter(thresholdDF, statistics == "Tajima_D")

par(mfrow = c(2,1))
# plot focusing chromosome
options(scipen = 1000)
plot(extractChr$window_start, extractChr$TajimaD, las = 1,
     xlab = "", ylab = "Tajima's D", ylim = c(-4, 4), col = "gray",
     xlim = c(2000000, 27000000), type = "n")
lines(extractChr$window_start, extractChr$TajimaD, col = "gray")
abline(tempThreshold$upperThreshold, 0, col = "gray", lty = "dotted")
abline(tempThreshold$lowerThreshold, 0, col = "gray", lty = "dotted")

# points the Tajima's D outlier
upperTD <- filter(extractChr, TajimaD > thresholdDF$upperThreshold[2])
lowerTD <- filter(extractChr, TajimaD < thresholdDF$lowerThreshold[2])
points(upperTD$window_start, upperTD$TajimaD, col = "gray30")
points(lowerTD$window_start, lowerTD$TajimaD, col = "gray30")

# Gene and centromere position on chromosome 1
# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(23530000, 25580000), c(-3.2, -3.2), lwd = 2)

# ------------------------------------------------------------------------------
# 2, import inter-sexual Fst dataset (calculated by pixy)
extractChr <- subset(sakStatDF, chromosome == 1)
extractChr <- dplyr::select(extractChr, chromosome, window_start, window_end, avg_wc_fst_sakFM) %>% distinct()
tempThreshold <- filter(thresholdDF, statistics == "avg_wc_fst_sakFM")

# plot focusing chromosome
plot(extractChr$window_start, extractChr$avg_wc_fst_sakFM, las = 1,
     xlab = "Positions in chromosome 1", ylab = "Inter-sexual Fst", ylim = c(-0.02, 0.10), col = "gray",
     xlim = c(2000000, 27000000), type = "n")
lines(extractChr$window_start, extractChr$avg_wc_fst_sakFM, col = "gray")
abline(tempThreshold$upperThreshold, 0, col = "gray", lty = "dotted")
abline(tempThreshold$lowerThreshold, 0, col = "gray", lty = "dotted")

# points the intersexual FST outlier
upperTD <- filter(extractChr, avg_wc_fst_sakFM > thresholdDF$upperThreshold[3])
points(upperTD$window_start, upperTD$avg_wc_fst_sakFM, col = "gray30")

# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(23530000, 25580000), c(-0.022, -0.022), lwd = 2)

dev.off()

# ------------------------------------------------------------------------------
# Female spawning latency, QTL region on chromosome 8
# import and modify Tajima's D, inter-sexual Fst files
# ------------------------------------------------------------------------------
# 1, import Tajima's D dataset (calculated by vcf-kit)
pdf(file = "Chr8_TajimaD_sexFst.pdf", width = 6.75, height = 6.75)

extractChr <- subset(sakStatDF, chromosome == 8)
extractChr <- dplyr::select(extractChr, chromosome, window_start, window_end, TajimaD) %>% distinct()
tempThreshold <- filter(thresholdDF, statistics == "Tajima_D")

par(mfrow = c(2,1))
# plot focusing chromosome
options(scipen = 1000)
plot(extractChr$window_start, extractChr$TajimaD, las = 1,
     xlab = "", ylab = "Tajima's D", ylim = c(-4, 4), col = "gray",
     xlim = c(2000000, 27000000), type = "n")
lines(extractChr$window_start, extractChr$TajimaD, col = "gray")
abline(tempThreshold$upperThreshold, 0, col = "gray", lty = "dotted")
abline(tempThreshold$lowerThreshold, 0, col = "gray", lty = "dotted")

# points the Tajima's D outlier
upperTD <- filter(extractChr, TajimaD > thresholdDF$upperThreshold[2])
lowerTD <- filter(extractChr, TajimaD < thresholdDF$lowerThreshold[2])
points(upperTD$window_start, upperTD$TajimaD, col = "gray30")
points(lowerTD$window_start, lowerTD$TajimaD, col = "gray30")

# Gene and centromere position on chromosome 8
# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(15190000, 17210000), c(-3, -3), lwd = 2)

# ------------------------------------------------------------------------------
# 2, import inter-sexual Fst dataset (calculated by pixy)
extractChr <- subset(sakStatDF, chromosome == 8)
extractChr <- dplyr::select(extractChr, chromosome, window_start, window_end, avg_wc_fst_sakFM) %>% distinct()
tempThreshold <- filter(thresholdDF, statistics == "avg_wc_fst_sakFM")

# plot focusing chromosome
plot(extractChr$window_start, extractChr$avg_wc_fst_sakFM, las = 1,
     xlab = "Positions in chromosome 8", ylab = "Inter-sexual Fst", ylim = c(-0.02, 0.10), col = "gray",
     xlim = c(2000000, 27000000), type = "n")
lines(extractChr$window_start, extractChr$avg_wc_fst_sakFM, col = "gray")
abline(tempThreshold$upperThreshold, 0, col = "gray", lty = "dotted")
abline(tempThreshold$lowerThreshold, 0, col = "gray", lty = "dotted")

# points the intersexual FST outlier
upperTD <- filter(extractChr, avg_wc_fst_sakFM > thresholdDF$upperThreshold[3])
points(upperTD$window_start, upperTD$avg_wc_fst_sakFM, col = "gray30")

# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(15190000, 17210000), c(-0.022, -0.022), lwd = 2)

dev.off()

# ------------------------------------------------------------------------------
# Spawning latency, QTL region on chromosome 13
# import and modify Tajima's D and intersexual Fst files
# 1, import Tajima's D dataset (calculated by vcf-kit)
# ------------------------------------------------------------------------------
pdf(file = "Chr13_TajimaD_sexFst.pdf", width = 6.75, height = 6.75)

extractChr <- subset(sakStatDF, chromosome == 13)
extractChr <- dplyr::select(extractChr, chromosome, window_start, window_end, TajimaD) %>% distinct()
tempThreshold <- filter(thresholdDF, statistics == "Tajima_D")

par(mfrow = c(2,1))
# plot focusing chromosome
options(scipen = 1000)
plot(extractChr$window_start, extractChr$TajimaD, las = 1,
     xlab = "", ylab = "Tajima's D", ylim = c(-4, 4), col = "gray",
     xlim = c(2000000, 27000000), type = "n")
lines(extractChr$window_start, extractChr$TajimaD, col = "gray")
abline(tempThreshold$upperThreshold, 0, col = "gray", lty = "dotted")
abline(tempThreshold$lowerThreshold, 0, col = "gray", lty = "dotted")

# points the Tajima's D outlier
upperTD <- filter(extractChr, TajimaD > thresholdDF$upperThreshold[2])
lowerTD <- filter(extractChr, TajimaD < thresholdDF$lowerThreshold[2])
points(upperTD$window_start, upperTD$TajimaD, col = "gray30")
points(lowerTD$window_start, lowerTD$TajimaD, col = "gray30")

# ------------------------------------------------------------------------------
# See also, QTL region information in Table S9, (Kawajiri et al. 2014, 2015)
lines(c(1600000, 8300000), c(-3, -3), col = "slategray2", lwd = 2)
text(x = mean(c(1600000, 8300000)), y = -3.6, labels = "Spawning latency", cex = 0.8)

# Gene and centromere position on chromosome 13
# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(27500000, 30850000), c(-3.3, -3.3), lwd = 2)

# ------------------------------------------------------------------------------
# 2, import pi Fst dataset (calculated by pixy)
extractChr <- subset(sakStatDF, chromosome == 13)
extractChr <- dplyr::select(extractChr, chromosome, window_start, window_end, avg_wc_fst_sakFM) %>% distinct()
tempThreshold <- filter(thresholdDF, statistics == "avg_wc_fst_sakFM")

# plot focusing chromosome
plot(extractChr$window_start, extractChr$avg_wc_fst_sakFM, las = 1,
     xlab = "Positions in chromosome 13", ylab = "Inter-sexual Fst", ylim = c(-0.02, 0.10), col = "gray",
     xlim = c(2000000, 27000000), type = "n")
lines(extractChr$window_start, extractChr$avg_wc_fst_sakFM, col = "gray")
abline(tempThreshold$upperThreshold, 0, col = "gray", lty = "dotted")
abline(tempThreshold$lowerThreshold, 0, col = "gray", lty = "dotted")

# points the intersexual FST outlier
upperTD <- filter(extractChr, avg_wc_fst_sakFM > thresholdDF$upperThreshold[3])
points(upperTD$window_start, upperTD$avg_wc_fst_sakFM, col = "gray30")

# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(27500000, 30850000), c(-0.0, -0.0), lwd = 2)

dev.off()
