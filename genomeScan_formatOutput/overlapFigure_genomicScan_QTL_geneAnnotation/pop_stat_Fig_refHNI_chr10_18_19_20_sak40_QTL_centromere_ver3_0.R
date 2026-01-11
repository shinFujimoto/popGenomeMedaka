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
# Anal fin length, QTL region on chromosome 10
# import and modify Tajima's D, inter-sexual Fst files
# ------------------------------------------------------------------------------
# 1, import Tajima's D dataset (calculated by vcf-kit)
pdf(file = "Chr10_TajimaD_sexFst.pdf", width = 6.75, height = 6.75)

extractChr <- subset(sakStatDF, chromosome == 10)
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
lines(c(4700000, 22340000), c(-3, -3), col = "slategray2", lwd = 2)
text(x = mean(c(4700000, 22340000)), y = -3.6, labels = "Anal fin length", cex = 0.8)

# Gene and centromere position on chromosome 18
# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(4530000, 6550000), c(-3.2, -3.2), lwd = 2)

# ------------------------------------------------------------------------------
# 2, import inter-sexual Fst dataset (calculated by pixy)
extractChr <- subset(sakStatDF, chromosome == 10)
extractChr <- dplyr::select(extractChr, chromosome, window_start, window_end, avg_wc_fst_sakFM) %>% distinct()
tempThreshold <- filter(thresholdDF, statistics == "avg_wc_fst_sakFM")

# plot focusing chromosome
plot(extractChr$window_start, extractChr$avg_wc_fst_sakFM, las = 1,
     xlab = "Positions in chromosome 10", ylab = "Inter-sexual Fst", ylim = c(-0.02, 0.10), col = "gray",
     xlim = c(2000000, 27000000), type = "n")
lines(extractChr$window_start, extractChr$avg_wc_fst_sakFM, col = "gray")
abline(tempThreshold$upperThreshold, 0, col = "gray", lty = "dotted")
abline(tempThreshold$lowerThreshold, 0, col = "gray", lty = "dotted")

# points the intersexual FST outlier
upperTD <- filter(extractChr, avg_wc_fst_sakFM > thresholdDF$upperThreshold[3])
points(upperTD$window_start, upperTD$avg_wc_fst_sakFM, col = "gray30")

# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(4530000, 6550000), c(-0.022, -0.022), lwd = 2)

dev.off()

# ------------------------------------------------------------------------------
# Female spawning latency, QTL region on chromosome 18
# import and modify Tajima's D, inter-sexual Fst files
# ------------------------------------------------------------------------------
# 1, import Tajima's D dataset (calculated by vcf-kit)
pdf(file = "Chr18_TajimaD_sexFst.pdf", width = 6.75, height = 6.75)

extractChr <- subset(sakStatDF, chromosome == 18)
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
# See also, QTL region information in Table S9, (Fujimoto et al. 2024, Mol. Ecol.)
lines(c(3010000, 12430000), c(-3, -3), col = "slategray2", lwd = 2)
text(x = mean(c(3010000, 12430000)), y = -3.6, labels = "Spawning latency", cex = 0.8)

# Gene and centromere position on chromosome 18
# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(16710000, 19070000), c(-3, -3), lwd = 2)

# ------------------------------------------------------------------------------
# 2, import inter-sexual Fst dataset (calculated by pixy)
extractChr <- subset(sakStatDF, chromosome == 18)
extractChr <- dplyr::select(extractChr, chromosome, window_start, window_end, avg_wc_fst_sakFM) %>% distinct()
tempThreshold <- filter(thresholdDF, statistics == "avg_wc_fst_sakFM")

# plot focusing chromosome
plot(extractChr$window_start, extractChr$avg_wc_fst_sakFM, las = 1,
     xlab = "Positions in chromosome 18", ylab = "Inter-sexual Fst", ylim = c(-0.02, 0.10), col = "gray",
     xlim = c(2000000, 27000000), type = "n")
lines(extractChr$window_start, extractChr$avg_wc_fst_sakFM, col = "gray")
abline(tempThreshold$upperThreshold, 0, col = "gray", lty = "dotted")
abline(tempThreshold$lowerThreshold, 0, col = "gray", lty = "dotted")

# points the intersexual FST outlier
upperTD <- filter(extractChr, avg_wc_fst_sakFM > thresholdDF$upperThreshold[3])
points(upperTD$window_start, upperTD$avg_wc_fst_sakFM, col = "gray30")

# ------------------------------------------------------------------------------
# Sex associated region on chromosome 18
# See also, Table S10 (Kitano et al. 2023, Am. Nat.)
lines(c(9980000, 25220000), c(0.08, 0.08), col = "darkorange", lwd = 2)
text(x = mean(c(9980000, 25220000)), y = 0.074, labels = "Minor sex association", cex = 0.8)	

# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(16710000, 19070000), c(-0.022, -0.022), lwd = 2)

# BMP 15
targetGene <- filter(extractChr, window_start %in% c(13000001, 15900001, 18300001))
points(targetGene$window_start, targetGene$avg_wc_fst_sakFM, col = "gray30")
text(x = targetGene$window_start, y = targetGene$avg_wc_fst_sakFM, labels = "bmp15", cex = 0.8)

dev.off()

# ------------------------------------------------------------------------------
# Anal fin length, QTL region on chromosome 19
# import and modify Tajima's D and Pi files
# 1, import Tajima's D dataset (calculated by vcf-kit)
# ------------------------------------------------------------------------------
pdf(file = "Chr19_TajimaD_sexFst.pdf", width = 6.75, height = 6.75)

extractChr <- subset(sakStatDF, chromosome == 19)
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
lines(c(1910000, 12100000), c(-3, -3), col = "slategray2", lwd = 2)
text(x = mean(c(0, 10100000)), y = -3.6, labels = "Sexual dimorphism (AFL)", cex = 0.8)

# papilally process
lines(c(5630000, 14580000), c(-3.2, -3.2), col = "slategray2", lwd = 2)
text(x = mean(c(9030000, 14580000)), y = -3.8, labels = "Sexual dimorphism(PP)", cex = 0.8)

# Gene and centromere position on chromosome 19
# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(13580000, 15595000), c(-3.3, -3.3), lwd = 2)

# ------------------------------------------------------------------------------
# 2, import pi Fst dataset (calculated by pixy)
extractChr <- subset(sakStatDF, chromosome == 19)
extractChr <- dplyr::select(extractChr, chromosome, window_start, window_end, avg_pi) %>% distinct()
tempThreshold <- filter(thresholdDF, statistics == "avg_pi")

# plot focusing chromosome
plot(extractChr$window_start, extractChr$avg_pi, las = 1,
     xlab = "Positions in chromosome 19", ylab = "Nucleotide diversity(ƒÎ)", ylim = c(0, 0.016), col = "gray",
     xlim = c(2000000, 27000000), type = "n")
lines(extractChr$window_start, extractChr$avg_pi, col = "gray")
abline(tempThreshold$upperThreshold, 0, col = "gray", lty = "dotted")
abline(tempThreshold$lowerThreshold, 0, col = "gray", lty = "dotted")

# points the Pi outlier
upperTD <- filter(extractChr, avg_pi > thresholdDF$upperThreshold[1])
points(upperTD$window_start, upperTD$avg_pi, col = "gray30")
lowerTD <- filter(extractChr, avg_pi < thresholdDF$lowerThreshold[1])
points(lowerTD$window_start, lowerTD$avg_pi, col = "gray30")

# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(13580000, 15595000), c(-0.0, -0.0), lwd = 2)

dev.off()


# ------------------------------------------------------------------------------
# Male courtship frequency, QTL region on chromosome 20
# import and modify Tajima's D and Pi files
# 1, import Tajima's D dataset (calculated by vcf-kit)
# ------------------------------------------------------------------------------
pdf(file = "Chr20_TajimaD_pi.pdf", width = 6.75, height = 6.75)

extractChr <- subset(sakStatDF, chromosome == 20)
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
points(upperTD$window_start, upperTD$TajimaD, col = "gray30")
lowerTD <- filter(extractChr, TajimaD < thresholdDF$lowerThreshold[2])
points(lowerTD$window_start, lowerTD$TajimaD, col = "gray30")

# ------------------------------------------------------------------------------
# See also, QTL region information in Table S9, (Kawajiri et al. 2014, 2015)
# ------------------------------------------------------------------------------
lines(c(13080000, 23440000), c(-3, -3), col = "slategray2", lwd = 2)
text(x = mean(c(13080000, 23440000)), y = -3.6, labels = "Male courtship frequency", cex = 0.8)

# Gene and centromere position on chromosome 19
# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(2920000, 7230000), c(-3.3, -3.3), lwd = 2)

# ------------------------------------------------------------------------------
# 2, import pi Fst dataset (calculated by pixy)
# ------------------------------------------------------------------------------
extractChr <- subset(sakStatDF, chromosome == 20)
extractChr <- dplyr::select(extractChr, chromosome, window_start, window_end, avg_pi) %>% distinct()
tempThreshold <- filter(thresholdDF, statistics == "avg_pi")

# plot focusing chromosome
plot(extractChr$window_start, extractChr$avg_pi, las = 1,
     xlab = "Positions in chromosome 20", ylab = "Nucleotide diversity(ƒÎ)", ylim = c(0, 0.016), col = "gray",
     xlim = c(2000000, 27000000), type = "n")
lines(extractChr$window_start, extractChr$avg_pi, col = "gray")
abline(tempThreshold$upperThreshold, 0, col = "gray", lty = "dotted")
abline(tempThreshold$lowerThreshold, 0, col = "gray", lty = "dotted")

# points the Pi outlier
upperTD <- filter(extractChr, avg_pi > thresholdDF$upperThreshold[1])
points(upperTD$window_start, upperTD$avg_pi, col = "gray30")
lowerTD <- filter(extractChr, avg_pi < thresholdDF$lowerThreshold[1])
points(lowerTD$window_start, lowerTD$avg_pi, col = "gray30")

# Centromere information: See also, Table Table S10 (Ansai et al. 2023, Mol. Ecol.)
lines(c(2920000, 7230000), c(0, 0), lwd = 2)

dev.off()



