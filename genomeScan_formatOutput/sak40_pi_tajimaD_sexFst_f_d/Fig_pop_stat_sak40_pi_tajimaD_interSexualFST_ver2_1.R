# ------------------------------------------------------------
# population_statistics_summary.R
# Tajima's D statistics
# 2021/05/13
# update: 2021/07/26
# Shingo Fujimoto
# ------------------------------------------------------------
library("tidyverse")
library("ggplot2")
source("../../scripts/population_statistics_summary.R") # self-made functions to make the 

# import population genetic statistics
popGenDF <- read.csv("O_sak40_pi_tajimaD_sexFst_f_d.csv")

# ------------------------------------------------------------------------------
# 1, import pi, Fst dataset (calculated by pixy)
# ------------------------------------------------------------------------------
# pi plot, O. sakaizumii (N = 40)
# Extract the nucleotide diversity (Pi) data and modified for figure plot function
piDF <- select(popGenDF, chromosome, window_start, window_end, avg_pi, pi_no_sites, count_diffs)
tempPopStatDF <- piDF

# Top and bottome Threshold information
N <- length(tempPopStatDF$avg_pi)
sPiThreshold99 <- sort(tempPopStatDF$avg_pi)[c(0.005*N, 0.995*N)] # 99%

# Manhattan plot of nucleotide diversity
g1 <- manhattanPlotPi(tempPopStatDF)
g1 <- g1 + geom_abline(intercept = sPiThreshold99[1], slope = 0, color = "gray")
g1 <- g1 + geom_abline(intercept = sPiThreshold99[2], slope = 0, color = "gray")
g1
g1 <- ggplotGrob(g1)

# ------------------------------------------------------------------------------
# 2, import Tajima's D dataset (calculated by vcf-kit)
# ------------------------------------------------------------------------------
tajimaDF <- select(popGenDF, chromosome, window_start, window_end, TajimaD, Tajima_no_sites, nSNPs)

# Top and bottome Threshold information
N <- length(tajimaDF$TajimaD)
threshold99 <- sort(tajimaDF$TajimaD)[c(0.005*N, 0.995*N)] # 99%
lowerD <- threshold99[1]
upperD <- threshold99[2]

# import the tajima's D window
g2 <- manhattanPlotTajima(tajimaDF)

# 3 SD in raw window distribution
g2 <- g2 + geom_hline(yintercept = upperD, color = "gray")
g2 <- g2 + geom_hline(yintercept = lowerD, color = "gray")
g2
g2 <- ggplotGrob(g2)

# ------------------------------------------------------------------------------
# 3, Inter-sexual Fst plot, O. sakaizumii(Male: N = 20, female: N = 20)
# ------------------------------------------------------------------------------
interSexualfstDF <- select(popGenDF, chromosome, window_start, window_end, avg_wc_fst_sakFM, fst_no_sites_sakFM)
names(interSexualfstDF)[4] <- "avg_wc_fst"

# Top threshold information
N <- length(interSexualfstDF$avg_wc_fst)
sFstThreshold99 <- sort(interSexualfstDF$avg_wc_fst)[c(0.005*N, 0.990*N)] # 99%

g3 <- manhattanPlotSexualFst(interSexualfstDF)
g3 <- g3 + geom_hline(yintercept = sFstThreshold99[2], color = "gray")
g3
g3 <- ggplotGrob(g3)

# write the genome scan figure
# ------------------------------------------------------------------------------
g <- rbind(g1, g2, g3, g1, g1,  size = "first")
g$widths = grid::unit.pmax(g1$widths, g2$widths, g3$widths, g1$widths, g1$widths)
# plot(g)

outputDIR <- "./"
ggsave(filename = paste0(outputDIR, "genomeScan_Pi_tajimaD_iFst_sak40samples.png"), g, width = 8, height = 8)

