# ------------------------------------------------------------
# population_statistics_summary.R
# Tajima's D statistics
# 2021/05/13
# update: 2021/07/26
# Shingo Fujimoto
# ------------------------------------------------------------
library("tidyverse")
library("ggplot2")
source("D:/GoogleDrive/lib/population_statistics_summary.R") # self-made functions to make the 

# ------------------------------------------------------------------------------
# 1, import pi, Fst dataset (calculated by pixy)
# ------------------------------------------------------------------------------
# pi plot, O. sakaizumii (N = 40)
setwd("//wsl$/Ubuntu/home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/genomescan/pixy_out/medakaWGS_sak40_pi")
piDF <- read.table("./pixy_pi.txt", header = TRUE, sep = "\t")
tempPopStatDF <- extractPopPi("sak")

N <- length(tempPopStatDF$avg_pi)
sPiThreshold99 <- sort(tempPopStatDF$avg_pi)[c(0.005*N, 0.995*N)] # 99%
g4 <- manhattanPlotPi(tempPopStatDF)
g4 <- g4 + geom_abline(intercept = sPiThreshold99[1], slope = 0, linetype = "dotted")
g4 <- g4 + geom_abline(intercept = sPiThreshold99[2], slope = 0, linetype = "dotted")
g4
g4 <- ggplotGrob(g4)

# ------------------------------------------------------------------------------
# 2, import Tajima's D dataset (calculated by vcf-kit)
# ------------------------------------------------------------------------------
setwd("//wsl$/Ubuntu/home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/genomescan/tajimaD/nmiss0")

# import the tajima's D window
tempPopName <- "sakaizumii_tajima"
sakStatDF <- modTajimaDF(tempPopName)
TajimaPlot <- manhattanPlotTajima(sakStatDF)
#sakStatDF <- subset(sakStatDF, nSNPs > 100)

# 3 SD in raw window distribution
N <- length(sakStatDF$TajimaD)
threshold99 <- sort(sakStatDF$TajimaD)[c(0.005*N, 0.995*N)] # 99%
lowerD <- threshold99[1]
upperD <- threshold99[2]
TajimaPlot <- TajimaPlot + geom_hline(yintercept = upperD, linetype = "dotted")
TajimaPlot <- TajimaPlot + geom_hline(yintercept = lowerD, linetype = "dotted")
TajimaPlot
g3 <- ggplotGrob(TajimaPlot)

# ------------------------------------------------------------------------------
# 3, Inter-sexual Fst plot, O. sakaizumii(Male: N = 20, female: N = 20)
# ------------------------------------------------------------------------------
setwd("//wsl$/Ubuntu/home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/genomescan/pixy_out/medakaWGS_sak40_intersfst")
interSexualfstDF <- read.table("./pixy_fst.txt", header = TRUE, sep = "\t")
tempPopStatDF <- extractPopIntersexualFst("sak")

tempPopStatDF <- subset(tempPopStatDF, fst_no_sites > 0)
N <- length(tempPopStatDF$avg_sex_fst)
sFstThreshold99 <- sort(tempPopStatDF$avg_sex_fst)[c(0.005*N, 0.990*N)] # 99%
g5 <- manhattanPlotSexualFst(tempPopStatDF)
g5 <- g5 + geom_abline(intercept = sFstThreshold99[2], slope = 0, linetype = "dotted")
g5 <- ggplotGrob(g5)

# write the genome scan figure
# ------------------------------------------------------------------------------
setwd("//wsl$/Ubuntu/home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/genomescan")
g <- rbind(g4, g3, g5, g4, g4,  size = "first")
g$widths = grid::unit.pmax(g4$widths, g3$widths, g5$widths, g1$widths, g2$widths)
# plot(g)

outputDIR <- "//wsl$/Ubuntu/home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/genomescan"
ggsave(filename = paste0(outputDIR, "/genomeScan_Pi_tajimaD_iFst_sak40samples_dxy_fst_saklat.png"), g, width = 8, height = 8)

