# -------------------------------------------------------------------
# pop_stat_Fig_sak8_HSOK_latWJ8EJ8_F_D
#
# plot the F_d value
# Date: 2025/02/04
# Author: Shingo Fujimoto
# -------------------------------------------------------------------
library("tidyverse")
library("ggplot2")
#setwd("D:/GoogleDrive/study/論文原稿/0_medakaPopGenomics/DatasetAndScript/genomescan")
setwd("C:/users/fujim/GoogleDrive/study/論文原稿/0_medakaPopGenomics/DatasetAndScript/genomescan")
source("formatting_output/scripts/population_statistics_summary.R") # self-made functions

# Extract more than 10000bp sequenced windows from merged tables
mergedDF <- read.csv("formatting_output/sak40_pi_tajimaD_sexFst_f_d/O_sak40_pi_tajimaD_sexFst_f_d.csv")
mergedDF$chromosome <- factor(mergedDF$chromosome)

# threshold of each pop-stat
thresholdDF <- read.csv("formatting_output/sak40_pi_tajimaD_sexFst_f_d/threshold0005.csv")

# ------------------------------------------------------------------------------
# 1, import F_d dataset (calculated by Dsuite)
# ------------------------------------------------------------------------------
# F_d plot (HSOK, Osak8, WJ8, EJ8)
Fd_DF <- select(mergedDF, chromosome, windowStart, windowEnd, f_d)
names(Fd_DF) <- c("chromosome", "window_start", "window_end", "fd")
        
# Exclude duplicated rows 
Fd_DF <- distinct(Fd_DF)

# Patterson's D has defined (BABA - ABBA), but the Dsuite returund (ABBA - BABA)
Fd_DF$fd <- -1 * Fd_DF$fd # Reverse the sign

# threshold of each pop-stat
tempThreshold <- filter(thresholdDF, statistics == "f_d")

g1 <- manhattanPlotfd(Fd_DF)
g1 <- g1 + geom_hline(yintercept = -1 * tempThreshold[,2], color = "gray")
g1 <- g1 + geom_hline(yintercept = -1 * tempThreshold[,3], color = "gray")
g2 <- g1
g2 <- ggplotGrob(g2)

g1 <- g1 + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), strip.text.x = element_blank())
g1 <- ggplotGrob(g1)
g1

# ------------------------------------------------------------------------------
# write the genomescan window  figure
# ------------------------------------------------------------------------------
g <- rbind(g1, g1, g2,  size = "first")
g$widths = grid::unit.pmax(g1$widths, g1$widths, g2$widths)
plot(g)

outputDIR <- paste0(getwd(), "/formatting_output/sak8_HSOK_latWJ8EJ8_F_d")
ggsave(filename = paste0(outputDIR, "/genomeScan_sak8_HSOK_latWJ8EJ8_F_d.png"), g, width = 6.75, height = 4.17)
ggsave(filename = paste0(outputDIR, "/genomeScan_sak8_HSOK_latWJ8EJ8_F_d.pdf"), g, width = 6.75, height = 4.17)

