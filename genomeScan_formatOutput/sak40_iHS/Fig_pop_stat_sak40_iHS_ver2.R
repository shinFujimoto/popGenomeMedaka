# ------------------------------------------------------------
# pop_stat_tableMerge_sak40_pi_tajimaD_sexFst
# Tajima's D statistics
# 2025
# 2025/04/14, Haplotype estimation: O. sakaizumii(N = 40) & O. latipes(N = 42) 
#             iHS score: O. sakaizumii(N = 40)
# 
# Author: Shingo Fujimoto
# ------------------------------------------------------------
library("tidyverse")
library("ggplot2")
library("MASS")
#setwd("D:/GoogleDrive/study/ò_ï∂å¥çe/0_medakaPopGenomics/DatasetAndScript/genomescan")
setwd("D:/GoogleDrive/study/ò_ï∂å¥çe/0_medakaPopGenomics/DatasetAndScript/deposit_Dryad_20251112-")
#setwd("C:/Users/fujim//GoogleDrive/study/ò_ï∂å¥çe/0_medakaPopGenomics/DatasetAndScript/genomescan")
#setwd("/mnt/d/GoogleDrive/study/ò_ï∂å¥çe/0_medakaPopGenomics/DatasetAndScript/genomescan")
source("scripts/population_statistics_summary.R") # self-made functions

# ------------------------------------------------------------------------------
# 1, import iHS dataset (calculated by selscan)
# ------------------------------------------------------------------------------
# iHS table
iHSDF <- read.table("Fig2-3_genomeScan_rawOutput/selscan_iHS/sak40_biallelic_snps_allchr_withpmap.ihs.out", header = FALSE)
names(iHSDF) <- c("locus", "physPOS", "1_freq", "ihh_1", "ihh_0", "ihs")

# Z-transform of iHS
z_trans_ihs <- (iHSDF$ihs - mean(iHSDF$ihs)) / sd(iHSDF$ihs)
hist(abs(z_trans_ihs))

# Normal distribution parameters in Z-transformed iHS
fitParam <- fitdistr(z_trans_ihs, densfun = "normal")
pvalue_ihs <- pnorm(z_trans_ihs, mean = fitParam$estimate[1], sd = fitParam$estimate[2])
# p_adjusted <- p.adjust(pvalue_ihs, method = "fdr") #
# hist(p_adjusted)

# Outlier of z_trans_ihs
N <- length(z_trans_ihs) # Total number of SNPs
iHS_threshold <- sort(z_trans_ihs)[c(0.0025*N, 0.9975*N)] # Thresholds of outlier (lower and upper 0.5%)

# Extract and output the outlier iHS, SNPs
chrNum <- as.numeric(factor(substr(iHSDF$locus, 1, 10)))
iHSDF <- cbind(iHSDF, z_trans_ihs, pvalue_ihs, chrNum)
outlier_iHSDF <- filter(iHSDF, z_trans_ihs < iHS_threshold[1] | z_trans_ihs > iHS_threshold[2])
outputFileName <- "Fig2-3_genomeScan_formatOutput/sak40_iHS/O_sak40_iHS_outlier_05percentile.csv"
write.csv(outlier_iHSDF, file = outputFileName, row.names = FALSE)

# Figure plot for iHS
# ggplot object
g1 <- manhattanPlotiHS(iHSDF)
g1 <- g1 + geom_hline(yintercept = max(abs(iHS_threshold)), color = "gray")
g1

# ------------------------------------------------------------------------------
# write the genomescan window  figure
# ------------------------------------------------------------------------------
g1 <- ggplotGrob(g1)
g <- rbind(g1, g1, g1, g1, g1,  size = "first")
g$widths = grid::unit.pmax(g1$widths, g1$widths, g1$widths, g1$widths, g1$widths)
plot(g)

outputDIR <- paste0(getwd(), "/Fig2-3_genomeScan_formatOutput/sak40_iHS")
ggsave(filename = paste0(outputDIR, "/genomeScan_sak40_lat42_iHS_05percentile.png"), g, width = 8, height = 8)
ggsave(filename = paste0(outputDIR, "/genomeScan_sak40_lat42_iHS_05percentile.pdf"), g, width = 8, height = 8)


# extract n SNPs in a chromosomal region
iHS_summaryDF <- function(x)
{
  nSNPsTable <- filter(iHSDF, chrNum == x$V1, between(physPOS, x$V2, x$V3))
  outlier_nSNPs_Table <- filter(outlier_iHSDF, chrNum == x$V1, between(physPOS, x$V2, x$V3))
  
  # Extract maximum iHS (z-transformed and absolute value)
  if(length(nSNPsTable$chrNum) == 0)
  {
    return(
      outputRow <- data.frame(chrNum = x$V1[1],
                              pysPOS_start = min(x$V2),
                              pysPOS_end = max(x$V3),
                              nSNPs = length(nSNPsTable$chrNum),
                              nSNPs_iHS_outlier =  length(outlier_nSNPs_Table$chrNum),
                              max_iHS_score = NA,
                              max_iHS_position = NA)
    )
  } else
  {
    if(length(outlier_nSNPs_Table$chrNum) == 0)
    {
      max_iHS_score <- max(abs(nSNPsTable$z_trans_ihs))
      max_iHS_position <- nSNPsTable$physPOS[abs(nSNPsTable$z_trans_ihs) == max_iHS_score]
    } else
    {
      max_iHS_score <- max(abs(outlier_nSNPs_Table$z_trans_ihs))
      max_iHS_position <- outlier_nSNPs_Table$physPOS[abs(outlier_nSNPs_Table$z_trans_ihs) == max_iHS_score]
    }
    
    outputRow <- data.frame(
                            chrNum = x$V1[1],
                            pysPOS_start = min(x$V2),
                            pysPOS_end = max(x$V3),
                            nSNPs = length(nSNPsTable$chrNum),
                            nSNPs_iHS_outlier =  length(outlier_nSNPs_Table$chrNum),
                            max_iHS_score = max_iHS_score,
                            max_iHS_position = max_iHS_position
    )

    return(outputRow)
  }
}

# ----------------------------------------------------------
# Calculate iHS snps in Fd outlier regions (Table. 2)
# ----------------------------------------------------------
fdRegion <- read.table("D:/GoogleDrive/study/ò_ï∂å¥çe/0_medakaPopGenomics/DatasetAndScript/genomescan/formatting_output/gene_Info_ensembl/sak8_HSOK_latWJ8EJ8_F_d/genomescan_windows_sak8_HSOK_latWJ8EJ8_F_d.bed")
as.data.frame(t(sapply(1:length(fdRegion$V1), function(x) { iHS_summaryDF(fdRegion[x,]) })))

# Concatenate some regions within 1M base pairs
tempChrRegion <- data.frame(V1 = 4, V2 = 26139070, V3 = 26809986)
iHS_summaryDF(tempChrRegion)
tempChrRegion <- data.frame(V1 = 19, V2 = 5386798, V3 = 6992578)
iHS_summaryDF(tempChrRegion)
tempChrRegion <- data.frame(V1 = 21, V2 = 3916687, V3 = 4201210)
iHS_summaryDF(tempChrRegion)
tempChrRegion <- data.frame(V1 = 21, V2 = 18275900, V3 = 18285814)
iHS_summaryDF(tempChrRegion)

# Calculate Pi, Tajimas's D, intersexual FST, outiler regions
tempChrRegion <- data.frame(V1 = 1, V2 = 16050001, V3 = 16650000)
iHS_summaryDF(tempChrRegion)
tempChrRegion <- data.frame(V1 = 8, V2 = 14350001, V3 = 14400000)
iHS_summaryDF(tempChrRegion)
tempChrRegion <- data.frame(V1 = 10, V2 = 13750001, V3 = 13900000)
iHS_summaryDF(tempChrRegion)
tempChrRegion <- data.frame(V1 = 10, V2 = 13850001, V3 = 13900000)
iHS_summaryDF(tempChrRegion)

tempChrRegion <- data.frame(V1 = 18, V2 = 4250001, V3 = 4300000)
iHS_summaryDF(tempChrRegion)
tempChrRegion <- data.frame(V1 = 18, V2 = 6850001, V3 = 6900000)
iHS_summaryDF(tempChrRegion)
tempChrRegion <- data.frame(V1 = 18, V2 = 6900001, V3 = 6950000)
iHS_summaryDF(tempChrRegion)
tempChrRegion <- data.frame(V1 = 18, V2 = 17800001, V3 = 17850000)
iHS_summaryDF(tempChrRegion)
tempChrRegion <- data.frame(V1 = 18, V2 = 26450001, V3 = 26500000)
iHS_summaryDF(tempChrRegion)

tempChrRegion <- data.frame(V1 = 19, V2 = 6750001, V3 = 6800000)
iHS_summaryDF(tempChrRegion)

tempChrRegion <- data.frame(V1 = 20, V2 = 21300001, V3 = 21350000)
iHS_summaryDF(tempChrRegion)


tempChrRegion <- data.frame(V1 = 23, V2 = 14885667, V3 = 16241702)
iHS_summaryDF(tempChrRegion)



# ----------------------------------------------------------

# 
candidateRegions <- read.table("clipboard", header = TRUE, sep = "\t")
for(ii in 1:length(candidateRegions$Chr))
{
  tempChrRegion <- data.frame(V1 = candidateRegions$Chr[ii],
                              V2 = candidateRegions$Start[ii],
                              V3 = candidateRegions$End..Mbp.[ii])
  print(iHS_summaryDF(tempChrRegion))
}

# plot the chromosomal region in chromosome 23
options(scipen = 10000)
tempChr <- filter(iHSDF, chrNum == 23, between(physPOS, 5000000, 18000000))
plot(tempChr$physPOS, tempChr$z_trans_ihs,
     las = 1, ylim = c(-4, 4), 
     xlab = "Physical position on chromosome 23 (bp)",
     ylab = "Z transformed iHS", type = "l")

points(tempChr$physPOS[tempChr$z_trans_ihs > iHS_threshold[2]],
       tempChr$z_trans_ihs[tempChr$z_trans_ihs > iHS_threshold[2]],
       pch = 16, col = "black")
points(tempChr$physPOS[tempChr$z_trans_ihs < iHS_threshold[1]],
       tempChr$z_trans_ihs[tempChr$z_trans_ihs < iHS_threshold[1]],
       pch = 16, col = "black")

abline(iHS_threshold[1], 0, lty = "dotted")
abline(iHS_threshold[2], 0, lty = "dotted")


# iHS outlier clustered region
lines(c(15752639, 16791743), c(-3.5, -3.5), col = "skyblue", lwd = 2)

# position of Leptin-b
lines(c(15675227, 15676114), c(-3.7, -3.7), lwd = 5)
text(c(15675227), c(-3.9), labels = "Leptin-B", cex = 0.75) 
