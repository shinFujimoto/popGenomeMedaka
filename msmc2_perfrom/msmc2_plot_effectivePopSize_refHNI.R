# ------------------------------------------------------------------
# MSMC2 plot 
# Date: 2023/06/07, 
# Update: 2024/05/13, 
# Update: 2024/06/27, 2ind2ind-8haplotype dataset 
# Author: Shingo Fujimoto
# ------------------------------------------------------------------
# setwd("D:/GoogleDrive/study/ò_ï∂å¥çe/0_medakaPopGenomics/DatasetAndScript/msmc2/msmc2_8ind16haplotype")
setwd("//wsl$/Ubuntu/home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/msmc2/msmc2input/4ind4ind_16haplotype")

# ---------------------------------------------------------------------------------
# Draw the effective population size
# ---------------------------------------------------------------------------------
options(scipen = 10) # plot parameter for axis

# Parameters to convert from generations to years
mu <- 3.5e-9
gen <- 1

# function to draw the effective population size 
epsLine <- function(x, lineCol = "black", lineType = 1) {
  tmpMSMC2table <- read.table(x, header = TRUE)
  # Exclude zero generation row
  tmpMSMC2table <- subset(tmpMSMC2table, (tmpMSMC2table$left_time_boundary != 0))
  logGen <- log10(tmpMSMC2table$left_time_boundary / mu) * gen # generations -> years
  logNe <- log10((1 / tmpMSMC2table$lambda) / (2 * mu)) # Effective population size
  
  lineCol <- rgb(col2rgb(lineCol)[1] / 255,
                 col2rgb(lineCol)[2] / 255,
                 col2rgb(lineCol)[3] / 255,
                 alpha = 0.4)
  lines(logGen, logNe,
        type="s", col = lineCol, lty = lineType, lwd = 2)
}

# Effective population size
pdf(file = "msmc2_Ne_4ind_6haplotype.pdf", width = 7.4, height = 10)
par(mfrow = c(3,2))
msmc2resultDF <- read.table("sakWK_4ind4ind/within1.msmc2.final.txt", header = TRUE)
msmc2resultDF <- subset(msmc2resultDF, (msmc2resultDF$left_time_boundary != 0))

plot(log10(msmc2resultDF$left_time_boundary / mu*gen),
     log10((1 / msmc2resultDF$lambda) / (2 * mu)),
     xaxt = "n", xlab = "Generations",
     xlim = c(2.8, 7), ylim = c(3, 7), las = 1,
     ylab = "Effective Population Size", type = "n",
     main = "")
axis(side = 1,
     at = c(1, 2, 3, 4, 5, 6, 7),
     labels = c(10, 100, 1000, 10000, 100000, 1000000, 10000000))

epsLine("sakWK_4ind4ind/within1.msmc2.final.txt", lineCol = "purple", lineType = 1) # O. sakaizumii, 
epsLine("sakWK_4ind4ind/within2.msmc2.final.txt", lineCol = "red", lineType = 1) # O. latipes, West Kyushu
epsLine("WJEJ_4ind4ind/within1.msmc2.final.txt", lineCol = "green", lineType = 1) # O. latipes, West Japan
epsLine("WJEJ_4ind4ind/within2.msmc2.final.txt", lineCol = "orange", lineType = 1) # O. latipes, East Japan

dev.off()

