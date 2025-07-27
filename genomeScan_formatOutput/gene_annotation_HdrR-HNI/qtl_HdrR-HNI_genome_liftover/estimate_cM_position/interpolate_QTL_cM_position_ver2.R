# ------------------------------------------------------------
# interpolate_QTL_cM_position.R
#
# interpolate the genetic marker position
# from centiMorgan to physical position in reference genome
# using generalized additive model
# 
# Author:Shingo Fujimoto
# Date: 2023/11/14
# update: 2024/03/22
# ------------------------------------------------------------
library("mgcv")
library("tidyverse")

# snp data convert to physical distance
setwd("D:/GoogleDrive/study/ò_ï∂å¥çe/0_behaviorQTL/ê}ï\/figureAndDatasets/submit_DataAndScript/gene_Info_ensembl/qtl_HdrR_genome_liftover/estimate_cM_position")

# data import
cM_POS_table <- read.csv("QTLmarker_cM_position.csv")
cM_POS_table$Strain <- as.factor(cM_POS_table$Strain)

# ------------------------------------------------------------------
# 1. Interpolate the centiMorgan in QTL analysis to the physical positions
# ------------------------------------------------------------------
interpolateResultDF <- data.frame()

for(ii in 1:24)
{
  tempChr <- ii
  tempSubset <- subset(cM_POS_table, CHROM_qtl == tempChr)
  
  # estimation curve between centiMorgan
  # and physical distance in HdrR reference (ASM223467v1)
  # -----------------------------------------------------------------
  pos <- tempSubset$POS_ASM223467v1
  centiMorgan <- tempSubset$cM
  tempSubset <- cbind(tempSubset, pos, centiMorgan)
  model <- gam(pos ~ s(centiMorgan, by = Strain, k = 4), data = tempSubset, family = gaussian)
  
  # Centimorgan for AFOM strain dataset
  HdrR_AFOM_samplex <- seq(trunc(min(tempSubset$cM[tempSubset$Strain == "AFOM"])),
                           ceiling(max(tempSubset$cM[tempSubset$Strain == "AFOM"])), 1)
  HdrR_AFOM_pred <- predict(model, newdata = list(centiMorgan = HdrR_AFOM_samplex,
                                                  Strain = rep("AFOM", length(HdrR_AFOM_samplex))))
  
  outputRows <- data.frame(Strain = rep("AFOM", length(HdrR_AFOM_samplex)),
                           reference = rep("HdrR_ASM223467v1", length(HdrR_AFOM_samplex)),
                           CHROM = rep(tempChr, length(HdrR_AFOM_samplex)),
                           cM = HdrR_AFOM_samplex,
                           POS_pred = HdrR_AFOM_pred
  )
  interpolateResultDF <- rbind(interpolateResultDF, outputRows)
  
  # Centimorgan for OFAM strain dataset
  HdrR_OFAM_samplex <- seq(trunc(min(tempSubset$cM[tempSubset$Strain == "OFAM"])),
                           ceiling(max(tempSubset$cM[tempSubset$Strain == "OFAM"])), 1)
  HdrR_OFAM_pred <- predict(model, newdata = list(centiMorgan = HdrR_OFAM_samplex,
                                                  Strain = rep("OFAM", length(HdrR_OFAM_samplex))))
  outputRows <- data.frame(Strain = rep("OFAM", length(HdrR_OFAM_samplex)),
                           reference = rep("HdrR_ASM223467v1", length(HdrR_OFAM_samplex)),
                           CHROM = rep(tempChr, length(HdrR_OFAM_samplex)),
                           cM = HdrR_OFAM_samplex,
                           POS_pred = HdrR_OFAM_pred
  )
  interpolateResultDF <- rbind(interpolateResultDF, outputRows)
}

write.csv(interpolateResultDF, file = "interporated_cM_physicalPosition.csv", row.names = FALSE)

# ------------------------------------------------------------
# 2. Estimation curve between centiMorgan and physical distance in reference genome
# ------------------------------------------------------------
pdf(width = 8, height =8, file = "interpolate_cM_pos_OFAM_AFOM.pdf")
par(mfrow = c(2,2))
options(scipen = 1000)
for(ii in 1:24)
{
  tempChr <- ii
  tempSubset <- subset(cM_POS_table, CHROM_qtl == tempChr)

  # estimation curve between centiMorgan and physical distance in HdrR reference
  pos <- tempSubset$POS_ASM223467v1
  centiMorgan <- tempSubset$cM
  tempSubset <- cbind(tempSubset, pos, centiMorgan)

  model <- gam(pos ~ s(centiMorgan, by = Strain, k = 4), data = tempSubset, family = gaussian)
  summary(model)

  # plot the curve
  plot(centiMorgan, pos,
       main = paste0("HdrR, Chr ", tempChr), xlab = "Centi morgan (cM)",
       ylab = "Physical position (bp)", las = 1, pch = as.numeric(tempSubset$Strain))

  samplex <- seq(min(tempSubset$cM[tempSubset$Strain == "AFOM"]),
                 max(tempSubset$cM[tempSubset$Strain == "AFOM"]), 0.1)
  predValue <- predict.gam(model,
                           newdata = list(centiMorgan = samplex, Strain = rep("AFOM", length(samplex))))
  lines(samplex, predValue, lwd =2)
  samplex <- seq(min(tempSubset$cM[tempSubset$Strain == "OFAM"]),
                 max(tempSubset$cM[tempSubset$Strain == "OFAM"]), 0.1)
  predValue <- predict.gam(model,
                           newdata = list(centiMorgan = samplex, Strain = rep("OFAM", length(samplex))))
  lines(samplex, predValue, lwd =2, lty = "dashed")
}
dev.off()

# ------------------------------------------------------------
# 3. Merge QTL regions, 95% CI(cM), and position information
# ------------------------------------------------------------
# data import for QTL regions
qtlregionTable <- read.csv("qtl_candidateRegion.csv")

merge95CI_resultDF <- data.frame()
for(ii in 1:length(qtlregionTable$Strain))
{
  tempRow <- qtlregionTable[ii,]
  
  HdrR_ASM223467v1_95CI_start <- filter(interpolateResultDF,
                          reference == "HdrR_ASM223467v1",
                          Strain == tempRow$Strain[1],
                          CHROM == tempRow$CHROM[1],
                          cM == tempRow$cM95CI_start[1]
  )[,5]
  
  HdrR_ASM223467v1_95CI_end <- filter(interpolateResultDF,
                        reference == "HdrR_ASM223467v1",
                        Strain == tempRow$Strain[1],
                        CHROM == tempRow$CHROM[1],
                        cM == tempRow$cM95CI_end[1]
  )[,5]
  print(ii)
  merge95CI_resultDF <- rbind(merge95CI_resultDF,
                              cbind(tempRow,
                                    HdrR_ASM223467v1_95CI_start,
                                    HdrR_ASM223467v1_95CI_end)
                              )
}

write.csv(merge95CI_resultDF, file = "QTL_95CI_interpolate_cMtoPOS.csv", row.names = FALSE)

# modify the BED like format 
write.table(data.frame(merge95CI_resultDF$CHROM,
                       round(merge95CI_resultDF$HdrR_ASM223467v1_95CI_start),
                       round(merge95CI_resultDF$HdrR_ASM223467v1_95CI_end),
                       merge95CI_resultDF$TraitName),
            file = "QTL_95CI_interpolate_cMtoPOS.bed", sep = "\t",
            row.names = FALSE, col.names = FALSE
            )
