# --------------------------------------------------
# convert_pos_cM.R
#
# physical position information convert to cM
# Author: Shinfo Fujimoto
# Date: 2024/02/09
# --------------------------------------------------
# input file: CHROM, position
fileName <- commandArgs(trailingOnly = TRUE)[1]
physicalPosDF <- read.table(fileName, header = FALSE, encoding="UTF-8")

posTocM <- physicalPosDF$V2 / 1000000
outputDF <- data.frame(chr = physicalPosDF$V1,
                       locusID = paste0(physicalPosDF$V1, "_", physicalPosDF$V2),
                       genetic_pos = posTocM,
                       physical_pos = physicalPosDF$V2
                       )

outputFileName <- paste0(gsub("\\.[0-9A-Za-z]+$", "", fileName), "_cM_pos.map")
write.table(outputDF, file = outputFileName,
            row.names = FALSE, col.names = FALSE, quote = FALSE,
            fileEncoding="UTF-8")



# imoport iHS output
iHS_DF <- read.table("sak40_biallelic_snps_CP020797.1_ihs.ihs.out", header = FALSE)

hist(iHS_DF$V6)
iHSthreshold_chr19 <- sort(iHS_DF$V6)[c(0.005*length(iHS_DF$V6), 0.995*length(iHS_DF$V6))]

plot(iHS_DF$V2, iHS_DF$V6, col = "gray", type = "l")
abline(iHSthreshold_chr19[1], 0, col = "gray", lty = "dotted")
abline(iHSthreshold_chr19[2], 0, col = "gray", lty = "dotted")

filter(iHS_DF, iHS_DF$V6 < iHSthreshold_chr19[1])
filter(iHS_DF, iHS_DF$V6 > iHSthreshold_chr19[2])

