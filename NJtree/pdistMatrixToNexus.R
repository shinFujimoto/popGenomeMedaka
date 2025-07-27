# ------------------------------------------------------------------
# pdistMatrixToNexus.R
#
# convert from vcf to genetic distance matrix
# 2022/06/15, 
# Author: Shingo Fujimoto
# ------------------------------------------------------------------
library("phangorn")

# get the p-distance matrix
dist.fn <- commandArgs(trailingOnly = TRUE)[1] # p-distance matrix file name from command line
pDisMatrix <- read.table(dist.fn, header = FALSE)
sampleNameList <- pDisMatrix[,1]
pDisMatrix <- pDisMatrix[, -1] # delete the species Name column

# get the file name list
#sample.fn <- commandArgs(trailingOnly = TRUE)[2] # sampleName list file from command line
#sampleNameList <- read.table(sample.fn, header = FALSE)[,1]

# output the distance matrix as the nexus format
outputDist <- as.matrix(pDisMatrix)
rownames(outputDist) <- sampleNameList
colnames(outputDist) <- sampleNameList

# make the output file name to delete the extension
outputFileName <- substr(dist.fn, 1, nchar(dist.fn) -5)
phangorn::writeDist(outputDist,
                    file = paste(getwd(), "/", outputFileName, ".nex", sep = ""),
                    format = "nexus")
