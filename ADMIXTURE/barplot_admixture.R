# ------------------------------------------------------------------
# barplot_admixture.R
#
# Date: 2025/04/09
# Author: Shingo Fujimoto
# ------------------------------------------------------------------
# import the sample names from the text file
sampleID <- read.table("./sampleNames.txt")[1] # Sample name list in ADMIXTURE result
sampleID <- as.vector(as.matrix(sampleID))

prefix <- "population_structure_ADMIXTURE"

par(mfrow = c(5,1))
admixtureQ <- read.table(paste0(prefix, ".2.Q"))
barplot(t(as.matrix(admixtureQ)), col = rainbow(length(admixtureQ)),
        xlab="Individual #", ylab="Ancestry (2)", border = NA,
        names.arg = sampleID, las = 2, cex.names = 0.6)

admixtureQ <- read.table(paste0(prefix, ".3.Q"))
barplot(t(as.matrix(admixtureQ)), col = rainbow(length(admixtureQ)),
        xlab="Individual #", ylab="Ancestry (3)", border = NA,
        names.arg = sampleID, las = 2, cex.names = 0.6)

admixtureQ <- read.table(paste0(prefix, ".4.Q"))
barplot(t(as.matrix(admixtureQ)), col = rainbow(length(admixtureQ)),
        xlab="Individual #", ylab="Ancestry (4)", border = NA,
        names.arg = sampleID, las = 2, cex.names = 0.6)

admixtureQ <- read.table(paste0(prefix, ".5.Q"))
barplot(t(as.matrix(admixtureQ)), col = rainbow(length(admixtureQ)),
        xlab="Individual #", ylab="Ancestry (5)", border = NA,
        names.arg = sampleID, las = 2, cex.names = 0.6)

admixtureQ <- read.table(paste0(prefix, ".6.Q"))
barplot(t(as.matrix(admixtureQ)), col = rainbow(length(admixtureQ)),
        xlab="Individual #", ylab="Ancestry (6)", border = NA,
        names.arg = sampleID, las = 2, cex.names = 0.6)
