# ------------------------------------------------------------------
# admixtools2_interCladeAdmix_dxyDirectionTestt.R
#
# admixture graph estimation using admixtools2
# 2025/04/13, 
# Author: Shingo Fujimoto
# ------------------------------------------------------------------
library(admixtools)
library(ggplot2)
library(tidyverse)
library(igraph)

# -------------------------------------------------------------------
snpDataPrefix <- "../genotype/medakaWGS_HSOK_sak_lat_33samples"
popList <- c("eKor", "wKyu", "westJ", "eastJ", "sak")
outDIRname <- paste0('../genotype/admixtools2_medakaWGS_', paste(popList, collapse = '_'))

# calculate only at the first time
# F2blocks <- extract_f2(snpDataPrefix, outdir = outDIRname,pops = popList, auto_only = FALSE, n_cores = 4, overwrite = TRUE)
# import the pre calculated F2
F2_blocks <- f2_from_precomp(outDIRname)

# 1. Calculate Outgroup F3 statistics
# -------------------------------------------------------------------
pop1 <- c("eKor") 
pop2 <- c("sak", "wKyu", "westJ") 
pop3 <- c("wKyu", "westJ", "eastJ") 
outgroupF3 <- f3(snpDataPrefix, pop1, pop2, pop3,
                 unique_only = TRUE, outgroupmode = TRUE, apply_corr = FALSE)
# Row order modified for figure output
outgroupF3_output <- outgroupF3[c(9, 6, 5, 3, 2, 1),] 
write.csv(outgroupF3_output, file = "../outgroupF3/outgroupF3.csv", row.names = FALSE)
