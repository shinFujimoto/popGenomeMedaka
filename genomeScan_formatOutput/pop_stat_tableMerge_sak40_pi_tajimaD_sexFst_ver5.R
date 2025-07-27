# ------------------------------------------------------------
# pop_stat_tableMerge_sak40_pi_tajimaD_sexFst
# Tajima's D statistics
# 2021/05/13
# update: 2021/07/26
# update: 2024/06/04
# update: 2024/10/08, iHS score
# update: 2025/01/24, O. sakaizumii(N = 40), pi, tajima's D, sex Fst 
# update: 2025/07/08, QQplot for intesexual FST
#
# Shingo Fujimoto
# ------------------------------------------------------------
library("tidyverse")
library("ggplot2")
#source("D:/GoogleDrive/lib/population_statistics_summary.R") # self-made functions to make the 
setwd("C:/Users/fujim/GoogleDrive/study/ò_ï∂å¥çe/0_medakaPopGenomics/DatasetAndScript/genomescan")
setwd("D:/GoogleDrive/study/ò_ï∂å¥çe/0_medakaPopGenomics/DatasetAndScript/genomescan")
source("formatting_output/scripts/population_statistics_summary.R") # self-made functions

# ------------------------------------------------------------------------------
# 1, import pi dataset (calculated by pixy)
# ------------------------------------------------------------------------------
# Merge pi table
piDF <- read.table("raw_output/pixy_out_pi_fst_dxy/medakaWGS_sak40_pi/pixy_pi.txt", header = TRUE, sep = "\t")
tempPopStatDF <- extractPopPi("sak")
statDF <- tempPopStatDF

# ------------------------------------------------------------------------------
# 2, import Tajima's D (calculated by vcf-kit)
# ------------------------------------------------------------------------------
# tajima D datasheet started: 0
tempPopStatDF <- modTajimaDF("raw_output/vcfkit_tajimaD/TajimaD_50000bp_sakaizumii_tajima.txt")
tempPopStatDF$window_start <- tempPopStatDF$window_start + 1
statDF <- merge(statDF, tempPopStatDF,
                by = c("chromosome", "window_start", "window_end"), all = TRUE)

# ------------------------------------------------------------------------------
# 3, import sex-FST between latF-latM and sakF-sakM dataset (calculated by pixy)
# ------------------------------------------------------------------------------
# Merge sakF-sakM fst table
fstDF <- read.table("raw_output/pixy_out_pi_fst_dxy/medakaWGS_sak40_intersfst/pixy_fst_20241025.txt", header = TRUE, sep = "\t")
tempPopStatDF <- extractPopFst("female", "male")
names(tempPopStatDF)[4:5] <- c("avg_wc_fst_sakFM", "fst_no_sites_sakFM")
statDF <- merge(statDF, tempPopStatDF,
                by = c("chromosome", "window_start", "window_end"), all = TRUE)

# ------------------------------------------------------------------------------
# 4, import fd dataset for HSOK-WK-WJ-EJ (calculated by Dsuite)
# ------------------------------------------------------------------------------
# fd: f_d_DF is added the window size (> 50000bp)
f_d_DF <- read.table("raw_output/Dsuite_f_d/eastJ_westJ_sak_localFstats__50_50.txt", header = TRUE, sep = "\t")

# convert chromosome ID -> chromosome number
modchromosomeID <- as.factor(f_d_DF$chr)
levels(modchromosomeID) <- seq(1, 24, 1)
f_d_DF$chr <- modchromosomeID

intervalStart <- (trunc(f_d_DF$windowStart / 50000) * 50000) + 1
intervalEnd <- (ceiling(f_d_DF$windowEnd / 50000) * 50000)
fd_modDF <- cbind(f_d_DF, intervalStart, intervalEnd)
fd_modDF <- dplyr::select(fd_modDF, chr, windowStart, windowEnd, f_d, intervalStart, intervalEnd)

statDF <- merge(statDF, fd_modDF,
                  by.x = c("chromosome", "window_start"),
                  by.y = c("chr", "intervalStart"),
                  all = TRUE)

# f_d interval is larger than 50000 bp 
currentInterval <- statDF[13,] # initial interval of fd window
for(ii in 14:length(statDF$chromosome))
{
  tempRow <- statDF[ii,]
  if(is.na(tempRow$windowStart))
  {
    # Add the interval information and fd value
    statDF$windowStart[ii] <- currentInterval$windowStart
    statDF$windowEnd[ii] <- currentInterval$windowEnd
    statDF$f_d[ii] <- currentInterval$f_d
    statDF$intervalEnd[ii] <- currentInterval$intervalEnd
  } else
  {
    # move the next interval
    currentInterval <- tempRow
  }
}

# ------------------------------------------------------------------------------
# 5, Exclude low sequence windows and output the CSV file
# ------------------------------------------------------------------------------
# less than 10000 bp windows was excluded, based on the (sak = 40 dataset)
outputStatDF <- subset(statDF, pi_no_sites > 10000)

# Exclude the duplicated windows of Fd records
outputStatDF <- dplyr::distinct(outputStatDF, chromosome, window_start, .keep_all = TRUE)

# output the merged table
outputFileName <- "formatting_output/sak40_pi_tajimaD_sexFst_f_d/O_sak40_pi_tajimaD_sexFst_f_d.csv"
write.csv(outputStatDF, file = outputFileName, row.names = FALSE)

# Threshold table
N <- length(outputStatDF$avg_pi)
thresoldDF <- data.frame()
# Average pi
thresoldDF <- rbind(thresoldDF,
                    data.frame( statistics = "avg_pi",
                                lowerThreshold = sort(outputStatDF$avg_pi)[c(0.005*N, 0.995*N)][1],
                                upperThreshold = sort(outputStatDF$avg_pi)[c(0.005*N, 0.995*N)][2]
                                )
                    )

# Average Tajima's D
thresoldDF <- rbind(thresoldDF,
                    data.frame( statistics = "Tajima_D",
                                lowerThreshold = sort(outputStatDF$TajimaD)[c(0.005*N, 0.995*N)][1],
                                upperThreshold = sort(outputStatDF$TajimaD)[c(0.005*N, 0.995*N)][2]
                    )
)

# Average wc fst sakFM
N <- length(outputStatDF$avg_pi)
thresoldDF <- rbind(thresoldDF,
                    data.frame( statistics = "avg_wc_fst_sakFM",
                                lowerThreshold = sort(outputStatDF$avg_wc_fst_sakFM)[c(0.005*N, 0.995*N)][1],
                                upperThreshold = sort(outputStatDF$avg_wc_fst_sakFM)[c(0.005*N, 0.995*N)][2]
                    )
)

# f_d
Fd_DF <- select(outputStatDF, chromosome, windowStart, windowEnd, f_d)
Fd_DF <- distinct(Fd_DF) # Exclude duplicated rows 

N <- length(Fd_DF$f_d)
thresoldDF <- rbind(thresoldDF,
                    data.frame( statistics = "f_d",
                                lowerThreshold = sort(Fd_DF$f_d)[c(0.005*N, 0.995*N)][1],
                                upperThreshold = sort(Fd_DF$f_d)[c(0.005*N, 0.995*N)][2]
                    )
)

# output the thresholdTable
outputFileName <- "formatting_output/sak40_pi_tajimaD_sexFst_f_d/threshold0005.csv"
write.csv(thresoldDF, file = outputFileName, row.names = FALSE)

