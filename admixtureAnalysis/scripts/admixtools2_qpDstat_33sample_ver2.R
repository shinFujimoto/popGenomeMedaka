# ------------------------------------------------------------------
# admixtools2_qpDstat_33sample.R
#
# admixture graph estimation using admixtools2
# 2025/04/12, 
# Author: Shingo Fujimoto
# ------------------------------------------------------------------
library(admixtools)
library(ggplot2)
library(tidyverse)
library(igraph)

# ABBA-BABA test (Patterson's D statistics)
# -------------------------------------------------------------------
# (EJ|WJ|WK), (Osak, HSOK)
snpDataPrefix <- "../genotype/medakaWGS_HSOK_sak_lat_33samples"

# Calculate D statistics
pop1 <- c("westJ", "eastJ", "eastJ", "eastJ", "eastJ", "westJ") 
pop2 <- c("wKyu", "wKyu", "westJ", "westJ", "wKyu", "wKyu") 
pop3 <- c("sak", "sak", "sak", "wKyu", "westJ", "eastJ") 
pop4 <- c("eKor") 
F4stats <- qpdstat(snpDataPrefix, pop1, pop2, pop3, pop4,
                   unique_only = TRUE, auto_only = FALSE, comb = FALSE, f4mode = FALSE)
# Extract and change the row order for figure output
F4stats_output <- F4stats[c(1, 2, 3, 4), ]
write.csv(F4stats_output, file = "../qpdstat/Dstat_EJWJWK_Osak_HSOK.csv", row.names = FALSE)

# -------------------------------------------------------------------
# (EJ, (WJ|WK)), ((Trg, Nig, Kwb, Ktg), HSOK)
snpDataPrefix <- "../genotype/EJ_WJWK_Osak-sep_HSOK/medakaWGS_HSOK_sak_lat_33samples"

# Calculate D statistics
pop1 <- c(rep("eastJ", 8)) 
pop2 <- c(rep("westJ", 4), rep("wKyu", 4)) 
pop3 <- c(rep(c("Trg", "Nig", "kwb", "ktg"), 2)) 
pop4 <- c(rep("eKor", 8))
F4stats <- qpdstat(snpDataPrefix, pop1, pop2, pop3, pop4,
                   unique_only = TRUE, auto_only = FALSE, comb = FALSE, f4mode = FALSE)
# Extract and change the row order for figure output
F4stats_output <- F4stats[c(4, 3, 2, 1), ]
write.csv(F4stats_output, file = "../qpdstat/Dstat_EJ_WJWK_Osak-sep_HSOK.csv", row.names = FALSE)

# -------------------------------------------------------------------
# (Iso, Csg, Mry, Ich, WJ|WK), (Osak, HSOK)
snpDataPrefix <- "../genotype/EJsep_WJWK_Osak_HSOK/medakaWGS_HSOK_sak_lat_33samples"

# Calculate D statistics
pop1 <- c(rep(c("iso", "Csg", "mry", "Ich"), 3)) 
pop2 <- c(rep("westJ", 4), rep("wKyu", 4), rep("westJ", 4)) 
pop3 <- c(rep("sak", 8), rep("wKyu", 4))
pop4 <- c(rep("eKor", 12))
F4stats <- qpdstat(snpDataPrefix, pop1, pop2, pop3, pop4,
                   unique_only = TRUE, auto_only = FALSE, comb = FALSE, f4mode = FALSE)
# Extract and change the row order for figure output
F4stats_output <- F4stats[c(1, 2, 3, 4), ]
write.csv(F4stats_output, file = "../qpdstat/Dstat_EJsep_WJWK_Osak_HSOK.csv", row.names = FALSE)

pop1 <- c(rep("iso", 9))
pop2 <- c(rep(c("Csg", "mry", "Ich"), 3))
pop3 <- c(rep("sak", 3), rep("wKyu", 3), rep("westJ", 3))
pop4 <- c(rep("eKor", 9))
F4stats <- qpdstat(snpDataPrefix, pop1, pop2, pop3, pop4,
                   unique_only = TRUE, auto_only = FALSE, comb = FALSE, f4mode = FALSE)
# Extract and change the row order for figure output
F4stats_output <- F4stats[c(1, 2, 3, 7, 8, 9), ]
write.csv(F4stats_output, file = "../qpdstat/Dstat_iso_EJsep_WJOsak_HSOK.csv", row.names = FALSE)

# -------------------------------------------------------------------
# (EJ, (Ehm, Eky, Myk, Miy, Tsm)), (Osak, HSOK)
snpDataPrefix <- "../genotype/EJ_WJsep_Osak_HSOK/medakaWGS_HSOK_sak_lat_33samples"

# Calculate D statistics
pop1 <- c(rep(c("Tsm", "Miy", "Ehm", "Eky", "Myk"), 3))
pop2 <- c(rep("eastJ", 5), rep("wKyu", 5), rep("eastJ", 5)) 
pop3 <- c(rep("sak", 10), rep("wKyu", 5)) 
pop4 <- c(rep("eKor", 15))
F4stats <- qpdstat(snpDataPrefix, pop1, pop2, pop3, pop4,
                   unique_only = TRUE, auto_only = FALSE, comb = FALSE, f4mode = FALSE)
# Extract and change the row order for figure output
F4stats_output <- F4stats[c(9, 8, 6, 7, 10), ]
write.csv(F4stats_output, file = "../qpdstat/Dstat_WJsep_WK_Osak_HSOK.csv", row.names = FALSE)

# -------------------------------------------------------------------
# (EJ, Gin,Izm,Nag,Hid), (Osak, HSOK)
snpDataPrefix <- "../genotype/EJ_WKsep_Osak_HSOK/medakaWGS_HSOK_sak_lat_33samples"

# Calculate D statistics
pop1 <- c(rep("eastJ", 4), rep("westJ", 4))
pop2 <- c(rep(c("Gin", "Izm", "Nag", "Hid"), 2))
pop3 <- c(rep("sak", 8))
pop4 <- c(rep("eKor", 8))
F4stats <- qpdstat(snpDataPrefix, pop1, pop2, pop3, pop4,
                   unique_only = TRUE, auto_only = FALSE, comb = FALSE, f4mode = FALSE)
# Extract and change the row order for figure output
F4stats_output <- F4stats[c(5, 6, 8, 7), ]
write.csv(F4stats_output, file = "../qpdstat/Dstat_WJ_WKsep_Osak_HSOK.csv", row.names = FALSE)
