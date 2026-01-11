# --------------------------------------------------
# permutationTestOverlapQTL_selectionScan.R
#
# Permutation test for the overlapping between
# QTL regions and selection scan windows
# Author: Shinfo Fujimoto
# Date: 2026/01/11
# --------------------------------------------------

#BiocManager::install("regioneR", force = TRUE)
library(regioneR)
library(GenomicRanges)

# --------------------------------------------------
# 1. Chrmosome length dummy data for permutation test
# dummy chromosome name: chr[number]
inputGenome <- toGRanges(data.frame(
  chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
          "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
          "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
          "chr19", "chr20", "chr21", "chr22", "chr23", "chr24"),
  start = c(rep(1, 24)),
  # Oryzias sakaizumii, HNI strain reference sequence length
  end = c(34611496, 21867142, 35737580, 31648822, 31424114, 29208093,
          32810049, 24535370, 30885934, 28873095, 25539739, 28246893,
          29950042, 28133707, 28784337, 30693862, 28809336, 27766529,
          22538199, 24186181, 29727532, 27485737, 22812639, 21356977)
  )
)

# --------------------------------------------------
# 2. definition of function for permutation test
# ntimes: number of trials in permutation (1000 times)
# randomize.function: "randomizeRegions" option will shufle using by the all chromosomes
# evaluate.function: "numOverlaps" return the overlapped number,
#                    "overlapTotalLength" was calculated the overlapped base pair length

overlapPermTest <- function(A, B)
{
  
  print("Performing the permutation test...") 

  # "overlapTotalLength" was calculated the overlapped base pair length
  overlapTotalLength <- function(A, B, ...) {
    ov <- intersect(A, B) # 重複部分を抽出
    return(sum(width(ov))) # 重複部分の合計長さを返す
  }
  
  # SetA: Target regions of QTLs
  # SetB: Target windows of selection scans
  pt <- permTest(A = SetA, ntimes = 1000,
                 randomize.function = randomizeRegions,
                 evaluate.function = overlapTotalLength,
                 genome = inputGenome,
                 B = SetB,
                 mask = NULL
  )
  
  # Show the calculation result
  print(pt)
  summary(pt)
  
  # Visualization of the distribution of permutation
  # and probability of the observed overlap
  plot(pt)
}  

# -------------------------------------------------------------
# SetA: Anal fin length QTL (chr 10)
# SetB: high genetic diversity
# -------------------------------------------------------------
# A: QTL region
SetA <- toGRanges(data.frame(
  chr = c("chr10"),
  start = c(470000),
  end = c(22340000)
)
)

#B: Selection scan windows (High SexFst, High Tajima's D)
SetB <- toGRanges(data.frame(
  chr = c("chr10", "chr18"),
  start = c(13850001, 26450001),
  end = c(13900000, 26500000)
)
)

overlapPermTest(SetA, SetB)
# -------------------------------------------------------------

# SetA: Female mate choice QTL (chr 18)
# SetB: high genetic diversity
# -------------------------------------------------------------
# A: QTL region
SetA <- toGRanges(data.frame(
  chr = c("chr18"),
  start = c(3010000),
  end = c(12430000)
)
)

#B: Selection scan windows (High Pi, High Tajima's D)
SetB <- toGRanges(data.frame(
  chr = c("chr8", "chr10", "chr10", "chr13",
          "chr18", "chr18", "chr19", "chr20"),
  start = c(14350001, 13750001, 13850001, 18550001,
            4250001, 6850001, 6750001, 21300001),
  end = c(14400000, 13900000, 13900000, 18600000,
          4300000, 6950000, 6800000, 21350000)
)
)

overlapPermTest(SetA, SetB)

# -------------------------------------------------------------
# SetA: Anal fin length QTL (chr 19)
# SetB: high genetic diversity
# -------------------------------------------------------------
# A: QTL region
SetA <- toGRanges(data.frame(
  chr = c("chr19"),
  start = c(1910000),
  end = c(12100000)
)
)

#B: Selection scan windows (High Pi, High Tajima's D)
SetB <- toGRanges(data.frame(
  chr = c("chr8", "chr10", "chr10", "chr13",
          "chr18", "chr18", "chr19", "chr20"),
  start = c(14350001, 13750001, 13850001, 18550001,
            4250001, 6850001, 6750001, 21300001),
  end = c(14400000, 13900000, 13900000, 18600000,
          4300000, 6950000, 6800000, 21350000)
)
)

overlapPermTest(SetA, SetB)
# -------------------------------------------------------------

# -------------------------------------------------------------
# SetA: Male courtship frequency (chr 20)
# SetB: high genetic diversity
# -------------------------------------------------------------
# A: QTL region
SetA <- toGRanges(data.frame(
  chr = c("chr20"),
  start = c(13080000),
  end = c(23440000)
)
)

#B: Selection scan windows (High Pi, High Tajima's D)
SetB <- toGRanges(data.frame(
  chr = c("chr8", "chr10", "chr10", "chr13",
          "chr18", "chr18", "chr19", "chr20"),
  start = c(14350001, 13750001, 13850001, 18550001,
            4250001, 6850001, 6750001, 21300001),
  end = c(14400000, 13900000, 13900000, 18600000,
          4300000, 6950000, 6800000, 21350000)
)
)

overlapPermTest(SetA, SetB)
# -------------------------------------------------------------
