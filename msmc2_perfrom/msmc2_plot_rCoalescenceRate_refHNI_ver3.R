# ------------------------------------------------------------------
# MSMC2 plot 
# Date: 2023/06/07, 
# Update: 2024/05/13, 
# Update: 2024/06/27, 2ind2ind-8haplotype dataset 
# Author: Shingo Fujimoto
# ------------------------------------------------------------------
# setwd("D:/GoogleDrive/study/ò_ï∂å¥çe/0_medakaPopGenomics/DatasetAndScript/msmc2/msmc2_8ind16haplotype")
setwd("//wsl$/Ubuntu/home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/msmc2/msmc2input/2ind2ind_8haplotype")

# ---------------------------------------------------------------
# 1. import the MSMC2 result of relative cross coalescent rate
# ---------------------------------------------------------------
inputPath <- "//wsl$/Ubuntu/home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/msmc2/msmc2input/2ind2ind_8haplotype"
txtFileNameList <- list.files(path = inputPath, recursive = TRUE) # MSMC2 combined result files
txtFileNameList <- txtFileNameList[grep(pattern = "combined", txtFileNameList)]

msmc2CombinedTable <- data.frame()

for(ii in 1:length(txtFileNameList))
{
  # obtain population and location combination
  fileName <- txtFileNameList[ii]
  splitfileName <- strsplit(txtFileNameList[ii], split = "-")
  pop1pop2Names <- splitfileName[[1]][1] # population name
  loc1loc2Name <- strsplit(splitfileName[[1]][2], split = "/")[[1]][1] # location name
  
  # modify the msmc2 combined result table
  tempDF <- read.table(txtFileNameList[ii], header=TRUE)
  tempDF <- cbind(tempDF, pop1pop2Names, loc1loc2Name, fileName)
  
  msmc2CombinedTable <- rbind(msmc2CombinedTable, tempDF)
}

# Exclude Inf value in the right time boundary
msmc2CombinedTable <- subset(msmc2CombinedTable,
                             msmc2CombinedTable$right_time_boundary != "Inf")
# Exclude zero value in the time_index
msmc2CombinedTable <- subset(msmc2CombinedTable,
                             msmc2CombinedTable$time_index != "0")

# ---------------------------------------------------------------------------------
# 2. Calculate the summary table of population diverergenet generations
#    based on the relative cross coalescent rate
# ---------------------------------------------------------------------------------
rCCRMinDivTimeDF <- data.frame()
mu <- 3.5e-9 # mutation rate
gen <- 1 # generation time 
for(xx in 1:length(txtFileNameList))
{
  crossPopDat <- subset(msmc2CombinedTable, fileName == txtFileNameList[xx])
  
  # Calculate the relative cross coalescent rate
  rCCR <- (2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11))
  crossPopDat <- cbind(crossPopDat, rCCR = rCCR)
  
  # Extract minimum generation time in rCCR > 0.5
  divTimeRow <- subset(crossPopDat, rCCR > 0.5)
  min_left_time_boundary <- min(divTimeRow$left_time_boundary, na.rm = TRUE)
  print(txtFileNameList[xx])
  print(min_left_time_boundary / mu * gen)
  
  # Output data frame
  minDivTimeRow <- subset(divTimeRow, left_time_boundary == min_left_time_boundary)
  minDivTimeRow <- cbind(minDivTimeRow, minDivTime = (min_left_time_boundary / mu * gen))
  rCCRMinDivTimeDF <- rbind(rCCRMinDivTimeDF, minDivTimeRow)
}
write.csv(rCCRMinDivTimeDF, file = "rCCRMinDivTimeDF.csv", row.names = FALSE)

# ---------------------------------------------------------------
# 3. import the MSMC2 result of effective populations size
# ---------------------------------------------------------------
inputPath <- "//wsl$/Ubuntu/home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/msmc2/msmc2input/2ind2ind_8haplotype"
txtFileNameList <- list.files(path = inputPath, recursive = TRUE) # MSMC2 combined result files
txtFileNameList <- c(txtFileNameList[grep(pattern = "within1.msmc2.final", txtFileNameList)],
                     txtFileNameList[grep(pattern = "within2.msmc2.final", txtFileNameList)])

# Extract and add the population combination information
msmc2NeTable <- data.frame()
for(ii in 1:length(txtFileNameList))
{
  # obtain population and location combination
  fileName <- txtFileNameList[ii]
  splitfileName <- strsplit(txtFileNameList[ii], split = "-")
  pop1pop2Names <- splitfileName[[1]][1] # population name
  loc1loc2Name <- strsplit(splitfileName[[1]][2], split = "/")[[1]][1] # location name
  popCombination <- gsub("O_sak", "Osak", pop1pop2Names)
  targetPop <- ifelse(length(grep(pattern = "within1", fileName)) == 0,
                      strsplit(popCombination, split = "_")[[1]][2],
                      strsplit(popCombination, split = "_")[[1]][1]
                      )
  # modify the msmc2 combined result table
  tempDF <- read.table(txtFileNameList[ii], header=TRUE)
  tempDF <- cbind(tempDF, pop1pop2Names, loc1loc2Name, fileName, targetPop)
   
  msmc2NeTable <- rbind(msmc2NeTable, tempDF)
}

# Exclude Inf value in the right time boundary
msmc2NeTable <- subset(msmc2NeTable,
                       msmc2NeTable$right_time_boundary != "Inf")
# Exclude zero value in the time_index
msmc2NeTable <- subset(msmc2NeTable,
                             msmc2NeTable$time_index != "0")

# Exclude pseudo replicate by the combination
msmc2NeTable <- rbind(msmc2NeTable[(msmc2NeTable$pop1pop2Names == "EJ_WJ" & msmc2NeTable$targetPop == "EJ"),],
                      msmc2NeTable[(msmc2NeTable$pop1pop2Names == "EJ_WK" & msmc2NeTable$targetPop == "WK"),],
                      msmc2NeTable[(msmc2NeTable$pop1pop2Names == "WJ_WK" & msmc2NeTable$targetPop == "WJ"),],
                      msmc2NeTable[(msmc2NeTable$pop1pop2Names == "O_sak_WK" & msmc2NeTable$targetPop == "Osak"),]
)

# ---------------------------------------------------------------------------------
# 2. Draw the relative coalescent rate line
# ---------------------------------------------------------------------------------
pdf(file = "msmc2_rcRate_Esize_2ind2ind_8haplotype.pdf", width = 7.4, height = 10)
options(scipen = 10)

par(mfrow = c(3,2))
mu <- 3.5e-9 # Mutation rate parameter
gen <- 1 # Generation year

# function to draw the relative coalescent rate line
rcrLine <- function(crossPopDat, lineCol = "black", lineType = 1) {
  lineCol <- rgb(col2rgb(lineCol)[1] / 255,
                 col2rgb(lineCol)[2] / 255,
                 col2rgb(lineCol)[3] / 255,
                 alpha = 0.3)
  lines(crossPopDat$left_time_boundary/mu*gen,
        2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11),
        type="s", col = lineCol, lty = lineType, lwd = 1.5)
  print(crossPopDat$left_time_boundary/mu*gen)
  print(2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11))
}

# Combination within O. latipes and O. sakaizumii vs. O. latipes(WK)
popCombination <- c("EJ_WJ", "EJ_WK", "WJ_WK") # Population combination
popCombColoer <- c("orange", "red", "lightgreen") # Line colors
popLineType <- c(1, 1, 1) # Line type

# x axis and y axis in the plot 
crossPopDat <- subset(msmc2CombinedTable, fileName == msmc2CombinedTable$fileName[1])
plot(crossPopDat$left_time_boundary / mu * gen,
     2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11),
     log="x", ylim = c(0, 1), las = 1, xaxt = "n", xlim = c(1000, 10000000),
     type="n", xlab="Generations", ylab="rCCR",
     main = "within O. latipes")
axis(side = 1, at = c(1000, 10000, 100000, 1000000, 10000000))

# rCCR line in each MSMC2 results
for(xx in 1:length(popCombination))
{
  # Figure outline
  tempCombination <- popCombination[xx]
  tempLineColor <- popCombColoer[xx]
  tempLineType <- popLineType[xx]

  # draw the line of each MSMC2 results
  tempCombinationDF <- subset(msmc2CombinedTable, pop1pop2Names == tempCombination)
  tempFilenamelist <- unique(tempCombinationDF$fileName)
  
  for(ii in 1:length(tempFilenamelist))
  { crossPopDat <- subset(tempCombinationDF, fileName == tempFilenamelist[ii])
    rcrLine(crossPopDat, lineCol = tempLineColor, lineType = tempLineType) }
}

# Combination within O. sakaizumii and O. sakaizumii vs. O. latipes (WK)
popCombination <- c("O_sak1_O_sak2", "O_sak_WK") # Population combination
popCombColoer <- c("skyblue", "purple") # Line colors
popLineType <- c(1, 1) # Line type

# x axis and y axis in the plot 
crossPopDat <- subset(msmc2CombinedTable, fileName == msmc2CombinedTable$fileName[1])
plot(crossPopDat$left_time_boundary / mu * gen,
     2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11),
     log="x", ylim = c(0, 1), las = 1, xaxt = "n", xlim = c(1000, 10000000),
     type="n", xlab="Generations", ylab="rCCR",
     main = "O. sakaizumii")
axis(side = 1, at = c(1000, 10000, 100000, 1000000, 10000000))

# rCCR line in each MSMC2 results
for(xx in 1:length(popCombination))
{
  # Figure outline
  tempCombination <- popCombination[xx]
  tempLineColor <- popCombColoer[xx]
  tempLineType <- popLineType[xx]
  
  # draw the line of each MSMC2 results
  tempCombinationDF <- subset(msmc2CombinedTable, pop1pop2Names == tempCombination)
  tempFilenamelist <- unique(tempCombinationDF$fileName)
  
  for(ii in 1:length(tempFilenamelist))
  { crossPopDat <- subset(tempCombinationDF, fileName == tempFilenamelist[ii])
  rcrLine(crossPopDat, lineCol = tempLineColor, lineType = tempLineType) }
}

# ---------------------------------------------------------------------------------
# 2. effective population size plot
# ---------------------------------------------------------------------------------
epsLine <- function(x, lineCol = "black", lineType = 1) {
  lineCol <- rgb(col2rgb(lineCol)[1] / 255,
                 col2rgb(lineCol)[2] / 255,
                 col2rgb(lineCol)[3] / 255,
                 alpha = 0.4)
  lines(x$left_time_boundary / mu * gen,
        (1 / x$lambda) / (2 * mu),
        type="s", col = lineCol, lty = lineType, lwd = 1.5)
}

# Effective population size
# Combination within O. latipes and O. sakziumii vs. WK
popCombination <- c("EJ_WJ", "WJ_WK", "EJ_WK") # Population combination
popCombColoer <- c("orange", "lightgreen", "red") # Line colors
popLineType <- c(1, 1, 1) # Line type

# x axis and y axis in the plot 
msmc2within1DF <- subset(msmc2NeTable, fileName == msmc2NeTable$fileName[1])
plot(msmc2within1DF$left_time_boundary / mu*gen, (1 / msmc2within1DF$lambda) / (2 * mu),
     log = "x", ylim = c(0, 6000000), xaxt = "n", xlim = c(1000, 1000000),
     type = "n", xlab = "Generations", ylab = "Population Size (Ne)",
     main = "")
axis(side = 1, at = c(1000, 10000, 100000, 1000000, 10000000))

# Ne line in each MSMC2 results
for(xx in 1:length(popCombination))
{
  # Figure outline
  tempCombination <- popCombination[xx]
  tempLineColor <- popCombColoer[xx]
  tempLineType <- popLineType[xx]
  
  # draw the line of each MSMC2 results
  tempCombinationDF <- subset(msmc2NeTable, pop1pop2Names == tempCombination)
  tempFilenamelist <- unique(tempCombinationDF$fileName)
  
  for(ii in 1:length(tempFilenamelist))
  { tempNeDF <- subset(tempCombinationDF, fileName == tempFilenamelist[ii])
  epsLine(tempNeDF, lineCol = tempLineColor, lineType = tempLineType) }
}

# within O. sakaizumii
popCombination <- c("Trg", "Nig", "Kwb", "Ktg") # Population combination
popCombColoer <- c("purple", "skyblue", "lightgray", "lightgray") # Line colors
popLineType <- c(1, 1, 2, 3) # Line type

# x axis and y axis in the plot 
msmc2within1DF <- subset(msmc2NeTable, fileName == msmc2NeTable$fileName[1])
plot(msmc2within1DF$left_time_boundary / mu*gen, (1 / msmc2within1DF$lambda) / (2 * mu),
     log = "x", ylim = c(0, 500000), xaxt = "n", xlim = c(1000, 1000000),
     type = "n", xlab = "Generations", ylab = "Population Size (Ne)",
     main = "")
axis(side = 1, at = c(1000, 10000, 100000, 1000000, 10000000))

# Ne line in each MSMC2 results
for(xx in 1:length(popCombination))
{
  # Figure outline
  tempCombination <- popCombination[xx]
  tempLineColor <- popCombColoer[xx]
  tempLineType <- popLineType[xx]
  
  # draw the line of each MSMC2 results
  tempCombinationDF <- msmc2NeTable[grep(pattern = tempCombination, msmc2NeTable$loc1loc2Name),]
  tempFilenamelist <- unique(tempCombinationDF$fileName)
  
  for(ii in 1:length(tempFilenamelist))
  { tempNeDF <- subset(tempCombinationDF, fileName == tempFilenamelist[ii])
  epsLine(tempNeDF, lineCol = tempLineColor, lineType = tempLineType) }
}

# Box plot for estimation of divergence times
rCCRMinDivTimeDF$pop1pop2Names <- factor(rCCRMinDivTimeDF$pop1pop2Names,
                                         levels = c("O_sak_WK", "WJ_WK", "EJ_WK",  "EJ_WJ", "O_sak1_O_sak2"))
boxplot(minDivTime ~ pop1pop2Names, data = rCCRMinDivTimeDF,
        xlab = "Populations", ylab = "Separation generation", log = "y", yaxt = "n",
        ylim = c(1000, 10000000), outline = FALSE)
axis(side = 2, at = c(1000, 10000, 100000, 1000000, 10000000),
     labels = c("10^3", "10^4", "10^5", "10^6", "10^7"), las = 1)


dev.off()
