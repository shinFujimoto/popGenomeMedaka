# ------------------------------------------------------------
# population_statistics_summary.R
# convert pixy and vcftools output of the genome scan
# Pi, Fst, and Tajima's D statistics
# 
# 2021/05/13
# update: 2021/07/26
# Shingo Fujimoto
# ------------------------------------------------------------
library("tidyverse")
library("ggplot2")

# import fst data calculated by vcftools
importWeirFst <- function(dataset)
{
        modDF <- data.frame(chromosome = dataset$CHROM,
                            window_start = dataset$BIN_START,
                            window_end = dataset$BIN_END,
                            avg_wc_fst = dataset$MEAN_FST,
                            weighted_wc_fst = dataset$WEIGHTED_FST,
                            fst_no_sites = dataset$N_VARIANTS
                            )
        modDF$chromosome <- as.factor(modDF$chromosome)
        levels(modDF$chromosome) <- seq(1, 24, 1)
        return(modDF)
}

# extract a population from pi file
extractPopPi <- function(x)
{
        popName <- x
        tempPopDF <- subset(piDF, pop == popName)
        
        modPiDF <- data.frame(chromosome = tempPopDF$chromosome,
                              window_start = tempPopDF$window_pos_1,
                              window_end = tempPopDF$window_pos_2,
                              avg_pi = tempPopDF$avg_pi,
                              pi_no_sites = tempPopDF$no_sites,
                              count_diffs = tempPopDF$count_diffs
        )
        modPiDF$chromosome <- as.factor(modPiDF$chromosome)
        levels(modPiDF$chromosome) <- seq(1, 24, 1)
        return(modPiDF)
}

# extract a population from inter-population Fst text file
extractPopFst <- function(popName1, popName2)
{
        tempFstDF <- subset(fstDF, pop1 == popName1 & pop2 == popName2)
        
        modFstDF <- data.frame(chromosome = tempFstDF$chromosome,
                               window_start = tempFstDF$window_pos_1,
                               window_end = tempFstDF$window_pos_2,
                               avg_wc_fst = as.numeric(tempFstDF$avg_wc_fst),
                               fst_no_sites = tempFstDF$no_snps
        )
        modFstDF$chromosome <- as.factor(modFstDF$chromosome)
        levels(modFstDF$chromosome) <- seq(1, 24, 1)
        return(modFstDF)
}

# extract a population from inter-population dxy text file
extractPopDxy <- function(popName1, popName2)
{
        tempDxyDF <- subset(dxyDF, pop1 == popName1 & pop2 == popName2)
        modDxyDF <- data.frame(chromosome = tempDxyDF$chromosome,
                               window_start = tempDxyDF$window_pos_1,
                               window_end = tempDxyDF$window_pos_2,
                               avg_dxy = tempDxyDF$avg_dxy,
                               dxy_no_sites = tempDxyDF$no_sites,
                               count_diffs = tempDxyDF$count_diffs
        )
        modDxyDF$chromosome <- as.factor(modDxyDF$chromosome)
        levels(modDxyDF$chromosome) <- seq(1, 24, 1)
        return(modDxyDF)
}

# Tajima's D
modTajimaDF <- function(x)
{
        
        tajimaDF <- read.table(x, header = TRUE, sep = "\t")
        tajimaDF$CHROM <- as.factor(tajimaDF$CHROM)
        levels(tajimaDF$CHROM) <- seq(1, 24, 1)
        return(data.frame(chromosome = tajimaDF$CHROM,
                   window_start = tajimaDF$BIN_START,
                   window_end = tajimaDF$BIN_END,
                   TajimaD = tajimaDF$TajimaD,
                   Tajima_no_sites = tajimaDF$N_Sites,
                   nSNPs = tajimaDF$N_SNPs)
        )
}

# Pi, Fst and Tajima's D data frame merged
modFstDxyFiles <- function(popName, dataset)
{
        tempDF <- subset(dataset, pop1 == popName)
        tempDF$chromosome <- as.factor(tempDF$chromosome)
        levels(tempDF$chromosome) <- seq(1, 24, 1)
        return(data.frame(chromosome = tempDF$chromosome,
                          window_start = tempDF$window_pos_1,
                          window_end = tempDF$window_pos_2,
                          dxy = tempDF$avg_dxy,
                          dxy_no_sites = tempDF$no_sites,
                          count_diffs = tempDF$count_diffs)
        )
}

# Pi, Fst and Tajima's D data frame merged
mergePopStatFiles <- function(popName)
{
        tempPiDF <- extractPopPi(popName)
        tempFstDF <- extractPopIntersexualFst(popName)
        tempTajimaDF <- modTajimaDF(popName)
        tempDF <- merge(tempPiDF, tempTajimaDF, by = c("chromosome", "window_end"), all.x = TRUE)
        tempDF <- merge(tempDF, tempFstDF, by = c("chromosome", "window_end"), all.x = TRUE)
        tempDF$chromosome <- as.factor(tempDF$chromosome)
        levels(tempDF$chromosome) <- seq(1, 24, 1)
        return(tempDF)
}

# draw the scatterplot between intersexual-Fst and Tajima's D
fstVsTajima <- function(dataset)
{
        # intersexual Fst vs. Tajima's D
        # Fst filtering by quantile (lower and upper 2.5%)
        lowerFst <- quantile(dataset$avg_wc_fst, probs = seq(0, 1, 0.005), na.rm = TRUE)[6]
        upperFst <- quantile(dataset$avg_wc_fst, probs = seq(0, 1, 0.005), na.rm = TRUE)[196]
        
        # Tajima's D filtering by quantile (lower and upper 2.5%)
        lowTajima <- quantile(dataset$TajimaD, probs = seq(0, 1, 0.005), na.rm = TRUE)[6]
        uppTajima <- quantile(dataset$TajimaD, probs = seq(0, 1, 0.005), na.rm = TRUE)[196]
        
        plot(dataset$TajimaD, dataset$avg_wc_fst,
             col = gray(level = 0.8), las = 1, ylim = c(min(dataset$avg_wc_fst, na.rm = TRUE), 1),
             xlab = "Tajima's D", ylab = "Average intersexual Fst")
        
        # Survival conflict region
        intersexualFstDF <- subset(dataset, dataset$avg_wc_fst > upperFst & dataset$TajimaD > uppTajima)
        intersexualFstDF
        points(dataset$TajimaD[dataset$avg_wc_fst > upperFst & dataset$TajimaD > uppTajima],
               dataset$avg_wc_fst[dataset$avg_wc_fst > upperFst & dataset$TajimaD > uppTajima],
               col = "red")
        
        # Balancing selection region
        balancingFstDF <- subset(dataset, dataset$avg_wc_fst < lowerFst & dataset$TajimaD > uppTajima)
        balancingFstDF
        points(dataset$TajimaD[dataset$avg_wc_fst < lowerFst & dataset$TajimaD > uppTajima],
               dataset$avg_wc_fst[dataset$avg_wc_fst < lowerFst & dataset$TajimaD > uppTajima],
               col = "blue")
        
        # purifing selection (No significant window)
        points(dataset$TajimaD[dataset$avg_wc_fst < lowerFst & dataset$TajimaD < lowTajima],
               dataset$avg_wc_fst[dataset$avg_wc_fst < lowerFst & dataset$TajimaD < lowTajima])
}

# make the manhattan plot
manhattanPlotTajima <- function(dataset)
{
        # even and odd number chromosomes were filled different colors
        chrom_color_group <- ifelse(as.numeric(dataset$chromosome) %% 2 == 0, "even", "odd")
        plotOut <- ggplot(dataset, aes(x = (window_start + window_end)/2,
                                             y = TajimaD, color = chrom_color_group))
        # ggplot object
        plotOut <- plotOut + geom_point(size = 0.5, alpha = 0.5, stroke = 0)
        plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
        plotOut <- plotOut + xlab("Chromsome")
        plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
        plotOut <- plotOut + ylab("Tajima's D (50000bp)")
        plotOut <- plotOut + scale_y_continuous(expand = c(0, 0), limits = c(-3,3))
        plotOut <- plotOut + facet_grid(. ~ chromosome, axis.labels = "margins",
                                        scales = "free_y", switch = "x", space = "free_x")
        plotOut <- plotOut + theme_classic()
        plotOut <- plotOut + theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   panel.spacing = unit(0.1, "cm"),
                                   strip.background = element_blank(),
                                   strip.placement = "outside",
                                   legend.position ="none")
        return(plotOut)
}

manhattanPlotNSites <- function(dataset)
{
  # even and odd number chromosomes were filled different colors
  chrom_color_group <- ifelse(as.numeric(dataset$chromosome) %% 2 == 0, "even", "odd")
  plotOut <- ggplot(dataset, aes(x = (window_start + window_end)/2,
                                 y = Tajima_no_sites, color = chrom_color_group))
  # ggplot object
  plotOut <- plotOut + geom_point(size = 0.5, alpha = 0.5, stroke = 0)
  plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
  plotOut <- plotOut + xlab("Chromsome")
  plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
  plotOut <- plotOut + ylab("Segregate sites (50000bp)")
  plotOut <- plotOut + scale_y_continuous(expand = c(0, 0))
  plotOut <- plotOut + facet_grid(. ~ chromosome,
                                  scales = "free_y", switch = "x", space = "free_x")
  plotOut <- plotOut + theme_classic()
  plotOut <- plotOut + theme(axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             panel.spacing = unit(0.1, "cm"),
                             strip.background = element_blank(),
                             strip.placement = "outside",
                             legend.position ="none")
  return(plotOut)
}

manhattanPlotAF <- function(dataset, yvalue = 0)
{
  # even and odd number chromosomes were filled different colors
  chrom_color_group <- ifelse(as.numeric(dataset$CHROM) %% 2 == 0, "even", "odd")
  plotOut <- ggplot(dataset, aes(x = startPos,
                                          y = yvalue, color = chrom_color_group))
  # ggplot object
  plotOut <- plotOut + geom_point(size = 1, alpha = 0.5, stroke = 0)
  plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
  plotOut <- plotOut + xlab("Chromsome")
  plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
  plotOut <- plotOut + ylab("Aleternative allele frequency")
  plotOut <- plotOut + scale_y_continuous(limits = c(0, 1))
  plotOut <- plotOut + facet_grid(. ~ CHROM,
                                  scales = "free_y", switch = "x", space = "free_x")
  plotOut <- plotOut + theme_classic()
  plotOut <- plotOut + theme(axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             panel.spacing = unit(0.1, "cm"),
                             strip.background = element_blank(),
                             strip.placement = "outside",
                             legend.position ="none")
  return(plotOut)
}


manhattanPlotPi <- function(dataset)
{
        # even and odd number chromosomes were filled different colors
        chrom_color_group <- ifelse(as.numeric(tempPopStatDF$chromosome) %% 2 == 0, "even", "odd")
        plotOut <- ggplot(tempPopStatDF, aes(x = (window_start + window_end)/2,
                                             y = avg_pi, color = chrom_color_group))
        # ggplot object
        plotOut <- plotOut + geom_point(size = 0.5, alpha = 0.5, stroke = 0)
        plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
        plotOut <- plotOut + xlab("Chromsome")
        plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
        plotOut <- plotOut + ylab("Pi (50000bp)")
        plotOut <- plotOut + scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
        plotOut <- plotOut + facet_grid(. ~ chromosome, axis.labels = "margins",
                                        scales = "free_y", switch = "x", space = "free_x")
        plotOut <- plotOut + theme_classic()
        plotOut <- plotOut + theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   panel.spacing = unit(0.1, "cm"),
                                   strip.background = element_blank(),
                                   strip.placement = "outside",
                                   legend.position ="none")
        return(plotOut)
}

manhattanPlotFst <- function(dataset)
{
        # even and odd number chromosomes were filled different colors
        chrom_color_group <- ifelse(as.numeric(dataset$chromosome) %% 2 == 0, "even", "odd")
        plotOut <- ggplot(dataset, aes(x = (window_start + window_end)/2,
                                             y = avg_wc_fst, color = chrom_color_group))
        # ggplot object
        plotOut <- plotOut + geom_point(size = 0.5, alpha = 0.5, stroke = 0)
        plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
        plotOut <- plotOut + xlab("Chromsome")
        plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
        plotOut <- plotOut + ylab("Fst (50000bp)")
        plotOut <- plotOut + scale_y_continuous(expand = c(0, 0), limits = c(0,1.0))
        plotOut <- plotOut + facet_grid(. ~ chromosome, axis.labels = "margins",
                                        scales = "free_y", switch = "x", space = "free_x")
        plotOut <- plotOut + theme_classic()
        plotOut <- plotOut + theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   panel.spacing = unit(0.1, "cm"),
                                   strip.background = element_blank(),
                                   strip.placement = "outside",
                                   legend.position ="none")
        return(plotOut)
}

manhattanPlotSexualFst <- function(dataset)
{
        # even and odd number chromosomes were filled different colors
        chrom_color_group <- ifelse(as.numeric(dataset$chromosome) %% 2 == 0, "even", "odd")
        plotOut <- ggplot(dataset, aes(x = (window_start + window_end)/2,
                                       y = avg_wc_fst, color = chrom_color_group))
        # ggplot object
        plotOut <- plotOut + geom_point(size = 0.5, alpha = 0.5, stroke = 0)
        plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
        plotOut <- plotOut + xlab("Chromsome")
        plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
        plotOut <- plotOut + ylab("Intersexual-Fst (50000bp)")
        plotOut <- plotOut + scale_y_continuous(expand = c(0, 0), limits = c(NA,0.1))
        plotOut <- plotOut + facet_grid(. ~ chromosome, axis.labels = "margins",
                                        scales = "free_y", switch = "x", space = "free_x")
        plotOut <- plotOut + theme_classic()
        plotOut <- plotOut + theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   panel.spacing = unit(0.1, "cm"),
                                   strip.background = element_blank(),
                                   strip.placement = "outside",
                                   legend.position ="none")
        return(plotOut)
}

manhattanPlotSexualFst2 <- function(dataset)
{
  # even and odd number chromosomes were filled different colors
  chrom_color_group <- ifelse(as.numeric(dataset$chromosome) %% 2 == 0, "even", "odd")
  plotOut <- ggplot(dataset, aes(x = (window_start + window_end)/2,
                                 y = avg_wc_fst, color = chrom_color_group))
  # ggplot object
  plotOut <- plotOut + geom_point(size = 0.5, alpha = 0.5, stroke = 0)
  plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
  plotOut <- plotOut + xlab("Chromsome")
  plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
  plotOut <- plotOut + ylab("Inter-sexual Fst (50000bp)")
  plotOut <- plotOut + scale_y_continuous(expand = c(0, 0), limits = c(-0.05,0.25))
  plotOut <- plotOut + facet_grid(. ~ chromosome, axis.labels = "margins",
                                  scales = "free_y", switch = "x", space = "free_x")
  plotOut <- plotOut + theme_classic()
  plotOut <- plotOut + theme(axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             panel.spacing = unit(0.1, "cm"),
                             strip.background = element_blank(),
                             strip.placement = "outside",
                             legend.position ="none")
  return(plotOut)
}


manhattanPlotDxy <- function(dataset)
{
        # even and odd number chromosomes were filled different colors
        chrom_color_group <- ifelse(as.numeric(dataset$chromosome) %% 2 == 0, "even", "odd")
        plotOut <- ggplot(dataset, aes(x = (window_start + window_end)/2,
                                             y = avg_dxy, color = chrom_color_group))
        # ggplot object
        plotOut <- plotOut + geom_point(size = 0.5, alpha = 0.5, stroke = 0)
        plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
        plotOut <- plotOut + xlab("Chromsome")
        plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
        plotOut <- plotOut + ylab("Dxy (50000bp)")
        plotOut <- plotOut + scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
        plotOut <- plotOut + facet_grid(. ~ chromosome,
                                        scales = "free_y", switch = "x", space = "free_x")
        plotOut <- plotOut + theme_classic()
        plotOut <- plotOut + theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   panel.spacing = unit(0.1, "cm"),
                                   strip.background = element_blank(),
                                   strip.placement = "outside",
                                   legend.position ="none")
        return(plotOut)
}

manhattanPlotfd <- function(dataset, upper = 2, lower = -2)
{
  # even and odd number chromosomes were filled different colors
  chrom_color_group <- ifelse(as.numeric(dataset$chromosome) %% 2 == 0, "even", "odd")
  plotOut <- ggplot(dataset, aes(x = (window_start + window_end)/2,
                                 y = fd, color = chrom_color_group))
  # ggplot object
  plotOut <- plotOut + geom_point(size = 0.75, alpha = 0.5, stroke = 0)
  #plotOut <- plotOut + geom_line()
  plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
  plotOut <- plotOut + xlab("Chromsome")
  plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
  plotOut <- plotOut + ylab("fd (50 SNVs)")
  plotOut <- plotOut + scale_y_continuous(expand = c(0, 0), limits = c(-1, 5))
  plotOut <- plotOut + facet_grid(. ~ chromosome,
                                  scales = "free_y", switch = "x", space = "free_x")
  plotOut <- plotOut + theme_classic()
  plotOut <- plotOut + theme(axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             panel.spacing = unit(0.1, "cm"),
                             strip.background = element_blank(),
                             strip.placement = "outside",
                             legend.position ="none")
  return(plotOut)
}

manhattanPlotfdM <- function(dataset)
{
  # even and odd number chromosomes were filled different colors
  chrom_color_group <- ifelse(as.numeric(dataset$chromosome) %% 2 == 0, "even", "odd")
  plotOut <- ggplot(dataset, aes(x = (window_start + window_end)/2,
                                 y = f_dM, color = chrom_color_group))
  # ggplot object
  plotOut <- plotOut + geom_point(size = 0.75, alpha = 0.5, stroke = 0)
  plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
  plotOut <- plotOut + xlab("Chromsome")
  plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
  plotOut <- plotOut + ylab("fdM (50 SNVs)")
  plotOut <- plotOut + scale_y_continuous(expand = c(0, 0), limits = c(-1,1))
  plotOut <- plotOut + facet_grid(. ~ chromosome,
                                  scales = "free_y", switch = "x", space = "free_x")
  plotOut <- plotOut + theme_classic()
  plotOut <- plotOut + theme(axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             panel.spacing = unit(0.1, "cm"),
                             strip.background = element_blank(),
                             strip.placement = "outside",
                             legend.position ="none")
  return(plotOut)
}

manhattanPlotiHS <- function(dataset)
{
  # iHS table modification
  chromsomeNum <- factor(as.numeric(factor(substr(dataset$locus, 1, 10))))
  dataset <- cbind(dataset, chromsomeNum)

  # Manhattan plot for absolute Z-transformed iHS
  dataset$ihs <- abs(dataset$z_trans_ihs) # Absolute Z transformed iHS
  chrom_color_group <- ifelse(as.numeric(dataset$chromsomeNum) %% 2 == 0, "even", "odd")
  plotOut <- ggplot(dataset, aes(x = physPOS,
                               y = ihs, color = chrom_color_group))
  # ggplot object
  plotOut <- plotOut + geom_point(size = 0.5, alpha = 0.5, stroke = 0)
  plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
  plotOut <- plotOut + xlab("Chromsome")
  plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
  plotOut <- plotOut + ylab("iHS")
  plotOut <- plotOut + scale_y_continuous(expand = c(0, 0))
  plotOut <- plotOut + facet_grid(. ~ chromsomeNum,
                                  scales = "free_y", switch = "x", space = "free_x")
  plotOut <- plotOut + theme_classic()
  plotOut <- plotOut + theme(axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             panel.spacing = unit(0.1, "cm"),
                             strip.background = element_blank(),
                             strip.placement = "outside",
                             legend.position ="none")

  return(plotOut)
}


linePlotRrate <- function(dataset)
{
  # even and odd number chromosomes were filled different colors
  chrom_color_group <- ifelse(as.numeric(dataset$chromosome) %% 2 == 0, "even", "odd")
  plotOut <- ggplot(dataset, aes(x = pos,
                                 y = rRate, color = chrom_color_group))
  # ggplot object
  plotOut <- plotOut + geom_line(size = 0.75, alpha = 0.5, stroke = 0)
  plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
  plotOut <- plotOut + xlab("Chromsome")
  plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
  plotOut <- plotOut + ylab("rRate (cM / Mbp)")
  plotOut <- plotOut + scale_y_continuous(expand = c(0, 0))
  plotOut <- plotOut + facet_grid(. ~ chromosome,
                                  scales = "free_y", switch = "x", space = "free_x")
  plotOut <- plotOut + theme_classic()
  plotOut <- plotOut + theme(axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             panel.spacing = unit(0.1, "cm"),
                             strip.background = element_blank(),
                             strip.placement = "outside",
                             legend.position ="none")
  return(plotOut)
}


manhattanPlotXPCLR <- function(dataset)
{
  # even and odd number chromosomes were filled different colors
  chrom_color_group <- ifelse(dataset$chromNum %% 2 == 0, "even", "odd")
  plotOut <- ggplot(dataset, aes(x = (start + stop)/2,
                                 y = xpclr_norm, color = chrom_color_group))
  # ggplot object
  plotOut <- plotOut + geom_point(size = 0.5, alpha = 0.5, stroke = 0)
  plotOut <- plotOut + scale_color_manual(values = c("grey50", "black"))
  plotOut <- plotOut + xlab("Chromsome")
  plotOut <- plotOut + scale_x_continuous(expand = c(0, 0))
  plotOut <- plotOut + ylab("XPCLR norm (50000bp)")
  plotOut <- plotOut + scale_y_continuous(expand = c(0, 0), limits = c(NA,NA))
  plotOut <- plotOut + facet_grid(. ~ chromNum,
                                  scales = "free_y", switch = "x", space = "free_x")
  plotOut <- plotOut + theme_classic()
  plotOut <- plotOut + theme(axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             panel.spacing = unit(0.1, "cm"),
                             strip.background = element_blank(),
                             strip.placement = "outside",
                             legend.position ="none")
  return(plotOut)
}
