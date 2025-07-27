# ----------------------------------------------------------------
# geneOrthologAnnotation_HdrR_HNI
# 
# Author: Fujimoto Shingo
# Date: 2024/10/16
# update: 2025/01/24, sak40_pi_tajimaD_sexFst_Fd
# ----------------------------------------------------------------
library("tidyverse")

setwd("D:/GoogleDrive/study/˜_•¶Œ´e/0_medakaPopGenomics/DatasetAndScript/genomescan/formatting_output/gene_Info_ensembl")
setwd("C:/users/fujim/GoogleDrive/study/˜_•¶Œ´e/0_medakaPopGenomics/DatasetAndScript/genomescan/formatting_output/gene_Info_ensembl")

# import the gene annotation table modified by Orthofinder one-to-one ortholog
hni_hdrR_geneDF <- read.csv("sak40_pi_tajimaD_sexFst_F_d/Table_S2_HNI_HdrR_gene_annotation_ortholog.csv")

# Delete the strand information in Ensembl gene ID of HNI reference
hni_hdrR_geneDF$ensemblGeneID_HNI <- str_sub(hni_hdrR_geneDF$ensemblGeneID_HNI,
                                             start = 1,
                                             end = (nchar(hni_hdrR_geneDF$ensemblGeneID_HNI[1]) - 2))
# Delete the strand information in Ensembl gene ID of HdrR reference
hni_hdrR_geneDF$ensemblGeneID_HdrR <- str_sub(hni_hdrR_geneDF$ensemblGeneID_HdrR,
                                              start = 1,
                                              end = (nchar(hni_hdrR_geneDF$ensemblGeneID_HdrR[1]) - 2))
# Modify the gene position information in HdrR reference
gene_start_HdrR <- str_split_i(hni_hdrR_geneDF$assembly_HdrR, pattern = ":", i = 2)
gene_end_HdrR <- str_split_i(hni_hdrR_geneDF$assembly_HdrR, pattern = ":", i = 3)
hni_hdrR_geneDF <- cbind(hni_hdrR_geneDF, gene_start_HdrR, gene_end_HdrR)

# Extract columns in HNI and HdrR reference genome information
hni_hdrR_geneDF <- dplyr::select(hni_hdrR_geneDF,
                          chrom_HNI,
                          gene_start_HNI,
                          gene_end_HNI,
                          gene_start_HdrR,
                          gene_end_HdrR,
                          ensemblGeneID_HNI,
                          ensemblGeneID_HdrR,
                          gene_symbol_HNI,
                          gene_symbol_HdrR,
                          description_HNI,
                          description_HdrR) %>% distinct(ensemblGeneID_HNI, ensemblGeneID_HdrR, .keep_all = TRUE)


# import the gene ID list of genome scan in HNI reference genome
genomeScanGenes <- read.table("sak40_lat42_iHS_50000bp/targetRegion_geneList_Oryzias_latipes_hni.ASM223471v1.109.txt", header = FALSE)
# Rename the column name
names(genomeScanGenes) <- c("chrom_HNI", "gene_start_HNI", "gene_end_HNI",
                            "ensemblGeneID_HNI", "gene_symbol_HNI", "popStat_phenotype_annotation")
# Exclude the duplicated Ensembl gene ID
genomeScanGenes <- distinct(genomeScanGenes, ensemblGeneID_HNI, .keep_all = TRUE)

# Merge the reference table and genome scan table
mergedGeneTable <- left_join(genomeScanGenes, hni_hdrR_geneDF,
                             by = join_by(chrom_HNI == chrom_HNI,
                                          ensemblGeneID_HNI == ensemblGeneID_HNI)
                             )
# Merge the gene symbol information 
gene_symbol_hni <- ifelse(mergedGeneTable$gene_symbol_HNI.x == "-1",
                          mergedGeneTable$gene_symbol_HNI.y,
                          mergedGeneTable$gene_symbol_HNI.x
                          )
mergedGeneTable <- cbind(mergedGeneTable, gene_symbol_hni)

# Extract the output column
mergedGeneTable <- dplyr::select(mergedGeneTable,
                                 chrom_HNI, gene_start_HNI.x, gene_end_HNI.x,
                                 ensemblGeneID_HNI, gene_symbol_hni, description_HNI,
                                 ensemblGeneID_HdrR, gene_start_HdrR, gene_end_HdrR, 
                                 gene_symbol_HdrR, description_HdrR,
                                 popStat_phenotype_annotation
)
names(mergedGeneTable)[c(2, 3)] <- c("gene_start_HNI", "gene_end_HNI")

# Output the csv file
write.csv(mergedGeneTable, file = "sak40_lat42_iHS_50000bp/genomeScanGenes_geneAnnotation.csv", row.names = FALSE)

