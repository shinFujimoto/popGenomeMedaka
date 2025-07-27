#!/bin/sh
#----------------------------------------------------------------
# populationStructure_supplymentaryAnalysis.sh
# Author：Shingo Fujimoto
# Date：2024/06/14
#----------------------------------------------------------------
# set the working directory
# These filteration conducted
# ----------------------------------------------------------------------------------------------
# 1-1, calculate p-distance matrix and nj tree, N = 45,
# with no missing and include invariant site dataset
# excluding individuals with multiple ancestry in ADMIXTURE, K=4(ykk, nck, oit, nss, tng)
# ----------------------------------------------------------------------------------------------
workDIR="/home/shifujimoto/medakaWGStmp/populationStat"
outputDIR="/home/shifujimoto/medakaWGStmp/populationStat/njtree45sample_nomissing_invariant"
mergeSampleList="sampleList_wgs_HSOK_sak_lat_45samples_njtree_invariant.txt"

cd $workDIR/vcf_populationStat

# convert vcf to p-distance matrix
cd $outputDIR
/home/shifujimoto/bin/VCF2Dis/bin/VCF2Dis \
  -InPut medakaWGS_HSOK_sak_lat_45samples_snpsQC_invariant.vcf.gz \
  -SubPop exclude_admixIndividual_nmg_nck_oit_nss_tng/sampleName.txt \
  -OutPut exclude_admixIndividual_nmg_nck_oit_nss_tng/medakaWGS_HSOK_sak_lat_38samples_snpsQC_invariant.dist
sed -i -e '1d' medakaWGS_HSOK_sak_lat_38samples_snpsQC_invariant.dist

# convert p-distance matrix to nexus format in local machine
cd /home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/njtree45sample_nomissing_invariant
Rscript /mnt/d/GoogleDrive/lib/pdistMatrixToNexus.R medakaWGS_HSOK_sak_lat_40samples_snpsQC_invariant.dist


# ----------------------------------------------------------------------------------------------
# vcftools
# Date: 2025/07/09
# ----------------------------------------------------------------------------------------------
cd /home/shifujimoto/medakaWGStmp/populationStat/genomescan/FST_popGenome_bp_resolution

zcat ../medakaWGS_sak40samples_snpsQC_exIndelMultialleic_snps.vcf.gz | \
sed 's/##FORMAT=<ID=AD,Number=,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">/##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">/' | bgzip -@ 24 -c > ./medakaWGS_sak40samples_snpsQC_exIndelMultialleic_snps.vcf.gz

vcftools --gzvcf ./medakaWGS_sak40samples_snpsQC_exIndelMultialleic_snps.vcf.gz \
         --weir-fst-pop pop1.txt \
         --weir-fst-pop pop2.txt \
         --out intersexualFST_sak40

