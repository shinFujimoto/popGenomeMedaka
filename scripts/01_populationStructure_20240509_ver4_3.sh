#!/bin/sh
#----------------------------------------------------------------
# populationStructure.sh
# Author：Shingo Fujimoto
# Date：2022/01/20
# update:2023/03/22, exclude not main results analyses
#----------------------------------------------------------------
# set the working directory
# These filteration conducted
# ---------------------------------------------------------------------------------------------------
# Exclude low depth and high depth positions
# ---------------------------------------------------------------------------------------------------
workDIR="/home/shifujimoto/medakaWGStmp/populationStat"

# Low quality positions in each VCF files were excluded
cd /home/shifujimoto/medakaWGStmp/vcf
ls -1 | grep "vcf.gz$" | sed "s/.vcf.gz//g" | xargs -n1 -P48 -I{} sh -c \
 "bcftools view {}.vcf.gz -i 'MIN(FORMAT/DP)>=10 & MAX(FORMAT/DP)<=200' | \
  bcftools view -i 'INFO/MQ>=40 & INFO/QD>=2.0 || GT=\"ref\"' | \
  sed -e \"5s/Number=R/Number=/g\" | bgzip -c > ../populationStat/vcf_qcfilter/{}_filtered.vcf.gz"

cd $workDIR/vcf_qcfilter
ls -1 | grep "filtered.vcf.gz$" | xargs -n1 -P48 -I{} sh -c "tabix -f ./{}"

# ----------------------------------------------------------------------------------------------
# 1-1, calculate p-distance matrix and nj tree, N = 45,
# with no missing and include invariant site dataset
# ----------------------------------------------------------------------------------------------
outputDIR="/home/shifujimoto/medakaWGStmp/populationStat/njtree45sample_nomissing_invariant"
mergeSampleList="/home/shifujimoto/medakaWGStmp/populationStat/njtree45sample_nomissing_invariant/sampleList_wgs_HSOK_sak_lat_45samples_njtree_invariant.txt"

cd $workDIR/vcf_qcfilter
# missing rate = 0%, include multi allelic and invariant site
cat chromosomeNameList.txt | xargs -n1 -P24 -I{} sh -c \
"bcftools merge -l $mergeSampleList -r {} | bcftools view -i 'N_MISSING=0'| bcftools view -e 'ALT=\"*\"' | bgzip -@ 4 -c > $outputDIR/tmp_{}_merged.vcf.gz"

cd $outputDIR
ls -v | grep "_merged.vcf.gz$" | grep "tmp" > tmpFileLsit.txt
bcftools concat -f tmpFileLsit.txt | bgzip -@ 12 -c > medakaWGS_HSOK_sak_lat_45samples_snpsQC_invariant.vcf.gz
tabix medakaWGS_HSOK_sak_lat_45samples_snpsQC_invariant.vcf.gz -f

# convert vcf to p-distance matrix
/home/shifujimoto/bin/VCF2Dis/bin/VCF2Dis \
  -InPut medakaWGS_HSOK_sak_lat_45samples_snpsQC_invariant.vcf.gz \
  -OutPut medakaWGS_HSOK_sak_lat_45samples_snpsQC_invariant.dist
sed -i -e '1d' medakaWGS_HSOK_sak_lat_45samples_snpsQC_invariant.dist

# convert p-distance matrix to nexus format in local machine
cd /home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/njtree45sample_nomissing_invariant
Rscript /mnt/d/GoogleDrive/lib/pdistMatrixToNexus.R medakaWGS_HSOK_sak_lat_45samples_snpsQC_invariant.dist

# ----------------------------------------------------------------------------------------------
# 1-2. ADMIXTURE dataset, N = 45 -> 44 (exclude HSOK),
# O. sak(Tsuruga, Niigata, Aomori, Kitatsugaru, N=2), O. latipes(all pop, N=1or2)
# ----------------------------------------------------------------------------------------------
cd /home/shifujimoto/medakaWGStmp/populationStat/ADMIXTURE
inputVCF="medakaWGS_HSOK_sak_lat_45samples_snpsQC_invariant.vcf.gz"

# Exclude HSOK strain, and extract the only biallelic sites 
bcftools view ../njtree45sample_nomissing_invariant/$inputVCF -s "^DRR002221_HSOK" \
 | bcftools view -m2 -M2 -v snps -c 1:minor \
 | bcftools annotate --rename-chrs chromosomeRename.txt | bgzip -c -@ 6 > medakaWGS_HSOK_sak_lat_44samples_snps_qcFilter_rename_exHSOK.vcf.gz
tabix medakaWGS_HSOK_sak_lat_44samples_snps_qcFilter_rename_exHSOK.vcf.gz

# convert PED file and perform ADMIXTURE
# considering LD and pruning and exclude HSOK
#  --set-all-var-ids: @_#_\$r_\$a, make the SNP_ID following, chromosome_positon_ref_alt
# ----------------------------------------------------------------------------------------------
inputVCF="medakaWGS_HSOK_sak_lat_44samples_snps_qcFilter_rename_exHSOK.vcf.gz"
inputPrefix=`basename $inputVCF`

plink2 --vcf $inputVCF --set-all-var-ids @_#_\$r_\$a --make-bed --out ./${inputPrefix%.vcf.gz} --allow-extra-chr --autosome-num 24
plink2 --bfile ./${inputPrefix%.vcf.gz} --indep-pairwise 50 10 0.1 --bad-ld --out ./${inputPrefix%.vcf.gz} --allow-extra-chr --autosome-num 24
plink2 --bfile ./${inputPrefix%.vcf.gz} --extract ./${inputPrefix%.vcf.gz}.prune.in --out ./${inputPrefix%.vcf.gz}.pruned --make-bed --allow-extra-chr --autosome-num 24

for K in 1 2 3 4 5 6 7 8 9 10; \
do
 admixture --cv  ./${inputPrefix%.vcf.gz}.pruned.bed $K | tee log${K}.out &
done
grep -h CV log*.out > CVerror.txt

# plot the admixture proportion in local machine
# cd /home/shifujimoto/processedSequence/medakaWGS/datasets/refHNI/output/populationStructure_ver3/ADMIXTURE
# Rscript ./barplot_admixture.R

# ----------------------------------------------------------------------------------------------
# 2-1. Outgroup F3, Dstat, and admixturegraph dataset, N = 45 -> select 33
# Population separated following NJtree and ADMIXTURE
# HSOK(N=1)
# O. sak(Tsuruga(2), Niigata(2), Aomori(2), Kitatsugaru(2), N=8),
# O. latipes, westKyu (Ginama(2), Izumi(2), Nagata(2), Hrd(2),  N=8)
# O. latipes, westJ (Myk(2), Ehm(2), Eky(2), Miy(1), tsm(1),  N=8)
# O. latipes, eastJ (Iso(2), Csg(2), Mry(2), Ich(2), N=8)
# ----------------------------------------------------------------------------------------------
cd /home/shifujimoto/medakaWGStmp/populationStat/admixtools

# excluded possible admixed individuals and adjust the sample number
inputVCF="medakaWGS_HSOK_sak_lat_45samples_snpsQC_invariant.vcf.gz"
bcftools view ../njtree45sample_nomissing_invariant/$inputVCF -S sampleName.txt | bcftools view -m2 -M2 -v snps -c 1:minor \
 | bcftools annotate --rename-chrs chromosomeRename.txt | bgzip -@ 12 -c > medakaWGS_HSOK_sak_lat_33samples_snps_qcFilter_rename.vcf.gz
tabix medakaWGS_HSOK_sak_lat_33samples_snps_qcFilter_rename.vcf.gz

# convert from vcf to plink BED format
plink2 --vcf medakaWGS_HSOK_sak_lat_33samples_snps_qcFilter_rename.vcf.gz --set-all-var-ids @_#_\$r_\$a --make-bed --out medakaWGS_HSOK_sak_lat_33samples --allow-extra-chr --autosome-num 24

# create an inputfile for convertf in admixtools, ver.7.0.2
# assignment of populations defined in ".ind" file
inputVCF="medakaWGS_HSOK_sak_lat_33samples"
echo "genotypename:    ${inputVCF}.bed" > par.PED.EIGENSTRAT
echo "snpname:         ${inputVCF}.bim" >> par.PED.EIGENSTRAT
echo "indivname:       ${inputVCF}.fam" >> par.PED.EIGENSTRAT
echo "outputformat:    EIGENSTRAT" >> par.PED.EIGENSTRAT
echo "genotypeoutname: ${inputVCF}.geno" >> par.PED.EIGENSTRAT
echo "snpoutname:      ${inputVCF}.snp" >> par.PED.EIGENSTRAT
echo "indivoutname:    ${inputVCF}.ind" >> par.PED.EIGENSTRAT
echo "familynames:     NO" >> par.PED.EIGENSTRAT

# convert BED to eigenstratgenotype format
convertf -p par.PED.EIGENSTRAT

# ---------------------------------------------------------------------------------------
# 3-1. f_d
# genome scan of f_d calculated by Dsuite
# 2024/07/19, 再計算
# ---------------------------------------------------------------------------------------
workDIR="/home/shifujimoto/medakaWGStmp/populationStat/genomescan/Dsuite_f_d"
inputVCF="/home/shifujimoto/medakaWGStmp/populationStat/njtree45sample_nomissing_invariant/medakaWGS_HSOK_sak_lat_45samples_snpsQC_invariant.vcf.gz"

cd $workDIR

# Extract only biallelic SNPs
cp /home/shifujimoto/medakaWGStmp/populationStat/admixtools/medakaWGS_HSOK_sak_lat_33samples_snps_qcFilter_rename.vcf.gz ./
gunzip medakaWGS_HSOK_sak_lat_33samples_snps_qcFilter_rename.vcf.gz

# Dinvestigate: calculation, Patterson'sD, f_dM and f_d statistics
Dsuite Dinvestigate ./medakaWGS_HSOK_sak_lat_33samples_snps_qcFilter_rename.vcf \
  medakaWGS_HSOK_sak_lat_33samples_ver2.ind \
  -w 50,50 \
  ./test_trios_eastJ_wKyu_sak.txt

