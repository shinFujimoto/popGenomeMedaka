# --------------------------------------------------
# GeneName_Exon_region_extraction.txt
# Date: 2023/11/20
#
# Update: 2024/06/28, HNI reference genome
# Update: 2025/01/24, HIGH TajimaD region extract
# Author: Shingo Fujimoto
# --------------------------------------------------
# set the working directory
workDIR="/mnt/d/GoogleDrive/study/論文原稿/0_medakaPopGenomics/DatasetAndScript/genomescan/formatting_output/gene_Info_ensembl"

# gtf file was retrieved from emsembl 109,
gtfFile="Oryzias_latipes_hni.ASM223471v1.109.gtf"

# extract genes based on the target gene list
# ---------------------------------------------------------------
cd $workDIR

# format conversion, gtf -> bed
awk -F'\t' '$3 == "exon" {print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}' $gtfFile | sort -k1,1 -k2,2n > Oryzias_latipes_hni.ASM223471v1.109.bed

# make the unique gene_id and gene_name list in all chromosome
gawk '{
  if (match($0, /gene_name "([^"]+)"/, arr)) {
   if (match($0, /gene_id "([^"]+)".+gene_name "([^"]+)"/, arr)) {
    print arr[1] "\t" arr[2];
   }
  } else {
   if (match($0, /gene_id "([^"]+)"/, arr)) {
    print arr[1] "\t" "-1";
    }
  }
}' $gtfFile | sort -u  > Oryzias_latipes_hni.ASM223471v1.109.geneName.txt &

# Extract the candidate regions from BED format file
# --------------------------------------------------------
# Output gene list text file
outputDIR="sak40_lat42_iHS_50000bp"
inputBED="genomescan_windows_sak40_lat42_iHS_50000bp.bed"
touch $outputDIR/targetRegion_geneList_Oryzias_latipes_hni.ASM223471v1.109.txt
mkdir $outputDIR/tmpDIR # confirm the empty in this directory

cat $outputDIR/$inputBED | while read -a line
do
 # chromosome names, start and end positions
 chrName=${line[0]} # chromosome name
 startPos=${line[1]} # start position
 endPos=${line[2]} # end position
 region_annotation=${line[3]} # Annotation information in the genetic region
 
 # extract target chromosome region
 bedtools intersect -a Oryzias_latipes_hni.ASM223471v1.109.bed -b <(echo -e "$chrName\t$startPos\t$endPos") > $outputDIR/tmpDIR/${chrName}_${startPos}_${endPos}_${region_annotation}.bed
 
 # make the unique gene_id and gene_name list
 gawk '{
   if (match($0, /gene_name "([^"]+)"/, arr)) {
    if (match($0, /gene_id "([^"]+)".+gene_name "([^"]+)"/, arr)) {
     print arr[1] "\t" arr[2];
    }
   } else {
    if (match($0, /gene_id "([^"]+)"/, arr)) {
     print arr[1] "\t" "-1";
     }
   }
 }' $outputDIR/tmpDIR/${chrName}_${startPos}_${endPos}_${region_annotation}.bed | sort -u  > $outputDIR/tmpDIR/genes_${chrName}_${startPos}_${endPos}_${region_annotation}.txt
 
  # extract start and end position in each gene
 cat $outputDIR/tmpDIR/genes_${chrName}_${startPos}_${endPos}_${region_annotation}.txt | while read -a gene
 do
  gene_id=${gene[0]} # gene id
  gene_name=${gene[1]} # gene name
  
  cat $outputDIR/tmpDIR/${chrName}_${startPos}_${endPos}_${region_annotation}.bed | grep $gene_id > tmp_gene.bed
  chrNum=`awk '{ print $1 }' tmp_gene.bed | head -n 1`
  gene_startPos=`awk '{ print $2+1 }' tmp_gene.bed | head -n 1` # start position
  gene_endPos=`awk '{ print $3+1 }' tmp_gene.bed | tail -n 1` # end position
  echo $chrNum $gene_startPos $gene_endPos $gene_id $gene_name $region_annotation >> $outputDIR/targetRegion_geneList_Oryzias_latipes_hni.ASM223471v1.109.txt
 done
done

# --------------------------------------------------------

