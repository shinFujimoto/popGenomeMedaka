#!/bin/sh
#----------------------------------------------------------------
# 01_fastqPreprocess.sh
# 
# Author：Shingo Fujimoto
# Date：2018/03/05
# update：2019/09/24
# update:2020/01/20, wildsamples 12 individuals
#----------------------------------------------------------------
# set the working directory
# Working directory
rawdata="/home/shifujimoto/pipelines/medakaWGS/datasets/refHdrR/fastq" # 生のfastqの置き場所
output="/home/shifujimoto/pipelines/medakaWGS/datasets/refHdrR/output" # 加工後ファイルの置き場
log="/home/shifujimoto/pipelines/medakaWGS/datasets/refHdrR/log" # logの置き場所
tmp="/home/shifujimoto/pipelines/medakaWGS/datasets/refHdrR/tmp" # 加工中のファイルの一時保管場所

# reference path
reference='/home/shifujimoto/pipelines/medakaWGS/datasets/refHdrR/reference/GCA_002234675.1_ASM223467v1_genomic.fna'
ref='/home/shifujimoto/pipelines/medakaWGS/datasets/refHdrR/reference/oryzias_index'
refDIR='/home/shifujimoto/pipelines/medakaWGS/datasets/refHdrR/reference'
gatkRef='/home/shifujimoto/pipelines/medakaWGS/datasets/refHdrR/reference/oryziasLatipes.dict'

# Software path
fastp='/home/shifujimoto/pipelines/library/fastp'
bwa="/usr/local/bin/bwa-0.7.16a/bwa"
picard='/usr/local/bin/picard-2.10.10'

# SNP call, filtering
gatk='/usr/local/bin/GenomeAnalysisTK-3.8-0-ge9d806836'
#----------------------------------------------------------------
cd $tmp

# sample name list 
ls -1 $rawdata | grep -v "_2" | sed -e 's/_1.fastq.gz//' > 10sampleName.txt

cat 10sampleName_20220603.txt | while read line; do
 sampleName=${line%_1.fastq.gz}
 $fastp -i $rawdata/${sampleName}_1.fastq.gz -I $rawdata/${sampleName}_2.fastq.gz \
  --detect_adapter_for_pe --cut_front \
  -h ${sampleName}.html -w 12 \
  -o ${sampleName}_R1_fastp.fastq.gz -O ${sampleName}_R2_fastp.fastq.gz

  # 2. Mapping on the reference genome
  bwa mem $ref -t 24 ./${sampleName}_R1_fastp.fastq.gz ./${sampleName}_R2_fastp.fastq.gz > ./${sampleName}_paired.sam
  samtools view -@ 6 -bS ./${sampleName}_paired.sam | samtools sort -@ 6 > ${sampleName}_sorted.bam
  
  java -jar $picard/picard.jar CleanSam INPUT=${sampleName}_sorted.bam OUTPUT=${sampleName}_clean.bam > out.log 2> $log/${sampleName}_picard_out.log
  java -jar $picard/picard.jar AddOrReplaceReadGroups I=${sampleName}_clean.bam \
	O=${sampleName}_sortrg.bam \
	SORT_ORDER=coordinate \
	RGID=${sampleName} \
	RGLB=medakaWGS \
	RGPL=ILLUMINA \
	RGPU=JN00006698 \
	RGSM=${sampleName} \
	CREATE_INDEX=True > out.log 2> $log/${sampleName}_picard_rg.out.log

# SNPのコールにはGATKを使う (Samuk et al. 2017)
# BWAは計算が高速だがindel周辺でマッピングエラーを生じやすい、その補正をgatkで行う
java -Xmx16g -jar $gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator \
 -nt 12 \
 -R $reference \
 -I ${sampleName}_sortrg.bam \
 -o ${sampleName}.intervals \
 -log $log/${sampleName}_gatkRTC.log

java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar -T IndelRealigner -R $reference \
 -I ${sampleName}_sortrg.bam \
 -targetIntervals ${sampleName}.intervals \
 -o ${sampleName}_realign.bam \
 -log $log/${sampleName}_gatkIndelRealigner.log

  # 3. genotype calling
 # g.vcf形式で書き出し
java -Xmx8g -jar $gatk/GenomeAnalysisTK.jar -T HaplotypeCaller \
  -R $reference \
  -nct 24 \
  -I ${sampleName}_realign.bam \
  -o ${sampleName}_realign.g.vcf \
  -log $log/{sampleName}_gatkHC.log \
  --emitRefConfidence GVCF \
  --variant_index_type LINEAR \
  --variant_index_parameter 128000

# コールしたgvcfからvcfを作成
java -Xmx18g -jar $gatk/GenomeAnalysisTK.jar -T GenotypeGVCFs \
  -R $reference \
  -nt 12 \
  --variant ${sampleName}_realign.g.vcf \
  -o ${sampleName}_rawvariants.vcf
  
done

  rm ${sampleName}_R1_fastp.fastq.gz
  rm ${sampleName}_R2_fastp.fastq.gz
  rm ${sampleName}_paired.sam
 
# -------------------------------------------------------------------------------------------------------
# 3. Base quality score recaliblation (BQSR)
# リファレンスゲノムと異なっていたSNP, indel の箇所を抽出後、ゲノム全体の base quality scoreをキャリブレーション
# 複数レーンに由来するシーケンスデータなどをコンカチして解析するときなどに必要
# -------------------------------------------------------------------------------------------------------
# SNPの抽出
java -jar $gatk/GenomeAnalysisTK.jar -T SelectVariants -R $reference \
  -nt 4 \
  -V ${sampleName}_rawvariants.vcf \
  -selectType SNP \
  -o ${sampleName}_snps.vcf

# SNPのフィルター作成
java -jar $gatk/GenomeAnalysisTK.jar -T VariantFiltration \
   -R $reference \
   -V ${sampleName}_snps.vcf \
   --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" \
   --filterName "basic_snp_filters" \
   -o ${sampleName}_snps_filltered.vcf

java -jar $gatk/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $reference \
   -V ${sampleName}_snps_filltered.vcf \
   -o ${sampleName}_output.vcf \
   -se 'SAMPLE.+basic_snp_filters'

# Indel の抽出
java -jar $gatk/GenomeAnalysisTK.jar -T SelectVariants \
  -nt 4 \
  -R $reference \
  -V ${sampleName}_rawvariants.vcf \
  -selectType INDEL \
  -o ${sampleName}_indels.vcf

# Indelのフィルター作成
java -jar $gatk/GenomeAnalysisTK.jar -T VariantFiltration \
   -R $reference \
   -V ${sampleName}_indels.vcf \
   --filterExpression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0" \
   --filterName "basic_indel_filters" \
   -o ${sampleName}_indels_filltered.vcf

# gatk, Base Quality Score Recalibration (BQSR) #1,
 # BQSRのテーブル作成
 java -jar $gatk/GenomeAnalysisTK.jar -T BaseRecalibrator \
  -R $reference \
  -I ${sampleName}_realign.bam \
  -knownSites ${sampleName}_snps_filltered.vcf \
  -knownSites ${sampleName}_indels_filltered.vcf \
  -o ${sampleName}_recal_data.table

# -------------------------------------------------------------------------------------------------------
# BQSRで較正後のBAMを作成
# -------------------------------------------------------------------------------------------------------
 java -jar $gatk/GenomeAnalysisTK.jar -T PrintReads \
 -R $reference \
 -I ${sampleName}_realign.bam \
 -BQSR ${sampleName}_recal_data.table \
 -o ${sampleName}_recal.bam

# 較正後のBAMを使ったg.vcfファイルを再計算
java -Xmx16g -jar $gatk/GenomeAnalysisTK.jar -T HaplotypeCaller \
   -nct 12 \
   -R $reference \
   -I ${sampleName}_recal.bam \
   -o ${sampleName}_recal.g.vcf \
   --emitRefConfidence GVCF \
   --variant_index_type LINEAR \
   --variant_index_parameter 128000 >out.log 2> $log/${sampleName}_gatkRecal.log

# includeNonVariantSitesは、リファレンスゲノムと同じ塩基も含めてすべて書き出す
# 設定し忘れると変異がなかったサイトの情報が欠落するので注意
java -Xmx18g -jar $gatk/GenomeAnalysisTK.jar -T GenotypeGVCFs \
  -nt 12 \
  -R $reference \
  --variant ${sampleName}_recal.g.vcf \
  -o ${sampleName}_recal.vcf \
  --max_alternate_alleles 2 \
  --includeNonVariantSites

