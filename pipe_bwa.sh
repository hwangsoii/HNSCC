#!/bin/bash
#$ -S /bin/bash

# Define Variables
PICARD=/opt/Yonsei/Picard/2.26.4/picard.jar
REFERENCE_GENOME=/data/resource/reference/human/UCSC/hg38/BWAIndex/genome.fa
# REFERENCE_GENOME=/data/resource/reference/human/UCSC/hg19/BWAIndex/genome.fa

DIR=$1
SubPath=$2
DataPath=$3
FASTQ1=$4
FASTQ2=$5
SAMPLE=$6
PROJECT=$7
PLATFORM=$8

#mkdir SubPath

cd SubPath

## convert FASTQ to uBAM and add read group information using FastqToSam
#
#java -Xmx8g -jar $PICARD FastqToSam \
#       FASTQ=$FASTQ1
#       FASTQ2=$FASTQ2
#       OUTPUT='unligned_reads_'${SAMPLE}'.bam'
#       READ_GROUP_NAME=YON69
#       SAMPLE_NAME=$SAMPLE
#       LIBRARY_NAME=
#
#export GATK PICARD JAVA REFERENCE_GENOME


## index reference genome

#bwa index -a bwtsw /data/resource/reference/human/NCBI/GRCh38/WholeGenomeFasta/genome.fa

##/data/resource/reference/human/NCBI/GRCh38_GATK/BWAIndex/genome.fa

## make sam from two fastq file

bwa mem -t 1 ${REFERENCE_GENOME} ${DataPath}/${FASTQ1} ${DataPath}/${FASTQ2} | java -jar -Xmx16g ${PICARD} SortSam \
I=/dev/stdin \
O=${SubPath}/${SAMPLE}'_sorted.bam' \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true


## add read group
java -Xmx16g -jar ${PICARD} AddOrReplaceReadGroups \
I=${SubPath}/${SAMPLE}_sorted.bam  \
O=${SubPath}/${SAMPLE}'_sort_rg.bam' \
RGLB=$PROJECT \
RGPL=$PLATFORM \
RGPU=$PLATFORM \
RGSM=${SAMPLE} \
SORT_ORDER=coordinate \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT
#TMP_DIR=${DIR}/'tempfile' ## temporary folder for large file

## MArkduplicate

java -jar -Xmx8g $PICARD MarkDuplicates \
I=${SubPath}/${SAMPLE}_sort_rg.bam \
O=${SubPath}/${SAMPLE}'_sort_rg_marked_duplicates.bam' \
M=${SubPath}/${SAMPLE}'_marked_dup_metrics.txt' \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT

# 1st base recalobratopm

gatk --java-options "-Xmx16g" BaseRecalibrator \
-I ${SubPath}/${SAMPLE}_sort_rg_marked_duplicates.bam \
-R ${REFERENCE_GENOME} \
--known-sites /data/public/GATK/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--known-sites /data/public/dbSNP/b155/GRCh38/GCF_000001405.39.re.common.vcf.gz \
--known-sites /data/public/GATK/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-O ${SubPath}/${SAMPLE}'_recal_data1.table'

# BQSR

gatk --java-options "-Xmx16g" ApplyBQSR \
-R ${REFERENCE_GENOME} \
-I ${SubPath}/${SAMPLE}_sort_rg_marked_duplicates.bam \
--bqsr-recal-file ${SubPath}/${SAMPLE}_recal_data1.table \
-O ${SubPath}/${SAMPLE}'_recal1.bam'

## 2nd pass
gatk --java-options "-Xmx16g" BaseRecalibrator \
-I ${SubPath}/${SAMPLE}_recal1.bam \
-R ${REFERENCE_GENOME} \
--known-sites /data/public/GATK/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--known-sites /data/public/dbSNP/b155/GRCh38/GCF_000001405.39.re.common.vcf.gz \
--known-sites /data/public/GATK/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-O ${SubPath}/${SAMPLE}'_recal_data2.table'



gatk --java-options "-Xmx16g" ApplyBQSR \
-R ${REFERENCE_GENOME} \
-I ${SubPath}/${SAMPLE}_recal1.bam \
--bqsr-recal-file ${SubPath}/${SAMPLE}_recal_data2.table \
-O ${SubPath}/${SAMPLE}'_recal2.bam'

## compare

gatk --java-options "-Xmx16g" AnalyzeCovariates \
-before ${SubPath}/${SAMPLE}_recal_data1.table \
-after ${SubPath}/${SAMPLE}_recal_data2.table \
-plots ${SubPath}/${SAMPLE}_AnalyzeCovariates.pdf




# --known-sites /data/public/GATK/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
# --known-sites /data/public/dbSNP/b155/GRCh37/GCF_000001405.25.re.common.vcf.gz \
# --known-sites /data/public/GATK/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \