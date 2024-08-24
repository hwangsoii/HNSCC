#!/bin/bash
#$ -S /bin/bash


# Define Variables
PICARD=/opt/Yonsei/picard/2.22.0/picard.jar
REF=/data/project/TRIUMPH/reference/hg38WGSFASTA/genome.fa #TRIUMPH
# REF=/data/project/TRIUMPH/tcga/download/reference/GRCh38.d1.vd1.fa #TCGA
# REFERENCE_GENOME=/data/resource/reference/human/UCSC/hg19/BWAIndex/genome.fa

DIR=$1
SubPath=$2
DataPath=$3
VCF=$4
SAMPLE=$5
# TUMOR=$6

filepath=/opt/Yonsei/vcf2maf/1.6.16/vcf2maf.pl

# perl ${filepath} \
vcf2maf.pl \
--input-vcf ${DataPath}/${VCF} \
--output-maf $SubPath/${SAMPLE}.vep.maf \
--vep-path /opt/Yonsei/ensembl-vep/104.3 \
--vep-data /data/public/VEP/104 \
--filter-vcf 0  \
--cache-version 104 \
--ref-fasta $REF \
--ncbi-build GRCh38 \
--tumor-id ${SAMPLE}'-FP' \
--normal-id ${SAMPLE}'-BL' \
--retain-info \
--retain-fmt 

# --custom-enst custome_file
###
##
#cutomfile

#ENSS1222
#ENSS11222
#like this
# --vep-data /data/public/VEP/104 \


#################newvcf2maf
# filepath=/opt/Yonsei/vcf2maf/1.6.21/vcf2maf.pl

# # perl ${filepath} \
# vcf2maf.pl \
# --input-vcf ${DataPath}/${VCF} \
# --output-maf $SubPath/${SAMPLE}.vep.maf \
# --vep-path /opt/Yonsei/ensembl-vep/104.3 \
# --vep-data /data/public/VEP/104 \
# --vep-custom dbNSFP,/home/goldpm1/tools/ANNOVAR/humandb/dbnsfp4.2a/dbNSFP4.0.zip,ALL  \
# --cache-version 104 \
# --ref-fasta $REF \
# --ncbi-build GRCh38 \
# --tumor-id ${SAMPLE}'-FP' \
# --normal-id ${SAMPLE}'-BL' \
# --retain-info \
# --retain-fmt
#