#!/bin/bash
#$ -cwd
#$ -S /bin/bash


VCF=$1
SAMPLE=$2
SubPath=$3
REF=/data/resource/reference/human/UCSC/hg19/BWAIndex/genome.fa
CHAIN=/data/project/TRIUMPH/download/liftover/hg38ToHg19.over.chain.gz

CrossMap.py  vcf  ${CHAIN}  ${VCF}  ${REF}  ${SubPath}/${SAMPLE}.hg19.vcf