#!/bin/bash
#$ -S /bin/bash

DIR=$1
FASTQ1=$2
FASTQ2=$3
DataPath=$4
SAMPLE=$5
SubPath=$6


fastp \
  -i ${DataPath}/${FASTQ1} -I ${DataPath}/${FASTQ2} \
  -o ${SubPath}/${SAMPLE}"_1.fastq.gz" -O ${SubPath}/${SAMPLE}"_2.fastq.gz" \
  -p -P 20 \
  -w 1 \
  --trim_poly_g \
  --trim_poly_x \
  --length_required 24 \
  -y --low_complexity_filter \
  --cut_front \
  --cut_tail \
  -h ${SubPath}/${SAMPLE}".html"