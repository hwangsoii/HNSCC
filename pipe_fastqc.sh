#!/bin/bash
#$ -cwd
#$ -S /bin/bash

# basic argument
DIR=$1
FASTQ1=$2
FASTQ2=$3
DataPath=$4

if [ ! -d $DIR/2.qc/fastqc ]
then
    	mkdir -p $DIR/2.qc/fastqc
fi
qcPath=$DIR/2.qc/fastqc

# print start time
date

fastqc -o $qcPath $DataPath/$FASTQ1
fastqc -o $qcPath $DataPath/$FASTQ2

# print end time
date