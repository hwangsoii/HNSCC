#!/bin/bash
#$ -S /bin/bash
#$ -cwd

# Define Variables
PICARD=/opt/Yonsei/picard/2.22.0/picard.jar
REFERENCE_GENOME=/data/resource/reference/human/UCSC/hg38/WholeGenomeFasta/genome.fa
REFERENCE_GENOME=/data/project/TRIUMPH/reference/hg38WGSFASTA/genome.fa
REF=/data/project/TRIUMPH/tcga/download/reference/GRCh38.d1.vd1.fa
# REFERENCE_GENOME=/data/resource/reference/human/UCSC/hg19/BWAIndex/genome.fa

# basic argument
# REF=/data/public/LongRanger/GRCh38_2.1.0/fasta/genome.fa
# INTERVAL=/home/hswbulls/project/MCD/Panel/0.85x_rsID_rescue_trial_1_Covered.refine.interval_list
INTERVAL=/data/project/TRIUMPH/bed/1612AHP-0021_KimHR_3033241_Cho_MutScape_V2_1_Regions.liftover.sorted.bed
INTERVAL=/data/project/TRIUMPH/tcga/download/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals_hg38.chradded.sorted.bed
refflat=/data/project/TRIUMPH/download/refene_UCSC/hg38/refFlat.txt
# INTERVAL=/data/project/TRIUMPH/bed/1612AHP-0021_KimHR_3033241_Cho_MutScape_V2_1_Regions.removechr.bed
DIR=$1
ID=$2
# REF=$3
TUMORPATH=$3
NORMALPATH=$4
NORMAL=$5
CNVkit_ANALYSIS_PATH=${DIR}/4.analysis/CNVkit
if [ ! -d $CNVkit_ANALYSIS_PATH ]
then
    	mkdir -p $CNVkit_ANALYSIS_PATH
fi
CNVkit_ANALYSIS_PATH=${DIR}/4.analysis/CNVkit/${ID}
if [ ! -d $CNVkit_ANALYSIS_PATH ]
then
    	mkdir -p $CNVkit_ANALYSIS_PATH
fi
STRELKA_ANALYSIS_PATH=${DIR}/4.analysis/strelka



# From baits and tumor/normal BAMs
cnvkit.py batch ${TUMORPATH} --normal ${NORMALPATH} \
--targets ${INTERVAL} --annotate ${refflat} \
--fasta ${REFERENCE_GENOME} \
--drop-low-coverage \
--output-reference ${CNVkit_ANALYSIS_PATH}/my_reference.cnn --output-dir ${CNVkit_ANALYSIS_PATH}/results/ \
--diagram --scatter



# --access /data/project/TRIUMPH/download/poor_mapability/hg38/access-5kb-mappable.hg38.bed
# ####with -p, parralel processing
# # cnvkit.py batch *.bam -r my_reference.cnn -p 8

# # Reusing a reference for additional samples
# cnvkit.py batch *Tumor.bam -r Reference.cnn -d results/

# # Reusing targets and antitargets to build a new reference, but no analysis
# cnvkit.py batch -n *Normal.bam --output-reference new_reference.cnn \
#     -t my_targets.bed -a my_antitargets.bed \
#     -f hg19.fasta -g data/access-5kb-mappable.hg19.bed


# cnvkit.py batch *.bam -r my_reference.cnn -p 8


# cnvkit.py target my_baits.bed --annotate refFlat.txt --split -o my_targets.bed


###############################no low coverage drop out#########################


# cnvkit.py batch ${TUMORPATH} --normal ${NORMALPATH} \
#  --targets ${INTERVAL} --annotate ${refflat} \
#  --fasta ${REFERENCE_GENOME} \
#  --output-reference ${CNVkit_ANALYSIS_PATH}/my_reference.cnn --output-dir ${CNVkit_ANALYSIS_PATH}/results/ \
#  --diagram --scatter
