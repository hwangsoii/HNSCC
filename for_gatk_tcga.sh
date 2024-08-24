#!/bin/bash
#$ -S /bin/bash
#$ -cwd

## mandatory argument
## 1. project name
PROJECT=TRIUMPH
DIR=/data/project/TRIUMPH/tcga

## 2. sequencing platform
PLATFORM=illumina

## 3. out path
SubPath=$DIR/4.analysis/GATK/DOC
if [ ! -d $SubPath ]
then
	mkdir -p $SubPath
fi


## 4. data path

DataPath=${DIR}/1.raw
DataPath=${DIR}/3.aligned
DataPath=/data/public/TCGA/HNSC_WES

## 6. percentage of reads
PCT=50
BAM_TD=$(ls $DataPath | egrep '*-TD.bam$' )
BAM_ND=$(ls $DataPath | egrep '*-ND.bam$' )
# BAM=$(ls $DataPath | egrep '*.bam$' )
BAM_FP_LIST=(${BAM_TD// / })
BAM_BL_LIST=(${BAM_ND// / })
SCRIPT=/data/project/TRIUMPH/script/TRIUMPH
TEMP=$DIR/temp/gatkdoc
#TEMP=$DIR/temp/pureCN
TEMP=$DIR/temp/afterCNV
#TEMP=$DIR/temp/CNVfilter
# TEMP=$DIR/temp/afterCNV_loh
#TEMP=$DIR/temp/CNV_no_drop
#TEMP=$DIR/temp/HC
# TEMP=$DIR/temp/mutect2
# TEMP=$DIR/temp/strelka2
# REF=/data/resource/reference/human/UCSC/hg19/BWAIndex/genome.fa
REF=/data/resource/reference/human/UCSC/hg38/BWAIndex/genome.fa
REF=/data/project/TRIUMPH/tcga/download/reference/GRCh38.d1.vd1.fa
INTERVAL=/data/project/TRIUMPH/tcga/download/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed
INTERVAL=/data/project/TRIUMPH/tcga/download/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals_hg38.chradded.bed
if [ ! -d $TEMP ]
then
    	mkdir -p $TEMP
fi


# MANTA_ANALYSIS_PATH=${DIR}/4.analysis/Manta
# STRELKA_ANALYSIS_PATH=${DIR}/4.analysis/strelka

# rm -r -f ${MANTA_ANALYSIS_PATH}
# rm -r -f ${STRELKA_ANALYSIS_PATH}

if [ ${#BAM_FP_LIST[@]} -eq ${#BAM_BL_LIST[@]} ]
then
for idx in ${!BAM_FP_LIST[@]};
do
     #echo ${BAM_FP_LIST[idx]}
	SAMPLE=${BAM_FP_LIST[idx]%-TD.bam*}
	# echo ${SAMPLE}
	TUMOR=${DataPath}/${SAMPLE}-TD.bam
	NORMAL=${DataPath}/${SAMPLE}-ND.bam
	a=$(samtools view -H $NORMAL | grep '^@RG' | sed 's/.*SM:\([^\\t]*\).*/\1/')
	normal_sample=$(echo $a | cut -d ':' -f 13 | cut -d ' ' -f 1)
	b=$(samtools view -H ${TUMOR} | grep '^@RG' | sed 's/.*SM:\([^\\t]*\).*/\1/')
	tumor_sample=$(echo $b | cut -d ':' -f 13 | cut -d ' ' -f 1)
	echo $TUMOR $NORMAL $SAMPLE
	echo $SAMPLE $idx
	echo $normal_sample $tumor_sample
	echo $DIR ${SAMPLE} ${REF} ${TUMOR} ${NORMAL}
	# echo $SAMPLE >> /data/project/TRIUMPH/update_sample_list.txt
	# # # # #########HaplotypeCaller_pair
	# sleep 0.1 | qsub -pe smp 2 -e $TEMP -o $TEMP -N 'HP'$SAMPLE"-FP" ${SCRIPT}/pipe_haplotypecaller_pair.sh $DIR ${SAMPLE}'-FP' ${TUMOR} ${NORMAL} $REF 
	# # # #########HaplotypeCaller
	# sleep 0.1 | qsub -pe smp 2 -e $TEMP -o $TEMP -N 'HP'$SAMPLE"-FP" ${SCRIPT}/pipe_haplotypecaller.sh $DIR ${SAMPLE}'-FP' ${TUMOR} $REF 
	# sleep 0.1 | qsub -pe smp 2 -e $TEMP -o $TEMP -N 'HP'$SAMPLE"-BL" ${SCRIPT}/pipe_haplotypecaller.sh $DIR ${SAMPLE}'-BL' ${NORMAL} $REF
	####MAKEPON
	# sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N 'MU'$SAMPLE ${SCRIPT}/pipe_PON.sh $DIR ${SAMPLE} ${REF} ${TUMOR} ${NORMAL} ${SAMPLE}'-BL'
	###########Mutect2
	# sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N 'MU'$SAMPLE ${SCRIPT}/pipe_mutect2_tcga.sh $DIR ${SAMPLE} ${REF} ${TUMOR} ${NORMAL} ${SAMPLE}'-ND' $normal_sample
	# ############strelka
	# sleep 0.1 | qsub -pe smp 10 -e $TEMP -o $TEMP -N 'St'$SAMPLE ${SCRIPT}/pipe_strelka_tcga.sh $DIR ${SAMPLE} ${REF} ${TUMOR} ${NORMAL}
    # # ##########CNVkit
	# TEMP=$DIR/temp/CNV
	# qsub -pe smp 6 -e $TEMP -o $TEMP -hold_jid 'St'$SAMPLE -N 'cnv'$SAMPLE ${SCRIPT}/pipe_cnvkit_tcga.sh $DIR ${SAMPLE} ${TUMOR} ${NORMAL} ${SAMPLE}'-ND'
	##qsub -pe smp 6 -e $TEMP -o $TEMP -N 'cnv_no_drop'$SAMPLE ${SCRIPT}/pipe_cnvkit.sh $DIR ${SAMPLE} ${TUMOR} ${NORMAL} ${SAMPLE}'-BL'
	# echo $TEMP 'cnv'$SAMPLE ${SCRIPT}/pipe_cnvkit.sh $DIR ${SAMPLE} ${TUMOR} ${NORMAL} ${SAMPLE}'-BL'
	# SAMPLE_LIST="$SAMPLE_LIST $SAMPLE"
	# SAMPLE1=${TUMOR_LIST[idx]%_recal2.bam*}
	# SAMPLE2=${NORMAL_LIST[idx]%_recal2.bam*}
    ## A=`echo ${LD_LIBRARY_PATH}`
    # echo ${BAM_LIST[idx]}
    # echo ${SAMPLE1} ${SAMPLE2}
	#echo ${idx}
	# if [${idx}-eq2]
	# then
	# 	break
	# fi
	# #############pureCN
	# # echo "${SAMPLE}"
	# # sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N "pureCN_"${SAMPLE} ${SCRIPT}/pipe_pureCN.sh ${SAMPLE} 
	# ############afterCNV
	# sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -hold_jid 'cnv'$SAMPLE -N "afterCNV_"${SAMPLE} ${SCRIPT}/pipe_afterCNV_tcga.sh ${DIR} ${SAMPLE}
	# sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N "CNVpureCN_"${SAMPLE} ${SCRIPT}/pipe_CNVpureCN.sh ${SAMPLE}
	##########add sample vcf2maf
	# sleep 0.1 | qsub -e $TEMP -o $TEMP -N 'v2m'$SAMPLE ${SCRIPT}/vcf2maf.sh ${DIR} ${SubPath} ${DataPath} ${VCF_LIST[idx]} ${SAMPLE}
	####################gatkdoc
	# sleep 0.1 | qsub -pe smp 2 -e $TEMP -o $TEMP -N 'gatkdoc'$SAMPLE ${SCRIPT}/pipe_gatkdoc.sh ${PROJECT} ${REF} ${TUMOR} ${INTERVAL} ${DIR} ${SAMPLE}'-FP' ${SubPath}
	# sleep 0.1 | qsub -pe smp 2 -e $TEMP -o $TEMP -N 'gatkdoc'$SAMPLE ${SCRIPT}/pipe_gatkdoc.sh ${PROJECT} ${REF} ${NORMAL} ${INTERVAL} ${DIR} ${SAMPLE}'-BL' ${SubPath}
	done
else
	echo "Not matched pair data"
fi


#TEMP=$DIR/temp/compensation
#if [ ! -d $TEMP ]
#then
#		mkdir -p $TEMP
#fi
# REF=/data/resource/reference/human/UCSC/hg19/BWAIndex/genome.fa
#REF=/data/resource/reference/human/UCSC/hg38/BWAIndex/genome.fa

#SAMPLE='27-S042'
# # # #########HaplotypeCaller
# qsub -pe smp 2 -e $TEMP -o $TEMP -N 'HP'$SAMPLE"-FP" ${SCRIPT}/pipe_haplotypecaller.sh $DIR ${SAMPLE}'-FP' ${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam $REF 
# qsub -pe smp 2 -e $TEMP -o $TEMP -N 'HP'$SAMPLE"-BL" ${SCRIPT}/pipe_haplotypecaller.sh $DIR ${SAMPLE}'-BL' ${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam $REF
# ##########CNVkit
# qsub -pe smp 6 -e $TEMP -o $TEMP -N 'cnv'$SAMPLE ${SCRIPT}/pipe_cnvkit.sh $DIR ${SAMPLE} ${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam ${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam ${SAMPLE}'-BL'
# # echo $TEMP 'cnv'$SAMPLE ${SCRIPT}/pipe_cnvkit.sh $DIR ${SAMPLE} ${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam ${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam ${SAMPLE}'-BL'
# # ############strelka
# qsub -pe smp 10 -e $TEMP -o $TEMP -N 'St'$SAMPLE ${SCRIPT}/pipe_strelka.sh $DIR ${SAMPLE} ${REF} ${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam ${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam

#############Mutect2
# sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N 'MU'$SAMPLE ${SCRIPT}/pipe_mutect2.sh $DIR ${SAMPLE} ${REF} ${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam ${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam ${SAMPLE}'-BL'

#SAMPLE='27-S046'
# # #########HaplotypeCaller
# qsub -pe smp 2 -e $TEMP -o $TEMP -N 'HP'$SAMPLE"-FP" ${SCRIPT}/pipe_haplotypecaller.sh $DIR ${SAMPLE}'-FP' ${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam $REF 
# qsub -pe smp 2 -e $TEMP -o $TEMP -N 'HP'$SAMPLE"-BL" ${SCRIPT}/pipe_haplotypecaller.sh $DIR ${SAMPLE}'-BL' ${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam $REF
# ##########CNVkit
# qsub -pe smp 6 -e $TEMP -o $TEMP -N 'cnv'$SAMPLE ${SCRIPT}/pipe_cnvkit.sh $DIR ${SAMPLE} ${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam ${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam ${SAMPLE}'-BL'
# echo $TEMP 'cnv'$SAMPLE ${SCRIPT}/pipe_cnvkit.sh $DIR ${SAMPLE} ${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam ${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam ${SAMPLE}'-BL'
# ############strelka
# qsub -pe smp 10 -e $TEMP -o $TEMP -N 'St'$SAMPLE ${SCRIPT}/pipe_strelka.sh $DIR ${SAMPLE} ${REF} ${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam ${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam

#############Mutect2
# sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N 'MU'$SAMPLE ${SCRIPT}/pipe_mutect2.sh $DIR ${SAMPLE} ${REF} ${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam ${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam ${SAMPLE}'-BL'

# TUMOR=$(find $DataPath -name '*_recal2.bam' | egrep 'FP'| grep -v 'bai')
# NORMAL=$(find $DataPath -name '*_recal2.bam' | egrep 'BL'| grep -v 'bai')
# #TUMOR=$(ls $DataPath | egrep '_duplicates.bam' | egrep 'TD'| grep -v 'bai')
# #NORMAL=$(ls $DataPath | egrep '_duplicates.bam' | egrep 'BD'| grep -v 'bai')
# #FASTQ1=$(ls $DataPath | egrep '_1.fastq')
# #FASTQ2=$(ls $DataPath | egrep '_2.fastq')

# #FASTQ1=$(ls $DataPath | egrep '_1.fastq'| egrep 'C57BL6')
# #FASTQ2=$(ls $DataPath | egrep '_2.fastq'| egrep 'C57BL6')

# TUMOR_LIST=(${TUMOR// / })
# NORMAL_LIST=(${NORMAL// / })

# #FASTQ1_LIST=(${FASTQ1// / })
# #FASTQ2_LIST=(${FASTQ2// / })

# if [ ${#TUMOR_LIST[@]} -eq ${#NORMAL_LIST[@]} ]
# then
# for idx in ${!TUMOR_LIST[@]};
# do
# 	#SAMPLE1=${TUMOR_LIST[idx]%_sort_rg_marked_duplicates.bam*}
# 	SAMPLE1=${TUMOR_LIST[idx]%_recal2.bam*}
# 	SAMPLE2=${NORMAL_LIST[idx]%_recal2.bam*}
# 	#SAMPLE2=${NORMAL_LIST[idx]%_sort_rg_marked_duplicates.bam*}
# 	SAMPLE1_LIST="${SAMPLE1_LIST} ${SAMPLE1}"
# 	SAMPLE2_LIST="${SAMPLE2_LIST} ${SAMPLE2}"
#     echo $idx
# 	echo "TUMOR "${TUMOR_LIST[idx]}
# 	echo "NORMAL "${NORMAL_LIST[idx]}
# 	echo "sample" ${SAMPLE1}
# 	echo "sample" ${SAMPLE2}
# 	# qsub -e ${DIR}/temp/final/mutect -o ${DIR}/temp/final/mutect -N 'hsw_MUTECT2_'${SAMPLE} ${DIR}/script/TRIUMPH/Mutect_prev.sh $DIR $SubPath $DataPath ${TUMOR_LIST[idx]} ${NORMAL_LIST[idx]} ${SAMPLE1} ${SAMPLE2}
# 	done
# else
# 	echo "Not matched pair data"
# fi

