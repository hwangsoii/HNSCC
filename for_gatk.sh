#!/bin/bash
#$ -S /bin/bash
#$ -cwd

## mandatory argument
## 1. project name
PROJECT=TRIUMPH
DIR=/data/project/TRIUMPH

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

## 6. percentage of reads
PCT=50
BAM_FP=$(ls $DataPath | egrep '*-FP$' )
BAM_BL=$(ls $DataPath | egrep '*-BL$' )
# BAM=$(ls $DataPath | egrep '*.bam$' )
BAM_FP_LIST=(${BAM_FP// / })
BAM_BL_LIST=(${BAM_BL// / })
SCRIPT=/data/project/TRIUMPH/script/TRIUMPH
# TEMP=$DIR/temp/gatkdoc
TEMP=$DIR/temp/haplotypecaller_pair
TEMP=$DIR/temp/msisensor_ori
#TEMP=$DIR/temp/pureCN
#TEMP=$DIR/temp/afterCNV
#TEMP=$DIR/temp/CNVfilter
# TEMP=$DIR/temp/afterCNV_loh
#TEMP=$DIR/temp/CNV_no_drop
#TEMP=$DIR/temp/HC
# REF=/data/resource/reference/human/UCSC/hg19/BWAIndex/genome.fa
REF=/data/resource/reference/human/UCSC/hg38/BWAIndex/genome.fa
INTERVAL=/data/project/TRIUMPH/bed/1612AHP-0021_KimHR_3033241_Cho_MutScape_V2_1_Regions.liftover.sorted.bed
if [ ! -d $TEMP ]
then
    	mkdir -p $TEMP
fi

OUTPATH=/data/project/TRIUMPH/4.analysis/msisensor
if [ ! -d $OUTPATH ]
then
    	mkdir -p $OUTPATH
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
	SAMPLE=${BAM_FP_LIST[idx]%-FP*}
	# echo ${SAMPLE}
	TUMOR=${DataPath}/${SAMPLE}'-FP'/${SAMPLE}'-FP'_recal2.bam
	NORMAL=${DataPath}/${SAMPLE}'-BL'/${SAMPLE}'-BL'_recal2.bam
	echo $TUMOR $NORMAL $SAMPLE
	echo $SAMPLE $idx
	# echo $SAMPLE >> /data/project/TRIUMPH/update_sample_list.txt
		# # #########HaplotypeCaller_pair
	# sleep 0.1 | qsub -pe smp 2 -e $TEMP -o $TEMP -N 'HP'$SAMPLE"-FP" ${SCRIPT}/pipe_haplotypecaller_pair.sh $DIR ${SAMPLE}'-FP' ${TUMOR} ${NORMAL} $REF 
	# # # #########HaplotypeCaller
	# sleep 0.1 | qsub -pe smp 2 -e $TEMP -o $TEMP -N 'HP'$SAMPLE"-FP" ${SCRIPT}/pipe_haplotypecaller.sh $DIR ${SAMPLE}'-FP' ${TUMOR} $REF 
	# sleep 0.1 | qsub -pe smp 2 -e $TEMP -o $TEMP -N 'HP'$SAMPLE"-BL" ${SCRIPT}/pipe_haplotypecaller.sh $DIR ${SAMPLE}'-BL' ${NORMAL} $REF
	####MAKEPON
	# sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N 'MU'$SAMPLE ${SCRIPT}/pipe_PON.sh $DIR ${SAMPLE} ${REF} ${TUMOR} ${NORMAL} ${SAMPLE}'-BL'
	############Mutect2
	# sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N 'MU'$SAMPLE ${SCRIPT}/pipe_mutect2.sh $DIR ${SAMPLE} ${REF} ${TUMOR} ${NORMAL} ${SAMPLE}'-BL'
	# ############strelka
	# # sleep 0.1 | qsub -pe smp 10 -e $TEMP -o $TEMP -N 'St'$SAMPLE ${SCRIPT}/pipe_strelka.sh $DIR ${SAMPLE} ${REF} ${TUMOR} ${NORMAL}
    # # ##########CNVkit
	# qsub -pe smp 6 -e $TEMP -o $TEMP  -N 'cnv'$SAMPLE ${SCRIPT}/pipe_cnvkit.sh $DIR ${SAMPLE} ${TUMOR} ${NORMAL} ${SAMPLE}'-BL'
	
	###MSIsensor
	MSI_bed=/data/project/TRIUMPH/download/msisensor-pro/hg38_loci.txt
	# qsub -pe smp 2 -e $TEMP -o $TEMP -N 'MSI'$SAMPLE ${SCRIPT}/msisensor/Pipe_msisensor_pro.sh 'YS' 'TN' ${SAMPLE} ${TUMOR} ${NORMAL} ${OUTPATH} ${INTERVAL} ${MSI_bed}
	qsub -pe smp 2 -e $TEMP -o $TEMP -N 'MSI'$SAMPLE ${SCRIPT}/msisensor/pipe_msisensor.sh 'YS' 'TN' ${SAMPLE} ${TUMOR} ${NORMAL} ${OUTPATH} ${INTERVAL} ${MSI_bed}
	# ##qsub -pe smp 6 -e $TEMP -o $TEMP -N 'cnv_no_drop'$SAMPLE ${SCRIPT}/pipe_cnvkit.sh $DIR ${SAMPLE} ${TUMOR} ${NORMAL} ${SAMPLE}'-BL'
	# # echo $TEMP 'cnv'$SAMPLE ${SCRIPT}/pipe_cnvkit.sh $DIR ${SAMPLE} ${TUMOR} ${NORMAL} ${SAMPLE}'-BL'
	# # SAMPLE_LIST="$SAMPLE_LIST $SAMPLE"
	# # SAMPLE1=${TUMOR_LIST[idx]%_recal2.bam*}
	# # SAMPLE2=${NORMAL_LIST[idx]%_recal2.bam*}
    # ## A=`echo ${LD_LIBRARY_PATH}`
    # # echo ${BAM_LIST[idx]}
    # # echo ${SAMPLE1} ${SAMPLE2}
	# #echo ${idx}
	# # if [${idx}-eq2]
	# # then
	# # 	break
	# # fi
	# #############pureCN
	# # echo "${SAMPLE}"
	# # sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N "pureCN_"${SAMPLE} ${SCRIPT}/pipe_pureCN.sh ${SAMPLE} 
	# ############afterCNV
	# sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -hold_jid 'cnv'$SAMPLE -N "afterCNV_"${SAMPLE} ${SCRIPT}/pipe_afterCNV.sh ${SAMPLE}
	# sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N "CNVpureCN_"${SAMPLE} ${SCRIPT}/pipe_CNVpureCN.sh ${SAMPLE}
	##########add sample vcf2maf
	# sleep 0.1 | qsub -e $TEMP -o $TEMP -N 'v2m'$SAMPLE ${SCRIPT}/vcf2maf.sh ${DIR} ${SubPath} ${DataPath} ${VCF_LIST[idx]} ${SAMPLE}
	# ####################gatkdoc
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

