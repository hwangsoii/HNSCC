#!/bin/bash
#$ S /bin/bash
#$ -cwd

## mandatory argument
# 1. project name
PROJECT=TRIUMPH
DIR=/data/project/TRIUMPH

# 2. sequencing platform
PLATFORM=illumina

# 3. out path
SubPath=$DIR/result/real_final

# 4. data path

DataPath=$DIR/1.raw/add_220623
DataPath=${DIR}/2.qc/add

SCRIPT=/home/hswbulls/script/TRIUMPH
SCRIPT=/data/project/TRIUMPH/script/TRIUMPH

FASTQ1=$(ls $DataPath | egrep '_1.fastq')
FASTQ2=$(ls $DataPath | egrep '_2.fastq')
FASTQ1=$(ls $DataPath | egrep '_1.fastq' | grep -v 'BL')
FASTQ2=$(ls $DataPath | egrep '_2.fastq'| grep -v 'BL')
FASTQ1=$(ls $DataPath | egrep '_1.fastq' | grep -v 'FP')
FASTQ2=$(ls $DataPath | egrep '_2.fastq'| grep -v 'FP')

#FASTQ1=$(ls $DataPath | egrep '_1.fastq'| egrep 'C57BL6')
#FASTQ2=$(ls $DataPath | egrep '_2.fastq'| egrep 'C57BL6')
TEMP=$DIR/temp
if [ ! -d $TEMP ]
then
    	mkdir -p $TEMP
fi
TEMP=$DIR/temp/bwa_20220623
if [ ! -d $TEMP ]
then
    	mkdir -p $TEMP
fi

FASTQ1_LIST=(${FASTQ1// / })
FASTQ2_LIST=(${FASTQ2// / })
# SubPath=$DIR/3.aligned
SubPath=$DIR/2.qc/add
if [ ! -d $SubPath ]
then
    	mkdir -p $SubPath
fi
if [ ${#FASTQ1_LIST[@]} -eq ${#FASTQ2_LIST[@]} ]
then
    for idx in ${!FASTQ1_LIST[@]};
    do
        SAMPLE=${FASTQ1_LIST[idx]%_1*}
        SAMPLE_LIST="$SAMPLE_LIST $SAMPLE"
        echo ${FASTQ1_LIST[idx]}
        echo ${FASTQ2_LIST[idx]}
        echo ${SAMPLE}  
        echo ${idx}
        # ##############qc
        # sleep 0.1 | qsub -e $TEMP -o $TEMP -N 'qc'$SAMPLE ${SCRIPT}/pipe_fastqc.sh ${DIR} ${FASTQ1_LIST[idx]} ${FASTQ2_LIST[idx]} ${DataPath}
        # #########fastp
        # SubPath=$DIR/2.qc/add
        # sleep 0.1 | qsub -e $TEMP -o $TEMP -N 'fastp'$SAMPLE ${SCRIPT}/pipe_fastp.sh ${DIR} ${FASTQ1_LIST[idx]} ${FASTQ2_LIST[idx]} ${DataPath} ${SAMPLE} ${SubPath}
        #########bwa
        SubPath=$DIR/3.aligned/add/
        if [ ! -d $SubPath ]
        then
                mkdir -p $SubPath
        fi
        SubPath=$DIR/3.aligned/add/${SAMPLE}
        if [ ! -d $SubPath ]
        then
                mkdir -p $SubPath
        fi
        sleep 0.1 | qsub -pe smp 2 -e $TEMP -o $TEMP -N 'bwaFP'$SAMPLE ${SCRIPT}/pipe_bwa.sh $DIR $SubPath $DataPath ${FASTQ1_LIST[idx]} ${FASTQ2_LIST[idx]} ${SAMPLE} ${PROJECT} ${PLATFORM}
	# #################HaplotypeCaller
        # REF=/data/resource/reference/human/UCSC/hg19/BWAIndex/genome.fa
        # AlignedPath=$DIR/3.aligned/${SAMPLE}
        # BAMPATH=${AlignedPath}/${SAMPLE}'_recal2.bam'
        # AnalysisPath=${DIR}/4.analysis/${SAMPLE}/

        done

else
        echo "Not matched pair data."
fi


# qsub -pe smp 2 -e $TEMP -o $TEMP -N 'bwaBL'$SAMPLE ${SCRIPT}/pipe_bwa.sh $DIR $SubPath $DataPath '27-S042-FP_1.fastq.gz' '27-S042-FP_2.fastq.gz' '27-S042' ${PROJECT} ${PLATFORM}
# qsub -pe smp 2 -e $TEMP -o $TEMP -N 'bwaBL'$SAMPLE ${SCRIPT}/pipe_bwa.sh $DIR $SubPath $DataPath '27-S046-FP_1.fastq.gz' '27-S046-FP_2.fastq.gz' '27-S046' ${PROJECT} ${PLATFORM}
                                                         
