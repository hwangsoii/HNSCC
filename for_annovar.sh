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
###HC
DataPath=${DIR}/4.analysis/GATK/HaplotypeCaller
# DataPath=/data/project/TRIUMPH/download/Korea1K
##Mutect2
# DataPath=${DIR}/4.analysis/GATK/Mutect2
###
DataPath=${DIR}/4.analysis/Korean_filter/HaplotypeCaller/integ
DataPath=${DIR}/4.analysis/Korean_filter/HaplotypeCaller/integ_all_220625_korea_filter
DataPath=${DIR}/4.analysis/Korean_filter/HaplotypeCaller/integ_all_220816_korea_filter
DataPath=${DIR}/4.analysis/Korean_filter/HaplotypeCaller/integ_all_230816_korea_filter

# DataPath=/data/project/TRIUMPH/4.analysis/Korean_filter/HaplotypeCaller/integ_all_hard_freq
# DataPath=${DIR}/4.analysis/vep
#DataPath=${DIR}/4.analysis/GATK/PON
# 6. percentage of reads
PCT=50
BAM_FP=$(ls $DataPath | egrep '*-FP$' )
BAM_BL=$(ls $DataPath | egrep '*-BL$' )
BAM=$(ls $DataPath | egrep '*.bam$' )
VCF=$(ls ${DataPath} | egrep '*_mutect2_filtered_ann.vcf$')
#####HC
VCF=$(ls ${DataPath} | egrep '*-FP.Haplot.Haplotype.raw.snps.filtered_snps.PASS.vcf$')
VCF=$(ls ${DataPath} | egrep '*-BL.Haplot.Haplotype.raw.snps.filtered_snps.PASS.vcf$')
VCF=$(ls ${DataPath} | egrep '*.Haplotype.raw.snps.filtered_snps.PASS.integ.vcf$')
VCF=$(ls ${DataPath} | egrep '.Haplotype.raw.snps.indels.All.PASS.integ.sort.vcf$')
VCF=$(ls ${DataPath} | egrep '.Haplotype.raw.snps.indels.All.PASS.integ.sort.vcf$')
# VCF=$(ls ${DataPath} | egrep '.normalized.vcf$')
###mutect2 all
# VCF=$(ls ${DataPath} | egrep '*.mutect2.somatic.confident.all.PASS.vcf$')

# VCF=$(ls ${DataPath} | egrep '*.KRGDB.Korea1K.vcf$')
# VCF=$(ls ${DataPath} | egrep '*.Haplotype.raw.snps.filtered_snps.PASS.integ.Korea_filter.vcf$')
VCF=$(ls ${DataPath} | egrep '*.Haplotype.raw.snps.filtered_snps.PASS.integ.KRGDB.Korea1K.vcf$')
VCF=$(ls ${DataPath} | egrep '*.All.PASS.integ.sort.Korea_filter.vcf$')
# VCF=$(ls ${DataPath} | egrep '*.All.PASS.integ.Korea_filter.vcf$')
VCF=$(ls ${DataPath} | egrep '*.All.PASS.integ.Korea_filter.vcf$')

# VCF=$(ls ${DataPath} )
# FOLD=$(find $DataPath -name '*pon.vcf.gz')
# BAM_FP_LIST=(${BAM_FP// / })
# BAM_BL_LIST=(${BAM_BL// / })
VCF_LIST=(${VCF// / })
SCRIPT=/data/project/TRIUMPH/script/TRIUMPH
SCRIPT_korea=/data/project/TRIUMPH/script/korea
# SCRIPT=/home/hswbulls/script/TRIUMPH
REF=/data/resource/reference/human/UCSC/hg19/BWAIndex/genome.fa
annotator='ANNOVAR'
annotator='vcf2maf'
call='germline'
# call='somatic'
# call='Korea1k'
TEMP=$DIR/temp/${annotator}
if [ ! -d $TEMP ]
then
    	mkdir -p $TEMP
fi
TEMP=$DIR/temp/${annotator}/${call}
if [ ! -d $TEMP ]
then
    	mkdir -p $TEMP
fi
TEMP=$DIR/temp/${annotator}/${call}/integ_230816
if [ ! -d $TEMP ]
then
    	mkdir -p $TEMP
fi



SubPath=/data/project/TRIUMPH/4.analysis/${annotator}
if [ ! -d $SubPath ]
then
    	mkdir -p $SubPath
fi
SubPath=/data/project/TRIUMPH/4.analysis/${annotator}/${call}
if [ ! -d $SubPath ]
then
    	mkdir -p $SubPath
fi
# SubPath=/data/project/TRIUMPH/4.analysis/${annotator}/${call}/all
# if [ ! -d $SubPath ]
# then
#     	mkdir -p $SubPath
# fi
SubPath=/data/project/TRIUMPH/4.analysis/${annotator}/${call}/integ_all
SubPath=/data/project/TRIUMPH/4.analysis/${annotator}/${call}/integ_all_220625
SubPath=/data/project/TRIUMPH/4.analysis/${annotator}/${call}/integ_all_220816
SubPath=/data/project/TRIUMPH/4.analysis/${annotator}/${call}/integ_all_230816
if [ ! -d $SubPath ]
then
    	mkdir -p $SubPath
fi

# if [ ${#BAM_FP_LIST[@]} -eq ${#BAM_BL_LIST[@]} ]
# then

echo $VCF
for idx in ${!VCF_LIST[@]};
do
    # echo ${BAM_FP_LIST[idx]}
	# SAMPLE=${VCF_LIST[idx]%.mutect2.somatic.confident.all.PASS.vcf*}
	# SAMPLE=${VCF_LIST[idx]%.mutect2.somatic.confident.SNV.PASS.vcf*}
	# SAMPLE=${VCF_LIST[idx]%-FP.Haplot.Haplotype.raw.snps.filtered_snps.PASS.vcf*}
	SAMPLE=${VCF_LIST[idx]%.Haplotype.raw.snps.indels.All.PASS.integ.sort.vcf*}
	SAMPLE=${VCF_LIST[idx]%.Haplotype.raw.snps.filtered_snps.PASS.integ.KRGDB.Korea1K.vcf*}
	SAMPLE=${VCF_LIST[idx]%-FP.All.PASS.integ.Korea_filter.vcf*}
	# SAMPLE=${VCF_LIST[idx]%.Haplotype.raw.snps.indel.PASS.integ.KRGDB.Korea1K.vcf*}
	# SAMPLE=${VCF_LIST[idx]%.recal.normalized.vcf*}
	# SAMPLE=${VCF_LIST[idx]%.All.PASS.integ.sort.Korea_filter.vcf*}
	#sort
	# cat ${DataPath}/${VCF_LIST[idx]} | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > ${DataPath}/${SAMPLE}.Haplotype.raw.snps.indels.All.PASS.integ.sort.vcf
	# SAMPLE=${VCF_LIST[idx]%.Haplotype.raw.snps.filtered_snps.PASS.integ.KRGDB.Korea1K.vcf*}
	# SAMPLE=${VCF_LIST[idx]%.Haplotype.raw.snps.All.normalized.PASS.integ.sort.KRGDB.Korea1K.vcf*}
	# SAMPLE=${VCF_LIST[idx]%.Haplotype.raw.snps.All.normalized.PASS.integ.sort.vcf*}
    echo ${VCF_LIST[idx]}
	echo ${SAMPLE} ${VCF_LIST[idx]}
	echo ${idx}
    # echo ${SAMPLE:7}
    # mv ${DataPath}/${VCF_LIST[idx]} ${DataPath}/${VCF_LIST[idx]}.mutect2.somatic.confident.vcf
	#echo  $DataPath/${VCF_LIST[idx]} >> ${DIR}/normals_for_pon_vcf.args
	#echo $SAMPLE >> /data/project/TRIUMPH/update_sample_list.txt
	#####MAKEPON
	#sleep 0.1 | qsub -pe smp 6 -e $TEMP -o $TEMP -N 'MU'$SAMPLE ${SCRIPT}/pipe_PON.sh $DIR ${SAMPLE} ${REF} ${TUMOR} ${NORMAL} ${SAMPLE}'-BL'
    ####VEP
    # sleep 0.1 | qsub -e $TEMP -o $TEMP -N 'vep'$SAMPLE ${SCRIPT}/VEP_pipeline.sh ${DIR} ${SubPath} ${DataPath} ${VCF_LIST[idx]} ${SAMPLE}
    #vcf2maf
    sleep 0.1 | qsub -e $TEMP -o $TEMP -N 'v2m'$SAMPLE ${SCRIPT}/vcf2maf.sh ${DIR} ${SubPath} ${DataPath} ${VCF_LIST[idx]} ${SAMPLE}
    ###annovar
    # sleep 0.1 | qsub -e $TEMP -o $TEMP -N 'ANV'$SAMPLE ${SCRIPT}/pipe_ANNOVAR.sh ${DIR} ${SubPath} ${DataPath} ${VCF_LIST[idx]} ${SAMPLE}
    ##parsing
    # sleep 0.1 | qsub -e $TEMP -o $TEMP -N 'ANV'$SAMPLE ${SCRIPT}/pass_finder.sh ${DIR} ${SubPath} ${DataPath} ${VCF_LIST[idx]} ${SAMPLE}
	# ###Korean_filter
	# TEMP=$DIR/temp/Korean/all_0208
	# if [ ! -d $TEMP ]
	# then
	# 		mkdir -p $TEMP
	# fi
	# SubPath=/data/project/TRIUMPH/4.analysis/Korean_filter
	# if [ ! -d $SubPath ]
	# then
	# 		mkdir -p $SubPath
	# fi
	# SubPath=/data/project/TRIUMPH/4.analysis/Korean_filter/HaplotypeCaller
	# if [ ! -d $SubPath ]
	# then
	# 		mkdir -p $SubPath
	# fi
	# SubPath=/data/project/TRIUMPH/4.analysis/Korean_filter/HaplotypeCaller/integ_all
	# if [ ! -d $SubPath ]
	# then
	# 		mkdir -p $SubPath
	# fi
	# sleep 0.1 | qsub -e $TEMP -o $TEMP -N 'Kor'$SAMPLE ${SCRIPT_korea}/pipe_annotate_koreanDBs.sh ${DIR} ${SCRIPT_korea} ${DataPath} ${VCF_LIST[idx]} ${SubPath}
done
# else
# 	echo "Not matched pair data"
# fi

#qsub -pe smp 6 -e $TEMP -o $TEMP -N 'PONMERGE' ${SCRIPT}/PON_merge.sh
