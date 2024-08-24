#!/bin/bash
#$ -S /bin/bash
#$ -cwd

DIR=$1
SAMPLE=$2

CNVkit_ANALYSIS_PATH=${DIR}/4.analysis/CNVkit/${SAMPLE}/results
VCF_DIR=/data/project/TRIUMPH/4.analysis/GATK/pureCN_Mu


#cnvkit.py genemetrics -y ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cnr" -s ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cns" | tail -n+2 | cut -f1 | sort > ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"_"segment-genes.txt

#cnvkit.py genemetrics -y ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cnr" | tail -n+2 | cut -f1 | sort > ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"_"ratio-genes.txt

#comm -12 ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"_"ratio-genes.txt \
# ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"_"segment-genes.txt > ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"_"trusted-genes.txt

#for gene in `cat ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"_"trusted-genes.txt`
#do
#	cnvkit.py scatter -s ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2".cn{s,r} -g $gene -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-${gene}-scatter.pdf
#done


#cnvkit.py metrics ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_recal2.targetcoverage.cnn ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_recal2.antitargetcoverage.cnn \
# ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_recal2.cnr \
# -s ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_recal2.cns  \
# > ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"_"metrics.txt

#cnvkit.py segmetrics ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cnr" -s ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2".cns --iqr 
#cnvkit.py segmetrics -s ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2".cn{s,r} --ci --pi

#cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cns" --drop-low-coverage --filter cn -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_cn.call.cns
#cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cns" --drop-low-coverage --filter ci -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_ci.call.cns
#cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cns" --drop-low-coverage --filter ampdel -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_ampdel.call.cns
#cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cns" --drop-low-coverage --filter ampdel -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_ampdel.call.cns
#cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cns" --drop-low-coverage --filter sem -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_sem.call.cns
#cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cns" --drop-low-coverage -m threshold -t=-1.1,-0.4,0.3,0.7 --filter cn -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_cn.call.cns

cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cns" --filter cn -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_cn.call.cns
cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cns" --filter ci -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_ci.call.cns
cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}"-FP_recal2.cns" --filter ampdel -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_ampdel.call.cns



#################################################################################################

# cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_recal2.cns --drop-low-coverage --vcf $VCF_DIR/${SAMPLE}.mutect2.somatic.confident.vcf -i ${SAMPLE}-FP -n ${SAMPLE}-BL -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_loh.call.cns

# cnvkit.py call ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_recal2.cns --drop-low-coverage --filter cn --vcf $VCF_DIR/${SAMPLE}.mutect2.somatic.confident.vcf -i ${SAMPLE}-FP -n ${SAMPLE}-BL -o ${CNVkit_ANALYSIS_PATH}/${SAMPLE}-FP_loh.cn.call.cns
