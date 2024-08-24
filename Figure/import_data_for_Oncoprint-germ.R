library(maftools)
library(readxl)
library(ComplexHeatmap)
library(colorspace)
library(RColorBrewer)
source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16//Figure/oncoprint_matrix_build.R")    ###original code
source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16//Figure/H15_inter_sample.R")
source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16//Figure/msi_reader.R")

###maf data
fileNames=Sys.glob("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/all_220625/*vep.maf")
fileNames_vaf=Sys.glob("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/all_220625_vaf/*vep.vaf.maf")
# fileNames=Sys.glob("/data/project/TRIUMPH/4.analysis/vcf2maf/somatic/all_220625/filter_maf/*vep.maf")
# fileNames_vaf=Sys.glob("/data/project/TRIUMPH/4.analysis/vcf2maf/somatic/all_220625_vaf/filter_maf/*vep.vaf.maf")
nonhnsccfp <- c("01-S005-FP", "21-S005-FP", "22-S011-FP", "22-S015-FP", "27-S024-FP", "28-S004-FP", "31-S013-FP", "46-S010-FP")
nonhnscc <- c("01-S005", "21-S005", "22-S011", "22-S015", "27-S024", "28-S004", "31-S013", "46-S010")
# nonhnscc_somatic<- paste0("/data/project/TRIUMPH/4.analysis/vcf2maf/somatic/all_220625/filter_maf/", nonhnscc, ".vep.maf")
nonhnscc_somatic<- paste0("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/all_220625/", nonhnscc, ".vep.maf")
nonhnscc_germ <- paste0("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/germline/integ_all_220625/", nonhnscc, ".vep.black.soft.maf")
nonhnscc_germ4 <- paste0("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/germline/integ_all_220625/", nonhnscc, ".vep.black.hard.maf")
# fileNames=Sys.glob("/data/project/TRIUMPH/4.analysis/vcf2maf/somatic/all_220625/*vep.maf")
fileNames_onlyhnscc <- fileNames[! fileNames %in% nonhnscc_somatic]
nonhnscc_somatic<- paste0("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/all_220625_vaf/", nonhnscc, ".vep.vaf.maf")
fileNames_onlyhnscc_vaf <- fileNames_vaf[! fileNames_vaf %in% nonhnscc_somatic]


###maf data
# fileNames=Sys.glob("/data/project/TRIUMPH/4.analysis/vcf2maf/somatic/all_220625/*vep.maf")
nonhnsccfp <- c("01-S005-FP", "21-S005-FP", "22-S011-FP", "22-S015-FP", "27-S024-FP", "28-S004-FP", "31-S013-FP", "46-S010-FP")
nonhnscc <- c("01-S005", "21-S005", "22-S011", "22-S015", "27-S024", "28-S004", "31-S013", "46-S010")
# nonhnscc_somatic<- paste0("/data/project/TRIUMPH/4.analysis/vcf2maf/somatic/all_220625/", nonhnscc, ".vep.maf")
nonhnscc_germ <- paste0("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/germline/integ_all_220816/", nonhnscc, ".vep.black.soft.maf")
nonhnscc_germ4 <- paste0("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/germline/integ_all_220816/", nonhnscc, ".vep.black.hard.maf")
# fileNames=Sys.glob("/data/project/TRIUMPH/4.analysis/vcf2maf/somatic/all_220625/*vep.maf")
# fileNames_onlyhnscc <- fileNames[! fileNames %in% nonhnscc_somatic]
# fileNames3=Sys.glob("/data/project/TRIUMPH/4.analysis/vcf2maf/germline/integ_all_220816/*.vep.black.soft.maf")
# fileNames3_onlyhnscc <- fileNames3[! fileNames3 %in% nonhnscc_germ]
nonhnsccfp <- c("01-S005-FP", "21-S005-FP", "22-S011-FP", "22-S015-FP", "27-S024-FP", "28-S004-FP", "31-S013-FP", "46-S010-FP")
nonhnscc <- c("01-S005", "21-S005", "22-S011", "22-S015", "27-S024", "28-S004", "31-S013", "46-S010")
fileNames4=Sys.glob("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/germline/integ_all_220816/*.vep.black.hard.maf")
# fileNames4=Sys.glob("/data/project/TRIUMPH/4.analysis/vcf2maf/germline/integ_all_230816/*.vep.black.hard.maf")
fileNames4_onlyhnscc <- fileNames4[! fileNames4 %in% nonhnscc_germ4]

###CNV data
cnvkit_0.4 <-  read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/CNV_filter_220626.txt", header = TRUE,  sep = "\t")
cnvkit_0.4_onlyhnscc <- cnvkit_0.4[!(cnvkit_0.4$sample %in% nonhnsccfp), ]


###Clinical data
total_clinical <- data.frame(read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/final_clinical_merged_221208_removal_nonhnscc_p16.txt", header = T, sep = '\t'))
total_clinical <- merge(total_clinical, result_df, by= 'Sample2')
total_clinical$MS

# ###add
# mRNA <- c('EGFR', '15-S034-FP', 'Amplified_expression')
# mRNA2 <- c('EGFR', '15-S038-FP', 'Amplified_expression')
# mRNA3 <- c('ERBB3', '27-S059-FP', 'Amplified_expression')
# mRNA4 <- c('FGFR1', '22-S008-FP', 'Amplified_expression')
# mRNA5 <- c('FGFR1', '15-S015-FP', 'Amplified_expression')
# doubt1 <- c('EGFR', '25-S020-FP', 'Doubtful_call')
# doubt2 <- c('PIK3R1', '22-S016-FP', 'Doubtful_call')
# doubt3 <- c('PIK3CA', '24-S022-FP', 'Doubtful_call')
# doubt4 <- c('CDKN2A', '31-S022-FP', 'Doubtful_call')
# doubt5 <- c('PIK3CA', '31-S006-FP', 'Doubtful_call')
# doubt6 <- c('FGFR3', '47-S013-FP', 'Doubtful_call')
# mild_amp <- c('EGFR', '47-S010-FP', 'Mild_amp')
# mild_amp2 <- c('CCND1', '28-S001-FP', 'Mild_amp')
# mild_amp3 <- c('EGFR', '31-S015-FP', 'Mild_amp')
# mild_amp5 <- c('FGFR1', '27-S075-FP', 'Mild_amp')
# mild_amp6 <- c('CDKN2A', '20-S009-FP', 'Mild_amp')
# tot_cnv <- rbind(cnvkit_0.4_onlyhnscc, mRNA, mRNA2, mRNA3, mRNA4, mRNA5, doubt1, doubt2,
#                  doubt3, doubt4, doubt5, doubt6, mild_amp, mild_amp2, mild_amp3, mild_amp5,
#                  mild_amp6)


#merge maf
All_cnvkit = merge_mafs(fileNames_onlyhnscc, clinicalData = total_clinical, cnTable = tot_cnv)
# getSampleSummary(All_cnvkit)
# getGeneSummary(All_cnvkit)
All_germ = merge_mafs(fileNames4_onlyhnscc, clinicalData = total_clinical)
# pdf("/data/project/TRIUMPH/Figure2/new_figure_manu/Figure_germ.pdf", width = 30, height = 17)
# oncoplot(All_germ, top=30)
# dev.off()
# All_germ@data

#matrix form
mat_tot <- oncoprint_matrix_build(All_cnvkit)
mat_tot <- mat_tot[,order(colnames(mat_tot))]
dim(mat_tot)
mat_tot = as.matrix(mat_tot)
mat_tot <- gsub(";Amp",'|Amp',mat_tot)
mat_tot <- gsub(";Del",'|Del',mat_tot)
mat_tot[grep(';',mat_tot)] <- paste0(mat_tot[grep(';',mat_tot)],";Multi_hit")
mat_tot <- gsub('\\|Amplification;Multi_hit', ';Multi_hit|Amplification', mat_tot)

###mat_tot_germ
mat_tot_germ <- oncoprint_matrix_build(All_germ)
mat_tot_germ  <- mat_tot_germ[,order(colnames(mat_tot_germ))]
mat_tot_germ = as.matrix(mat_tot_germ)
mat_tot_germ <- gsub(";Amp",'|Amp',mat_tot_germ)
mat_tot_germ <- gsub(";Del",'|Del',mat_tot_germ)
mat_tot_germ[grep(';',mat_tot_germ)] <- paste0(mat_tot_germ[grep(';',mat_tot_germ)],";Multi_hit")
mat_tot_germ <- gsub('\\|Amplification;Multi_hit', ';Multi_hit|Amplification', mat_tot_germ)

###clinical data of patient in mat_tot
inter_clinical <- subset(total_clinical, total_clinical$Tumor_Sample_Barcode%in%colnames(mat_tot));inter_clinical <- inter_clinical[order(inter_clinical$Tumor_Sample_Barcode),]


###importing oncoplot top30 gene
setwd("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/")
oncoplot(maf = All_cnvkit, writeMatrix = T, top = 30);top30_genes <- rownames(read.table(file = "onco_matrix.txt", sep = "\t"))
mat30 = mat_tot[top30_genes, ]


######matrix sorting function
sample_order <- function(mat){
  val_list = numeric()
  for(i in 1:length(mat[1,])){
    num = seq(1:30)[mat[,i] != ""]
    value = 0
    for(j in num){
      value = value + 1/(10^j)
    }
    val_list = c(val_list, value)
  }
  o = order(val_list, decreasing = TRUE)
  return(o)
}

cat("merged maf: All_cnvkit\nClinical data: total_clinical\nCNV data: cnvkit_0.4_onlyhnscc\nMatrix: mat_tot\nTop30 Matrix: mat30\nMatrix sorting func: sample_order")
