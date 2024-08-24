library(maftools)
library(readxl)
# BiocManager::install("ComplexHeatmap", force = T)
library(ComplexHeatmap)
library(colorspace)
library(RColorBrewer)
library(stringr)
library(gplots)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ConsensusClusterPlus", force = T)
library(ConsensusClusterPlus)
library(ggrepel)
library(ggplot2)
# BiocManager::install("qvalue", force = T)
library(qvalue)
library(ggsignif)
library(survival)
library(survminer)
library(readxl)
library(dplyr)
library(RColorBrewer)
# source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/oncoprint_matrix_build.R")
# source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/randomForest_source.r")


####read file
fileNames=Sys.glob("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/all_220625/*vep.maf")
fileNames_vaf=Sys.glob("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/all_220625_vaf/*vep.vaf.maf")

###CNV data
cnvkit_0.4 <-  read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/CNV_filter_220626.txt", header = TRUE,  sep = "\t")



H15_16_data_1 <- as.data.frame(read_excel("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/HN15-16_20220816.xlsx", sheet = 1, col_names = T))
# new_clinical <- read.table("/data/project/TRIUMPH/raw_data/total_clinical_220816", header = TRUE,  sep = "\t")
# H15_16_data <- read.table("/data/project/TRIUMPH/raw_data/new_hnscc_clinical_combine_220816", header = TRUE,  sep = "\t")
H15_16_data <- read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/total_clinical_221122.txt", header = TRUE,  sep = "\t")
# H15_16_data <- read.table("/data/project/TRIUMPH/raw_data/total_clinical_221127_p16.txt", header = TRUE,  sep = "\t")
new_clinical <- as.data.frame(read_excel("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/clinicaldata-2022-0816.xlsx",
                                         sheet = "Sheet1",
                                         col_types = c("numeric", "text", "text", "text", "text", "text", 
                                                       "text", "text", "text", "date", "text", 
                                                       "text", "text", "text", "text", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "date", "text", "text", "text", "text", 
                                                       "text", "text", "text", "date", "text", 
                                                       "date", "text", "text", "date", "date", 
                                                       "numeric", "text")))

new_clinical$HPV
new_clinical[7]
H15_16_data['AGE1']
merged_clinical <- apply(new_clinical, 1, function(x) c(x,H15_16_data[which(H15_16_data[,1]%in%x[7]),]))
merged_clinical <- data.frame(do.call(rbind, merged_clinical))
# merged_clinical <- unlist(merged_clinical)
merged_clinical[merged_clinical=="character(0)"] <- "NA"
merged_clinical[merged_clinical=="integer(0)"] <- "NA"
merged_clinical[merged_clinical=="numeric(0)"] <- "NA"
merged_clinical[is.na(merged_clinical)] <- "NA"
merged_clinical$final_HPV <- unlist(ifelse(merged_clinical$HPV_final!="NA", merged_clinical$HPV_final, merged_clinical$HPV))
merged_clinical$final_HPV[merged_clinical$Patient_ID == '2018-072631'] <- "Negative"
merged_clinical$final_site <- unlist(ifelse(merged_clinical$Primary_Site_final!="NA",merged_clinical$Primary_Site_final, merged_clinical$Primary_Site))
merged_clinical$final_smoking <- unlist(ifelse(merged_clinical$Smoking_final !="NA",merged_clinical$Smoking_final, merged_clinical$Smoking))

merged_clinical$final_smoking
merged_clinical$final_smoking[merged_clinical$final_smoking == ' 0'] <- "Never smoker"
merged_clinical$final_smoking[merged_clinical$final_smoking == ' 1'] <- "Former smoker"
merged_clinical$final_smoking[merged_clinical$final_smoking == ' 2'] <- "Current smoker"
merged_clinical$final_smoking[merged_clinical$final_smoking == ' 3'] <- "NA"
merged_clinical$final_age <- unlist(ifelse(merged_clinical$AGE1 != "NA", merged_clinical$AGE1, merged_clinical$Age_diagnosis))
#as.numeric(unlist(merged_clinical$days_to_last_fu))
merged_clinical$days_to_last_fu <- as.numeric((merged_clinical$days_to_last_fu))
merged_clinical$days_to_last_fu <- unlist(ifelse(merged_clinical$days_to_last_fu !='NA', merged_clinical$days_to_last_fu/365.25*12, "NA"))
merged_clinical$final_ostime <- unlist(ifelse(merged_clinical$OS_TIME1!="NA", merged_clinical$OS_TIME1, merged_clinical$days_to_last_fu))
merged_clinical$final_death <- unlist(ifelse(merged_clinical$DEATH!="NA", merged_clinical$DEATH, merged_clinical$Death))

merged_clinical$final_arm <- ifelse(!is.na(merged_clinical$TRTP1), substr(merged_clinical$TRTP1, 1,5), substr(merged_clinical$ARM, 1,5))
merged_clinical$Tumor_Sample_Barcode <- unlist(merged_clinical$Tumor_Sample_Barcode)
merged_clinical$Tumor_Sample_Barcode <- paste0(merged_clinical$Tumor_Sample_Barcode,'-FP')




# 221208 making file
merged_clinical$HPV_final[merged_clinical$Patient_ID == '2018-072631'] <- "Negative"
merged_clinical$final_HPV[merged_clinical$Patient_ID == '2018-072631'] <- "Negative"
merged_clinical2 <- apply(merged_clinical, 2, as.character)
merged_clinical2 <- merged_clinical2[,-c(1,2,3,4)]
write.table(merged_clinical2, file = "/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/final_clinical_merged_221128.txt", quote = F, sep = "\t", row.names = F, na = "NA")
nonhnsccfp <- c("01-S005-FP", "21-S005-FP", "22-S011-FP", "22-S015-FP", "27-S024-FP", "28-S004-FP", "31-S013-FP", "46-S010-FP")
total_clinical <- data.frame(read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/final_clinical_merged_221128.txt", header = T, sep = '\t'))
total_clinical_onlyhnscc <-  total_clinical[!(total_clinical$Tumor_Sample_Barcode %in% nonhnsccfp),]
total_clinical_onlyhnscc$final_arm <- ifelse(!is.na(total_clinical_onlyhnscc$TRTP1), substr(total_clinical_onlyhnscc$TRTP1, 1,5), substr(total_clinical_onlyhnscc$ARM, 1,5))
###############################
H15_inter <- merged_clinical[merged_clinical$Patient_ID %in% H15_16_data_1$Patient_ID,]

