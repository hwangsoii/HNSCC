library(maftools)
library(readxl)
library(ComplexHeatmap)
library(colorspace)
library(RColorBrewer)
source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/oncoprint_matrix_build.R")    ###original code
source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/msi_reader.R")

###maf data
fileNames=Sys.glob("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/all_220625/*vep.maf")
fileNames_vaf=Sys.glob("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/all_220625_vaf/*vep.vaf.maf")


###CNV data
cnvkit_0.4 <-  read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/CNV_filter_220626.txt", header = TRUE,  sep = "\t")
# cnvkit_0.4_onlyhnscc <- cnvkit_0.4[!(cnvkit_0.4$sample %in% nonhnsccfp), ]

library(stringr)
library(gplots)
library(ConsensusClusterPlus)
library(ggrepel)
library(ggplot2)
library(qvalue)
library(ggsignif)
library(survival)
library(survminer)
library(readxl)
library(dplyr)
library(RColorBrewer)

H15_16_data_1 <- as.data.frame(read_excel("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/HN15-16_20220816.xlsx", sheet = 1, col_names = T))
# new_clinical <- read.table("/data/project/TRIUMPH/raw_data/total_clinical_220816", header = TRUE,  sep = "\t")
# H15_16_data <- read.table("/data/project/TRIUMPH/raw_data/new_hnscc_clinical_combine_220816", header = TRUE,  sep = "\t")
# H15_16_data <- read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/total_clinical_221122.txt", header = TRUE,  sep = "\t")
H15_16_data <- read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/total_clinical_221127_p16.txt", header = TRUE,  sep = "\t")
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
write.table(merged_clinical2, file = "/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/final_clinical_merged_221128_p16.txt", quote = F, sep = "\t", row.names = F, na = "NA")
nonhnsccfp <- c("01-S005-FP", "21-S005-FP", "22-S011-FP", "22-S015-FP", "27-S024-FP", "28-S004-FP", "31-S013-FP", "46-S010-FP")
total_clinical <- data.frame(read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/final_clinical_merged_221128_p16.txt", header = T, sep = '\t'))
total_clinical_onlyhnscc <-  total_clinical[!(total_clinical$Tumor_Sample_Barcode %in% nonhnsccfp),]
total_clinical_onlyhnscc$final_arm <- ifelse(!is.na(total_clinical_onlyhnscc$TRTP1), substr(total_clinical_onlyhnscc$TRTP1, 1,5), substr(total_clinical_onlyhnscc$ARM, 1,5))
write.table(total_clinical_onlyhnscc, file = "/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/final_clinical_merged_221208_removal_nonhnscc_p16.txt", quote = F, sep = "\t", row.names = F, na = "NA")
###############################
H15_inter <- merged_clinical[merged_clinical$Patient_ID %in% H15_16_data_1$Patient_ID,]
H15_inter$Sample2 <- unlist(H15_inter$Sample2)
H15_inter <- merge(H15_inter, result_df, by= 'Sample2')

select_barcode <- as.vector(unlist(H15_inter$Tumor_Sample_Barcode))
new_cnv <- subset(cnvkit_0.4, sample %in% select_barcode)


H15_16_cnv <-  new_cnv

#######HPV final change arm4 221127
H15_inter$HPV_final[H15_inter$Patient_ID == '2018-072631'] <- "Negative"


#merge maf
# fileNames
# select_barcode
# substr(select_barcode, 1, 7)
H15_16_somatic<- paste0("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/all_220625/", substr(select_barcode, 1, 7), ".vep.maf")
H15_16_somatic_vaf<- paste0("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/maf/all_220625_vaf/", substr(select_barcode, 1, 7), ".vep.vaf.maf")
fileNames_H15_16 <- fileNames[ fileNames %in% H15_16_somatic]
fileNames_H15_16_vaf <- fileNames_vaf[ fileNames_vaf %in% H15_16_somatic_vaf]

All_cnvkit = merge_mafs(fileNames_H15_16, clinicalData = H15_inter, cnTable = H15_16_cnv)

getGeneSummary(All_cnvkit)


# Initialize an empty vector to store the extracted values
extracted_values <- character()

# Loop through each file path and extract the desired substring
for (file_path in fileNames_H15_16) {
  match <- regmatches(file_path, regexpr("[0-9]+-[A-Z0-9]+", file_path))
  extracted_values <- c(extracted_values, match)
}

# Print the extracted values
print(extracted_values)

# Add "-FP" to each extracted value
modified_names <- paste0(extracted_values, "-FP")
other_names <- getSampleSummary(All_cnvkit)[,1]
# Find which modified names are not in the other vector
setdiff(modified_names,getSampleSummary(All_cnvkit)$Tumor_Sample_Barcode)

All_cnvkit@data


#matrix form
mat_tot <- oncoprint_matrix_build(All_cnvkit)
# write.table(mat_tot, file = "H15_matrix.txt",sep = '\t', row.names = TRUE)

mat_tot <- mat_tot[,order(colnames(mat_tot))]
dim(mat_tot)
mat_tot = as.matrix(mat_tot)
mat_tot <- gsub(";Amp",'|Amp',mat_tot)
mat_tot <- gsub(";Del",'|Del',mat_tot)
mat_tot[grep(';',mat_tot)] <- paste0(mat_tot[grep(';',mat_tot)],";Multi_hit")
mat_tot <- gsub('\\|Amplification;Multi_hit', ';Multi_hit|Amplification', mat_tot)
###clinical data of patient in mat_tot
inter_clinical <- subset(H15_inter, H15_inter$Tumor_Sample_Barcode%in%colnames(mat_tot))
inter_clinical <- inter_clinical[order(unlist(inter_clinical$Tumor_Sample_Barcode)),]


###importing oncoplot top30 gene
setwd("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/")
oncoplot(maf = All_cnvkit, writeMatrix = T, top = 30);top30_genes <- rownames(read.table(file = "onco_matrix.txt", sep = "\t"))
file.remove('onco_matrix.txt')
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


mutation_type <- unique(unlist(strsplit(unlist(strsplit(unlist(as.list(mat_tot)),";")),'\\|')))
mutation_type <- mutation_type[mutation_type!=""]
mutation_color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#F8C8DC",
                    "#C6CE3E", "#5EDF82", "#00E4D0", "#00D7FF", "#C7BAFF", "#808B96")
barplot(rep(1,length(mutation_type)), col = mutation_color, names.arg = names(mutation_color), cex.names = 0.5)
names(mutation_color) <- mutation_type[c(3,1,8,2,4,5,7,9,6,11,12,13,15,14,10)]



height = 0.97
width = 0.95


###pathway
genes_arm1 <- c("PIK3CA", "PIK3CB", "PIK3C2A", "PIK3C3", "PIK3CD", "PIK3R1", "PIK3R2","PIK3R3",  "PTEN") ##PI3K
genes_arm2 <- c("EGFR", "ERBB2", "ERBB3","ERBB4", "MET", "PDGFRA") #EGFR
genes_arm3 <- c("FGFR1", "FGFR2","FGFR3", "FGFR4", "KIT", "IGF1R") ##FGFR
genes_arm4 <- c("CDKN2A", "CDKN2B", "CDKN2C", "CDKN1B", "CCND1", "CCND2", "CCND3","CCNE1", "RB1") ##Cell cycle
genes_arm5 <- c("AKT1","AKT2", "AKT3","TSC1","TSC2","MTOR","RICTOR","RPTOR","PPP2R1A","STK11") ##MTOR
genes_arm6 <- c("NOTCH1","NOTCH2", "FBXW7", "CREBBP", "EP300","KDM5A") #NOTCH
genes_arm7 <- c("KRAS","HRAS", "NRAS", "RIT1","ARAF","BRAF", "RAF1","RAC1","MAPK1","MAP2K1", "MAP2K2" ) #RASRAF
genes_arm8 <- c("TP53", "ATM", "CHEK2", "RPS6KA1") #TP53
genes_arm9 <- c("MYC", "MYCN", "MYCL1") #MYC
genes_arm10 <- c( "FAT1","FAT2", "FAT3","FAT4", "NF2")#HIPPO
genes_arm11 <- c("KEAP1", "CUL3", "NFE2L2") #Nrf2
genes_arm12 <- c("TGFBR2", "CTNNB1","APC") #wnt

gene_li <- c(genes_arm1, genes_arm2, genes_arm3, genes_arm4, genes_arm5, genes_arm6, genes_arm7, genes_arm8, genes_arm9, genes_arm10, genes_arm11, genes_arm12)





Pathways <- c('PI3Kinase', 'EGFR', 'FGFR', 'Cellcycle', 'MTOR', 'NOTCH', 'RASRAF', 'TP53','MYC', 'HIPPO', 'Nrf2', 'WNT')

H15_inter<- H15_inter[-c(1:4)]
data_H15 <- All_cnvkit@data

All_cnvkit = merge_mafs(fileNames_H15_16, clinicalData = H15_inter, cnTable = new_cnv)
All_cnvkit_nocnv = merge_mafs(fileNames_H15_16, clinicalData = H15_inter)
gene_li <- c(genes_arm1, genes_arm2, genes_arm3, genes_arm4, genes_arm5, genes_arm6, genes_arm7, genes_arm8, genes_arm9, genes_arm10, genes_arm11, genes_arm12)
for (gen in gene_li){
  print(gen)
  # print(gene)
  assign(paste0('somatic_', gen, '_barcode'), unlist(genesToBarcodes(maf = All_cnvkit_nocnv, genes = gen, justNames = TRUE)))
  H15_inter[get('gen')] <- H15_inter$Tumor_Sample_Barcode %in% get(paste0('somatic_', gen, '_barcode')) 
  assign(paste0('somatic_', gen, '_barcode_total'), unlist(genesToBarcodes(maf = All_cnvkit, genes = gen, justNames = TRUE)))
  H15_inter[paste0(gen, '_total')] <- H15_inter$Tumor_Sample_Barcode %in% get(paste0('somatic_', gen, '_barcode_total')) 
  tmp_df <- new_cnv %>% filter(gene == gen)
  assign(paste0('cnv_', gen, '_barcode'), tmp_df$sample)
  H15_inter[paste0(gen, '_cnv')] <- H15_inter$Tumor_Sample_Barcode %in% get(paste0('cnv_', gen, '_barcode'))
}




# Create an empty matrix to store the count for each sample and each pathway
pathway_counts <- matrix(0, nrow = nrow(H15_inter), ncol = 12)

for (i in c(1:12)){
  gene_arm <- paste0('genes_arm', i)
  # H15_inter[paste0('genes_arm', i)] <- rowSums(H15_inter[, paste0(get(gene_arm),'')] == TRUE) > 0
  H15_inter[paste0('genes_arm', i)] <- rowSums(H15_inter[, paste0(get(gene_arm),'_total')] == TRUE) > 0
  # Create a new column "n_mut_pathway" that counts the number of pathways with TRUE mutations
  # H15_inter[paste0('n_mut_', gene_arm)] <- rowSums(H15_inter[, paste0(get(gene_arm),'')] == TRUE)
  H15_inter[paste0('n_mut_', gene_arm)] <- rowSums(H15_inter[, paste0(get(gene_arm),'_total')] == TRUE)
  
  # Loop through each sample and count pathways with at least one mutation
  for (j in 1:nrow(H15_inter)){
    if (H15_inter[j, paste0('genes_arm', i)]){
      pathway_counts[j, i] <- 1
    }
  }
}

# Define the number of gene arms (pathways) and samples
num_gene_arms <- 12
num_samples <- nrow(H15_inter)  # Replace with the actual number of samples

H15_inter$Tumor_Sample_Barcode
# Create an empty matrix to store the mutation status
mutation_matrix <- matrix(0, nrow = num_gene_arms, ncol = num_samples)
# Set the column names of the mutation_matrix
colnames(mutation_matrix) <- H15_inter$Tumor_Sample_Barcode
# Assign row names to the mutation_matrix
rownames(mutation_matrix) <- paste0("gene_arm", 1:num_gene_arms) 
# Loop through gene arms

for (i in 1:num_gene_arms) {
  gene_arm_col <- paste0('genes_arm', i)
  
  # Loop through samples and populate the matrix
  for (j in 1:num_samples) {
    if (H15_inter[j, paste0('n_mut_', gene_arm_col)] > 0) {
      mutation_matrix[i, j] <- 1
      # mutation_matrix[i, j] <- 'mut'
      # mutation_matrix[i, j] <- Pathways[i]
      # mutation_matrix[i, j] <- "Mut"
    }
  }
}

