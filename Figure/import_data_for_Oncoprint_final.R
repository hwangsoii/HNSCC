library(maftools)          ###P
library(readxl)            ###P
library(ComplexHeatmap)    ###P
library(colorspace)        ###N
library(RColorBrewer)      ###N
library(gridtext)          ###P
source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/oncoprint_matrix_build.R")    ###original code
source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/H15_inter_sample.R")
source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/msi_reader.R")


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

###CNV data
cnvkit_0.4 <-  read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/CNV_filter_220626.txt", header = TRUE,  sep = "\t")
cnvkit_0.4_onlyhnscc <- cnvkit_0.4[!(cnvkit_0.4$sample %in% nonhnsccfp), ]

###Clinical data
total_clinical <- data.frame(read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/final_clinical_merged_221208_removal_nonhnscc_p16.txt", header = T, sep = '\t'))
total_clinical <- merge(total_clinical, result_df, by= 'Sample2')
total_clinical$MS


tot_cnv <- cnvkit_0.4_onlyhnscc

#merge maf
All_cnvkit = merge_mafs(fileNames_onlyhnscc, clinicalData = total_clinical, cnTable = tot_cnv)
All_cnvkit_vaf = merge_mafs(fileNames_onlyhnscc_vaf, clinicalData = total_clinical, cnTable = tot_cnv)
All_cnvkit_nocnv = merge_mafs(fileNames_onlyhnscc, clinicalData = total_clinical)
b <- getSampleSummary(All_cnvkit)
a <- getGeneSummary(All_cnvkit)

# data = All_cnvkit@data
# write.table(data, file = "tot_data.txt", sep='\t')

#matrix form
mat_tot <- oncoprint_matrix_build(All_cnvkit)
write.table(mat_tot, file = "tot_matrix.txt", sep='\t', row.names = TRUE)
# mat_tot <- mat_tot[,order(colnames(mat_tot))]
dim(mat_tot)
mat_tot = as.matrix(mat_tot)
mat_tot <- gsub(";Amp",'|Amp',mat_tot)
mat_tot <- gsub(";Del",'|Del',mat_tot)
mat_tot[grep(';',mat_tot)] <- paste0(mat_tot[grep(';',mat_tot)],";Multi_hit")
mat_tot <- gsub('\\|Amplification;Multi_hit', ';Multi_hit|Amplification', mat_tot)
mat_tot <- gsub('\\|Deletion;Multi_hit', ';Multi_hit|Deletion', mat_tot)
mat_tot <- gsub("\\|Amp",';Amp',mat_tot)
mat_tot <- gsub("\\|Del",';Del',mat_tot)



# total_clinical$cnv <-ifelse(total_clinical$Tumor_Sample_Barcode %in% cnv_sample_barcode,"Yes", "No")
# clinical_trial_sample <- unlist(H15_inter$Tumor_Sample_Barcode.1)
# total_clinical$clinical_trial <- ifelse(total_clinical$Tumor_Sample_Barcode %in% clinical_trial_sample, "YES", "NO")

# top30_genes <- c(top30_genes, 'CCND1')
# for (gene in top30_genes){
#   assign(paste0('somatic_', gene, '_barcode'), unlist(genesToBarcodes(maf = All_cnvkit, genes = gene, justNames = TRUE)))
#   total_clinical[paste0(get('gene'),'_total')] <- total_clinical$Tumor_Sample_Barcode %in% get(paste0('somatic_', gene, '_barcode')) 
# }
# for (gen in top30_genes){
#   assign(paste0('somatic_', gen, '_barcode'), unlist(genesToBarcodes(maf = All_cnvkit_nocnv, genes = gen, justNames = TRUE)))
#   total_clinical[get('gen')] <- total_clinical$Tumor_Sample_Barcode %in% get(paste0('somatic_', gen, '_barcode')) 
#   tmp_df <- as.data.frame(tot_cnv) %>% filter(gene == gen)
#   assign(paste0('cnv_', gen, '_barcode'), tmp_df$sample)
#   total_clinical[paste0(gen, '_cnv')] <- total_clinical$Tumor_Sample_Barcode %in% get(paste0('cnv_', gen, '_barcode'))
# }

# for (gene in top20_genes){
#   print(gene)
#   
#   count_true <- sum(total_clinical[,paste0(gene,'_total')] == TRUE)
#   print(count_true)
#   print(round(count_true/419,2))
# }

###clinical data of patient in mat_tot
inter_clinical <- subset(total_clinical, total_clinical$Tumor_Sample_Barcode%in%colnames(mat_tot));inter_clinical <- inter_clinical[order(inter_clinical$Tumor_Sample_Barcode),]


###importing oncoplot top30 gene
setwd("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/")
pdf("temp.pdf")
oncoplot(maf = All_cnvkit, writeMatrix = T, top = 30);top30_genes <- rownames(read.table(file = "onco_matrix.txt", sep = "\t"))
# top30_genes <- as.data.frame(getGeneSummary(All_cnvkit)[1:30])[,1]
dev.off()
file.remove('onco_matrix.txt')
file.remove('temp.pdf')
mat30 = mat_tot[top30_genes[1:20], ]
top20_genes <- top30_genes[1:20]



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

cat("merged maf: All_cnvkit\nClinical data: total_clinical\nCNV data: cnvkit_0.4_onlyhnscc\nMatrix: mat_tot\nTop20 Matrix: mat30\nMatrix sorting func: sample_order")

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

# H15_inter<- H15_inter[-c(1:4)]
All_cnvkit = merge_mafs(fileNames_onlyhnscc, clinicalData = total_clinical, cnTable = tot_cnv)
All_cnvkit_nocnv = merge_mafs(fileNames_onlyhnscc, clinicalData = total_clinical)
gene_li <- c(genes_arm1, genes_arm2, genes_arm3, genes_arm4, genes_arm5, genes_arm6, genes_arm7, genes_arm8, genes_arm9, genes_arm10, genes_arm11, genes_arm12)
for (gen in gene_li){
  print(gen)
  # print(gene)
  assign(paste0('somatic_', gen, '_barcode'), unlist(genesToBarcodes(maf = All_cnvkit_nocnv, genes = gen, justNames = TRUE)))
  total_clinical[get('gen')] <- total_clinical$Tumor_Sample_Barcode %in% get(paste0('somatic_', gen, '_barcode')) 
  assign(paste0('somatic_', gen, '_barcode_total'), unlist(genesToBarcodes(maf = All_cnvkit, genes = gen, justNames = TRUE)))
  total_clinical[paste0(gen, '_total')] <- total_clinical$Tumor_Sample_Barcode %in% get(paste0('somatic_', gen, '_barcode_total')) 
  tmp_df <- cnvkit_0.4_onlyhnscc %>% filter(gene == gen)
  assign(paste0('cnv_', gen, '_barcode'), tmp_df$sample)
  total_clinical[paste0(gen, '_cnv')] <- total_clinical$Tumor_Sample_Barcode %in% get(paste0('cnv_', gen, '_barcode'))
}

# summary(total_clinical$TP53_total)
# total_clinical$Tumor_Sample_Barcode

# Create an empty matrix to store the count for each sample and each pathway
pathway_counts <- matrix(0, nrow = nrow(total_clinical), ncol = 12)

for (i in c(1:12)){
  gene_arm <- paste0('genes_arm', i)
  # total_clinical[paste0('genes_arm', i)] <- rowSums(total_clinical[, paste0(get(gene_arm),'')] == TRUE) > 0
  total_clinical[paste0('genes_arm', i)] <- rowSums(total_clinical[, paste0(get(gene_arm),'_total')] == TRUE) > 0
  # Create a new column "n_mut_pathway" that counts the number of pathways with TRUE mutations
  # total_clinical[paste0('n_mut_', gene_arm)] <- rowSums(total_clinical[, paste0(get(gene_arm),'')] == TRUE)
  total_clinical[paste0('n_mut_', gene_arm)] <- rowSums(total_clinical[, paste0(get(gene_arm),'_total')] == TRUE)
  
  # Loop through each sample and count pathways with at least one mutation
  for (j in 1:nrow(total_clinical)){
    if (total_clinical[j, paste0('genes_arm', i)]){
      pathway_counts[j, i] <- 1
    }
  }
}

# Define the number of gene arms (pathways) and samples
num_gene_arms <- 12
num_samples <- nrow(total_clinical)  # Replace with the actual number of samples

total_clinical$Tumor_Sample_Barcode
# Create an empty matrix to store the mutation status
mutation_matrix <- matrix(0, nrow = num_gene_arms, ncol = num_samples)
# Set the column names of the mutation_matrix
colnames(mutation_matrix) <- total_clinical$Tumor_Sample_Barcode
# Assign row names to the mutation_matrix
rownames(mutation_matrix) <- paste0("gene_arm", 1:num_gene_arms) 
# Loop through gene arms

for (i in 1:num_gene_arms) {
  gene_arm_col <- paste0('genes_arm', i)
  
  # Loop through samples and populate the matrix
  for (j in 1:num_samples) {
    if (total_clinical[j, paste0('n_mut_', gene_arm_col)] > 0) {
      mutation_matrix[i, j] <- 1
      # mutation_matrix[i, j] <- 'mut'
      # mutation_matrix[i, j] <- Pathways[i]
      # mutation_matrix[i, j] <- "Mut"
    }
  }
}

# Assuming total_clinical is your dataframe
write.csv(total_clinical, file = "/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/total_clinical.csv", row.names = FALSE)


