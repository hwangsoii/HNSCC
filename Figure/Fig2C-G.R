source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_for_H1516.R")

library(survival)
library(survminer)
# if (!requireNamespace("openxlsx", quietly = TRUE)) {
#   install.packages("openxlsx")
# }
library(openxlsx)
library(finalfit)
library(dplyr)
library(survival)

select_barcode <- as.vector(unlist(H15_inter$Tumor_Sample_Barcode))
new_cnv <- subset(cnvkit_0.4, sample %in% select_barcode)
H15_inter <- data.frame(lapply(H15_inter, function(x) unlist(x)))


genes_arm1 <- c("PIK3CA", "PIK3CB", "PIK3C2A", "PIK3C3", "PIK3CD", "PIK3R1", "PIK3R2","PIK3R3",  "PTEN") ##PI3K
genes_arm2 <- c("EGFR", "ERBB2", "ERBB3","ERBB4", "MET", "PDGFRA") #EGFR
genes_arm3 <- c("FGFR1", "FGFR2","FGFR3", "FGFR4", "KIT", "IGF1R") ##FGFR
genes_arm4 <- c("CDKN2A", "CDKN2B", "CDKN2C", "CDKN1B", "CCND1", "CCND2", "CCND3","CCNE1", "RB1") ##Cell cycle
genes_arm5 <- c("AKT1","AKT2", "AKT3","TSC1","TSC2","MTOR","RICTOR","RPTOR","PPP2R1A","STK11") ##MTOR
genes_arm6 <- c("NOTCH1","NOTCH2", "FBXW7", "CREBBP", "EP300","KDM5A") #NOTCH
# genes_arm6 <- c("NOTCH2", "FBXW7", "CREBBP", "EP300","KDM5A") #NOTCH
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
H15_inter$Tumor_Sample_Barcode <- unlist(H15_inter$Tumor_Sample_Barcode.1)
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


#######################################################
library(survival)
library(survminer)
library(lubridate)
# library(ranger)
library(ggplot2)
library(dplyr)
# df_tmp <- H15_inter[H15_inter$final_arm == "Arm 1",]
# 
# 
# H15_inter <- as.data.frame(H15_inter)
# H15_inter$의뢰.기관명 <- unlist(H15_inter$의뢰.기관명)
# H15_inter$의뢰.기관명.1 <- unlist(H15_inter$의뢰.기관명.1)
# 
# typeof(H15_inter$final_death)
H15_inter$BOR_final2 <- ifelse(H15_inter$BOR_final=="PR", "Responser", "Nonresponder")
H15_inter$final_ostime <- as.numeric(H15_inter$final_ostime)
H15_inter$OS_TIME1 <- as.numeric(H15_inter$OS_TIME1)
H15_inter$PFS_TIME1 <- as.numeric(H15_inter$PFS_TIME1)
H15_inter$final_death <- as.numeric(H15_inter$final_death)
H15_inter$DEATH <- as.numeric(H15_inter$DEATH)
H15_inter$PFS_EVNT1 <- as.numeric(H15_inter$PFS_EVNT1)
H15_inter[is.na(H15_inter) ==TRUE] 
H15_inter[H15_inter ==""] <- NA
H15_inter[H15_inter =="NA"] <- NA



custom_colors <-  c('#ED553B','#253F7F')

##ARM1

##2C
clinical_arm<- subset(H15_inter, TRTP1 == "Arm 1 (BYL719)")

clinical_arm_spe <- clinical_arm
clinical_arm1_PIK3CA <- clinical_arm_spe[clinical_arm_spe$PIK3CA_total == TRUE,]
clinical_arm1_PIK3CA
clinical_arm1_PIK3CA <- clinical_arm1_PIK3CA[!(clinical_arm1_PIK3CA$PIK3CA_cnv == TRUE & clinical_arm1_PIK3CA$PIK3CA == TRUE),]

pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2C.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ PIK3CA, data=clinical_arm1_PIK3CA)
# ggsurvplot(km_AG_fit, title = "PIK3CA only ARM1",data =clinical_arm1_PIK3CA, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
C2 <- ggsurvplot(km_AG_fit, pval = TRUE, 
           pval.size = 8,
           data=clinical_arm1_PIK3CA,
           palette = custom_colors,
           # title = "Kaplan-Meier Survival Plot",
           xlab = "Time (months)",
           ylab = "Survival Probability",
           legend.title = "Groups",
           risk.table = F,
           risk.table.col = "strata",
           risk.table.y.text = FALSE,
           break.time.by = 10,
           linetype = 1,
           font.x = c(24),
           font.y = c(24),
           font.tickslab = c(22),
           # xscale = "y_m",
           font.legend=12,
           legend.labs = c('With PIK3CA CNV','With PIK3CA SNV/Indel'))  # To display time in years
C2
dev.off()

res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ final_HPV + PIK3CA , data=clinical_arm1_PIK3CA)
summary(res.cox)


dependent_os  <- "Surv(OS_TIME1, DEATH)"

explanatory   <- c("final_HPV", "PIK3CA")
clinical_arm1_PIK3CA %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  dplyr::rename("Overall survival" = label) %>%
  dplyr::rename(" " = levels) %>%
  dplyr::rename("  " = all) -> t1
t1

##2D
clinical_arm<- subset(H15_inter, TRTP1 == "Arm 1 (BYL719)")

clinical_arm_spe <- clinical_arm
tmp_mat <- clinical_arm_spe[clinical_arm_spe$PIK3CA_total == TRUE,]
# tmp_mat <- clinical_arm_spe[clinical_arm_spe$AKT1 == TRUE | clinical_arm_spe,]

con_mat <- tmp_mat[(tmp_mat$PIK3CB_total == TRUE|tmp_mat$PIK3C2A_total == TRUE|tmp_mat$PIK3C3_total == TRUE|tmp_mat$PIK3CD_total == TRUE|tmp_mat$PIK3R1_total == TRUE|tmp_mat$PIK3R2_total == TRUE|tmp_mat$PIK3R3_total == TRUE|tmp_mat$PTEN_total == TRUE),]
tmp_mat$more_PIK3 <- ifelse(tmp_mat$Tumor_Sample_Barcode %in% con_mat$Tumor_Sample_Barcode, TRUE, FALSE)

# Define custom colors for different groups or strata
custom_colors <- c("#990000", "#000099")
custom_colors <- c('#ED553B','#253F7F')


# Customize the Kaplan-Meier survival plot with your custom color palette
pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2D.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ more_PIK3, data=tmp_mat)
D2 <- ggsurvplot(km_AG_fit, pval = TRUE, 
           pval.size = 8,
           palette = c('#ED553B','#253F7F'),
           # title = "Kaplan-Meier Survival Plot",
           xlab = "Time (months)",
           ylab = "Survival Probability",
           legend.title = "Groups",
           risk.table = F,
           risk.table.col = "strata",
           risk.table.y.text = FALSE,
           break.time.by = 10,
           linetype = 1,
           font.x = c(24),
           font.y = c(24),
           font.tickslab = c(22),
           # xscale = "y_m",
           font.legend=12,
           legend.labs = c('More PI3K','No more PI3K'))  # To display time in years
D2
dev.off()
res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ final_HPV + more_PIK3 , data=tmp_mat)
summary(res.cox)


explanatory   <- c("final_HPV", "more_PIK3")
tmp_mat %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  dplyr::rename("Overall survival" = label) %>%
  dplyr::rename(" " = levels) %>%
  dplyr::rename("  " = all) -> t2
t2



clinical_arm <- H15_inter[H15_inter$TRTP1 == "Arm 1 (BYL719)",]
custom_colors <-  c('#ED553B','#253F7F')
pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2E.pdf", height = 9, width = 6)
clinical_arm[genes_arm5 == TRUE,]
clinical_arm$genes_arm5
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ genes_arm5, data=clinical_arm)
# ggsurvplot(km_AG_fit, title = "ARM1",data =clinical_arm, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
E2 <- ggsurvplot(km_AG_fit, pval = TRUE, 
           pval.size = 8,
           data=clinical_arm,
           palette = custom_colors,
           # title = "Kaplan-Meier Survival Plot",
           xlab = "Time (months)",
           ylab = "Survival Probability",
           legend.title = "Groups",
           risk.table = F,
           risk.table.col = "strata",
           risk.table.y.text = FALSE,
           break.time.by = 10,
           linetype = 1,
           font.x = c(24),
           font.y = c(24),
           font.tickslab = c(22),
           # xscale = "y_m",
           font.legend=12,
           legend.labs = c('Without MTOR pathway mut','With MTOR pathway mut'))  # To display time in years

E2
dev.off()
res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ final_HPV + genes_arm5 , data=clinical_arm)
summary(res.cox)

dependent_os  <- "Surv(OS_TIME1, DEATH)"

explanatory   <- c("final_HPV", "genes_arm5")
clinical_arm %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  dplyr::rename("Overall survival" = label) %>%
  dplyr::rename(" " = levels) %>%
  dplyr::rename("  " = all) -> t3
t3


pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2F.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ NOTCH1_total, data=clinical_arm)
F2 <- ggsurvplot(km_AG_fit, pval = TRUE, 
           pval.size = 8,
           data=clinical_arm,
           palette = custom_colors,
           # title = "Kaplan-Meier Survival Plot",
           xlab = "Time (months)",
           ylab = "Survival Probability",
           legend.title = "Groups",
           risk.table = F,
           risk.table.col = "strata",
           risk.table.y.text = FALSE,
           break.time.by = 10,
           linetype = 1,
           font.x = c(24),
           font.y = c(24),
           font.tickslab = c(22),
           # xscale = "y_m",
           font.legend=12,
           legend.labs = c('Without NOTCH1 mut','With NOTCH1 mut'))  # To display time in years

F2
dev.off()
res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ final_HPV + NOTCH1_total , data=clinical_arm)
summary(res.cox)

dependent_os  <- "Surv(OS_TIME1, DEATH)"

explanatory   <- c("final_HPV", "NOTCH1_total")
clinical_arm %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  dplyr::rename("Overall survival" = label) %>%
  dplyr::rename(" " = levels) %>%
  dplyr::rename("  " = all) -> t4
t4


pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2G.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ MYC_total, data=clinical_arm)
G2 <- ggsurvplot(km_AG_fit, pval = TRUE, 
           pval.size = 8,
           data=clinical_arm,
           palette = custom_colors,
           # title = "Kaplan-Meier Survival Plot",
           xlab = "Time (months)",
           ylab = "Survival Probability",
           legend.title = "Groups",
           risk.table = F,
           risk.table.col = "strata",
           risk.table.y.text = FALSE,
           break.time.by = 10,
           linetype = 1,
           font.x = c(24),
           font.y = c(24),
           font.tickslab = c(22),
           # xscale = "y_m",
           font.legend=12,
           legend.labs = c('Without MYC pathway mut','With MYC pathway mut'))  # To display time in years

G2
dev.off()
res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ final_HPV + MYC_total , data=clinical_arm)
summary(res.cox)

dependent_os  <- "Surv(OS_TIME1, DEATH)"

explanatory   <- c("final_HPV", "MYC_total")
explanatory   <- c("final_HPV", "MYC_total", "NOTCH1_total")
clinical_arm %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  dplyr::rename("Overall survival" = label) %>%
  dplyr::rename(" " = levels) %>%
  dplyr::rename("  " = all) -> t5
t5

library(survminer)
library(cowplot)

# ... Your code for creating the survival plots ...

# Extract just the plot component from each ggsurvplot object
B2_plot <- B2$plot
C2_plot <- C2$plot
D2_plot <- D2$plot
E2_plot <- E2$plot
F2_plot <- F2$plot
G2_plot <- G2$plot

B2$plot_env
# Now, combine these plots
pdf("/data/project/TRIUMPH/Figure_Combined_Horizontal.pdf")
combined_plot <- plot_grid(B2_plot, C2_plot, D2_plot, E2_plot, F2_plot, G2_plot, nrow = 1, align = 'h', axis = 'tb')
print(combined_plot)
dev.off()

###Fig2 추가부분
res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ final_HPV + PIK3CA , data=clinical_arm)
res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ ERBB2_cnv , data=clinical_arm)
res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ ERBB2_total , data=clinical_arm)
res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ MET_total , data=clinical_arm)
res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ genes_arm7 , data=clinical_arm)
summary(res.cox)

t1
t2
t3
t4
t5


