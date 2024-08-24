source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim/Figure/import_data_for_Oncoprint_for_H1516.R")



#########add check



select_barcode <- as.vector(unlist(H15_inter$Tumor_Sample_Barcode))
new_cnv <- subset(cnvkit_0.4, sample %in% select_barcode)
H15_inter <- data.frame(lapply(H15_inter, function(x) unlist(x)))


#######################################################
library(survival)
library(survminer)
library(lubridate)
# library(ranger)
library(ggplot2)
library(dplyr)

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

H15_inter$TRTP1
(H15_inter$TRTP1 == "Arm 4 (Abemaciclib)") == (H15_inter$ARM4 == " 1")
((H15_inter$TRTP1 == "Arm 5 (Durvalumab +/- Tremelimumab)") |H15_inter$CTRTP1  =="Arm 5 (Durvalumab +/- Tremelimumab)") == (H15_inter$ARM5 == " 1")
H15_inter$CTRTP1

# H15_inter$final_arm1 <- ifelse((H15_inter$PIK3CA_total == TRUE|H15_inter$PIK3CB_total == TRUE|H15_inter$PIK3C2A_total == TRUE|H15_inter$PIK3C3_total == TRUE|H15_inter$PIK3CD_total == TRUE|H15_inter$PIK3R1_total == TRUE|H15_inter$PTEN_total == TRUE),TRUE,FALSE)
# H15_inter$final_arm3 <- ifelse((H15_inter$FGFR1_total == TRUE|H15_inter$FGFR2_total == TRUE|H15_inter$FGFR3_total == TRUE|H15_inter$FGFR4_total == TRUE),TRUE,FALSE)
# H15_inter$MAPK <- ifelse((H15_inter$KRAS_total == TRUE|H15_inter$HRAS_total == TRUE|H15_inter$RAF1_total == TRUE|H15_inter$BRAF_total == TRUE|H15_inter$RASSF1_total == TRUE),TRUE,FALSE)





#####
clinical_arm<- subset(H15_inter, TRTP1 == "Arm 4 (Abemaciclib)" )
clinical_arm_spe <- clinical_arm

(clinical_arm$CDKN2A == TRUE) &(clinical_arm$CDKN2A_cnv == TRUE)
summary(clinical_arm$CDKN2A)
summary(clinical_arm$CDKN2A_cnv)



# pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/ARM4_CDKN2Acnv.pdf")
pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2M.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ CDKN2A_cnv, data=clinical_arm)
# ggsurvplot(km_AG_fit, title = "ARM4",data =clinical_arm, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
ggsurvplot(km_AG_fit, pval = TRUE, 
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
           legend.labs = c('Without CDKN2A del','With CDKN2A del'))  # To display time in years
dev.off()
dependent_os  <- "Surv(OS_TIME1, DEATH)"

explanatory   <- c("final_HPV", "CDKN2A_cnv")
clinical_arm$final_HPV



# pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/ARM4_CDKN2A.pdf", height = 9, width = 6)
pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2N.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ CDKN2A, data=clinical_arm)
# ggsurvplot(km_AG_fit, title = "ARM4",data =clinical_arm, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
ggsurvplot(km_AG_fit, pval = TRUE, 
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
           legend.labs = c('Without CDKN2A SNV/Indel','With CDKN2A SNV/Indel'))  # To display time in years
dev.off()

pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2O.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ CCND1_cnv, data=clinical_arm)
# ggsurvplot(km_AG_fit, title = "ARM4",data =clinical_arm, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
ggsurvplot(km_AG_fit, pval = TRUE, 
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
           legend.labs = c('Without CCND1 cnv','With CCND1 cnv'))  # To display time in years
dev.off()
pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2P.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ genes_arm1, data=clinical_arm)
# ggsurvplot(km_AG_fit, title = "ARM4",data =clinical_arm, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
ggsurvplot(km_AG_fit, pval = TRUE, 
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
           legend.labs = c('Without PI3K pathway mut','With PI3K pathway mut'))  # To display time in years
dev.off()
pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2Q.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ genes_arm2, data=clinical_arm)
# ggsurvplot(km_AG_fit, title = "ARM4",data =clinical_arm, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
ggsurvplot(km_AG_fit, pval = TRUE, 
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
           legend.labs = c('Without EGFR pathway mut','With EGFR pathway mut'))  # To display time in years
dev.off()
pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2R.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ genes_arm7, data=clinical_arm)
# ggsurvplot(km_AG_fit, title = "ARM4",data =clinical_arm, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
ggsurvplot(km_AG_fit, pval = TRUE, 
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
           legend.labs = c('Without RAS-RAF pathway mut','With RAS-RAF pathway mut'))  # To display time in years
dev.off()
