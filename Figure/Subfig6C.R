source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim/Figure/import_data_for_Oncoprint_for_H1516.R")
#########add check


#######################################################
library(survival)
library(survminer)
library(lubridate)
# library(ranger)
library(ggplot2)
library(dplyr)
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



clinical_tp53 <- subset(H15_inter, TP53_total == TRUE)
clinical_PIK3CA <- subset(H15_inter, PIK3CA_cnv == TRUE)
clinical_PIK3CA <- subset(H15_inter, PIK3CA_total == TRUE)
clinical_PIK3CA <- subset(H15_inter, PIK3CA == TRUE)
clinical_CCND1cnv <- subset(H15_inter, CCND1_cnv == TRUE)
clinical_TP53_CCND1cnv <- subset(clinical_tp53, CCND1_cnv == TRUE)
clinical_not_tp53 <- subset(H15_inter, TP53 == FALSE)
clinical_HPV_H15 <- subset(H15_inter, final_HPV == 'Positive')
clinical_HPV_H15_TP53nega <- subset(H15_inter, final_HPV == 'Positive') %>% filter(TP53 ==FALSE)
clinical_tp53_arm <- subset(clinical_tp53, TRTP1 == "Arm 5 (Durvalumab +/- Tremelimumab)")
clinical_arm<- subset(H15_inter, TRTP1 == "Arm 5 (Durvalumab +/- Tremelimumab)")
clinical_tp53$TRTP1
clinical_arm[clinical_arm$CCND1_cnv ==TRUE,]$OS_TIME1




pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Subfig6C.pdf", height = 9, width = 9)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ TP53_total, data=clinical_arm)
custom_colors <-  c('#ED553B','#253F7F')
# ggsurvplot(km_AG_fit, title = "ARM1",data =clinical_arm, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
C6 <- ggsurvplot(km_AG_fit, pval = TRUE, 
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
                 font.legend=20,
                 legend.labs = c('Without TP53 mut','With TP53 mut')) 
# Modify the plot using ggplot2's theme function
C6$plot <- C6$plot + 
theme(legend.title = element_text(size = 20),  # Adjust legend title size
      legend.text = element_text(size = 20))   # Adjust legend text size

C6
dev.off()
