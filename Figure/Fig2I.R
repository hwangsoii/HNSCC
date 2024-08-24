source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_for_H1516.R")


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











#####
# clinical_arm<- subset(H15_inter, ARM2 == ' 1')
clinical_arm<- subset(H15_inter, TRTP1 == "Arm 2 (Poziotinib)")
clinical_arm_spe <- clinical_arm
clinical_arm$EGFR

summary(as.factor(clinical_arm$EGFR_total))
pdf("/data/project/TRIUMPH/Figure_swkim/Supplementary fig4E.pdf")
pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/ARM2_EGFR.pdf", height = 9, width = 6)
pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig2I.pdf", height = 9, width = 6)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ EGFR_total, data=clinical_arm)
# km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ PDGFRA_total, data=clinical_arm)
# ggsurvplot(km_AG_fit, title = "ARM2",data =clinical_arm, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
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
           legend.labs = c('Without EGFR mut','With EGFR mut'))  # To display time in years
dev.off()
dependent_os  <- "Surv(OS_TIME1, DEATH)"

explanatory   <- c("final_HPV", "EGFR_total")

clinical_arm %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  dplyr::rename("Overall survival" = label) %>%
  dplyr::rename(" " = levels) %>%
  dplyr::rename("  " = all) -> t1
t1
