source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim/Figure/import_data_for_Oncoprint_for_H1516.R")


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
H15_inter[is.na(H15_inter$CTRTP1),]$CTRTP1 <- FALSE
H15_inter$final_arm5 <- ifelse((H15_inter$TRTP1 =="Arm 5 (Durvalumab +/- Tremelimumab)" |H15_inter$CTRTP1=="Arm 5 (Durvalumab +/- Tremelimumab)"), TRUE, FALSE)
H15_inter$arm5 <- ifelse(H15_inter$TRTP1 =="Arm 5 (Durvalumab +/- Tremelimumab)" ,TRUE,FALSE)

# Inspect the structure of your dataframe
str(H15_inter)



# H15_inter$final_arm5 == (H15_inter$ARM5 == " 1")
clinical_HPV_H15 <- subset(H15_inter, final_HPV == 'Positive')
# pdf("/data/project/TRIUMPH/tmp_graph_OStime1/HPV_ARM5.pdf")
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ arm5, data=clinical_HPV_H15)
custom_colors <-  c('#ED553B','#253F7F')

pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정3/Figure/Fig4F.pdf", height = 6, width = 8)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ arm5, data=clinical_HPV_H15)
# ggsurvplot(km_AG_fit, title = "ARM1",data =clinical_arm, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
F4 <- ggsurvplot(km_AG_fit, pval = TRUE, 
                 pval.size = 8,
                 data=clinical_HPV_H15,
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
                 legend.labs = c('Without immunotherapy','With immunotherapy'))  # To display time in years
F4$plot <- F4$plot + 
  theme(legend.title = element_text(size = 20),  # Adjust legend title size
        legend.text = element_text(size = 20))   # Adjust legend text size

F4
dev.off()
