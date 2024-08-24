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


pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정4/Figure/supplementary fig3C.pdf", height = 6, width = 8)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ EGFR_cnv, data=H15_inter)
subC3 <- ggsurvplot(km_AG_fit, pval = TRUE, 
                 pval.size = 8,
                 data=H15_inter,
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
                 legend.labs = c('Without EGFR cnv','With EGFR cnv'))  # To display time in years
subC3$plot <- subC3$plot + 
  theme(legend.title = element_text(size = 20),  # Adjust legend title size
        legend.text = element_text(size = 20))   # Adjust legend text size
subC3
dev.off()


km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ PIK3CA_cnv, data=H15_inter)
pdf("G:/내 드라이브/작업/HNSCC/수정2/추가수정4/Figure/supplementary fig3D.pdf", height = 6, width = 8)
km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ PIK3CA_cnv, data=H15_inter)
subD3 <- ggsurvplot(km_AG_fit, pval = TRUE, 
                 pval.size = 8,
                 data=H15_inter,
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
                 legend.labs = c('Without PIK3CA cnv','With PIK3CA cnv'))  # To display time in years
subD3$plot <- subD3$plot + 
  theme(legend.title = element_text(size = 20),  # Adjust legend title size
        legend.text = element_text(size = 20))   # Adjust legend text size
subD3
dev.off()
