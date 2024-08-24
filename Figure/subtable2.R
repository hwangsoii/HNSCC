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

dependent_os  <- "Surv(OS_TIME1, DEATH)"
# # dependent_dss <- "Surv(time, status_dss)"
# # dependent_crr <- "Surv(time, status_crr)"
# explanatory   <- c("TP53_total","EGFR_cnv","MYC_total", "CDKN2A", "NOTCH1_total")
# explanatory   <- c("TP53_total", "CDKN2A_total", "NOTCH1_total")
explanatory   <- c("TP53_total", "CDKN2A_total", "NOTCH1_total")
library(finalfit)
library(dplyr)
library(survival)
H15_inter %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  dplyr::rename("Overall survival" = label) %>%
  dplyr::rename(" " = levels) %>%
  dplyr::rename("  " = all) -> t1
t1
t1 %>%
  kableExtra::kbl(caption = "Recreating booktabs style table") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")
