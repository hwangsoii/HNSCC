source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim/Figure/import_data_for_Oncoprint_for_H1516.R")


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

# pdf("/data/project/TRIUMPH/tmp_graph_OStime1/H15_TP53.pdf")
# km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ TP53_total, data=H15_inter)
# ggsurvplot(km_AG_fit, title = "Overall",data =H15_inter, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
# dev.off()
# 
# pdf("/data/project/TRIUMPH/tmp_graph_OStime1/H15_NOTCH1.pdf")
# km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ NOTCH1_total, data=H15_inter)
# ggsurvplot(km_AG_fit, title = "Overall",data =H15_inter, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
# dev.off()
# 
# pdf("/data/project/TRIUMPH/tmp_graph_OStime1/H15_CDKN2A.pdf")
# km_AG_fit <- survfit(Surv(OS_TIME1, DEATH) ~ CDKN2A_total, data=H15_inter)
# ggsurvplot(km_AG_fit, title = "Overall",data =H15_inter, pval = TRUE,conf.int = TRUE, risk.table = TRUE)
# dev.off()
# 
# ####COX
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ TP53+TRTP1+cnv+final_age , data=H15_inter)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ TP53+final_site+cnv , data=H15_inter)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ TP53+final_site+cnv , data=H15_inter)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ cnv , data=H15_inter)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ TP53 , data=H15_inter)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ NOTCH1 , data=H15_inter)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ MYC_total , data=H15_inter)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ CDKN2A_total , data=H15_inter)
res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ TP53_total + CDKN2A_total + NOTCH1_total , data=H15_inter)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ TP53_total + CDKN2A + NOTCH1_total , data=H15_inter)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ TP53_total + CDKN2A_cnv  + CDKN2A + NOTCH1_total , data=H15_inter)
# H15_inter$final_age <- as.numeric(H15_inter$final_age)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ TP53_total + CDKN2A + final_age+ TRTP1 + final_smoking  + NOTCH1_total , data=H15_inter)
# res.cox <- coxph(Surv(OS_TIME1, DEATH) ~ TP53_total + CDKN2A + EGFR_cnv  + NOTCH1_total + MYC_total , data=H15_inter)
summary(res.cox)
dependent_os  <- "Surv(OS_TIME1, DEATH)"
# # dependent_dss <- "Surv(time, status_dss)"
# # dependent_crr <- "Surv(time, status_crr)"
# explanatory   <- c("TP53_total","EGFR_cnv","MYC_total", "CDKN2A", "NOTCH1_total")
# explanatory   <- c("TP53_total", "CDKN2A_total", "NOTCH1_total")
explanatory   <- c("TP53_total", "CDKN2A_total", "NOTCH1_total")
explanatory   <- c("TP53_total", "CDKN2A_total", "NOTCH1_total", "final_HPV")
library(finalfit)
library(dplyr)
library(survival)
H15_inter %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  dplyr::rename("Overall survival" = label) %>%
  dplyr::rename(" " = levels) %>%
  dplyr::rename("  " = all) -> t1
t1
H15_inter_HPV_posi <- H15_inter[(H15_inter$final_HPV=="Positive")&!(is.na(H15_inter$final_HPV)),]
H15_inter_HPV_posi$final_HPV
H15_inter_HPV_posi %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  dplyr::rename("Overall survival" = label) %>%
  dplyr::rename(" " = levels) %>%
  dplyr::rename("  " = all) -> t2
t2
H15_inter_HPV_nega <- H15_inter[(H15_inter$final_HPV=="Negative")&!(is.na(H15_inter$final_HPV)),]
H15_inter_HPV_nega %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  dplyr::rename("Overall survival" = label) %>%
  dplyr::rename(" " = levels) %>%
  dplyr::rename("  " = all) -> t3
t3
# 
# knitr::kable(t1, align=c("l", "l", "r", "r", "r"))
# setwd("/data/project/TRIUMPH/tmp_graph_OStime1")
# save(t1, dependent_os, explanatory, 
#      file = here::here(".", "out.rda"))
#      
# pdf("/data/project/TRIUMPH/Figure_swkim/Overall_OR_plot.pdf")
# H15_inter %>%
#   hr_plot(dependent_os, explanatory)
# dev.off()

# 
# pdf("/data/project/TRIUMPH/Figure_swkim/proportinal_hazard.pdf")
# H15_inter %>% 
#   coxphmulti(dependent_os, explanatory) %>% 
#   cox.zph() %>% 
#   {zph_result <<- .} %>% 
#   plot(var=3)
# dev.off()
cox.zph(coxph(formula = Surv(OS_TIME1, DEATH) ~ TP53_total + CDKN2A_total + NOTCH1_total, data = H15_inter))
cox.zph(coxph(formula = Surv(OS_TIME1, DEATH) ~ CDKN2A_total, data = H15_inter))
cox.zph(coxph(formula = Surv(OS_TIME1, DEATH) ~ NOTCH1_total, data = H15_inter))
