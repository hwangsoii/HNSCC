source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")
library(ggsci)
# All_finalcnvkit_nocnv= merge_mafs(fileNames_onlyhnscc, clinicalData = final_clinical)
# clinical <- final_clinical
# clinical$cnv <-ifelse(clinical$Tumor_Sample_Barcode %in% cnv_sample_barcode,"Yes", "No")
cnv_sample_barcode <- cnvkit_0.4_onlyhnscc$sample
total_clinical$cnv <-ifelse(total_clinical$Tumor_Sample_Barcode %in% cnv_sample_barcode,"Yes", "No")
clinical_practice <- total_clinical


vet2 <- mutate(clinical_practice, TP63_amp = ifelse((Tumor_Sample_Barcode %in% TP63_amp_barcode),"Yes", "No"))
vet_test <- vet2 %>% filter(final_site =='oropharynx'|final_site =='Hypopharynx'|final_site =='Lip or oral cavity'|final_site =='Larynx')
vet_test$final_site
##TP53 and cnv Fig3D


TP53_cnv <- table(clinical_practice$TP53, clinical_practice$cnv)
rownames(TP53_cnv) <- c('without TP53', 'with TP53 mutation')
colnames(TP53_cnv) <- c('without cnv', 'with cnv amplification')
# TP53_cnv %>%
#   kbl() %>%
#   kable_paper("hover", full_width=F)
fisher.test(TP53_cnv)
####Fig3D
g <- ggplot(vet_test, aes(x = TP53, fill = cnv)) + geom_bar(position = "dodge")+ theme_bw() +
  scale_fill_grey(start = 0.9, end = 0.4)+
  labs(x= 'TP53',y= 'CNV', fill = 'CNV') +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
g
######
g <- ggplot(vet_test, aes(x = TP53, fill = cnv)) + geom_bar(position = "fill")+ theme_bw() +
  scale_fill_grey(start = 0.9, end = 0.4)+
  labs(x= 'TP53',y= 'CNV', fill = 'CNV') +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
g
vet_test$final_HPV
################교수님 요청 HPV
vet_test_hpv <- vet_test[!is.na(vet_test$final_HPV),]
g <- ggplot(vet_test_hpv, aes(x = final_HPV, fill = TP53)) + geom_bar(position = "dodge")+ theme_bw() +
  scale_fill_grey(start = 0.9, end = 0.4)+
  labs(x= 'p16',y= 'TP53', fill = 'TP53') +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
g
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig5A_2.pdf", height = 10, width = 22)   ##########wilcox & Benjamini and Hochberg
g
dev.off()

g <- ggplot(vet_test_hpv, aes(x = final_HPV, fill = cnv)) + geom_bar(position = "dodge")+ theme_bw() +
  scale_fill_grey(start = 0.9, end = 0.4)+
  labs(x= 'p16',y= 'CNV', fill = 'CNV') +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
g
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig5B_2.pdf", height = 10, width = 22)   ##########wilcox & Benjamini and Hochberg
g
dev.off()
HPV_TP53 <- table(vet_test_hpv$final_HPV, vet_test_hpv$TP53)
fisher.test(HPV_TP53) #p-value = 2.096e-14
HPV_cnv <- table(vet_test_hpv$final_HPV, vet_test_hpv$cnv)
fisher.test(HPV_cnv) #p-value = 0.3968
######################################################
g <- ggplot(vet_test, aes(x = TP53, fill = cnv)) + geom_bar(position = "dodge")+ theme_bw()  +
  scale_fill_grey(start = 0.9, end = 0.4)+
  labs(x= 'TP53',y= 'CNV', fill = 'CNV') +
  
  # geom_signif(
  #   comparisons = list(c("FALSE", "TRUE")),
  #   map_signif_level = TRUE, y_position = 180,
  #   annotations="***"
  # )+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
g
g_tp53_cnv <- g
vet_test %>% 
  group_by(TP53) %>% 
  count(cnv) %>% 
  mutate(Freq = n) %>%
  ggplot() + 
  aes(x = TP53, Freq, fill = cnv) +
  scale_fill_grey(start = 0.9, end = 0.4)+
  labs(x= 'TP53',y= 'CNV', fill = 'CNV') +
  geom_bar(position = "dodge", stat='identity') +
  geom_signif(
    comparisons = list(c("FALSE", "TRUE")),
    map_signif_level = TRUE, y_position = 220,
    annotations="***"
  )+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

fisher.test(table(vet_test$TP53, vet_test$cnv))

vet_test %>% 
  group_by(TP53) %>% 
  count(cnv) %>% 
  mutate(Freq = n) %>%
  ggplot() + 
  aes(x = TP53, Freq, fill = cnv) +
  scale_fill_grey(start = 0.9, end = 0.4)+
  labs(x= 'TP53',y= 'CNV', fill = 'CNV') +
  geom_bar(position = "fill", stat='identity') +
  # geom_signif(
  #   comparisons = list(c("FALSE", "TRUE")),
  #   map_signif_level = TRUE, y_position = 220,
  #   annotations="***"
  # )+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


#############
vet_test$final_site == 'Maxillary sinus'
vet_test [is.na(vet_test$final_smoking),]$final_smoking <- "NA"
g <- ggplot(vet_test, aes(x = final_site, fill = final_smoking)) + geom_bar(position = "fill")+ theme_bw()  +
  # scale_fill_grey(start = 0.4, end = 0.9)+
  scale_fill_startrek()+
  labs(x= 'Primary Site',y= 'Smoking', fill = 'Smoking') +
  scale_fill_igv()+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
g_site_smoking <- g
g
#################################fig 3E
clinical_smoking <- clinical_practice[!(is.na(clinical_practice$final_smoking)),]
# clinical_final_smoking <- clinical_smoking %>% 
#   mutate(Smoking_status = ifelse((final_smoking == 'Current smoker')|(final_smoking == 'Former smoker'),'Smoking experience','Never smoking'))
clinical_final_smoking <- clinical_smoking %>% 
  mutate(Smoking_status = ifelse((final_smoking == 'Current smoker')|(final_smoking == 'Former smoker'),'Smoking experience','Never smoking'))
g <- ggplot(clinical_final_smoking, aes(x = Smoking_status, fill = cnv)) + geom_bar(position = "dodge") + ylab("The number of patients")+ theme_bw() +
  scale_x_discrete(limits=c('Smoking experience','Never smoking')) +
  # scale_fill_grey(start = 0.8, end = 0.4)+
  scale_fill_nejm(alpha = 0.75)+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
gsmoking <- g
g
Smoking_cnv2 <- table(clinical_final_smoking$Smoking_status, clinical_final_smoking$cnv)
fisher.test(Smoking_cnv2)
#########subfig3A
#####################################################################smoking
# clinical$final_smoking[is.na(clinical$final_smoking )] <- "NA"
# g <- ggplot(clinical, aes(x = final_smoking, fill = TP53)) + 
#   # scale_fill_grey(start = 0.9, end = 0.4)+
#   scale_fill_d3() +
#   labs(x= 'Smoking status',y= 'TP53', fill = 'TP53') +
#   geom_bar(position = "dodge") + ylab("The number of patients")+ theme_bw() +
#   scale_x_discrete(limits=c("Current smoker", "Former smoker", "Never smoker", "NA")) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) 
# g
g <- ggplot(clinical_final_smoking, aes(x = Smoking_status, fill = TP53)) + geom_bar(position = "dodge") + ylab("The number of patients")+ theme_bw() +
  scale_x_discrete(limits=c('Smoking experience','Never smoking')) +
  # scale_fill_grey(start = 0.8, end = 0.4)+
  scale_fill_nejm(alpha = 0.75)+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
g
######################3DEF

library(ggpubr)
ggarrange(g_tp53_cnv, g_site_smoking, gsmoking, ncol=3, nrow=1, legend = "top", align ="hv")
getwd()
setwd("D:/황신원/프로젝트/HNSCC/figure_1124/정리/figure_tmp/")
ggsave("def.pdf", width =6.5, height =5.1 )

