source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/nanostring_import_data.R")

cnvkit_fil <- cnvkit_0.4_onlyhnscc[cnvkit_0.4_onlyhnscc[,2]%in%formed,]
cnvkit_fil <- cnvkit_fil[cnvkit_fil[,1]%in%gene_list,]
cnvkit_fil <- unique(cnvkit_fil)
cnvkit_fil$amp_data <- apply(cnvkit_fil, 1, function(x) fil_FC[x[1], substr(x[2], 1, 7)])

cnvkit_fil[cnvkit_fil$gene == "ERBB2",]
FGFR_amp <- cnvkit_fil[cnvkit_fil$gene == "FGFR1",]
EGFR_amp <- cnvkit_fil[cnvkit_fil$gene == "EGFR",]
PIK3CA_amp <- cnvkit_fil[cnvkit_fil$gene == "PIK3CA",]
EGFRvsPIK3CA <- rbind(cnvkit_fil[cnvkit_fil$gene == "EGFR",], cnvkit_fil[cnvkit_fil$gene == "PIK3CA",])




for_chi_EGFR <- for_chi("EGFR")[[1]]
chisq.test(for_chi_EGFR, correct = T)  ### 2.2e-16
p1 <- fisher.test(for_chi_EGFR) ## <2.2e-16 lowest value
t1 <- for_chi("EGFR")[[2]]

for_chi_ERBB2 <- for_chi("ERBB2")[[1]]
chisq.test(for_chi_ERBB2, correct = T) ###  < 2.2e-16 
p2 <- fisher.test(for_chi_ERBB2) ##2.323e-09
t2 <- for_chi("ERBB2")[[2]]

for_chi_FGFR1 <- for_chi("FGFR1")[[1]]
chisq.test(for_chi_FGFR1, correct = T) ###no Amplification patient over logFC  > 1  ###0.7922
p3 <- fisher.test(for_chi_FGFR1) ##1
t3 <- for_chi("FGFR1")[[2]]

for_chi_FGFR2 <- for_chi("FGFR2")[[1]]
chisq.test(for_chi_FGFR2, correct = T) ###no NA patient over logFC  > 1  ###< 2.2e-16
p4 <- fisher.test(for_chi_FGFR2) ##0.0001583
t4 <- for_chi("FGFR2")[[2]]

for_chi_FGFR3 <- for_chi("FGFR3")[[1]]
chisq.test(for_chi_FGFR3, correct = T) ###only two Amplification patient  ###2.47e-06
p5 <- fisher.test(for_chi_FGFR3) ##0.01116
t5 <- for_chi("FGFR3")[[2]]

for_chi_PIK3CA <- for_chi("PIK3CA")[[1]]
chisq.test(for_chi_PIK3CA, correct = T)  ###0.0002246
p6 <- fisher.test(for_chi_PIK3CA) ##0.00109
t6 <- for_chi("PIK3CA")[[2]]

EGFRvsPIK3CA_chi <- rbind(for_chi_EGFR[1,],for_chi_PIK3CA[1,])


###############################2Fold
for_bar <- EGFRvsPIK3CA[EGFRvsPIK3CA[,3]=="Amplification",]
for_bar$FC <- ifelse(for_bar$amp_data >=1, "over2Fold", "under2Fold")

fig9 <- ggplot(data = for_bar, aes(x=gene,y=..count..,fill=FC), stat = "count")+
  geom_bar()+theme_bw() +
  geom_signif(aes(x = 1,xend = 2, y=60, yend=60, annotation = fisher.test(EGFRvsPIK3CA_chi)$p.value),
              stat = "identity", inherit.aes=F, manual = T, tip_length = 1)+
  # scale_fill_grey(start = grey1, end = grey2)+
  geom_text(aes(label=..count..), stat = "count", position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")



pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig3B.pdf", height = 10, width = 15)
fig9
dev.off()

# ################################4Fold
# for_bar <- EGFRvsPIK3CA[EGFRvsPIK3CA[,3]=="Amplification",]
# for_bar$FC <- ifelse(for_bar$amp_data >=2, "over4Fold", "under4Fold")
# 
# FC_4 <- data.frame(c(18,53),c(24,3))
# 
# ggplot(data = for_bar, aes(x=gene,y=..count..,fill=FC), stat = "count")+
#   geom_bar()+
#   geom_signif(aes(x = 1,xend = 2, y=60, yend=60, annotation = fisher.test(FC_4)$p.value),
#               stat = "identity", inherit.aes=F, manual = T, tip_length = 1)+
#   scale_fill_grey(start = grey1, end = grey2)+
#   scale_color_grey(start = grey3, end = grey4)+
#   geom_text(aes(label=..count..), stat = "count", position = position_stack(vjust = 0.5))