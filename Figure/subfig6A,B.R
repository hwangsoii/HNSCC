library(readxl)    ###P
# library(circlize)  ###N
library(ComplexHeatmap)  ###P
library('ggrepel')
library(dplyr)
library(tidyverse)
grey1 <- 0.9
grey2 <- 0.4
grey3 <- 0.5
grey4 <- 0.0
source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")
#################################################################################
#####################################reboot data#################################
#################################################################################
##############reboot_data = nsolver reboot normalization data (not reboot2)
Norm_Flag = c("16-S003","27-S018","31-S009","27-S026","45-S003","27-S035","27-S047",
              "15-S003","25-S003","26-S003","45-S004","25-S005","46-S011","08-S008",
              "25-S006","14-S007","26-S004","28-S007","25-S026","08-S014","24-S038",
              "25-S031","27-S068","34-S005","27-S084","27-S085","25-S039","25-S040",
              "45-S001","20-S011","15-S049","47-S011","45-S001","15-S050","47-S011")
# "45-S001","15-S050","47-S011"
nonhnscc <- c("01-S005", "21-S005", "22-S011", "22-S015", "27-S024", "28-S004", "31-S013", "46-S010")
reboot_data <- read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/reboot.txt", header = T, sep = "\t")
reboot_data.df <- as.data.frame(reboot_data)
reboot_coln <- colnames(reboot_data.df)
reboot_coln.modi1 <- gsub("_..\\.RCC","",reboot_coln)
reboot_coln.modi2 <- gsub("\\.+0","-S0",reboot_coln.modi1)
reboot_coln.modi3 <- str_sub(reboot_coln.modi2,-7,-1)
reboot_coln.modi4 <- c(reboot_coln.modi2[1:5], gsub("\\.","-", reboot_coln.modi3[6:length(reboot_coln.modi2)]))
colnames(reboot_data.df) <- toupper(reboot_coln.modi4)
reboot_data.df[is.na(reboot_data.df)] <- ""
reboot_data.df <- reboot_data.df[,reboot_data.df[1,]!="FLAG"]
reboot_data.df <- reboot_data.df[reboot_data.df[,2]=="Endogenous",]
reboot_rown <- reboot_data.df[,1]
reboot_data.df <- reboot_data.df[,-c(1,2)]
reboot_data.df <- reboot_data.df[!(colnames(reboot_data.df)%in%Norm_Flag)]
reboot_data.df <- reboot_data.df[!(colnames(reboot_data.df)%in%nonhnscc)]
reboot_coln <- colnames(reboot_data.df)
reboot_data.df <- reboot_data.df[,order(names(reboot_data.df))]
rownames(reboot_data.df) <- reboot_rown
reboot_data.m <- matrix(as.numeric(as.matrix(reboot_data.df)), ncol = ncol(reboot_data.df))
rownames(reboot_data.m) <- reboot_rown
colnames(reboot_data.m) <- sort(reboot_coln)
log2reboot_data <- log2(reboot_data.m)
reboot_data.z <- t(scale(t(log2reboot_data)))
samples <- colnames(reboot_data.z)
cnvkit_0.4 <-  read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/CNV_filter_220626.txt", header = TRUE,  sep = "\t")
clinical <- read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/final_clinical_merged_221208_removal_nonhnscc_p16.txt", header = T, sep = '\t')
genes_arm1 = c("PIK3CA", "PIK3CB", "PIK3C2A", "PIK3C3", "PIK3CD", "PIK3R1", "PTEN")
genes_arm2 = c("EGFR", "ERBB2", "ERBB3","ERBB4")
genes_arm3= c("FGFR1", "FGFR2","FGFR3", "FGFR4")
genes_arm4 = c("CDKN2A", "CCND1", "CCND2", "CCND3", "CDKN2C")
genes_arm5 = c("HRAS", "JAK2", "BRAF", "KRAS", "CCNE1")
cnvkit_0.4_onlyhnscc <- cnvkit_0.4[!(cnvkit_0.4$sample %in% nonhnscc), ]
All_finalcnvkit_nocnv= merge_mafs(fileNames_onlyhnscc, clinicalData = total_clinical)
All_finalcnvkit= merge_mafs(fileNames_onlyhnscc, clinicalData = total_clinical, cnTable = cnvkit_0.4_onlyhnscc)
gene_li <- c(genes_arm1, genes_arm2, genes_arm3, genes_arm4, genes_arm5, 'TP53', 'FAT1', 'KMT2D')

as.data.frame(tot_cnv)
tot_cnv %>% filter(gene == gen)

for (gen in gene_li){
  print(gen)
  # print(gene)
  assign(paste0('somatic_', gen, '_barcode'), unlist(genesToBarcodes(maf = All_finalcnvkit_nocnv, genes = gen, justNames = TRUE)))
  clinical[get('gen')] <- clinical$Tumor_Sample_Barcode %in% get(paste0('somatic_', gen, '_barcode')) 
  assign(paste0('somatic_', gen, '_barcode_total'), unlist(genesToBarcodes(maf = All_finalcnvkit, genes = gen, justNames = TRUE)))
  clinical[paste0(gen, '_total')] <- clinical$Tumor_Sample_Barcode %in% get(paste0('somatic_', gen, '_barcode_total')) 
  tmp_df <- tot_cnv %>% filter(gene == gen)
  assign(paste0('cnv_', gen, '_barcode'), tmp_df$sample)
  clinical[paste0(gen, '_cnv')] <- clinical$Tumor_Sample_Barcode %in% get(paste0('cnv_', gen, '_barcode'))
}


cnvkit_0.4_onlyhnscc <- cnvkit_0.4[!(cnvkit_0.4$sample %in% nonhnscc), ]
patients <- clinical[apply(clinical, 1, function(x) x["Sample2"]%in%samples),]

with_clinical <- reboot_data.m[,patients$Sample2]
log2with_clinical <- log2(with_clinical)
with_clinical.z <- t(scale(t(log2with_clinical)))

# wssplot(with_clinical.z)
km <- kmeans(with_clinical.z, centers=2)
with_clinical.z.km <- cbind(with_clinical.z, km$cluster)
o <- order(with_clinical.z.km[,length(with_clinical.z.km[1,])])
with_clinical.z.km.sort <- with_clinical.z.km[o, 1:(length(with_clinical.z.km[1,])-1)]

with_clinical.z.km.sort.modi <- ifelse(with_clinical.z.km.sort > 4, 4, with_clinical.z.km.sort)
with_clinical.z.km.sort.modi <- ifelse(with_clinical.z.km.sort.modi < -4, -4, with_clinical.z.km.sort.modi)


### RNA expression scaling


normal_maxmin <- function(x){
  return((x-min(x))/(max(x)-min(x)))
}

patients %>% group_by(final_site) %>% count(final_site)
want_site <- patients %>%
  filter(final_site =='oropharynx'|final_site =='Hypopharynx'|final_site =='Lip or oral cavity'|final_site =='Larynx')
want_site_clinical <- reboot_data.m[,want_site$Sample2]

with_clinical.maxmin <- t(apply(want_site_clinical, 1,normal_maxmin))

###matrix build

Immune_Profile <- read_excel("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/LBL-10043-08_nCounter_PanCancer_Immune_Profiling_Panel_Gene_List.xlsx",
                             sheet = "Annotations", col_names = T, skip = 1)
Immune_Profile.df <- as.data.frame(Immune_Profile)
rownames(Immune_Profile.df) <- Immune_Profile.df[,"Gene Name"]

Immune_Profile.df_used <- Immune_Profile.df[rownames(with_clinical),]
Immune_Profile.df_used <- Immune_Profile.df_used[!is.na(Immune_Profile.df_used$`Gene Name`),]

Category <- unique(unlist(strsplit(Immune_Profile.df_used$`Immune Response Category`,", ")))
Category <- Category[!is.na(Category)]

GenesByCategory = data.frame()

for(i in Category){
  assign(i , Immune_Profile.df_used$`Gene Name`[grep(pattern = i, Immune_Profile.df_used$`Immune Response Category`)])
  Names = Immune_Profile.df_used$`Gene Name`[grep(pattern = i, Immune_Profile.df_used$`Immune Response Category`)]
  GenesByCategory =  rbind(GenesByCategory, paste(Names, collapse = ","))
}

colnames(GenesByCategory) <- "Genes"
rownames(GenesByCategory) <- Category

###plot build

geom.mean_m <- numeric()

for(i in Category){
  genes = unlist(strsplit(GenesByCategory[i,], ","))
  if (length(genes) > 1){
    geom.mean = apply(with_clinical.maxmin[genes,], 2, function(x) exp(mean(log(x))))
  }
  else{
    geom.mean = with_clinical.maxmin[genes,]
  }
  geom.mean_m <- rbind(geom.mean_m, geom.mean)
}
geom.mean_df <- as.data.frame(geom.mean_m)
rownames(geom.mean_df) <- Category


merged_df <- as.data.frame(t(rbind(geom.mean_m,want_site$final_HPV, want_site$final_site)))

colnames(merged_df) <- c(rownames(geom.mean_m),"HPV","Site")

library(ggpubr)



T_cell_function = c("CD2","CD27","CD274","CD38","CD3E","CD3G","CD80","CD86","CD8A","CTLA4","CXCL10","CXCL9","CXCR5","IDO1","IFNG","IL18","IRF1","LAG3","LCK","TIGIT")
T_cell_P = c("CD3E","CD3G","CD8A","IFNG","IL18","LAG3")
T_cell_N = c("TIGIT","CTLA4","CD274")
Cytotoxicity_nonHLA = c("GZMB","GZMK")
Cytotoxicity = c("GZMB","GZMK", "HLA-B", "HLA-C")
rownames(with_clinical)
Chemokines = Chemokines



TP_geom.mean = apply(with_clinical.maxmin[T_cell_P,], 2, function(x) exp(mean(log(x))))
shapiro.test(TP_geom.mean)
TN_geom.mean = apply(with_clinical.maxmin[T_cell_N,], 2, function(x) exp(mean(log(x))))
shapiro.test(TN_geom.mean)
Cyto_geom.mean = apply(with_clinical.maxmin[Cytotoxicity_nonHLA,], 2, function(x) exp(mean(log(x))))
Cyto_geom.mean = apply(with_clinical.maxmin[Cytotoxicity,], 2, function(x) exp(mean(log(x))))
shapiro.test(Cyto_geom.mean)
Chemo_geom.mean = apply(with_clinical.maxmin[Chemokines,], 2, function(x) exp(mean(log(x))))
shapiro.test(Chemo_geom.mean)
Tfun_geom.mean = apply(with_clinical.maxmin[T_cell_function,], 2, function(x) exp(mean(log(x))))
shapiro.test(Tfun_geom.mean)

Select_merged_Site <- as.data.frame(t(rbind(TP_geom.mean,TN_geom.mean,Cyto_geom.mean,Chemo_geom.mean,Tfun_geom.mean,
                                            want_site$final_site)))
Select_merged_HPV <- as.data.frame(t(rbind(Cyto_geom.mean,Chemo_geom.mean,Tfun_geom.mean, want_site$final_HPV)))
Select_merged_TP53_forgraph <- as.data.frame(t(rbind(Cyto_geom.mean,Chemo_geom.mean,Tfun_geom.mean,want_site$final_HPV,want_site$TP53)))
Select_merged_TP53_forgraph <- Select_merged_TP53_forgraph[-c(4)]
Select_merged_CCND1cnv_forgraph <- as.data.frame(t(rbind(Cyto_geom.mean,Chemo_geom.mean,Tfun_geom.mean, want_site$CCND1_cnv)))
colnames(Select_merged_CCND1cnv_forgraph) <- c("Cytotoxicity_nonHLA","Chemokines","T_cell_function","CCND1_cnv")
Select_merged_CCND1cnv_forgraph$CCND1_cnv <- as.character(Select_merged_CCND1cnv_forgraph$CCND1_cnv)
Select_merged <- as.data.frame(t(rbind(TP_geom.mean,TN_geom.mean,Cyto_geom.mean,Chemo_geom.mean,Tfun_geom.mean,
                                       want_site$final_site,want_site$final_HPV, want_site$TP53, want_site$PIK3CA, want_site$EGFR_cnv, want_site$CDKN2A,want_site$CCND1_cnv)))
colnames(Select_merged) <- c("T_cell_P","T_cell_N","Cytotoxicity_nonHLA","Chemokines","T_cell_function", "Site","HPV", "TP53", "PIK3CA", "EGFR_cnv", "CDKN2A", "CCND1_cnv")
colnames(Select_merged_Site) <- c("T_cell_P","T_cell_N","Cytotoxicity_nonHLA","Chemokines","T_cell_function","Site")
colnames(Select_merged_HPV) <- c("Cytotoxicity_nonHLA","Chemokines","T_cell_function","HPV")
colnames(Select_merged_TP53_forgraph) <- c("Cytotoxicity_nonHLA","Chemokines","T_cell_function","TP53")

Select_merged_Site$T_cell_function <- as.numeric(Select_merged_Site$T_cell_function)
Select_merged_Site %>% group_by(Site)

colors_in=c("#4faf4b4D", "#94509a4D","#e41d1f4D","#377fb94D")


###subfig7A (이중에서 하나)
boxdraw <-  function(Anno){
  trimed = cbind(Select_merged_Site[,Anno], Select_merged_Site[,"Site"])
  colnames(trimed) = c("Annotation", "Final")
  trimed = na.omit(trimed)
  plot <- ggplot(as.data.frame(trimed))+
    geom_boxplot(aes(x = Final, y = as.numeric(Annotation)), fill = colors_in)+
    stat_compare_means(aes(x = Final, y = as.numeric(Annotation)), method = 'kruskal.test',
                       label.x = 4) +
    scale_fill_manual(values = colors_in) +
    theme_classic()
  
  return(plot)
}
# setwd("/data/project/TRIUMPH/Figure2/new_figure_manu/")
plot <- boxdraw("T_cell_function")
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig6A.pdf", width = 15, height = 10)   ####Lip or oral cavity_Hypopharynx p.value  0.03037633
print(plot)
dev.off()


###subfig7B
colors_in=c("#4faf4b4D", "#94509a4D","#e41d1f4D","#377fb94D", "#ee4036","e08f5e")
Select_merged_Site

colors_in=c("#4faf4b4D", "#94509a4D","#e41d1f4D","#377fb94D")



Anno <- "T_cell_function"
colors_in=c("#4faf4b4D", "#94509a4D","#e41d1f4D","#377fb94D")
boxdraw2 <-  function(Anno){
  trimed = cbind(Select_merged_TP53_forgraph[,Anno], Select_merged_TP53_forgraph[,"TP53"])
  colnames(trimed) = c("Annotation", "Final")
  trimed = na.omit(trimed)
  plot <- ggplot(as.data.frame(trimed))+
    geom_boxplot(aes(x = Final, y = as.numeric(Annotation)), fill = c("#a6cee3","#fb9a99"))+
    stat_compare_means(aes(x = Final, y = as.numeric(Annotation)), method = 'wilcox',
                       label.x = 2) +
    scale_fill_manual(values = colors_in) +
    theme_classic() +
    geom_signif(aes(x = Final, y = as.numeric(Annotation)),
                comparisons = list(c("FALSE", "TRUE")),
                map_signif_level = TRUE)
  
  return(plot)
}


plot <- boxdraw2("T_cell_function")
plot
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/subfig6B_1.pdf", width = 15, height = 10)   ####Lip or oral cavity_Hypopharynx p.value  0.03037633
print(plot)
dev.off()


as.character(Select_merged_CCND1cnv_forgraph$CCND1_cnv)
boxdraw2 <-  function(Anno){
  trimed = as.data.frame(cbind(Select_merged_CCND1cnv_forgraph[,Anno], as.character(Select_merged_CCND1cnv_forgraph$CCND1_cnv)))
  colnames(trimed) = c("Annotation", "Final")
  trimed = na.omit(trimed)
  plot <- ggplot(as.data.frame(trimed))+
    geom_boxplot(aes(x = Final, y = as.numeric(Annotation)), fill = c("#a6cee3","#fb9a99"))+
    stat_compare_means(aes(x = Final, y = as.numeric(Annotation)), method = 'wilcox',
                       label.x = 2) +
    scale_fill_manual(values = colors_in) +
    theme_classic() +
    geom_signif(aes(x = Final, y = as.numeric(Annotation)),
                comparisons = list(c("FALSE", "TRUE")),
                map_signif_level = TRUE)

  return(plot)
}

plot <- boxdraw2("T_cell_function")
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig6B_2.pdf", width = 15, height = 10)   ####Lip or oral cavity_Hypopharynx p.value  0.03037633
print(plot)
dev.off()






Select_merged_TP53 <-  Select_merged %>% filter(TP53==TRUE)
boxdraw2 <-  function(Anno, target){

  trimed = cbind(Select_merged_TP53[,Anno], Select_merged_TP53[,target])
  colnames(trimed) = c("Annotation", "Final")
  trimed = na.omit(trimed)
  plot <- ggplot(as.data.frame(trimed))+
    geom_boxplot(aes(x = Final, y = as.numeric(Annotation)), fill = c("#a6cee3","#fb9a99"))+
    stat_compare_means(aes(x = Final, y = as.numeric(Annotation)), method = 'wilcox',
                       label.x = 2) +
    ggtitle(paste0(Anno,"_",target,"TP53_posi"))+
    scale_fill_manual(values = colors_in) +
    theme_classic() +
    geom_signif(aes(x = Final, y = as.numeric(Annotation)),
                comparisons = list(c("FALSE", "TRUE")),
                map_signif_level = TRUE)

}

plot <- boxdraw2("T_cell_function", "CCND1_cnv")
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig6B_3.pdf", width = 15, height = 10)   ####Lip or oral cavity_Hypopharynx p.value  0.03037633
print(plot)
dev.off()
