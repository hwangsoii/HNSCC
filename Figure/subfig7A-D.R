library(readxl)    ###P
# library(circlize)  ###N
library(ComplexHeatmap)  ###P
library('ggrepel')
library(tibble)

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
cnvkit_0.4_onlyhnscc <- cnvkit_0.4[!(cnvkit_0.4$sample %in% nonhnscc), ]
patients <- clinical[apply(clinical, 1, function(x) x["Sample2"]%in%samples),]

with_clinical <- reboot_data.m[,patients$Sample2]


#################################################################################
#################################################################################
#################################################################################
### RNA expression scaling


normal_maxmin <- function(x){
  return((x-min(x))/(max(x)-min(x)))
}

patients %>% group_by(final_site) %>% count(final_site)
want_site <- patients %>%
  filter(final_site =='oropharynx'|final_site =='Hypopharynx'|final_site =='Lip or oral cavity'|final_site =='Larynx')
want_site_clinical <- reboot_data.m[,want_site$Sample2]

with_clinical.maxmin <- t(apply(want_site_clinical, 1,normal_maxmin))

patients$Sample2

with_clinical.maxmin.t <- t(with_clinical.maxmin)
rownames(with_clinical.maxmin.t)
new_with_clinical.maxmin.t <- as.data.frame(with_clinical.maxmin.t) %>% rownames_to_column(var = "Sample2")

new_merged <- merge(patients, new_with_clinical.maxmin.t, by ='Sample2')




##############
colors_in=c("#4faf4b4D", "#94509a4D","#e41d1f4D","#377fb94D")

boxdraw2 <-  function(Anno){
  trimed = cbind(new_merged[,Anno], new_merged[,"final_HPV"])
  colnames(trimed) = c("Expression", "HPV_status")
  trimed = na.omit(trimed)
  plot <- ggplot(as.data.frame(trimed))+
    geom_boxplot(aes(x = HPV_status, y = as.numeric(Expression)), fill = c("#a6cee3","#fb9a99"))+
    ggtitle(Anno)+
    stat_compare_means(aes(x = HPV_status, y = as.numeric(Expression)), method = 'wilcox',
                       label.x = 2) +
    scale_fill_manual(values = colors_in) +
    theme_classic() +
    geom_signif(aes(x = HPV_status, y = as.numeric(Expression)),
                comparisons = list(c("Negative", "Positive")),
                map_signif_level = TRUE)
  
  
  return(plot)
}

# 
# plot <- boxdraw2("PDCD1")
# pdf("PDCD1_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
# plot <- boxdraw2("CD274")
# pdf("CD274_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
# plot <- boxdraw2("CD276")
# pdf("CD276_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
# plot <- boxdraw2("LAG3")
# pdf("LAG3_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
# plot <- boxdraw2("CTLA4")
# pdf("CTLA4_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
# plot <- boxdraw2("CD80")
# pdf("CD80_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
# plot <- boxdraw2("CTLA4")
# pdf("CTLA4_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
plot <- boxdraw2("IFNG")
plot
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig7A.pdf", width = 15, height = 10)
print(plot)
dev.off()
# plot <- boxdraw2("TIGIT")
# pdf("TIGIT_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
# plot <- boxdraw2("CD86")
# pdf("CD86_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
# plot <- boxdraw2("CXCL9")
# pdf("CXCL9_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
# plot <- boxdraw2("CXCL10")
# pdf("CXCL10_HPV.pdf", width = 15, height = 10)
# print(plot)
# dev.off()
plot <- boxdraw2("GZMB")
plot
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig7C.pdf", width = 15, height = 10)
print(plot)
dev.off()
plot <- boxdraw2("GZMK")
plot
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig7D.pdf", width = 15, height = 10)
print(plot)
dev.off()

# positive_values <- subset(new_merged, final_HPV == "Positive")$PDCD1
# negative_values <- subset(new_merged, final_HPV == "Negative")$PDCD1
# na_values <- subset(new_merged, is.na(final_HPV))$PDCD1
# 
# mean(positive_values)
# mean(negative_values)
# 
# wilcox.test(positive_values, negative_values)
# 
# positive_values <- subset(new_merged, final_HPV == "Positive")$CD274
# negative_values <- subset(new_merged, final_HPV == "Negative")$CD274
# na_values <- subset(new_merged, is.na(final_HPV))$CD274
# mean(positive_values)
# mean(negative_values)
# wilcox.test(positive_values, negative_values)
# 
# positive_values <- subset(new_merged, final_HPV == "Positive")$CD276
# negative_values <- subset(new_merged, final_HPV == "Negative")$CD276
# na_values <- subset(new_merged, is.na(final_HPV))$CD276
# mean(positive_values)
# mean(negative_values)
# wilcox.test(positive_values, negative_values)
# 
# 
# positive_values <- subset(new_merged, final_HPV == "Positive")$CD276
# negative_values <- subset(new_merged, final_HPV == "Negative")$CD276
# na_values <- subset(new_merged, is.na(final_HPV))$CD276
# mean(positive_values)
# mean(negative_values)
# wilcox.test(positive_values, negative_values)
# 
# positive_values <- subset(new_merged, final_HPV == "Positive")$LAG3
# negative_values <- subset(new_merged, final_HPV == "Negative")$LAG3
# na_values <- subset(new_merged, is.na(final_HPV))$LAG3
# mean(positive_values)
# mean(negative_values)
# wilcox.test(positive_values, negative_values)
# 
# new_merged$LAG3
# 
# 
# positive_values <- subset(new_merged, final_HPV == "Positive")$CTLA4
# negative_values <- subset(new_merged, final_HPV == "Negative")$CTLA4
# na_values <- subset(new_merged, is.na(final_HPV))$CTLA4
# mean(positive_values)
# mean(negative_values)
# wilcox.test(positive_values, negative_values)
# 
# positive_values <- subset(new_merged, final_HPV == "Positive")$CD80
# negative_values <- subset(new_merged, final_HPV == "Negative")$CD80
# na_values <- subset(new_merged, is.na(final_HPV))$CD80
# mean(positive_values)
# mean(negative_values)
# wilcox.test(positive_values, negative_values)
# 
# positive_values <- subset(new_merged, final_HPV == "Positive")$`TIMD4 (TIM-4)`
# negative_values <- subset(new_merged, final_HPV == "Negative")$`TIMD4 (TIM-4)`
# na_values <- subset(new_merged, is.na(final_HPV))$`TIMD4 (TIM-4)`
# mean(positive_values)
# mean(negative_values)
# wilcox.test(positive_values, negative_values)
# 
# positive_values <- subset(new_merged, final_HPV == "Positive")$CD86
# negative_values <- subset(new_merged, final_HPV == "Negative")$CD86
# na_values <- subset(new_merged, is.na(final_HPV))$CD86
# mean(positive_values)
# mean(negative_values)
# wilcox.test(positive_values, negative_values)
# 
# 
# positive_values <- subset(new_merged, final_HPV == "Positive")$IFNG
# negative_values <- subset(new_merged, final_HPV == "Negative")$IFNG
# na_values <- subset(new_merged, is.na(final_HPV))$IFNG
# mean(positive_values)
# mean(negative_values)
# wilcox.test(positive_values, negative_values)
# 
# positive_values <- subset(new_merged, final_HPV == "Positive")$TIGIT
# negative_values <- subset(new_merged, final_HPV == "Negative")$TIGIT
# na_values <- subset(new_merged, is.na(final_HPV))$TIGIT
# mean(positive_values)
# mean(negative_values)
# wilcox.test(positive_values, negative_values)
# 
# new_merged$CD2

# ###matrix build
# 
# Immune_Profile <- read_excel("/data/project/TRIUMPH/raw_data/LBL-10043-08_nCounter_PanCancer_Immune_Profiling_Panel_Gene_List.xlsx",
#                              sheet = "Annotations", col_names = T, skip = 1)
# Immune_Profile.df <- as.data.frame(Immune_Profile)
# rownames(Immune_Profile.df) <- Immune_Profile.df[,"Gene Name"]
# 
# Immune_Profile.df_used <- Immune_Profile.df[rownames(with_clinical),]
# Immune_Profile.df_used <- Immune_Profile.df_used[!is.na(Immune_Profile.df_used$`Gene Name`),]
# 
# Category <- unique(unlist(strsplit(Immune_Profile.df_used$`Immune Response Category`,", ")))
# Category <- Category[!is.na(Category)]
# 
# GenesByCategory = data.frame()
# 
# for(i in Category){
#   assign(i , Immune_Profile.df_used$`Gene Name`[grep(pattern = i, Immune_Profile.df_used$`Immune Response Category`)])
#   Names = Immune_Profile.df_used$`Gene Name`[grep(pattern = i, Immune_Profile.df_used$`Immune Response Category`)]
#   GenesByCategory =  rbind(GenesByCategory, paste(Names, collapse = ","))
# }
# 
# colnames(GenesByCategory) <- "Genes"
# rownames(GenesByCategory) <- Category
# 
# ###plot build
# 
# geom.mean_m <- numeric()
# 
# for(i in Category){
#   genes = unlist(strsplit(GenesByCategory[i,], ","))
#   if (length(genes) > 1){
#     geom.mean = apply(with_clinical.maxmin[genes,], 2, function(x) exp(mean(log(x))))
#   }
#   else{
#     geom.mean = with_clinical.maxmin[genes,]
#   }
#   geom.mean_m <- rbind(geom.mean_m, geom.mean)
# }
# geom.mean_df <- as.data.frame(geom.mean_m)
# rownames(geom.mean_df) <- Category
# 
# 
# merged_df <- as.data.frame(t(rbind(geom.mean_m,want_site$final_HPV, want_site$final_site)))
# 
# colnames(merged_df) <- c(rownames(geom.mean_m),"HPV","Site")
# 
# library(ggpubr)
# 
# for(i in c(1:15)){
#   for(j in c(16:17)){
#     trimed = cbind(merged_df[,i], merged_df[,j])
#     colnames(trimed) = c("Category", "Final")
#     trimed = na.omit(trimed)
#     plot <- ggplot(as.data.frame(trimed))+
#       geom_boxplot(aes(x = Final, y = as.numeric(Category)))+
#       stat_compare_means(aes(x = Final, y = as.numeric(Category)), method = 'kruskal.test')
#     pdf(paste0(Category[i],"_",colnames(merged_df)[j],".pdf"), width = 15, height = 10)
#     print(plot)
#     dev.off()
#   }
# }






#################################################################################
#################################################################################
#################################################################################


normal_maxmin <- function(x){
  return((x-min(x))/(max(x)-min(x)))
}

patients %>% group_by(final_site) %>% count(final_site)
want_site <- patients %>%
  filter(final_site =='oropharynx'|final_site =='Hypopharynx'|final_site =='Lip or oral cavity'|final_site =='Larynx')
want_site_clinical <- reboot_data.m[,want_site$Sample2]

with_clinical.maxmin <- t(apply(want_site_clinical, 1,normal_maxmin))



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

# Select_merged <- as.data.frame(t(rbind(TP_geom.mean,TN_geom.mean,Cyto_geom.mean,Chemo_geom.mean,Tfun_geom.mean,
                                       # want_site$final_site,want_site$final_HPV, want_site$TP53, want_site$PIK3CA, want_site$EGFR_cnv, want_site$CDKN2A,want_site$CCND1_cnv)))
Select_merged_Site <- as.data.frame(t(rbind(TP_geom.mean,TN_geom.mean,Cyto_geom.mean,Chemo_geom.mean,Tfun_geom.mean,
                                            want_site$final_site)))
Select_merged_HPV <- as.data.frame(t(rbind(Cyto_geom.mean,Chemo_geom.mean,Tfun_geom.mean, want_site$final_HPV)))
Select_merged_TP53_forgraph <- as.data.frame(t(rbind(Cyto_geom.mean,Chemo_geom.mean,Tfun_geom.mean,want_site$final_HPV,want_site$TP53)))
# Select_merged_TP53_forgraph <- Select_merged_TP53_forgraph[-c(4)]
# colnames(Select_merged) <- c("T_cell_P","T_cell_N","Cytotoxicity_nonHLA","Chemokines","T_cell_function", "Site","HPV", "TP53", "PIK3CA", "EGFR_cnv", "CDKN2A", "CCND1_cnv")
colnames(Select_merged_Site) <- c("T_cell_P","T_cell_N","Cytotoxicity_nonHLA","Chemokines","T_cell_function","Site")
colnames(Select_merged_HPV) <- c("Cytotoxicity_nonHLA","Chemokines","T_cell_function","HPV")
colnames(Select_merged_TP53_forgraph) <- c("Cytotoxicity_nonHLA","Chemokines","T_cell_function","TP53")

Select_merged_Site$T_cell_function <- as.numeric(Select_merged_Site$T_cell_function)
Select_merged_Site %>% group_by(Site)


###subfig8B 
colors_in=c("#377fb94D","#e41d1f4D")
boxdraw <-  function(Anno){
  trimed = cbind(Select_merged_HPV[,Anno], Select_merged_HPV[,"HPV"])
  colnames(trimed) = c("Annotation", "Final")
  trimed = na.omit(trimed)
  plot <- ggplot(as.data.frame(trimed))+
    geom_boxplot(aes(x = Final, y = as.numeric(Annotation)), fill = colors_in)+
    stat_compare_means(aes(x = Final, y = as.numeric(Annotation)), method = 'wilcox.test',
                       label.x = 2) +
    scale_fill_manual(values = colors_in) +
    theme_classic()
  
  return(plot)
}

# setwd("/data/project/TRIUMPH/Figure2/new_figure_manu/")
plot <- boxdraw("T_cell_function")
# print(plot)
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig7B.pdf", width = 15, height = 10)   ####Lip or oral cavity_Hypopharynx p.value  0.03037633
print(plot)
dev.off()
print(plot)
# setwd("/data/project/TRIUMPH/Figure/for/")
# for(i in c(1:3)){
#   for(j in c(4:5)){
#     trimed = cbind(Select_merged[,i], Select_merged[,j])
#     colnames(trimed) = c("Annotation", "Final")
#     trimed = na.omit(trimed)
#     plot <- ggplot(as.data.frame(trimed))+
#       geom_boxplot(aes(x = Final, y = as.numeric(Annotation)))+
#       stat_compare_means(aes(x = Final, y = as.numeric(Annotation)), method = 'kruskal.test')
#     pdf(paste0(colnames(Select_merged)[i],"_",colnames(Select_merged)[j],".pdf"), width = 15, height = 10)
#     print(plot)
#     dev.off()
#   }
# }

# with_clinical.maxmin