library(readxl)  ###P
# library(circlize)  ###N
library(ComplexHeatmap)  ###P
library('ggrepel')
library(stringr)


ht_opt$message = FALSE

# # 
# grey1 <- 0.9
# grey2 <- 0.4
# grey3 <- 0.5
# grey4 <- 0.0


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
reboot_data <- read.table("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/NanoString/reboot.txt", header = T, sep = "\t")
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
#################################Normalized counts###############################
#################################################################################
log2with_clinical <- log2(with_clinical)
with_clinical.z <- t(scale(t(log2with_clinical)))

# wssplot(with_clinical.z)
km <- kmeans(with_clinical.z, centers=2)
with_clinical.z.km <- cbind(with_clinical.z, km$cluster)
o <- order(with_clinical.z.km[,length(with_clinical.z.km[1,])])
with_clinical.z.km.sort <- with_clinical.z.km[o, 1:(length(with_clinical.z.km[1,])-1)]

with_clinical.z.km.sort.modi <- ifelse(with_clinical.z.km.sort > 4, 4, with_clinical.z.km.sort)
with_clinical.z.km.sort.modi <- ifelse(with_clinical.z.km.sort.modi < -4, -4, with_clinical.z.km.sort.modi)

normal_mean <- apply(reboot_data.m, 1, function(x) mean(as.numeric(x[c("SS0912168","SS1007561","SS1138427")])))

fil_FC <- log2((reboot_data.m[,c(1:(length(reboot_data.m[1,])-3))])/normal_mean)

