source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_for_H1516.R")

H15_inter$Tumor_Sample_Barcode <- unlist(H15_inter$Tumor_Sample_Barcode.1)

All_cnvkit = merge_mafs(fileNames_H15_16, clinicalData = H15_inter, cnTable = H15_16_cnv)

H15_inter #from import oncorpint H1516

data_H15 <- All_cnvkit@data


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
H15_inter$final_ostime <- as.numeric(H15_inter$final_ostime)
H15_inter$final_death <- as.numeric(H15_inter$final_death)
# H15_inter[is.na(H15_inter) ==TRUE] 
H15_inter[H15_inter ==""] <- NA
H15_inter[H15_inter =="NA"] <- NA


df_final <- data.frame()
for (i in c(1:5)){
  arm = paste0('Arm ', i)
  print(arm)
  df_tmp <- H15_inter[H15_inter$final_arm == arm,]
  newtmp <- data.frame()
  for (n in c(1:nrow(df_tmp))){
    tmp <- df_tmp[n,]
    other_gene_primary <- c()
    for (gene in get(paste0('genes_arm', i))){
      print(gene)
      if (!('Main_mutation' %in% colnames(tmp))){
        if (tmp[get('gene')] == TRUE){
          tmp['Main_mutation'] <- gene
          tmp['Variant_Classification'] <-paste(data_H15$Variant_Classification[(data_H15$Tumor_Sample_Barcode == tmp$Tumor_Sample_Barcode)&(data_H15$Hugo_Symbol == get('gene'))], collapse = ';') 
          tmp['HGVSp_Short'] <-paste(data_H15$HGVSp_Short[(data_H15$Tumor_Sample_Barcode == tmp$Tumor_Sample_Barcode)&(data_H15$Hugo_Symbol == get('gene'))], collapse = ';') 
        }  
      } else {
        if (tmp[get('gene')] == TRUE){
          other_gene_primary <- c(other_gene_primary, gene)
        }  
      }
    }
    tmp['Other_primary_alteration'] <- paste(other_gene_primary,collapse=";")
    gene_other <- gene_li[!(gene_li %in% get(paste0('genes_arm', i)))]
    others_gene <- data.frame()
    for (gene in gene_other){
      if (tmp[get('gene')] == TRUE){
        others_gene <- c(others_gene, gene)
      }        
    }
    tmp['Other_not_primary'] <- paste(others_gene, collapse = ';')
    
    if (!('Main_mutation' %in% colnames(tmp))){
      tmp['Main_mutation'] <- NA
      tmp['Variant_Classification'] <- NA
      tmp['HGVSp_Short'] <- NA
    } 
    newtmp <- rbind(newtmp, tmp)
  }
  df_final <- rbind(df_final,newtmp)
}

H15_inter$BOR_final2 <- ifelse(H15_inter$BOR_final=="PR", "Responder", "Nonresponder")
df_final$BOR_final2 <- ifelse(df_final$BOR_final=="PR", "Responder", "Nonresponder")
df_final$BOR_final3 <- ifelse(df_final$BOR_final=="PD", "Nonresponder", "Responder")
df_final$HGVSp_Short[df_final$Tumor_Sample_Barcode == '31-S010-FP' & df_final$Main_mutation == "PIK3CA" & df_final$final_arm == 'Arm 1']
df_final$HGVSp_Short[df_final$Main_mutation == "PIK3CA"& df_final$TRTP1 == "Arm 1 (BYL719)"]
df_final$HGVSp_Short[df_final$Main_mutation == "EGFR"& df_final$final_arm == 'Arm 2']
df_final$HGVSp_Short[df_final$Main_mutation == "FGFR1"& df_final$final_arm == 'Arm 3']
df_final$HGVSp_Short[df_final$Main_mutation == "FGFR2"& df_final$final_arm == 'Arm 3']
df_final$HGVSp_Short[df_final$Main_mutation == "CDKN2A"& df_final$final_arm == 'Arm 4']

data_H15[data_H15$Tumor_Sample_Barcode == '31-S010-FP']
####add CDKN2A and EGFR
data_CDKN2A <- data_H15[data_H15$Hugo_Symbol=='CDKN2A']
data_EGFR <- data_H15[data_H15$Hugo_Symbol=='EGFR']
data_PIK3CA <- data_H15[data_H15$Hugo_Symbol=='PIK3CA']
data_ERBB2 <- data_H15[data_H15$Hugo_Symbol=='ERBB2']
data_TP53 <- data_H15[data_H15$Hugo_Symbol=='TP53']

data_CDKN2A$Variant_Classification

data_CDKN2A$Tumor_Sample_Barcode

df$Tumor_Sample_Barcode
# c(df_final$Tumor_Sample_Barcode, df_final$TRTP1, df_final$CTRTP1)


library(dplyr)


library(dplyr)


library(dplyr)


data_CDKN2A$CLIN_SIG

# Perform a left join using dplyr, and create a new column with the variant_classification values
merged_df <- left_join(df, data_CDKN2A %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(Variant_Classification_CDKN2A = paste(Variant_Classification, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(Variant_Classification_CDKN2A = ifelse(is.na(Variant_Classification_CDKN2A), "NA", Variant_Classification_CDKN2A))

# Perform a left join using dplyr, and create a new column with the HGVSp_Short values
merged_df <- left_join(merged_df, data_CDKN2A %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(HGVSp_Short_CDKN2A = paste(HGVSp_Short, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(HGVSp_Short_CDKN2A = ifelse(is.na(HGVSp_Short_CDKN2A), "NA", HGVSp_Short_CDKN2A))

# Perform a left join using dplyr, and create a new column with the HGVSc values
merged_df <- left_join(merged_df, data_CDKN2A %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(HGVSc_CDKN2A = paste(HGVSc, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(HGVSc_CDKN2A = ifelse(is.na(HGVSc_CDKN2A), "NA", HGVSc_CDKN2A))

# Perform a left join using dplyr, and create a new column with the IMPACT values
merged_df <- left_join(merged_df, data_CDKN2A %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(IMPACT_CDKN2A = paste(IMPACT, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(IMPACT_CDKN2A = ifelse(is.na(IMPACT_CDKN2A), "NA", IMPACT_CDKN2A))
# Perform a left join using dplyr, and create a new column with the gnomAD_AF values
merged_df <- left_join(merged_df, data_CDKN2A %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(gnomAD_AF_CDKN2A = paste(gnomAD_AF, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(gnomAD_AF_CDKN2A = ifelse(is.na(gnomAD_AF_CDKN2A), "NA", gnomAD_AF_CDKN2A))
# Perform a left join using dplyr, and create a new column with the i_TumorVAF_WU values
# merged_df <- left_join(merged_df, data_CDKN2A %>% 
#                          group_by(Tumor_Sample_Barcode) %>% 
#                          summarize(i_TumorVAF_WU_CDKN2A = paste(i_TumorVAF_WU, collapse = ";")), 
#                        by = "Tumor_Sample_Barcode") %>%
#   mutate(i_TumorVAF_WU_CDKN2A = ifelse(is.na(i_TumorVAF_WU_CDKN2A), "NA", i_TumorVAF_WU_CDKN2A))
# Perform a left join using dplyr, and create a new column with the CLIN_SIG values
merged_df <- left_join(merged_df, data_CDKN2A %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(CLIN_SIG_CDKN2A = paste(CLIN_SIG, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(CLIN_SIG_CDKN2A = ifelse(is.na(CLIN_SIG_CDKN2A), "NA", CLIN_SIG_CDKN2A))
data_H15[data_H15$Tumor_Sample_Barcode == '24-S023-FP' & data_H15$Hugo_Symbol == 'FGFR1',]$HGVSp_Short
data_H15[data_H15$Tumor_Sample_Barcode == '24-S023-FP' & data_H15$Hugo_Symbol == 'FGFR1',]$HGVSc
data_H15[data_H15$Tumor_Sample_Barcode == '24-S023-FP' & data_H15$Hugo_Symbol == 'FGFR1',]$Start_Position


############EGFR
merged_df <- left_join(merged_df, data_EGFR %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(Variant_Classification_EGFR = paste(Variant_Classification, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(Variant_Classification_EGFR = ifelse(is.na(Variant_Classification_EGFR), "NA", Variant_Classification_EGFR))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_EGFR %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(HGVSp_Short_EGFR = paste(HGVSp_Short, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(HGVSp_Short_EGFR = ifelse(is.na(HGVSp_Short_EGFR), "NA", HGVSp_Short_EGFR))


# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_EGFR %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(HGVSc_EGFR = paste(HGVSc, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(HGVSc_EGFR = ifelse(is.na(HGVSc_EGFR), "NA", HGVSc_EGFR))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_EGFR %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(IMPACT_EGFR = paste(IMPACT, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(IMPACT_EGFR = ifelse(is.na(IMPACT_EGFR), "NA", IMPACT_EGFR))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_EGFR %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(gnomAD_AF_EGFR = paste(gnomAD_AF, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(gnomAD_AF_EGFR = ifelse(is.na(gnomAD_AF_EGFR), "NA", gnomAD_AF_EGFR))
# # Perform a left join using dplyr, and create a new column with the Variant_Classification values
# merged_df <- left_join(merged_df, data_EGFR %>% 
#                          group_by(Tumor_Sample_Barcode) %>% 
#                          summarize(i_TumorVAF_WU_EGFR = paste(i_TumorVAF_WU, collapse = ";")), 
#                        by = "Tumor_Sample_Barcode") %>%
#   mutate(i_TumorVAF_WU_EGFR = ifelse(is.na(i_TumorVAF_WU_EGFR), "NA", i_TumorVAF_WU_EGFR))
# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_EGFR %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(CLIN_SIG_EGFR = paste(CLIN_SIG, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(CLIN_SIG_EGFR = ifelse(is.na(CLIN_SIG_EGFR), "NA", CLIN_SIG_EGFR))
############ERBB2
merged_df <- left_join(merged_df, data_ERBB2 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(Variant_Classification_ERBB2 = paste(Variant_Classification, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(Variant_Classification_ERBB2 = ifelse(is.na(Variant_Classification_ERBB2), "NA", Variant_Classification_ERBB2))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_ERBB2 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(HGVSp_Short_ERBB2 = paste(HGVSp_Short, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(HGVSp_Short_ERBB2 = ifelse(is.na(HGVSp_Short_ERBB2), "NA", HGVSp_Short_ERBB2))


# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_ERBB2 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(HGVSc_ERBB2 = paste(HGVSc, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(HGVSc_ERBB2 = ifelse(is.na(HGVSc_ERBB2), "NA", HGVSc_ERBB2))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_ERBB2 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(IMPACT_ERBB2 = paste(IMPACT, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(IMPACT_ERBB2 = ifelse(is.na(IMPACT_ERBB2), "NA", IMPACT_ERBB2))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_ERBB2 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(gnomAD_AF_ERBB2 = paste(gnomAD_AF, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(gnomAD_AF_ERBB2 = ifelse(is.na(gnomAD_AF_ERBB2), "NA", gnomAD_AF_ERBB2))
# # Perform a left join using dplyr, and create a new column with the Variant_Classification values
# merged_df <- left_join(merged_df, data_ERBB2 %>% 
#                          group_by(Tumor_Sample_Barcode) %>% 
#                          summarize(i_TumorVAF_WU_ERBB2 = paste(i_TumorVAF_WU, collapse = ";")), 
#                        by = "Tumor_Sample_Barcode") %>%
#   mutate(i_TumorVAF_WU_ERBB2 = ifelse(is.na(i_TumorVAF_WU_ERBB2), "NA", i_TumorVAF_WU_ERBB2))
# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_ERBB2 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(CLIN_SIG_ERBB2 = paste(CLIN_SIG, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(CLIN_SIG_ERBB2 = ifelse(is.na(CLIN_SIG_ERBB2), "NA", CLIN_SIG_ERBB2))

############TP53
merged_df <- left_join(merged_df, data_TP53 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(Variant_Classification_TP53 = paste(Variant_Classification, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(Variant_Classification_TP53 = ifelse(is.na(Variant_Classification_TP53), "NA", Variant_Classification_TP53))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_TP53 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(HGVSp_Short_TP53 = paste(HGVSp_Short, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(HGVSp_Short_TP53 = ifelse(is.na(HGVSp_Short_TP53), "NA", HGVSp_Short_TP53))


# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_TP53 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(HGVSc_TP53 = paste(HGVSc, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(HGVSc_TP53 = ifelse(is.na(HGVSc_TP53), "NA", HGVSc_TP53))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_TP53 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(IMPACT_TP53 = paste(IMPACT, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(IMPACT_TP53 = ifelse(is.na(IMPACT_TP53), "NA", IMPACT_TP53))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_TP53 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(gnomAD_AF_TP53 = paste(gnomAD_AF, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(gnomAD_AF_TP53 = ifelse(is.na(gnomAD_AF_TP53), "NA", gnomAD_AF_TP53))
# # Perform a left join using dplyr, and create a new column with the Variant_Classification values
# merged_df <- left_join(merged_df, data_TP53 %>% 
#                          group_by(Tumor_Sample_Barcode) %>% 
#                          summarize(i_TumorVAF_WU_TP53 = paste(i_TumorVAF_WU, collapse = ";")), 
#                        by = "Tumor_Sample_Barcode") %>%
#   mutate(i_TumorVAF_WU_TP53 = ifelse(is.na(i_TumorVAF_WU_TP53), "NA", i_TumorVAF_WU_TP53))
# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_TP53 %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(CLIN_SIG_TP53 = paste(CLIN_SIG, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(CLIN_SIG_TP53 = ifelse(is.na(CLIN_SIG_TP53), "NA", CLIN_SIG_TP53))
############PIK3CA
merged_df <- left_join(merged_df, data_PIK3CA %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(Variant_Classification_PIK3CA = paste(Variant_Classification, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(Variant_Classification_PIK3CA = ifelse(is.na(Variant_Classification_PIK3CA), "NA", Variant_Classification_PIK3CA))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_PIK3CA %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(HGVSp_Short_PIK3CA = paste(HGVSp_Short, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(HGVSp_Short_PIK3CA = ifelse(is.na(HGVSp_Short_PIK3CA), "NA", HGVSp_Short_PIK3CA))


# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_PIK3CA %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(HGVSc_PIK3CA = paste(HGVSc, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(HGVSc_PIK3CA = ifelse(is.na(HGVSc_PIK3CA), "NA", HGVSc_PIK3CA))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_PIK3CA %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(IMPACT_PIK3CA = paste(IMPACT, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(IMPACT_PIK3CA = ifelse(is.na(IMPACT_PIK3CA), "NA", IMPACT_PIK3CA))

# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_PIK3CA %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(gnomAD_AF_PIK3CA = paste(gnomAD_AF, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(gnomAD_AF_PIK3CA = ifelse(is.na(gnomAD_AF_PIK3CA), "NA", gnomAD_AF_PIK3CA))
# # Perform a left join using dplyr, and create a new column with the Variant_Classification values
# merged_df <- left_join(merged_df, data_PIK3CA %>% 
#                          group_by(Tumor_Sample_Barcode) %>% 
#                          summarize(i_TumorVAF_WU_PIK3CA = paste(i_TumorVAF_WU, collapse = ";")), 
#                        by = "Tumor_Sample_Barcode") %>%
#   mutate(i_TumorVAF_WU_PIK3CA = ifelse(is.na(i_TumorVAF_WU_PIK3CA), "NA", i_TumorVAF_WU_PIK3CA))
# Perform a left join using dplyr, and create a new column with the Variant_Classification values
merged_df <- left_join(merged_df, data_PIK3CA %>% 
                         group_by(Tumor_Sample_Barcode) %>% 
                         summarize(CLIN_SIG_PIK3CA = paste(CLIN_SIG, collapse = ";")), 
                       by = "Tumor_Sample_Barcode") %>%
  mutate(CLIN_SIG_PIK3CA = ifelse(is.na(CLIN_SIG_PIK3CA), "NA", CLIN_SIG_PIK3CA))

merged_df$Variant_Classification_TP53
merged_df$HGVSp_Short_TP53


merged_df$HGVSc_EGFR
merged_df$CLIN_SIG_CDKN2A
merged_df$CLIN_SIG_EGFR
merged_df$CLIN_SIG_ERBB2
# merged_df$i_TumorVAF_WU_ERBB2
#######이거 돌려줘야  merge 문제안생김
merged_df <- merged_df %>%
  mutate(CLIN_SIG_CDKN2A = gsub(",", ";", CLIN_SIG_CDKN2A))
merged_df <- merged_df %>%
  mutate(CLIN_SIG_EGFR = gsub(",", ";", CLIN_SIG_EGFR))
merged_df <- merged_df %>%
  mutate(CLIN_SIG_ERBB2 = gsub(",", ";", CLIN_SIG_ERBB2))
merged_df <- merged_df %>%
  mutate(CLIN_SIG_TP53 = gsub(",", ";", CLIN_SIG_TP53))






# Load the tidyr library
library(tidyr)

df_final_split_variant <- df_final %>%
  tidyr::separate(HGVSp_Short, into = c("HGVSp_Short"), sep = ",", remove = FALSE) %>%
  tidyr::spread(HGVSp_Short, HGVSp_Short)

df_final$PIK3CA
df_final$Main_mutation


df_final$HGVSp_Short[df_final$Main_mutation == "PIK3CA"& df_final$TRTP1 == "Arm 1 (BYL719)"]


df_final$PIK3CA_diff <- ifelse(df_final$Main_mutation == "PIK3CA", df_final$HGVSp_Short, NA)



# The resulting data frame

merged_df$HGVSp_Short_PIK3CA

df_final$HGVSp_Short

df_final2 <- df_final %>%
  tidyr::separate_rows(HGVSp_Short, sep = ";") %>%
  mutate(original_index = row_number() - 1) %>%
  group_by(original_index) %>%
  mutate(HGVSp_Short_num = row_number()) %>%
  ungroup() %>%
  dplyr::select(-original_index)

# df_final3 <- merged_df %>%
#   tidyr::separate_rows(HGVSp_Short_PIK3CA, sep = ";") %>%
#   mutate(original_index = row_number() - 1) %>%
#   group_by(original_index) %>%
#   mutate(HGVSp_Short_PIK3CA_num = row_number()) %>%
#   ungroup() %>%
#   dplyr::select(-original_index)


library(dplyr)

###BORfinal 종류
H15_inter$BOR_final2 <- ifelse(H15_inter$BOR_final=="PR", "Responder", "Nonresponder")
df_final$BOR_final2 <- ifelse(df_final$BOR_final=="PR", "Responder", "Nonresponder")
df_final$BOR_final3 <- ifelse(df_final$BOR_final=="PD", "Nonresponder", "Responder")

df_final2$Tumor_Sample_Barcode

df_final2$PIK3CA_diff <- ifelse(df_final2$Main_mutation == "PIK3CA", df_final2$HGVSp_Short, NA)
# df_final3$PIK3CA_diff <- ifelse(df_final3$Main_mutation == "PIK3CA", df_final3$HGVSp_Short_PIK3CA, NA)
df_final2$BOR_final <- unlist(df_final2$BOR_final)
table(df_final2$BOR_final, df_final2$PIK3CA_diff)
table(df_final2$BOR_final3, df_final2$PIK3CA_diff)
fisher.test(table(df_final2$BOR_final, df_final2$PIK3CA_diff))
df_frequentPIK3CA <- subset(df_final2, PIK3CA_diff %in% c("p.E542K", "p.E545K", "p.H1047R"))
# df_frequentPIK3CA <- subset(df_final3, PIK3CA_diff %in% c("p.E542K", "p.E545K", "p.H1047R"))
table(df_frequentPIK3CA$BOR_final, df_frequentPIK3CA$PIK3CA_diff)
table(df_frequentPIK3CA$BOR_final2, df_frequentPIK3CA$PIK3CA_diff)
table(df_frequentPIK3CA$BOR_final3, df_frequentPIK3CA$PIK3CA_diff)
fisher.test(table(df_frequentPIK3CA$BOR_final, df_frequentPIK3CA$PIK3CA_diff))
fisher.test(table(df_frequentPIK3CA$BOR_final3, df_frequentPIK3CA$PIK3CA_diff))
fisher.test(table(df_frequentPIK3CA$BOR_final3, df_frequentPIK3CA$PIK3CA_diff))


library(RColorBrewer)
# Choose an academic-friendly color palette (e.g., "Set1")
my_colors <- brewer.pal(3, "Set1")  # 3 colors for 3 mutation types

pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Figure2B.pdf", height = 6, width = 9)

##portrait (5,5)
# Assuming df_frequentPIK3CA is your data frame
# Define colors for the different statuses with the specified color codes
my_colors <- c("PD" = "#E41A1C", "PR" = "#377EB8", "SD" = "#4DAF4A")

# Load necessary libraries
library(ggplot2)

# Ensure that 'BOR_final' is a factor with the correct levels
df_frequentPIK3CA$BOR_final <- factor(df_frequentPIK3CA$BOR_final)


# Adjust the factor levels of BOR_final to ensure the desired order
df_frequentPIK3CA$BOR_final <- factor(df_frequentPIK3CA$BOR_final, levels = c("PR", "SD", "PD"))

# Define colors for the different statuses with the specified color codes
my_colors <- c("PD" = "#E41A1C", "PR" = "#377EB8", "SD" = "#4DAF4A")

# ggplot code with bars and legend in the specified order, including vacant column for "SD"
B2 <- ggplot(data = subset(df_frequentPIK3CA, !is.na(BOR_final)), aes(x = PIK3CA_diff)) +
  geom_bar(aes(fill = BOR_final), position = position_dodge(preserve = "single"), width = 0.7) +
  labs(x = "Mutation type", y = "Number of Patients", fill = "Best overall response") +
  scale_fill_manual(values = my_colors) +  # Assuming my_colors is predefined
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20), # Increase x-axis text size
    axis.text.y = element_text(angle = 0, hjust = 1, size = 20),   # Increase y-axis text size
    axis.title.x = element_text(size = 22), # Increase x-axis title size
    axis.title.y = element_text(size = 22),
    axis.ticks = element_line(color = "black"), # Add ticks on axes
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), # Remove panel border
    axis.line.x = element_line(color = "black"), # Add x-axis line
    axis.line.y = element_line(color = "black"), # Add y-axis line
    legend.title = element_text(size=20),
    legend.text = element_text(size=20),
    legend.key.height= unit(1.2, 'cm'),
    legend.key.width= unit(1.2, 'cm')
  )

B2

dev.off()
