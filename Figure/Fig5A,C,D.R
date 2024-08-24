source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")




#######color change####### 
mutation_type <- unique(unlist(strsplit(unlist(as.list(mat_tot)),";")))
#mutation_type <- unique(unlist(strsplit(unlist(strsplit(unlist(as.list(mat_tot)),";")),'\\|')))
mutation_type <- mutation_type[mutation_type!=""]
# mutation_color <- as.vector(qualitative_hcl(n=7, h = c(0,270), c = 80, l=80))
# mutation_color2 <- brewer.pal(7, "Paired")
# mutation_color <- c(mutation_color2, mutation_color)
mutation_color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#F8C8DC",
                    "#C6CE3E", "#5EDF82", "#00E4D0", "#00D7FF", "#C7BAFF", "#808B96")
barplot(rep(1,length(mutation_type)), col = mutation_color, names.arg = names(mutation_color), cex.names = 0.5)
names(mutation_color) <- mutation_type[c(3,1,8,2,4,5,7,9,6,11,12,13,15,14,10)]

height = 0.98
width = 0.98

alter_fun = list(background = function(x, y, w, h) grid.rect(x, y, w, h, 
                                                             gp = gpar(fill = NA)),
                 Nonsense_Mutation=  function(x, y, w, h) grid.rect(x, y, w*width, h*height,
                                                                    gp = gpar(fill = mutation_color["Nonsense_Mutation"], col = NA)),
                 Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w*width, h*height,
                                                                    gp = gpar(fill = mutation_color["Missense_Mutation"], col = NA)),
                 Amplification = function(x, y, w, h) grid.rect(x, y, w*width, h*0.3,
                                                                gp = gpar(fill = mutation_color["Amplification"], col = NA)),
                 Frame_Shift_Ins = function(x, y, w, h) grid.rect(x, y, w*width, h*height, 
                                                                  gp = gpar(fill = mutation_color["Frame_Shift_Ins"], col = NA)),
                 In_Frame_Del = function(x, y, w, h) grid.rect(x, y, w*width, h*height,
                                                               gp = gpar(fill = mutation_color["In_Frame_Del"], col = NA)),
                 Splice_Site =function(x, y, w, h) grid.rect(x, y, w*width, h*height,
                                                             gp = gpar(fill = mutation_color["Splice_Site"], col = NA)),
                 Deletion =function(x, y, w, h) grid.rect(x, y, w*width, h*0.3,
                                                          gp = gpar(fill = mutation_color["Deletion"], col = NA)),
                 Frame_Shift_Del = function(x, y, w, h) grid.rect(x, y, w*width, h*height,
                                                                  gp = gpar(fill = mutation_color["Frame_Shift_Del"], col = NA)),
                 Translation_Start_Site = function(x, y, w, h) grid.rect(x, y, w*width, h*height,
                                                                         gp = gpar(fill = mutation_color["Translation_Start_Site"], col = NA)),
                 In_Frame_Ins = function(x, y, w, h) grid.rect(x, y, w*width, h*height, 
                                                               gp = gpar(fill = mutation_color["In_Frame_Ins"], col = NA)),
                 Doubtful_call = function(x, y, w, h) grid.rect(x, y, w*width, h*0.3, 
                                                                gp = gpar(fill = mutation_color["Doubtful_call"], col = NA)),
                 Nonstop_Mutation = function(x, y, w, h) grid.rect(x, y, w*width, h*height, 
                                                                   gp = gpar(fill = mutation_color["Nonstop_Mutation"], col = NA)),
                 Amplified_expression = function(x, y, w, h) grid.rect(x, y, w*width, h*0.3, 
                                                                       gp = gpar(fill = mutation_color["Amplified_expression"], col = NA)),
                 Mild_amp = function(x, y, w, h) grid.rect(x, y, w*width, h*0.3,
                                                           gp = gpar(fill = mutation_color["Mild_amp"], col = NA)),
                 Multi_hit = function(x, y, w, h) grid.rect(x, y, w*width, h*height,
                                                            gp = gpar(fill = mutation_color["Multi_hit"], col = NA)))

###test color
test_alter_fun(alter_fun)
# clinical_final_smoking <- clinical_smoking %>% 
#   mutate(Smoking_status = ifelse((final_smoking == 'Current smoker')|(final_smoking == 'Former smoker'),'Smoking experience','Never smoking'))
clinical_smoking <- clinical[!(is.na(total_clinical$final_smoking)),]
clinical_final_smoking <- total_clinical %>% 
  mutate(Smoking_status = ifelse((final_smoking == 'Current smoker')|(final_smoking == 'Former smoker'),'Smoking experience','Never smoking'))

total_clinical$Tumor_Sample_Barcode
total_clinical <- total_clinical %>% 
  mutate(Smoking_status = ifelse((final_smoking == 'Current smoker')|(final_smoking == 'Former smoker'),'Smoking experience','Never smoking'))
clinical_age <- total_clinical %>% filter(!is.na(final_age))
clinical_age_40_below <- clinical_age %>% filter(final_age <=40)
clinical_oral <-  total_clinical %>% filter(final_site == "Lip or oral cavity") 
clinical_oral$Tumor_Sample_Barcode
clinical_age_40_below_oral <- clinical_age_40_below %>% filter(final_site == "Lip or oral cavity")
clinical_age_40_below_oral$Tumor_Sample_Barcode
clinical_age_40_below_oral <- clinical_age_40_below_oral[!(is.na(clinical_age_40_below_oral$Smoking_status)),]
clinical_oral <- clinical_oral[!(is.na(clinical_oral$Smoking_status)),]
# table(clinical_oral$Smoking_status, clinical_age_40_below_oral$Smoking_status)
summary(as.factor(clinical_oral$Smoking_status))
summary(as.factor(clinical_oral$final_HPV))
summary(as.factor(clinical_age_40_below_oral$Smoking_status))



young <- clinical_age_40_below$Tumor_Sample_Barcode
clinical_oral$low_age <- ifelse(clinical_oral$final_age<=40, "young", "old")
table(clinical_oral$Smoking_status, clinical_oral$low_age)
###fig for smoking
library(ggplot2)
clinical_oral$low_age <- as.factor(clinical_oral$low_age)
clinical_oral$Smoking_status <- as.factor(clinical_oral$Smoking_status)
g <- ggplot(clinical_oral, aes(x = low_age, fill = Smoking_status), na.rm = TRUE) + geom_bar(position = "dodge") + ylab("The number of patients")+ theme_bw() +
  scale_x_discrete(limits=c('young','old')) +
  # scale_fill_grey(start = 0.8, end = 0.4)+
  scale_fill_nejm(alpha = 0.75)+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
# pdf("/data/project/TRIUMPH/temp/smoking_young.pdf", width = 6, height = 6)
# pdf("/data/project/TRIUMPH/Figure_swkim/Fig5C.pdf", width = 6, height = 6)
g
# dev.off()
# fisher.test(table(clinical_oral$Smoking_status, clinical_oral$low_age))


# Display the plot


g2 <- ggplot(clinical_oral, aes(x = low_age, fill = final_HPV)) + geom_bar(position = "dodge") + ylab("The number of patients")+ theme_bw() +
  scale_x_discrete(limits=c('young','old')) +
  # scale_fill_grey(start = 0.8, end = 0.4)+
  scale_fill_nejm(alpha = 0.75)+
  scale_y_continuous(breaks=seq(0, 50, 10), limits=c(0, 50)) +  # Adjusting y-axis
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text( hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Fig5E.pdf", width = 10, height = 10)
g2
dev.off()
clinical_oral[clinical_oral$low_age == 'young',]$final_HPV

clinical_age_40_below_oral$final_HPV

clinical_age_40_below_oral$final_smoking
clinical_age_40_below_oral
maf_young_oral = subsetMaf(maf = All_cnvkit,tsb = clinical_age_40_below_oral$Tumor_Sample_Barcode)
maf_young_oral = subsetMaf(maf = All_cnvkit_vaf,tsb = clinical_age_40_below_oral$Tumor_Sample_Barcode)
data_vaf = All_cnvkit_vaf@data
mean(data_vaf[data_vaf$Hugo_Symbol == 'TP53']$i_TumorVAF_WU, na.rm = T)

data_oral = maf_young_oral@data
data_oral[data_oral$Hugo_Symbol == 'TP53']$i_TumorVAF_WU
mean(data_oral[data_oral$Hugo_Symbol == 'TP53']$i_TumorVAF_WU)
# oropharynx_clinical <- inter_clinical[inter_clinical$final_site == "oropharynx",]



maf_young = subsetMaf(maf = All_cnvkit,tsb = young)
oncoplot(maf=maf_young, top=10, clinicalFeatures = c('final_HPV'))
# pdf("temp_oral.pdf", width = 20, height = 17)
# oncoplot(maf=maf_young_oral, top=30, clinicalFeatures = c('final_HPV', 'final_smoking'))
# dev.off()

setwd("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result")
oncoplot(maf = maf_young, top =10, writeMatrix = TRUE);young_top_30 = rownames(read.table(file = "/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/onco_matrix.txt", sep = "\t"))
file.remove('onco_matrix.txt')
getwd()

young_top_30
young <- intersect(young, colnames(mat_tot))
mat_young <- mat_tot[,young]
mat_young.top30 <- mat_young[young_top_30,]
dim(mat_young.top30)


ord <- sample_order(mat_young.top30)
# ord %in% colnames(mat_young.top30)
clinical_age_40_below <- clinical_age_40_below %>% filter((Tumor_Sample_Barcode %in% young))
clinical_age_40_below.sort <- clinical_age_40_below[ord,]

sample_order3 <- function(mat){
  clinical_age_40_below.site.sort <- clinical_age_40_below.sort[order(clinical_age_40_below.sort$final_site),]
  o = clinical_age_40_below.site.sort$Tumor_Sample_Barcode
  return(o)
}

o3 <- sample_order3(mat_young.top30)

HPV_col <- c("Positive"="#ff0000","Negative"="#0738e5")
site_col <- c("Hypopharynx"="#3b4992", "Larynx" = "#ee0000", "Lip or oral cavity" = "#008b45",
              "Maxillary sinus" = "#ff6600", "Nasal cavity" = "#ffff33", "oropharynx" = "#631879")
Smoking_col <- c("Current smoker" = "#003865", "Former smoker" = "#D61C4E", "Never smoker" = "#3330E4")

# setwd("/data/project/TRIUMPH/Figure2/new_figure_manu/")
# pdf("Figure6_age40_below2.pdf", width = 20, height = 17)
# pdf("/data/project/TRIUMPH/Figure_swkim/Figure5A.pdf", width = 20, height = 17)
size = 1
fig6_age40 <- oncoPrint(mat_young.top30, alter_fun = alter_fun, col = mutation_color, remove_empty_columns = F,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(
            height = unit(1, "cm")),
            HPV = clinical_age_40_below$final_HPV,
            Smoking = clinical_age_40_below$final_smoking,
            Site = clinical_age_40_below$final_site,
            annotation_name_side = 'left',
            col = list(HPV = HPV_col,
                       Site = site_col,
                       Smoking = Smoking_col),
            simple_anno_size = unit(1, "cm"),
            show_legend = F, 
            annotation_name_gp = gpar(fontsize = 30)),
          right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(
            width = unit(4, "cm"))),
          row_names_side = "left", pct_side =  "right", column_order = o3,
          show_heatmap_legend = F,
          column_title_gp = gpar(fontsize = 30),
          row_names_gp = gpar(fontsize = 30))

lgd <- list(Legend(labels = names(alter_fun)[2:16], title = "Alteration", labels_gp = gpar(fontsize = 24), 
                   title_gp = gpar(fontsize = 24), grid_height = unit(size,"cm"), grid_width = unit(size,"cm"),
                   graphics = alter_fun[2:16]),
            Legend(labels = c("Positive", "Negative"), title = "p16", labels_gp = gpar(fontsize = 24), 
                   title_gp = gpar(fontsize = 24), legend_gp = gpar(fill = HPV_col), grid_height = unit(size, "cm"),
                   grid_width = unit(size, "cm")),
            Legend(labels = names(site_col), title = "Site", labels_gp = gpar(fontsize = 24), 
                   title_gp = gpar(fontsize = 24), legend_gp = gpar(fill = site_col), 
                   grid_height = unit(size, "cm"), grid_width = unit(size, "cm")),
            Legend(labels = names(Smoking_col), title = "Smoking", labels_gp = gpar(fontsize = 24), 
                   title_gp = gpar(fontsize = 24), legend_gp = gpar(fill = Smoking_col), 
                   grid_height = unit(size, "cm"), grid_width = unit(size, "cm"))
)
draw(fig6_age40, padding = unit(c(2, 2, 7, 2), "mm"), heatmap_legend_list = lgd, heatmap_legend_side = "right")
# dev.off()
