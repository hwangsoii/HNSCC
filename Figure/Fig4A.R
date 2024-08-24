source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")


#######color change####### 
mutation_type <- unique(unlist(strsplit(unlist(strsplit(unlist(as.list(mat_tot)),";")),'\\|')))
mutation_type <- mutation_type[mutation_type!=""]
# mutation_color <- as.vector(qualitative_hcl(n=7, h = c(0,270), c = 80, l=80))
# mutation_color2 <- brewer.pal(7, "Paired")
# mutation_color <- c(mutation_color2, mutation_color)
mutation_color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#F8C8DC",
                    "#C6CE3E", "#5EDF82", "#00E4D0", "#00D7FF", "#C7BAFF", "#808B96")
mutation_color <- c("#91b5c9", "#1F78B4", "#a7d281", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#F8C8DC",
                    "#C6CE3E", "#5EDF82", "#00E4D0", "#00D7FF", "#C7BAFF", "#808B96")
# barplot(rep(1,length(mutation_type)), col = mutation_color, names.arg = names(mutation_color), cex.names = 0.5)
names(mutation_color) <- mutation_type[c(3,1,8,2,4,5,7,9,6,11,12,13,15,14,10)]


height = 0.97
width = 0.95
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
                 # Doubtful_call = function(x, y, w, h) grid.rect(x, y, w*width, h*0.3, 
                 #                                                gp = gpar(fill = mutation_color["Doubtful_call"], col = NA)),
                 Nonstop_Mutation = function(x, y, w, h) grid.rect(x, y, w*width, h*height, 
                                                                   gp = gpar(fill = mutation_color["Nonstop_Mutation"], col = NA)),
                 # Amplified_expression = function(x, y, w, h) grid.rect(x, y, w*width, h*0.3, 
                 #                                                       gp = gpar(fill = mutation_color["Amplified_expression"], col = NA)),
                 # Mild_amp = function(x, y, w, h) grid.rect(x, y, w*width, h*0.3,
                 #                                           gp = gpar(fill = mutation_color["Mild_amp"], col = NA)),
                 Multi_hit = function(x, y, w, h) grid.rect(x, y, w*width, h*height,
                                                            gp = gpar(fill = mutation_color["Multi_hit"], col = NA)))

###test color
# test_alter_fun(alter_fun)




oropharynx_clinical <- inter_clinical[inter_clinical$final_site == "oropharynx",]


oropharynx_barcode <- oropharynx_clinical$Tumor_Sample_Barcode
maf_oropharynx = subsetMaf(maf = All_cnvkit,tsb = oropharynx_barcode)
oncoplot(maf = maf_oropharynx, top =30, writeMatrix = TRUE);oro_top_30 = rownames(read.table(file = "onco_matrix.txt", sep = "\t"))
file.remove('onco_matrix.txt')
dev.off()

mat_oro <- mat_tot[,oropharynx_clinical$Tumor_Sample_Barcode]

mat_oro.top30 <- mat_oro[oro_top_30[1:20],]

ord <- sample_order(mat_oro.top30)

oropharynx_clinical.sort <- oropharynx_clinical[ord,]

sample_order3 <- function(mat){
  oropharynx_clinical.HPV.sort <- oropharynx_clinical.sort[order(oropharynx_clinical.sort$final_HPV),]
  o = oropharynx_clinical.HPV.sort$Tumor_Sample_Barcode
  return(o)
}


HPV_col <- c("Positive"="#ff0000","Negative"="#0738e5")
site_col <- c("Hypopharynx"="#3b4992", "Larynx" = "#ee0000", "Lip or oral cavity" = "#008b45",
              "Maxillary sinus" = "#ff6600", "Nasal cavity" = "#ffff33", "oropharynx" = "#631879")
Smoking_col <- c("Current smoker" = "#003865", "Former smoker" = "#D61C4E", "Never smoker" = "#3330E4")
o3 <- sample_order3(mat_oro.top30)

size = 1
fig6 <- oncoPrint(mat_oro.top30, alter_fun = alter_fun, col = mutation_color, remove_empty_columns = F,
                  top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(
                    height = unit(1, "cm")),
                    Smoking = oropharynx_clinical$final_smoking,
                    p16 = oropharynx_clinical$final_HPV,
                    col = list(p16 = HPV_col,
                               Smoking = Smoking_col),
                    annotation_name_side = "left",
                    annotation_name_gp = gpar(fontsize = 24),
                    simple_anno_size = unit(1, "cm"),
                    show_legend = F),
                  right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(
                    width = unit(4, "cm"))),
                  row_names_side = "left", pct_side =  "right",
                  pct_gp = gpar(fontsize = 24),
                  column_split = oropharynx_clinical$final_HPV,
                  column_gap = unit(3, "mm"),
                  column_order = o3,
                  column_title_gp = gpar(fontsize = 24),
                  show_heatmap_legend =F,
                  row_names_gp = gpar(fontsize = 24))

lgd <- list(Legend(labels = names(alter_fun)[2:length(alter_fun)], title = "Alteration", labels_gp = gpar(fontsize = 24),
                   title_gp = gpar(fontsize = 24), grid_height = unit(size,"cm"),
                   grid_width = unit(size,"cm"), graphics = alter_fun[2:length(alter_fun)]),
            Legend(labels = c("Positive", "Negative"), title = "p16", labels_gp = gpar(fontsize = 24),
                   title_gp = gpar(fontsize = 24), legend_gp = gpar(fill = HPV_col),
                   grid_height = unit(size, "cm"), grid_width = unit(size, "cm")),
            Legend(labels = names(site_col), title = "Site", labels_gp = gpar(fontsize = 24),
                   title_gp = gpar(fontsize = 24), legend_gp = gpar(fill = site_col),
                   grid_height = unit(size, "cm"), grid_width = unit(size, "cm"))
)

# setwd("/data/project/TRIUMPH/Figure2/new_figure_manu/")
# pdf("Figure6_alter5.pdf", width = 35, height = 17)
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Fig4A.pdf", width = 35, height = 17)
draw(fig6, padding = unit(c(2, 2, 7, 400), "mm"), heatmap_legend_list = lgd, heatmap_legend_side = "right")
dev.off()
