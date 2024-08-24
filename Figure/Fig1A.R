source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16//Figure/import_data_for_Oncoprint_final.R")
BiocManager::install("qvalue")


#######color change####### 
mutation_type <- unique(unlist(strsplit(unlist(strsplit(unlist(as.list(mat_tot)),";")),'\\|')))
mutation_type <- mutation_type[mutation_type!=""]
mutation_color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#F8C8DC",
                    "#C6CE3E", "#5EDF82", "#00E4D0", "#00D7FF", "#C7BAFF", "#808B96")
# barplot(rep(1,length(mutation_type)), col = mutation_color, names.arg = names(mutation_color), cex.names = 0.5)
names(mutation_color) <- mutation_type[c(3,1,8,2,4,5,7,9,6,11,12,13,15,14,10)]

height = 0.9
width = 0.8
# background = function(x, y, w, h) grid.rect(x, y, w, h, 
#                                             gp = gpar(fill = NA)),
#background = function(...) NULL,
alter_fun = list(background = function(x, y, w, h) grid.rect(x, y, w*1, h*1, 
                                                             gp = gpar(fill = "#f2f2f2", col = NA)),
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

#gene sorting
o <- sample_order(mat30)


HPV_col <- c("Positive"="#ff0000","Negative"="#0738e5")
site_col <- c("Hypopharynx"="#3b4992", "Larynx" = "#ee0000", "Lip or oral cavity" = "#008b45",
              "Maxillary sinus" = "#ff6600", "Nasal cavity" = "#ffff33", "oropharynx" = "#631879")
Smoking_col <- c("Current smoker" = "#003865", "Former smoker" = "#D61C4E", "Never smoker" = "#3330E4")



size = 1
fig2 <- oncoPrint(mat30, alter_fun = alter_fun, col = mutation_color, remove_empty_columns = F,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(
            height = unit(1, "cm")),
            Smoking = inter_clinical$final_smoking,
            p16 = inter_clinical$final_HPV, 
            Site = inter_clinical$final_site,
            col = list(p16 = HPV_col,
                       Site = site_col,
                       Smoking = Smoking_col),
            show_legend = F,
            simple_anno_size = unit(1, "cm"),
            annotation_name_gp = gpar(fontsize = 24)
            ),
          right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(
            width = unit(4, "cm"))),
          row_names_side = "left", pct_side =  "right",
          pct_gp = gpar(fontsize = 24),
          row_order = top30_genes[1:20],
          column_order = o,
          show_heatmap_legend =F,
          column_title_gp = gpar(fontsize = 24),
          row_names_gp = gpar(fontsize = 24)
          )
lgd <- list(Legend(labels = names(alter_fun)[2:length(alter_fun)], title = "Alteration", labels_gp = gpar(fontsize = 20),
                   title_gp = gpar(fontsize = 20), grid_height = unit(size,"cm"),
                   grid_width = unit(size,"cm"), graphics = alter_fun[2:length(alter_fun)]),
            Legend(labels = c("Positive", "Negative"), title = "p16", labels_gp = gpar(fontsize = 20),
                   title_gp = gpar(fontsize = 20), legend_gp = gpar(fill = HPV_col),
                   grid_height = unit(size, "cm"), grid_width = unit(size, "cm")),
            Legend(labels = names(site_col), title = "Site", labels_gp = gpar(fontsize = 20),
                   title_gp = gpar(fontsize = 20), legend_gp = gpar(fill = site_col),
                   grid_height = unit(size, "cm"), grid_width = unit(size, "cm")),
            Legend(labels = names(Smoking_col), title = "Smoking", labels_gp = gpar(fontsize = 20),
                   title_gp = gpar(fontsize = 20), legend_gp = gpar(fill = Smoking_col), 
                   grid_height = unit(size, "cm"), grid_width = unit(size, "cm"))
)

length(alter_fun)

set.seed(9533)  #for random annotation color

pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Fig1A.pdf", width = 35, height = 15)
draw(fig2, padding = unit(c(2, 8, 7, 2), "mm"), heatmap_legend_list = lgd, heatmap_legend_side = "right")
dev.off()

