source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")


mutation_type <- unique(unlist(strsplit(unlist(strsplit(unlist(as.list(mat_tot)),";")),'\\|')))
mutation_type <- mutation_type[mutation_type!=""]

unwant <- mutation_type[(mutation_type!="Amplification" & mutation_type!='Deletion')]
mat_cnv <- mat_tot
for(i in unwant){
  mat_cnv <- gsub(i, "", mat_cnv)
  mat_cnv <- gsub(";","", mat_cnv)
  mat_cnv <- gsub("\\|","", mat_cnv)
}

# pal_ucscgb("default", alpha = 0.79)(26)[1]
# pal_ucscgb("default", alpha = 0.79)(26)[15]

Site_color <- c("Amplification"="#e64b35","Deletion"="#4dbbd5",
                "Hypopharynx"="#3b4992", "Larynx" = "#ee0000", "Lip or oral cavity" = "#008b45",
                "Maxillary sinus" = "#ff6600", "Nasal cavity" = "#ffff33", "oropharynx" = "#631879")

height = 0.97
width = 0.9


alter_fun_CNV = list(background = function(x, y, w, h) grid.rect(x, y, w, h, 
                                                                 gp = gpar(fill = "#f2f2f2", col = NA)),
                     Amplification = function(x, y, w, h) grid.rect(x, y, w*width, h*height, 
                                                                    gp = gpar(fill = Site_color["Amplification"], col = NA)),
                     Deletion =function(x, y, w, h) grid.rect(x, y, w*width, h*height,
                                                              gp = gpar(fill = Site_color["Deletion"], col = NA))
)


top_30 = order(apply(mat_cnv, 1, function(x) sum(x != "")), decreasing = TRUE)[1:30]
mat_cnv_30 = as.data.frame(mat_cnv[top_30,])


site_info <- data.frame(inter_clinical$Tumor_Sample_Barcode, inter_clinical$final_site)
val_list = cumprod(rep(0.1,30))
# arm_info$count2 <- apply(arm_info, 1,
#                          function(x) mat_arm.sort[gene_list[x[2]][[1]],x[1]]!="")
site_info$count2 <- apply(site_info, 1,
                          function(x) sum(val_list[mat_cnv_30[,x[[1]]]!=""]))

site_info.sort <- site_info[order(site_info$count2, decreasing = T),]
site_info.sort2 <- site_info.sort[order(site_info.sort[,2]),]


mat_cnv_30.sort = mat_cnv_30[,site_info.sort2[,1]]
size = 1

fig3 <- oncoPrint(mat_cnv_30.sort, alter_fun = alter_fun_CNV, col = Site_color, remove_empty_columns = F,
                  top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(
                    height = unit(1, "cm")),
                    `Anatomic site` = site_info.sort2$inter_clinical.final_site,
                    col = list(`Anatomic site` = Site_color), show_annotation_name = T,
                    annotation_name_side = "left", 
                    annotation_name_gp = gpar(fontsize = 24),
                    simple_anno_size = unit(1, "cm"),
                    show_legend = F),
                  right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(
                    width = unit(4, "cm"))),
                  row_names_side = "left", pct_side =  "right",
                  pct_gp = gpar(fontsize = 24),
                  column_order = site_info.sort2[,1],
                  show_heatmap_legend =F,
                  column_title_gp = gpar(fontsize = 24),
                  row_names_gp = gpar(fontsize = 24)
)

lgd <- list(Legend(labels = c("Amplification", "Deletion"), title = "Alteration",
                   labels_gp = gpar(fontsize = 24), grid_height = unit(size,"cm"),
                   grid_width = unit(size,"cm"), graphics = alter_fun_CNV[2:3],
                   title_gp = gpar(fontsize = 24)),
            Legend(labels = names(Site_color)[3:8], title = "Site", labels_gp = gpar(fontsize = 24),
                   legend_gp = gpar(fill = Site_color[3:8]), grid_height = unit(size, "cm"),
                   grid_width = unit(size, "cm"), title_gp = gpar(fontsize = 24))
)

# setwd("/data/project/TRIUMPH/Figure2/new_figure_manu/")

pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig4B.pdf",width = 35, height = 17)
draw(fig3, padding = unit(c(2, 2, 7, 2), "mm"), heatmap_legend_list = lgd, heatmap_legend_side = "right")
dev.off()
