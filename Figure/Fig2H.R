source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")

data <- All_cnvkit@data
data[data$Hugo_Symbol == 'MET']
data[data$Hugo_Symbol == 'IGF1R']
'NOTCH3' %in% data$Hugo_Symbol
'NRAS' %in% data$Hugo_Symbol
unique(data$Hugo_Symbol)[order(unique(data$Hugo_Symbol))]


mutation_type <- unique(unlist(strsplit(unlist(strsplit(unlist(as.list(mat_tot)),";")),'\\|')))
mutation_type <- mutation_type[mutation_type!=""]
# mutation_color <- as.vector(qualitative_hcl(n=7, h = c(0,270), c = 80, l=80))
# mutation_color2 <- brewer.pal(8, "Paired")
# mutation_color <- c(mutation_color2, mutation_color,"#000000")
mutation_color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#F8C8DC",
                    "#C6CE3E", "#5EDF82", "#00E4D0", "#00D7FF", "#C7BAFF", "#808B96")
barplot(rep(1,length(mutation_type)), col = mutation_color, names.arg = names(mutation_color), cex.names = 0.5)
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

test_alter_fun(alter_fun)


# test_alter_fun(alter_fun)
# Step 1: Create a new matrix to duplicate samples for each gene group they have mutations in
mutation_matrix <- mat_tot # create mutation matrix

genes_arm1 <- c("PIK3CA", "PIK3CB", "PIK3C2A", "PIK3C3", "PIK3CD", "PIK3R1", "PIK3R2","PIK3R3",  "PTEN", "NOTCH1", "MYC") ##PI3K
# genes_arm1 <- c("PIK3CA", "PIK3CB", "PIK3C2A", "PIK3C3", "PIK3CD", "PIK3R1", "PIK3R2","PIK3R3",  "PTEN") ##PI3K
genes_arm2 <- c("EGFR", "ERBB2", "ERBB3","ERBB4", "MET", "PDGFRA") #EGFR
genes_arm3 <- c("FGFR1", "FGFR2","FGFR3", "FGFR4", "KIT", "IGF1R") ##FGFR
genes_arm4 <- c("CDKN2A", "CDKN2B", "CDKN2C", "CDKN1B", "CCND1", "CCND2", "CCND3","CCNE1") ##Cell cycle
genes_arm5 <- c("AKT1","AKT2", "AKT3","TSC1","TSC2","MTOR","RICTOR","RPTOR","PPP2R1A","STK11") ##MTOR
genes_arm6 <- c("NOTCH1","NOTCH2", "FBXW7", "CREBBP", "EP300","KDM5A") #NOTCH patwhay
genes_arm7 <- c("KRAS","HRAS", "NRAS", "RIT1","ARAF","BRAF", "RAF1","RAC1","MAPK1","MAP2K1", "MAP2K2" ) #RAS-RAF pathway
genes_arm8 <- c("TP53", "ATM", "CHEK2", "RPS6KA1") ## p53 pathway
genes_arm9 <- c("MYC", "MYCN", "MYCL1") ##  Myc pathway
genes_arm10 <- c( "FAT1","FAT2", "FAT3","FAT4", "NF2") ##Hippo pathway
genes_arm11 <- c("KEAP1", "CUL3", "NFE2L2") ##Nrf2 pathway
genes_arm12 <- c("TGFBR2","APC", "CTNNB1") ##WNT TGFB pathway





gene_groups <- list(genes_arm1, genes_arm2, genes_arm3, genes_arm4, genes_arm5, genes_arm6, genes_arm7, genes_arm8, genes_arm9, genes_arm10, genes_arm11, genes_arm12)
genes_arm <- c(genes_arm1, genes_arm2, genes_arm3, genes_arm4, genes_arm5, genes_arm6, genes_arm7, genes_arm7, genes_arm8, genes_arm9, genes_arm10, genes_arm11, genes_arm12)
# Step 2: Create sub-matrices for each gene group
sub_matrices <- list()


# # Step 1: Get the gene group of interest (gene_arm1)
# arm1_matrix <- find_gmat(genes_arm1)
# arm2_matrix <- find_gmat(genes_arm2)
# arm3_matrix <- find_gmat(genes_arm3)
# arm4_matrix <- find_gmat(genes_arm4)
# arm5_matrix <- find_gmat(genes_arm5)

# matrix_dup <- cbind(arm1_matrix, arm2_matrix, arm3_matrix, arm4_matrix, arm5_matrix)

sort(rownames(mat_tot))

find_gmat <- function(g_arm){
  # Step 2: Find the sample indices that have mutations in the gene group
  gene_group_indices <- match(g_arm, rownames(mat_tot))
  sample_indices <- which(colSums(mat_tot[gene_group_indices, ] != "") > 0)
  
  # Step 3: Create a sub-matrix of the original mat_tot for samples with mutations in the gene group
  sub_matrix <- mat_tot[, sample_indices]
  return(sub_matrix)
}


find_gmat2 <- function(g_arm){
  # Step 2: Find the sample indices that have mutations in the gene group
  sub_matrix <- NULL
  for (gen in g_arm) {
    gene_group_indices <- match(gen, rownames(mat_tot))
    sample_indices <- which(mat_tot[gene_group_indices, ] != "")
    sub_matrix <- cbind(sub_matrix, mat_tot[, sample_indices])
    
  }
  
  # Step 3: Create a sub-matrix of the original mat_tot for samples with mutations in the gene group
  return(sub_matrix)
}



# tm <- mat_tot[gene_group_indices, ]
# Step 1: Get the gene group of interest (gene_arm1)
arm1_matrix <- find_gmat(genes_arm1)
# arm2_matrix <- find_gmat(genes_arm2)
# arm3_matrix <- find_gmat(genes_arm3)
# arm4_matrix <- find_gmat(genes_arm4)
# arm5_matrix <- find_gmat(genes_arm5)

# mat_tmp <- arm1_matrix[genes_arm,]
# mat_tmp <- prioritize_genes(genes_arm1, mat_tmp)
# for (i in rev(seq_along(genes_arm1))){
#   # mat_tmp <- mat_tmp[,order(mat_tmp[i,] != "")]
#   print(i)
#   samples_with_mutation <- which(mat_tmp[i,] != "")
#   mat_tmp <- cbind(mat_tmp[, samples_with_mutation], mat_tmp[, -samples_with_mutation])
# }
# 
# mat_tmp <- arm5_matrix[genes_arm,]
# mat_tmp <- prioritize_genes(genes_arm5, mat_tmp)
# print(dim(mat_tmp))
# for (i in rev(seq_along(genes_arm5))){
#   # mat_tmp <- mat_tmp[,order(mat_tmp[i,] != "")]
#   print(i)
#   samples_with_mutation <- which(mat_tmp[genes_arm5[i],] != "")
#   mat_tmp <- cbind(mat_tmp[, samples_with_mutation], mat_tmp[, -samples_with_mutation])
#   print(dim(mat_tmp))
# }
# samples_with_mutation <- which(mat_tmp[1,] != "")
# samples_with_mutation <- which(mat_tmp[5,] != "")
# dim(mat_tmp[, samples_with_mutation])
# dim(mat_tmp[, -samples_with_mutation])
# 
# mat_tmp <- cbind(mat_tmp[, samples_with_mutation], mat_tmp[, -samples_with_mutation])
# 
# mat_tmp <- prioritize_genes(genes_arm1, arm1_matrix)
# # 
# mat_tmp2 <- prioritize_genes(genes_arm2, arm2_matrix)
# mat_tmp2 <- prioritize_genes(genes_arm2, arm2_matrix)
# mat_tmp2 <- prioritize_genes(genes_arm2, arm2_matrix)
# mat_tmp2 <- prioritize_genes(genes_arm5, arm5_matrix)
# # 
# matrix_sorted <- mat_tmp[, order(mat_tmp[1,] == "")]
# dim(matrix_dup)

prioritize_genes <- function(gene_order, mutation_matrix) {
  # Subset the mutation matrix based on gene order
  tmp_matrix <- mutation_matrix[genes_arm, ]
  
  # gene_sums <- colSums(mutation_matrix[gene_order,])
  
  # Get the sum of mutations for each sample in the sub-matrix
  sum_mutations <- colSums(mutation_matrix[gene_order,] != "")
  
  # Get the indices of samples with at least one mutation
  # samples_with_mutation <- which(sum_mutations > 0)
  
  # # Get the number of mutations for each gene group
  # num_mutations_group <- rowSums(tmp_matrix != "")
  # 
  # # Order the gene groups based on the number of mutations
  group_order <- order(sum_mutations, decreasing = TRUE)
  
  # Rearrange the rows of the sub-matrix based on group order
  tmp_matrix <- tmp_matrix[,group_order ]
  
  # Move the samples with at least one mutation to the front of the sub-matrix
  # tmp_matrix <- cbind(tmp_matrix[, samples_with_mutation], tmp_matrix[, -samples_with_mutation])
  
  # # Rearrange the columns of the sub-matrix based on gene order
  # new_tmp_matrix <- tmp_matrix[, gene_order]
  
  # Return the rearranged sub-matrix
  return(tmp_matrix)
}
# 
# myfunc <- function(){
#   mat_sorted <- NULL
#   # mat_tmp_list <- list()
#   # ncol(mat_sorted)
#   # ncol(find_gmat(genes_arm1)[genes_arm,])
#   # ncol(mat_tmp)
#   for (group in gene_groups){
#     # print(group)
#     mat_tmp <- find_gmat(group)[genes_arm,]
#     # mat_tmp <- prioritize_genes(group, find_gmat(group)[genes_arm,])
#     # mat_tmp <- find_gmat(group)[genes_arm,]
#     # print(dim(mat_tmp))
#     mat_tmp <- prioritize_genes(group, mat_tmp)
#     for (i in rev(seq_along(group))){
#       # mat_tmp <- mat_tmp[,order(mat_tmp[i,] != "")]
#       samples_with_mutation <- which(mat_tmp[group[i],] != "")
#       mat_tmp <- cbind(mat_tmp[, samples_with_mutation], mat_tmp[, -samples_with_mutation])
#       # mat_tmp <- mat_tmp[,order(mat_tmp[i,] != "",sum_mutations, decreasing = TRUE)]
#     }
#     # sum_mutations <- colSums(mat_tmp[group,] != "")
#     # print(dim(mat_tmp))
#     # print(dim(mat_sorted))
#     mat_sorted  <- cbind(mat_sorted, mat_tmp)
#     # print(dim(mat_sorted))
#     # mat_tmp_list[[length(mat_tmp_list) + 1]] <- mat_tmp[, colSums(mat_tmp != "") > 0] # Only add columns with at least one non-empty value
#   }
#   return(mat_sorted)
# }
# 
# mat_sorted <- NULL
# # mat_tmp_list <- list()
# # ncol(mat_sorted)
# # ncol(find_gmat(genes_arm1)[genes_arm,])
# # ncol(mat_tmp)
# for (group in gene_groups){
#   # print(group)
#   mat_tmp <- find_gmat(group)[genes_arm,]
#   # mat_tmp <- prioritize_genes(group, find_gmat(group)[genes_arm,])
#   # mat_tmp <- find_gmat(group)[genes_arm,]
#   # print(dim(mat_tmp))
#   mat_tmp <- prioritize_genes(group, mat_tmp)
#   for (i in rev(seq_along(group))){
#     # mat_tmp <- mat_tmp[,order(mat_tmp[i,] != "")]
#     samples_with_mutation <- which(mat_tmp[group[i],] != "")
#     mat_tmp <- cbind(mat_tmp[, samples_with_mutation], mat_tmp[, -samples_with_mutation])
#     # mat_tmp <- mat_tmp[,order(mat_tmp[i,] != "",sum_mutations, decreasing = TRUE)]
#   }
#   # sum_mutations <- colSums(mat_tmp[group,] != "")
#   # print(dim(mat_tmp))
#   # print(dim(mat_sorted))
#   mat_sorted  <- cbind(mat_sorted, mat_tmp)
#   # print(dim(mat_sorted))
#   # mat_tmp_list[[length(mat_tmp_list) + 1]] <- mat_tmp[, colSums(mat_tmp != "") > 0] # Only add columns with at least one non-empty value
# }
# # mat_sorted <- do.call(cbind, mat_tmp_list)


# mat_sorted <- NULL
genes_arm <- genes_arm2
# genes_arm_2 <- genes_arm2_2
mat_tmp <- find_gmat(genes_arm)[genes_arm,]
# mat_tmp <- find_gmat2(genes_arm)[genes_arm,]
mat_tmp <- prioritize_genes(genes_arm, mat_tmp)

for (i in rev(seq_along(genes_arm))){
  # mat_tmp <- mat_tmp[,order(mat_tmp[i,] != "")]
  samples_with_mutation <- which(mat_tmp[genes_arm[i],] != "")
  mat_tmp <- cbind(mat_tmp[, samples_with_mutation], mat_tmp[, -samples_with_mutation])
  # mat_tmp <- mat_tmp[,order(mat_tmp[i,] != "",sum_mutations, decreasing = TRUE)]
}
# sum_mutations <- colSums(mat_tmp[group,] != "")
# print(dim(mat_tmp))
# print(dim(mat_sorted))
# mat_sorted  <- cbind(mat_sorted, mat_tmp)
# ############################duplicated
# for (i in rev(seq_along(genes_arm))){
#   # mat_tmp <- mat_tmp[,order(mat_tmp[i,] != "")]
#   
#   samples_with_mutation <- which(mat_tmp[genes_arm[i],] != "")
#   mat_tmp <- cbind(mat_tmp[, samples_with_mutation], mat_tmp[, -samples_with_mutation])
#   # mat_tmp <- mat_tmp[,order(mat_tmp[i,] != "",sum_mutations, decreasing = TRUE)]
# }
# colnames(mat_tmp[,mat_tmp[genes_arm1[5],]!=""])
# 
# samples_with_mutation <- which(mat_tmp[genes_arm[5],] != "")
# dup_indices <- which(duplicated(colnames(mat_tmp), fromLast = FALSE))
# (samples_with_mutation %in% dup_indices)
# 
# first_indices <- unique(dup_indices[, 1], byrow = FALSE)
# 
# colnames(mat_sorted[,duplicated(colnames(mat_sorted))])
# which((mat_tmp[genes_arm1[i],] != "")&duplicated(colnames(mat_tmp)))
# 
# 
# tt <-  which(mat_tmp[genes_arm[i],] != "")
# unique(names(tt))
# #########################################################



# which(duplicated(colnames(mat_sorted))) %in% which(unique(mat_sorted[,duplicated(colnames(mat_sorted))],fromLast = TRUE))


mat_sorted <- mat_tmp
vacant_mat <- matrix("", nrow=7-nrow(mat_sorted), ncol = ncol(mat_sorted))
rownames(vacant_mat) <- rep("vacant", 7-nrow(mat_sorted))
mat_sorted <- rbind(mat_sorted, vacant_mat)



# 
# 
# 
# arm_info[2]
# gene_list['Arm 2,4']
# math = cumprod(rep(0.1,7))
# 
# arm_info$count <- apply(arm_info, 1,
#                         function(x) sum(mat_arm.sort[gene_list[x[2]][[1]],x[1]]!=""))
# arm_info$count2 <- apply(arm_info, 1,
#                          function(x) mat_arm.sort[gene_list[x[2]][[1]],x[1]]!="")
# arm_info$count3 <- apply(arm_info, 1,
#                          function(x) sum(math[1:length(x[4][[1]])][x[4][[1]]]))
# arm_info.sort <- arm_info[order(arm_info$count3, decreasing = T),]
# arm_info.sort2 <- arm_info.sort[order(arm_info.sort[,2]),]
# mat_arm.sort2 <- mat_arm.sort[,c(arm_info.sort2[,1])]
# level <- c(arm_info.sort2[,1])
# level <- c(colnames(mat_sorted))

# 
# Arm_dup_split <- c(rep('Arm 1', ncol(arm1_matrix)), rep('Arm 2', ncol(arm2_matrix)), rep('Arm 3', ncol(arm3_matrix)),rep('Arm 4', ncol(arm4_matrix)),rep('Arm 5', ncol(arm5_matrix)))
# Arm_dup_split <- c(rep(1, ncol(arm1_matrix)), rep(2, ncol(arm2_matrix)), rep(3, ncol(arm3_matrix)),rep(4, ncol(arm4_matrix)),rep(5, ncol(arm5_matrix)))
# length(Arm_dup_split)
# ncol(mat_sorted)
# ncol(arm1_matrix)+ncol(arm2_matrix) +ncol(arm3_matrix) + ncol(arm4_matrix) + ncol(arm5_matrix)
# colnames(mat_arm.sort2) == arm_info.sort2[,1]
# 
# column_gap_v <- rep(0.1,ncol(mat_sorted)-1)
# column_gap_v[ncol(arm1_matrix)] <- 3
# column_gap_v[ncol(arm2_matrix)] <- 3
# column_gap_v[ncol(arm3_matrix)] <- 3
# column_gap_v[ncol(arm4_matrix)] <- 3
# column_gap_v[ncol(arm5_matrix)] <- 3



inter_clinical.arm <- inter_clinical

sample_names <- inter_clinical.arm[["Tumor_Sample_Barcode"]]
# Use match to get the row indices in mat_sorted that match the sample names
col_indices <- match(colnames(mat_sorted),inter_clinical.arm$Tumor_Sample_Barcode)

# Subset mat_sorted using the row indices
inter_clinical.arm <- inter_clinical.arm[col_indices, ]

colnames(mat_sorted[,duplicated(colnames(mat_sorted))])
mat_sorted[,duplicated(colnames(mat_sorted))]
colnames(mat_sorted[,duplicated(colnames(mat_sorted))])


# inter_clinical.arm$Tumor_Sample_Barcode <- factor(inter_clinical.arm$Tumor_Sample_Barcode, levels = level)
# inter_clinical.arm.sort <- inter_clinical.arm[order(inter_clinical.arm$Tumor_Sample_Barcode,decreasing = FALSE),]

# inter_clinical.arm$Tumor_Sample_Barcode <- factor(inter_clinical.arm$Tumor_Sample_Barcode, levels = level)
# inter_clinical.arm.sort <- inter_clinical.arm[order(inter_clinical.arm$Tumor_Sample_Barcode,decreasing = FALSE),]
H15_inter$TRTP1 <- unlist(H15_inter$TRTP1)
clinical_arms <- H15_inter[H15_inter$TRTP1=="Arm 2 (Poziotinib)",]
ARMS_barcode <- unlist(clinical_arms$Tumor_Sample_Barcode.1)
inter_clinical.arm$TRIUMPH <- ifelse(inter_clinical.arm$Tumor_Sample_Barcode %in% ARMS_barcode, "Yes", "No")
HPV_df <- inter_clinical.arm$final_HPV
site_df <- inter_clinical.arm$final_site
TRIUMPH_df <- inter_clinical.arm$TRIUMPH
# HPV_df <- inter_clinical.arm.sort$final_HPV
# site_df <- inter_clinical.arm.sort$final_site
HPV_col <- c("Positive"="#ff0000","Negative"="#0738e5")
site_col <- c("Hypopharynx"="#3b4992", "Larynx" = "#ee0000", "Lip or oral cavity" = "#008b45",
              "Maxillary sinus" = "#ff6600", "Nasal cavity" = "#ffff33", "oropharynx" = "#631879")
TRIUMPH_col <-  c("Yes" = "#336699", "No" = "#CCCCCC") 


# oncoPrint(mat_arm.sort2, alter_fun = alter_fun, col = mutation_color, remove_empty_columns = F,
#           top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(
#             height = unit(1, "cm")),
#             HPV = HPV_df,
#             Site = site_df,
#             Arm = inter_clinical.arm.sort$final_arm
#           ),
#           right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(
#             width = unit(4, "cm"))),
#           column_order = arm_info.sort2[,1],
#           row_names_side = "left", pct_side =  "right",
#           row_order = rownames(mat_arm.sort),
#           row_split = c(rep("Arm1",7), rep("Arm2",4), rep("Arm3",4), rep("Arm4",5), rep("Arm5",5)),
#           row_title = c('PIK3CA pathway', 'EGFR pathway', 'FGFR pathway', 'Cell cycle pathway', 'Others'),
#           row_gap = unit(3, "mm"),
#           column_split = arm_info.sort2[,2],
#           column_gap = unit(3, "mm"),
#           show_column_names = F)
########################################
# create a vector to keep track of duplicate sample names

# mat_sorted <- myfunc()
# dup_count <- rep(1, ncol(mat_sorted))
# duplicated(colnames(mat_sorted))
# # loop through each column of the matrix
# for (i in seq_along(colnames(mat_sorted))) {
#   # check if the column name is a duplicate
#   if (duplicated(colnames(mat_sorted))[i]==T) {
#     # get the index of the column name
#     index <- which(colnames(mat_sorted) == colnames(mat_sorted)[i])
#     # rename the duplicate column names with unique names
#     print(i)
#     print(index)
#     colnames(mat_sorted)[index] <- paste0(colnames(mat_sorted)[index], "_dup", dup_count[index])
#     # update the dup_count vector for the duplicate column names
#     dup_count[index] <- dup_count[index] + 1
#   }
# }
# 
# 
# which(colnames(mat_sorted) == colnames(mat_sorted)[517])
# dup_count[517]
# 


# nrow(mat_sorted)

nrow(mat_sorted)
sum(mutation_matrix['EGFR',] != "")
ncol(mat_sorted)

sum(mutation_matrix['EGFR',] != "")/ncol(mat_sorted)
sum(mutation_matrix['EGFR',] != "")/419
ncol(mat_sorted)/419


###
ratio_finder <-  function(gene1, gene2, mat){
  tmp_mat <- mat[,mat[gene1,] != '']
  con_mat <- tmp_mat[,tmp_mat[gene2, ] != '']
  ratio = ncol(con_mat)/ ncol(tmp_mat)
  ans <-  c(ncol(con_mat), ncol(tmp_mat), ratio)
  return(ans)
}
#### arm2 ratio
ratio_finder('EGFR', 'KRAS', find_gmat(genes_arm))

tmp_mat['EGFR',]
tmp_mat <- find_gmat(genes_arm)[,find_gmat(genes_arm)['EGFR',] != '']
con_mat <- tmp_mat[,tmp_mat['KRAS', ] != '']
tmp_mat <- mat_sorted[,mat_sorted['EGFR',] != '']
con_mat <- tmp_mat[,(tmp_mat['ERBB2', ] != ''|tmp_mat['ERBB3', ] != ''|tmp_mat['ERBB4', ] != '')]
genes_arm1
print(ncol(con_mat))
print(ncol(tmp_mat))

ncol(con_mat)/ ncol(tmp_mat)





###########################
fig1_2 <- oncoPrint(mat_sorted, alter_fun = alter_fun, col = mutation_color, remove_empty_columns = F,
                    top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(
                      height = unit(1, "cm")),
                      p16 = HPV_df,
                      Site = site_df,
                      TRIUMPH = TRIUMPH_df,
                      col = list(p16 = HPV_col,
                                 Site = site_col, TRIUMPH = TRIUMPH_col),
                      simple_anno_size = unit(1, "cm"),
                      show_legend = F,
                      annotation_name_gp = gpar(fontsize = 14)
                    ),
                    right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(
                      width = unit(3, "cm")),
                      annotation_name_gp = gpar(fontsize = 20)),
                    column_order = colnames(mat_sorted),
                    row_names_side = "left", pct_side =  "right",
                    pct_gp = gpar(fontsize = 25),
                    row_order = rownames(mat_sorted),
                    # row_split = c(rep("Arm1",7), rep("Arm2",4), rep("Arm3",4), rep("Arm4",5), rep("Arm5",5)),
                    row_title = c('PI3kinase pathway'),
                    # row_title = c('PIK3CA pathway', 'EGFR pathway', 'FGFR pathway', 'Cell cycle pathway', 'Others'),
                    row_names_gp = gpar(fontsize = 25),
                    # row_gap = unit(3, "mm"),
                    # column_split = arm_info.sort2[,2],
                    # column_gap = unit(3, "mm"),
                    # show_column_names = T,
                    column_title = c("PI3Kinase pathway inhibitor"),
                    column_title_gp = gpar(fontsize = 40),
                    show_heatmap_legend =F
)
size = 1
lgd <- list(Legend(labels = names(alter_fun)[2:length(alter_fun)], title = "Alteration", labels_gp = gpar(fontsize = 14),  title_gp = gpar(fontsize = 14),
                   grid_height = unit(size,"cm"), grid_width = unit(size,"cm"), graphics = alter_fun[2:length(alter_fun)]),
            Legend(labels = c("Positive", "Negative"), title = "p16", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14),
                   legend_gp = gpar(fill = HPV_col), grid_height = unit(size, "cm"), grid_width = unit(size, "cm")),
            Legend(labels = names(site_col), title = "Site", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14),
                   legend_gp = gpar(fill = site_col), grid_height = unit(size, "cm"), grid_width = unit(size, "cm")),
            Legend(labels = c("Yes", "No"), title = "TRIUMPH", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14),
                   legend_gp = gpar(fill = TRIUMPH_col), grid_height = unit(size, "cm"), grid_width = unit(size, "cm")))


# setwd("/data/project/TRIUMPH/Figure2/")
# setwd("/data/project/TRIUMPH/Figure_dup/")
# setwd("/data/project/TRIUMPH/Figure_swkim//")
# pdf("Figure1_dup2.pdf", width = 35, height = 15)
# pdf("Figure1_dup_gap2.pdf", width = 35, height = 35)
# pdf("Figure1_dup_gap3.pdf", width = 500, height = 200)
# pdf("Figure1_ARM2.pdf", width = 35, height = 20)
# pdf("Figure1_ARM2_dup.pdf", width = 35, height = 20)
# pdf("Figure1_ARM2_2.pdf", width = 35, height = 7)
# pdf("/data/project/TRIUMPH/Figure4/figure_Arm/Figure1_ARM2.pdf", width = 35, height = 6.5)
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Figure2H.pdf", width = 35, height = 10)
# draw(fig1, padding = unit(c(2, 2, 10, 2), "mm"), heatmap_legend_list = lgd, heatmap_legend_side = "right")
draw(fig1_2, padding = unit(c(2, 2, 10, 2), "mm"), heatmap_legend_list = lgd, heatmap_legend_side = "right")
dev.off()

###missing patients check
check <- arm_info.sort2[arm_info.sort2[,3]==0,]
check_sample <- check[check[,2]!="Arm 5",1]





