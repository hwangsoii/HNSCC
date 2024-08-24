source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")
# #####signature
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)


library('NMF')

# install.packages("pheatmap")
library('pheatmap')

###(4.1, 2)
pheatmap::pheatmap(mat = HNSCC.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")


### (4.1, 2)
maftools::plotSignatures(nmfRes = HNSCC.sig, title_size = 1.2, sig_db = "SBS")

install.packages("reticulate")
library("reticulate")
use_python("path_to_your_python3")
py_config()
library("devtools")
install_github("AlexandrovLab/SigProfilerMatrixGeneratorR")


#####24.05.14
# Load necessary libraries
library(ggplot2)
library(pheatmap)
library(reshape2)



# Melt the signatures data for ggplot
signatures_melted <- melt(signatures)
names(signatures_melted) <- c("Sample", "Signature", "Exposure")

# Plot signatures using ggplot2
ggplot(signatures_melted, aes(x = Signature, y = Exposure, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Signature Exposures per Sample", y = "Relative Exposure", x = "Signature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


####
# Load necessary library
library(ggplot2)

# Read the activities data
activities <- read.csv("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/COSMIC_SBS96_Activities_tcga.txt", sep = "\t", header = TRUE, row.names = 1)



library(pheatmap)
# pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Fig1B.pdf", width = 4.25, height = 1.81)
# Plotting a heatmap of the activities
pheatmap::pheatmap(as.matrix(activities),
                   main = "Activities of COSMIC Signatures",
                   scale = "row",
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   color = colorRampPalette(c("blue", "white", "red"))(100))
# dev.off()


# Load necessary libraries
library(proxy)  # for cosine similarity calculation
library(reshape2)  # for data manipulation
library(pheatmap)  # for heatmap visualization

# Example data loading
my_signatures <- read.csv("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/SBS96_De-Novo_Signatures_tcga.txt", sep ='\t')
cosmic_signatures <- read.csv("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/COSMIC_v3.4_SBS_GRCh38.txt", sep='\t')

# Assuming my_signature and cosmic_signature are your matrices
my_signature_transposed <- t(my_signatures)
cosmic_signature_transposed <- t(cosmic_signatures)

rownames(my_signatures) <- my_signatures[,1]
my_signatures <- my_signatures[,-1]


# Setting the first column as row names
colnames(my_signature_transposed) <- my_signature_transposed[1,]
colnames(cosmic_signature_transposed) <- cosmic_signature_transposed[1, ]

# Removing the first column
my_signature_transposed <- my_signature_transposed[-1,]
cosmic_signature_transposed <- cosmic_signature_transposed[-1 ,]



# Calculate cosine similarities
cosine_sim_matrix <- as.matrix(proxy::simil(my_signature_transposed, cosmic_signature_transposed, method = "cosine"))

# Check the result
cosine_sim_matrix
# Plotting a heatmap of the cosine similarities
# pheatmap::pheatmap(cosine_sim_matrix,
#                    main = "Cosine Similarity against COSMIC Signatures",
#                    color = colorRampPalette(c("blue", "white", "red"))(100),
#                    cluster_rows = FALSE,
#                    cluster_cols = FALSE)
# pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Fig1B_tcga.pdf", width = 4.25, height = 1.81)
pheatmap::pheatmap(mat = cosine_sim_matrix, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
# dev.off()

### (4.1, 2)
# maftools::plotSignatures(nmfRes = HNSCC.sig, title_size = 1.2, sig_db = "SBS")

# Finding the maximum similarity for each of your signatures and corresponding COSMIC signature
max_indices <- max.col((cosine_sim_matrix), ties.method = "first")  # transposing to find row-wise maxima in original
max_values <- cosine_sim_matrix[cbind(max_indices, 1:ncol(cosine_sim_matrix))]

# Assuming cosine_sim_matrix is already loaded and has proper row and column names

# Find indices of the maximum values in each row
max_indices <- max.col(cosine_sim_matrix, ties.method = "first")

# Extract these maximum values
max_values <- cosine_sim_matrix[cbind(1:nrow(cosine_sim_matrix), max_indices)]

# Retrieve the names of the COSMIC signatures that each of your signatures best matches
best_match_cosmic_sigs <- colnames(cosine_sim_matrix)[max_indices]

# Create a data frame for easier handling and visualization
best_matches_df <- data.frame(
  MySignature = rownames(cosine_sim_matrix),
  BestMatchCOSMIC = best_match_cosmic_sigs,
  CosineSimilarity = max_values
)

# Print out the dataframe to see the results
print(best_matches_df)

plotCustomSignatures <- function(signature_matrix, best_fit, max_values, title_size = 1.3, axis_lwd = 2, font_size = 1.2, yaxisLim = 0.3) {
  # Check that the number of columns in signature_matrix matches the number of rows in best_fit
  if (ncol(signature_matrix) != nrow(best_fit)) {
    stop("The number of columns in signature_matrix must match the number of rows in best_fit")
  }
  
  # Plot settings
  nsigs <- ncol(signature_matrix)
  color <- c('coral4', 'lightcyan4', 'deeppink3', 'lightsalmon1', 'forestgreen', 'cornflowerblue')
  colors <- c('coral4', 'lightcyan4', 'deeppink3', 'lightsalmon1', 'forestgreen', 'cornflowerblue')
  colors <- rep(colors, each=16)
  
  # Adjusting oma and mar to ensure titles are not cut off
  par(mfrow = c(nsigs, 1), oma = c(5, 4, 4, 0) + 0.1, mar = c(0, 0, 3, 0) + 0.1, las=1, tcl=-.25, font.main=4, xpd = NA)
  
  for (i in 1:nsigs) {
    etiology <- paste0(best_fit$BestMatch[i], " \n Aetiology: ", best_fit$Aetiology[i])
    cosine_sim <- round(max_values[i], 3)  # Round the cosine similarity value for display
    title_text <- paste0(etiology, "\n Cosine Similarity: ", cosine_sim)
    
    d <- as.matrix(signature_matrix[, i])
    
    if (is.na(yaxisLim)) {
      bh <- ceiling(max(d, na.rm = TRUE) * 10) / 10  # Bar height
    } else {
      bh <- yaxisLim
    }
    
    barplot(d, xaxt = "n", yaxt = "n", col = colors, beside = TRUE, ylim = c(-0.1, bh),
            cex.main = 1, border = NA, font.axis = 2, font.lab = 2,
            adj = 0.25)
    
    title(main = title_text, cex.main = title_size, line = 0, font.main = 3)
    
    axis(side = 2, at = seq(0, bh, 0.1),
         pos = -2, las = 2, lwd = axis_lwd, hadj = 1.1,
         font = 1, cex.axis = font_size)
    
    rect(xleft = seq(0, 96, 16), ybottom = -0.05, xright = 96, ytop = -0.02, col = color, border = 'gray70')
    # rect(xleft = seq(0, 160, 32), ybottom = -0.05, xright = seq(8, 96, 16), ytop = 0.05, col = colors[1:6], border = 'gray70')
    
    if (i == nsigs) {
      text(labels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
           y = rep(-0.15, 6), x = seq(8, 96, 16), cex = font_size,
           font = 1, font.lab = 2, pos = 3)
    }
  }
}
##8.1,4

best_match_cosmic_sigs
etiology_vector <- c("APOBEC Cytidine Deaminase (C>T)", "Defective DNA mismatch repair", "Unknown : Clock-like Signature", "Tobacco smoking")
best_fit <- data.frame(
  BestMatch = best_match_cosmic_sigs,
  Aetiology = etiology_vector  # Assuming you have this vector; replace with actual data
)

plotCustomSignatures(my_signatures, best_fit, max_values)

plotCustomSignatures <- function(signature_matrix, best_fit, max_values, title_size = 1.3, axis_lwd = 2, font_size = 1.2, yaxisLim = 0.3) {
  # Check that the number of columns in signature_matrix matches the number of rows in best_fit
  if (ncol(signature_matrix) != nrow(best_fit)) {
    stop("The number of columns in signature_matrix must match the number of rows in best_fit")
  }
  
  # Plot settings
  nsigs <- ncol(signature_matrix)
  colors <- c('coral4', 'lightcyan4', 'deeppink3', 'lightsalmon1', 'forestgreen', 'cornflowerblue')
  colors <- rep(colors, each = 16)
  
  # Adjusting oma and mar to ensure titles and annotations are not cut off
  par(mfrow = c(nsigs, 1), oma = c(7, 4, 4, 2) + 0.1, mar = c(0, 0, 4, 2) + 0.1, las = 1, tcl = -0.25, font.main = 4, xpd = NA)
  
  for (i in 1:nsigs) {
    etiology <- paste0(best_fit$BestMatch[i], " \n Aetiology: ", best_fit$Aetiology[i])
    cosine_sim <- round(max_values[i], 3)  # Round the cosine similarity value for display
    title_text <- paste0(etiology, "\n Cosine Similarity: ", cosine_sim)
    
    d <- as.matrix(signature_matrix[, i])
    
    if (is.na(yaxisLim)) {
      bh <- ceiling(max(d, na.rm = TRUE) * 10) / 10  # Bar height
    } else {
      bh <- yaxisLim
    }
    
    barplot(d, xaxt = "n", yaxt = "n", col = colors, beside = TRUE, ylim = c(-0.1, bh),
            cex.main = 1, border = NA, font.axis = 2, font.lab = 2,
            adj = 0.25)
    
    title(main = title_text, cex.main = title_size, line = 0, font.main = 3)
    
    axis(side = 2, at = seq(0, bh, 0.1),
         pos = -2, las = 2, lwd = axis_lwd, hadj = 1.1,
         font = 1, cex.axis = font_size)
  }
  
  # Add the lower bar with proper alignment and colors
  par(new = TRUE, oma = c(2, 4, 0, 2) + 0.1, mar = c(0, 0, 0, 0) + 0.1, xpd = NA)
  plot(0, 0, type = "n", xlim = c(0, 192), ylim = c(-0.1, 0.1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  # rect(xleft = seq(0, 160, 32), ybottom = -0.05, xright = seq(32, 192, 32), ytop = 0.05, col = colors[1:6], border = 'gray70')
  text(labels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
       y = rep(-0.05, 6), x = seq(16, 192, 32), cex = font_size,
       font = 1, font.lab = 2, pos = 3)
}



# Use the function to plot the data
plotCustomSignatures(my_signatures, best_fit, max_values)



