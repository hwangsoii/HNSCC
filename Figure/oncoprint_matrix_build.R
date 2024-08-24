oncoprint_matrix_build <- function(maf_data){
  samples <- unique(maf_data@data$Tumor_Sample_Barcode)
  genes <- unique(maf_data@data$Hugo_Symbol)
  mat <- matrix("",length(genes),length(samples))
  colnames(mat) <- samples; rownames(mat) <- genes; mat
  pb <- txtProgressBar(min = 0, max = nrow(mat), style = 3)
  for(i in 1:nrow(mat)) {
    curGene <- rownames(mat)[i]
    setTxtProgressBar(pb, i)
    for(j in 1:ncol(mat)) {
      curSample <- colnames(mat)[j]
      if(length(intersect(maf_data@data$Tumor_Sample_Barcode, curSample))==1){
        mat1 <- maf_data@data[maf_data@data$Tumor_Sample_Barcode == curSample,]
        if(length(intersect(mat1$Hugo_Symbol, curGene))==1){
          mat3 <- mat1[mat1$Hugo_Symbol == curGene,]
          mat[curGene,curSample] <- paste(as.character(mat3$Variant_Classification), collapse = ";")
        }
      }
    }
  }
  close(pb)
  return(mat)
}
