source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")

library(maftools)

####lollipop plot
#p53
dat_mat <- All_finalcnvkit@data
lollipopPlot(
  maf = All_finalcnvkit, gene = 'TP53', AACol = 'HGVSp_Short', showMutationRate = FALSE,labelPos = 'all'
)
#pik3ca
dat_mat <- All_finalcnvkit@data
lollipopPlot(
  maf = All_finalcnvkit, gene = 'PIK3CA', AACol = 'HGVSp_Short', showMutationRate = FALSE,labelPos = 'all'
)
#EGFR
dat_mat <- All_finalcnvkit@data
lollipopPlot(
  maf = All_finalcnvkit, gene = 'EGFR', AACol = 'HGVSp_Short',domainLabelSize = 0.5,
  axisTextSize = c(1, 1), showMutationRate = FALSE,labelPos = 'all'
)

#EGFR
dat_mat <- All_finalcnvkit@data
lollipopPlot(
  maf = All_finalcnvkit, gene = 'CDKN2A', AACol = 'HGVSp_Short',domainLabelSize = 0.5,
  axisTextSize = c(1, 1), showMutationRate = FALSE,labelPos = 'all'
)
