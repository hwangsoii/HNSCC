source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")

install.packages("VennDiagram")
library(VennDiagram)
###pik3ca in oro
total_clinical$cnv <-ifelse(total_clinical$Tumor_Sample_Barcode %in% cnv_sample_barcode,"Yes", "No")
clinical_practice <- total_clinical

clinical_oropharynx <- clinical_practice %>% filter(final_site =='oropharynx')
clinical_oralcavity <- clinical_practice %>% filter(final_site =='Lip or oral cavity')
clinical_larynx <- clinical_practice %>% filter(final_site =='Larynx')
clinical_hypopharynx <- clinical_practice %>% filter(final_site =='Hypopharynx')

g <- ggplot(clinical_oropharynx, aes(x = TP53, fill = PIK3CA)) + geom_bar(position = "fill")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
g
g <- ggplot(clinical_oropharynx, aes(x = TP53, fill = PIK3CA)) + geom_bar(position = "dodge")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
g

TP53_pik3ca_oro <- table(clinical_oropharynx$TP53, clinical_oropharynx$PIK3CA)
TP53_pik3ca_oro
fisher.test(TP53_pik3ca_oro)
colnames(TP53_pik3ca_oro) <- c('No PIK3CA mutation', 'with PIK3CA mutation')
rownames(TP53_pik3ca_oro) <- c('No TP53 mutation', 'with TP53 mutation')
TP53_pik3ca_oro %>%
  kbl() %>%
  kable_paper("hover", full_width=F)
chisq.test(TP53_pik3ca_oro)
###tp53 pik3ca except oro
clinical_exceptoro <- clinical_practice %>% filter(final_site !='oropharynx')

TP53_pik3ca_exceptoro <- table(clinical_exceptoro$TP53, clinical_exceptoro$PIK3CA)
colnames(TP53_pik3ca_exceptoro) <- c('No PIK3CA mutation', 'with PIK3CA mutation')
rownames(TP53_pik3ca_exceptoro) <- c('No TP53 mutation', 'with TP53 mutation')
TP53_pik3ca_exceptoro
fisher.test(TP53_pik3ca_exceptoro)
(TP53_pik3ca_exceptoro) %>%
  kbl() %>%
  kable_paper("hover", full_width=F)

g <- ggplot(clinical_exceptoro, aes(x = TP53, fill = PIK3CA)) + geom_bar(position = "fill")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
g
g <- ggplot(clinical_exceptoro, aes(x = TP53, fill = PIK3CA)) + geom_bar(position = "dodge")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

g


########revise as vendiagram

# load the VennDiagram package
library(VennDiagram)

# subset the data frames based on TP53 and PIK3CA columns
oropharynx_TP53 <- clinical_oropharynx[clinical_oropharynx$TP53 == TRUE, ]
oropharynx_PIK3CA <- clinical_oropharynx[clinical_oropharynx$PIK3CA == TRUE, ]
exceptoro_TP53 <- clinical_exceptoro[clinical_exceptoro$TP53 == TRUE, ]
exceptoro_PIK3CA <- clinical_exceptoro[clinical_exceptoro$PIK3CA == TRUE, ]

# create the Venn diagrams
venn_oropharynx <- venn.diagram(
  x = list(TP53 = rownames(oropharynx_TP53), PIK3CA = rownames(oropharynx_PIK3CA)),
  filename = NULL,
  fill = c("dodgerblue", "tomato"),
  alpha = c(0.5, 0.5),
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  cat.dist = c(0.06, 0.06),
  cat.col = c("dodgerblue", "tomato"),
  cat.fontface = "bold",
  cex = 1.5,
  fontfamily = "sans",
  main = "Venn Diagram of TP53 and PIK3CA in clinical_oropharynx"
)

venn_exceptoro <- venn.diagram(
  x = list(TP53 = rownames(exceptoro_TP53), PIK3CA = rownames(exceptoro_PIK3CA)),
  filename = NULL,
  fill = c("darkseagreen", "goldenrod"),
  alpha = c(0.5, 0.5),
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  cat.dist = c(0.06, 0.06),
  cat.col = c("darkseagreen", "goldenrod"),
  cat.fontface = "bold",
  cex = 1.5,
  fontfamily = "sans",
  main = "Venn Diagram of TP53 and PIK3CA in clinical_exceptoro"
)

# display the diagrams
grid.newpage()
grid.draw(venn_oropharynx)
grid.newpage()
grid.draw(venn_exceptoro)


# define custom colors
my_colors <- c("dodgerblue", "tomato", "darkseagreen", "goldenrod")

# define custom shapes
my_shapes <- c("circle", "square", "triangle-up", "triangle-down")

# define custom font sizes
my_font_sizes <- c(3, 2)

# create the first Venn diagram
venn_oropharynx <- venn.diagram(
  x = list(TP53 = rownames(oropharynx_TP53), PIK3CA = rownames(oropharynx_PIK3CA)),
  filename = NULL,
  fill = my_colors[1:2],
  alpha = c(0.7, 0.7),
  cat.cex = my_font_sizes[1],
  cat.pos = c(-150, 150),
  cat.dist = c(0.1, 0.1),
  cat.col = my_colors[1:2],
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cex = my_font_sizes[2],
  fontfamily = "serif",
  main = "Clinical Oropharynx",
  main.fontfamily = "serif",
  main.cex = my_font_sizes[1],
  main.col = "black",
  margin = 0.1,
  lwd = 3,
  col = my_colors[1:2],
  fill.alpha = 0.7,
  shape = my_shapes[1:2]
)

# create the second Venn diagram
venn_exceptoro <- venn.diagram(
  x = list(TP53 = rownames(exceptoro_TP53), PIK3CA = rownames(exceptoro_PIK3CA)),
  filename = NULL,
  fill = my_colors[3:4],
  alpha = c(0.7, 0.7),
  cat.cex = my_font_sizes[1],
  cat.pos = c(-150, 150),
  cat.dist = c(0.1, 0.1),
  cat.col = my_colors[3:4],
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cex = my_font_sizes[2],
  fontfamily = "serif",
  main = "Clinical Exceptoro",
  main.fontfamily = "serif",
  main.cex = my_font_sizes[1],
  main.col = "black",
  margin = 0.1,
  lwd = 3,
  col = my_colors[3:4],
  fill.alpha = 0.7,
  shape = my_shapes[3:4]
)

# display the diagrams
grid.newpage()
grid.draw(venn_oropharynx)
grid.newpage()
grid.draw(venn_exceptoro)


