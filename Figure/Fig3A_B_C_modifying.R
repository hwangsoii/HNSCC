library("gridExtra")   ###P
library(ggsci)   ###P
library(ggplot2) ###P
library(ggsignif)###P
library(dplyr)   ###P

library(fmsb)    ###N
library(reshape2)###N
library(tidyverse)###P
library(viridis)###N

source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")

# install.packages(devtools)



cnv_sample_barcode <- cnvkit_0.4_onlyhnscc$sample
df_tsb <- total_clinical$Tumor_Sample_Barcode
df_tsb <- as.data.frame(df_tsb)
getwd()
amp_sample_barcode <- 

# Step 1: Create amp_barcode and del_barcode based on amp_del values
amp_barcode <- cnvkit_0.4_onlyhnscc$sample[cnvkit_0.4_onlyhnscc$amp_del == "Amplification"]
del_barcode <- cnvkit_0.4_onlyhnscc$sample[cnvkit_0.4_onlyhnscc$amp_del == "Deletion"]

# Step 2: Update total_clinical dataframe to include whether each sample is an amplification or deletion
total_clinical$cnv_amp <- ifelse(total_clinical$Tumor_Sample_Barcode %in% amp_barcode, "Yes", "No")
total_clinical$cnv_del <- ifelse(total_clinical$Tumor_Sample_Barcode %in% del_barcode, "Yes", "No")

total_clinical$cnv <-ifelse(total_clinical$Tumor_Sample_Barcode %in% cnv_sample_barcode,"Yes", "No")

for (gene in top30_genes){
  assign(paste0('somatic_', gene, '_barcode'), unlist(genesToBarcodes(maf = All_cnvkit, genes = gene, justNames = TRUE)))
  total_clinical[paste0(get('gene'),'_total')] <- total_clinical$Tumor_Sample_Barcode %in% get(paste0('somatic_', gene, '_barcode')) 
}
for (gene in top30_genes){
  assign(paste0('somatic_', gene, '_barcode'), unlist(genesToBarcodes(maf = All_cnvkit_nocnv, genes = gene, justNames = TRUE)))
  total_clinical[get('gene')] <- total_clinical$Tumor_Sample_Barcode %in% get(paste0('somatic_', gene, '_barcode')) 
}


df_top_gene <- data.frame()
for (gene in top30_genes){
  tmp <- total_clinical  %>% group_by(final_site, total_clinical[get('gene')]) %>%
    summarise(Counts = n(),.groups="drop") %>% mutate(freq = Counts/sum(Counts))
  newtmp <- cbind(rep(gene, times = nrow(tmp)), tmp[1:4])
  names(newtmp) <- names(df_top_gene)
  df_top_gene <-  rbind(df_top_gene,  newtmp)
  names(df_top_gene) <- c('Gene','Site', 'Type', 'Count', 'Ratio')
  
}



 
# ##################


df_top_gene_onlytrue <- df_top_gene %>% filter(Type == TRUE)
df_top_gene_onlytrue <- df_top_gene_onlytrue %>% 
  filter (Site =='oropharynx'|Site =='Hypopharynx'|Site =='Lip or oral cavity'|Site =='Larynx')


g <- ggplot(df_top_gene_onlytrue, aes(x=Gene, y=Ratio)) + geom_bar(stat = 'identity', position = "dodge")+
  scale_x_discrete(limits = top30_genes)+ scale_y_continuous(expand = c(0,0), limits = c(0,1), 
                                                             breaks = seq(0,1, by = 0.25))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),axis.text.x=element_text(size=rel(1), angle=45)) 

g + geom_signif(stat = "identity",
                data = data.frame(x = c(0.7, 1.7, 2.7),
                                  xend = c(1.3, 2.3, 3.3),
                                  y = c(0.8, 0.9, 0.5),
                                  annotation = c("NS", "***", "*")),
                aes(x = x,
                    xend = xend,
                    y = y,
                    yend = y,
                    annotation = annotation))

g

df_top_gene2 <- df_top_gene %>%
  filter (Site =='oropharynx'|Site =='Hypopharynx'|Site =='Lip or oral cavity'|Site =='Larynx')
df_top_gene2 %>% 
  group_by(Site) %>% 
  summarise(mean(Count))
threshold <- 0.05 /6
genesss_mut <- c()
for (gene in top30_genes){
  sites <- c("oropharynx", "Hypopharynx", "Lip or oral cavity", "Larynx")
  for (organ in sites){
    df_gene <- df_top_gene2 %>% filter ((gene == Gene) & (organ == Site))
    vector1 <- c(df_gene$Count[1], df_gene$Count[2])
    vector1[is.na(vector1)] <- 0
    sites <- sites [! sites %in% organ]
    for (organ2 in sites){
      df_gene2 <- df_top_gene2 %>% filter ((gene == Gene) & (organ2 == Site))
      vector2 <- c(df_gene2$Count[1], df_gene2$Count[2])
      vector2[is.na(vector2)] <- 0
      res <- fisher.test(cbind(vector1,vector2))
      if (res$p.value < threshold){
        genesss_mut <- c(genesss_mut, gene)
      }
    }
  }
}




CNV_top30 <- c("CTTN" ,"CCND1", "AMER1" , "SOX2",   "SOX2-OT",  "KLHL6",    "KDM6A",    "TP63",     "KDM5C",    "HDAC6",    "ARAF",     "MIR6895",  "MIR6894", 
               "AR",       "PIK3CA",   "GATA1",    "ATRX",     "PRKCI",    "PAK3",     "EGFR",     "EGFR-AS1", "MYC",      "CDKN2A",   "FGFR1",    "NSD3",     "ATR",     
               "TRIO",     "EPHB4",    "SLC12A9",  "IL7R" )
for (cnvgene in CNV_top30){
  tmp_df <- cnvkit_0.4_onlyhnscc %>% filter(gene == cnvgene)
  assign(paste0('cnv_', cnvgene, '_barcode'), tmp_df$sample)
  # print(setdiff(cnvkit_0.4_onlyhnscc[ cnvkit_0.4_onlyhnscc$gene == gene,]$sample, get(paste0('somatic_', gene, '_barcode'))))
  total_clinical[paste0(cnvgene, '_cnv')] <- total_clinical$Tumor_Sample_Barcode %in% get(paste0('cnv_', cnvgene, '_barcode'))
}


df_top_cnv <- data.frame()
for (cnvgene in CNV_top30){
  tmp <- total_clinical  %>% 
    group_by(final_site, total_clinical[paste0(cnvgene, '_cnv')]) %>% summarise(Counts = n(),.groups="drop") %>% mutate(cnvfreq = Counts/sum(Counts))
  newtmp <- cbind(rep(cnvgene, times = nrow(tmp)), tmp[1:4])
  names(newtmp) <-  c('Gene','Site', 'Type', 'Count', 'Ratio')
  df_top_cnv <-  rbind(df_top_cnv,  newtmp)
  names(df_top_cnv) <- c('Gene','Site', 'Type', 'Count', 'Ratio')
}



vet2 <- total_clinical
bar_size <- 0.2
text_size <- 2
a = ggsci::pal_aaas()
cols <- c("Hypopharynx" = a(4)[1], "Larynx" = a(4)[2], "Lip or oral cavity" = a(4)[3], "oropharynx" = a(4)[4])


#################################################################################
grid_func <- function(x, y){
  vet_filtered <- vet2 %>%
    filter(final_site =='oropharynx'|final_site =='Hypopharynx'|final_site =='Lip or oral cavity'|final_site =='Larynx')

  vet_filtered$want <- vet_filtered[,x]
  
  count_plot <- vet_filtered %>% 
    group_by(final_site) %>% 
    count(want) %>% 
    mutate(Freq = n) %>%
    ggplot() + 
    aes(x = final_site, Freq, fill = want) +
    labs(x= 'Primary Site',y= y, fill = y) +
    scale_fill_grey(start = 0.9, end = 0.4)+
    geom_bar(position = "dodge", stat='identity') + 
    geom_signif(
      comparisons = list(c("Hypopharynx", 'Lip or oral cavity')),
      map_signif_level = TRUE, y_position = 140,
      size=bar_size, # This doesn't work!
      textsize=text_size,
      annotations="***"
    )+
    geom_signif(
      comparisons = list(c("Lip or oral cavity", "oropharynx")),
      map_signif_level = TRUE, y_position = 150,
      size=bar_size, # This doesn't work!
      textsize=text_size,
      annotations="*"
    ) +
    geom_signif(
      comparisons = list(c("Lip or oral cavity", "Larynx")),
      map_signif_level = TRUE, y_position = 160,
      size=bar_size, # This doesn't work!
      textsize=text_size,
      annotations="***"
    )+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),legend.title=element_text(size=8),     legend.text=element_text(size=6),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  
  
  ratio_plot <- vet_filtered %>% 
    group_by(final_site) %>% 
    count(want) %>% 
    mutate(Freq = n) %>%
    ggplot() + 
    # scale_fill_npg()+
    scale_fill_grey(start = 0.9, end = 0.4)+
    aes(x = final_site, Freq, fill = want) +
    labs(x= 'Primary Site',y= y, fill = y) +
    geom_bar(position = "fill", stat='identity') + 
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),legend.title=element_text(size=8),     legend.text=element_text(size=6),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  
  return(list(count_plot, ratio_plot))
}
print("grid_func() return: list(count_plot,ratio_plot)")
#################################################################################

PIK3CA_plot <- grid_func("PIK3CA_cnv","PIK3CA_amp")
EGFR_plot <- grid_func("EGFR_cnv","EGFR_amp")
FGFR1_plot <- grid_func("FGFR1_cnv","FGFR1_amp")
CCND1_plot <- grid_func("CCND1_cnv","CCND1_amp")
CDKN2A_plot <- grid_func("CDKN2A_cnv","CDKN2A_amp")

# setwd("/data/project/TRIUMPH/Figure2/new_figure_manu/")
# pdf("Count&ratio.pdf", width = 20, height = 25)
# pdf("Fig3A.pdf", width = 20, height = 25)
# grid.arrange(PIK3CA_plot[[1]], PIK3CA_plot[[2]],
#              EGFR_plot[[1]], EGFR_plot[[2]],
#              FGFR1_plot[[1]], FGFR1_plot[[2]],
#              CCND1_plot[[1]], CCND1_plot[[2]],
#              CDKN2A_plot[[1]], CDKN2A_plot[[2]], ncol = 2)
# dev.off()
# 


vet_filtered <- vet2 %>%
  filter(final_site =='oropharynx'|final_site =='Hypopharynx'|final_site =='Lip or oral cavity'|final_site =='Larynx')

site_num <- as.data.frame(vet_filtered %>% group_by(final_site) %>% count(final_site) )

PIK3CA_site_num <- vet_filtered %>% group_by(final_site) %>% count(PIK3CA_cnv) 
PIK3CA_site_num <- as.data.frame(PIK3CA_site_num[PIK3CA_site_num$PIK3CA_cnv==TRUE,c(1,3)])

EGFR_site_num <- vet_filtered %>% group_by(final_site) %>% count(EGFR_cnv)
EGFR_site_num <- as.data.frame(EGFR_site_num[EGFR_site_num$EGFR_cnv==TRUE,c(1,3)])

FGFR1_site_num <- vet_filtered %>% group_by(final_site) %>% count(FGFR1_cnv)
FGFR1_site_num <- as.data.frame(FGFR1_site_num[FGFR1_site_num$FGFR1_cnv==TRUE,c(1,3)])

CCND1_site_num <- vet_filtered %>% group_by(final_site) %>% count(CCND1_cnv)
CCND1_site_num <- as.data.frame(CCND1_site_num[CCND1_site_num$CCND1_cnv==TRUE,c(1,3)])

CDKN2A_site_num <- vet_filtered %>% group_by(final_site) %>% count(CDKN2A_cnv)
CDKN2A_site_num <- as.data.frame(CDKN2A_site_num[CDKN2A_site_num$CDKN2A_cnv==TRUE,c(1,3)])

# rchart_df <- cbind(PIK3CA_site_num[,2]/site_num[,2],
#                    EGFR_site_num[,2]/site_num[,2],
#                    FGFR1_site_num[,2]/site_num[,2],
#                    CCND1_site_num[,2]/site_num[,2],
#                    CDKN2A_site_num[,2]/site_num[,2])

rchart_df <- cbind(PIK3CA_site_num[,2]/site_num[,2],
                   EGFR_site_num[,2]/site_num[,2],
                   FGFR1_site_num[,2]/site_num[,2],
                   CCND1_site_num[,2]/site_num[,2])

rownames(rchart_df) <- PIK3CA_site_num[,1]
# colnames(rchart_df) <- c("PIK3CA_amp" , "EGFR_amp" , "FGFR1_amp" , "CCND1_amp", "CDKN2A_amp")
colnames(rchart_df) <- c("PIK3CA_amp" , "EGFR_amp" , "FGFR1_amp" , "CCND1_amp")


rchart_df2 <- data.frame(PIK3CA_site_num[,1],
                    PIK3CA_site_num[,2]/site_num[,2],
                    EGFR_site_num[,2]/site_num[,2],
                    FGFR1_site_num[,2]/site_num[,2],
                    CCND1_site_num[,2]/site_num[,2])
rownames(rchart_df2) <- PIK3CA_site_num[,1]
colnames(rchart_df2) <- c("Site","PIK3CA_amp" , "EGFR_amp" , "FGFR1_amp" , "CCND1_amp")

rchart_df
rchart_df2
rchart_df.melt.origin <- melt(rchart_df)
# rchart_df <- as.data.frame(rbind(rep(0.5,5) , rep(0,5) , rchart_df))
rchart_df <- as.data.frame(rbind(rep(0.5,4) , rep(0,4) , rchart_df))

###ratcio check
mean(as.numeric(rchart_df['Hypopharynx',]))

colors_border=c("#4faf4b", "#94509a","#e41d1f","#377fb9")
colors_in=c("#4faf4b4D", "#94509a4D","#e41d1f4D","#377fb94D")

# pdf("radarchart.pdf", width = 20, height = 20)
# pdf("Fig3A.pdf", width = 20, height = 25)
# radarchart(rchart_df  , axistype=1 ,
#            pcol=colors_border , pfcol=colors_in , plwd=6 , plty=1,
#            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.5,0.1), cglwd=0.1,
#            vlcex=1.5, seg = 5
# )
# legend(x=0.7, y=1, legend = rownames(rchart_df[-c(1,2),]), 
#        bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3)
# dev.off()


# setwd("/data/project/TRIUMPH/Figure2/new_figure_manu/")
# pdf("radarchart2.pdf", width = 20, height = 20)
pdf("/data/project/TRIUMPH/Figure_swkim/Fig3B.pdf", width = 20, height = 25)
ggradar(rchart_df2, values.radar = c("0","0.25","0.5"),
        grid.min = 0, grid.mid = 0.25, grid.max = 0.5,
        group.line.width = 3, 
        group.point.size = 7,
        group.colours = colors_border,
        background.circle.colour = "white",
        gridline.max.colour = "#000000",
        gridline.mid.colour = "#0000004C",
        gridline.min.colour = "#0000004C",
        grid.line.width = 1,
        legend.position = "bottom",
        fill = T,
        fill.alpha = 0.3,
        legend.text.size = 20,
        axis.label.size = 10)
dev.off()

rchart_df.melt <- rchart_df.melt.origin
rchart_df.melt$Name <- paste(rchart_df.melt$Var1,rchart_df.melt$Var2, sep = "-")

empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(rchart_df.melt$Var2), ncol(rchart_df.melt)) )
colnames(to_add) <- colnames(rchart_df.melt)
to_add$Var2 <- rep(levels(rchart_df.melt$Var2), each=empty_bar)
rchart_df.melt <- rbind(rchart_df.melt, to_add)
rchart_df.melt <- rchart_df.melt %>% arrange(Var2)
rchart_df.melt$id <-seq(1, nrow(rchart_df.melt))


label_data <- rchart_df.melt
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

p <- ggplot(rchart_df.melt, aes(x=as.factor(id), y=value, fill=Var2)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-0.1,0.55) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+0.05, label=Var1, hjust=hjust),
            color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE) 

# pdf("test.pdf", width = 20, height = 20)
# p
# dev.off()


PIK3CA_site_TF <- vet_filtered %>% group_by(final_site) %>% count(PIK3CA_cnv) %>%
  mutate(PIK3CA_cnv = ifelse(PIK3CA_cnv,"PIK3CA_amp","PIK3CA"))
colnames(PIK3CA_site_TF)<-c("site","Amp","count")

EGFR_site_TF <- vet_filtered %>% group_by(final_site) %>% count(EGFR_cnv) %>%
  mutate(EGFR_cnv = ifelse(EGFR_cnv,"EGFR_amp","EGFR"))
colnames(EGFR_site_TF)<-c("site","Amp","count")

FGFR1_site_TF <- vet_filtered %>% group_by(final_site) %>% count(FGFR1_cnv) %>%
  mutate(FGFR1_cnv = ifelse(FGFR1_cnv,"FGFR1_amp","FGFR1"))
colnames(FGFR1_site_TF)<-c("site","Amp","count")

CCND1_site_TF <- vet_filtered %>% group_by(final_site) %>% count(CCND1_cnv) %>%
  mutate(CCND1_cnv = ifelse(CCND1_cnv,"CCND1_amp","CCND1"))
colnames(CCND1_site_TF)<-c("site","Amp","count")

cir_bar <- as.data.frame(rbind(PIK3CA_site_TF, EGFR_site_TF, FGFR1_site_TF, CCND1_site_TF))
cir_bar$Name <- paste(cir_bar$site,cir_bar$Amp, sep = "_")


empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(as.factor(cir_bar$Amp)), ncol(cir_bar)))
colnames(to_add) <- colnames(cir_bar)
to_add$Amp <- rep(levels(as.factor(cir_bar$Amp)), each=empty_bar)
cir_bar <- rbind(cir_bar, to_add)
cir_bar <- cir_bar %>% arrange(Amp)
cir_bar$id <-seq(1, nrow(cir_bar))


label_data <- cir_bar
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

p <- ggplot(cir_bar, aes(x=as.factor(id), y=count, fill=Amp)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-50,150) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  # coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=count+10, label=site, hjust=hjust),
            color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 



cir_bar2 <- as.data.frame(rbind(PIK3CA_site_TF, EGFR_site_TF, FGFR1_site_TF, CCND1_site_TF))
cir_bar2$Name <- paste(cir_bar2$site,cir_bar2$Amp, sep = "_")
cir_bar2$Amp <- gsub("_amp","",cir_bar2$Amp)


empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(as.factor(cir_bar2$Amp)), ncol(cir_bar2)))
colnames(to_add) <- colnames(cir_bar2)
to_add$Amp <- rep(levels(as.factor(cir_bar2$Amp)), each=empty_bar)
cir_bar2 <- rbind(cir_bar2, to_add)
cir_bar2 <- cir_bar2 %>% arrange(Amp)
cir_bar2$id <- rep(seq(1,nrow(cir_bar2)/2), each= 2)


label_data <- cir_bar2 %>% group_by(id, site) %>% summarize(tot=sum(count))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

base_data <- cir_bar2 %>% 
  group_by(Amp) %>% 
  summarize(start=min(id), end=max(id) - 2) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
base_data$angle <- c(-30,60,-30,60)
base_data$id <- apply(base_data,1,function(x) mean(c(as.numeric(x[2]),as.numeric(x[3]))))

grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:(nrow(grid_data)-1))] + 1
grid_data$start <- grid_data$start - 1
grid_data$start[1] <- max(cir_bar2$id)


manual_color <- brewer.pal(8,"Paired")[c(rep(1:2,4),rep(3:4,4),rep(5:6,4),rep(7:8,4))]
manual_color <- rep(brewer.pal(8,"Paired"),4)


p2 <- ggplot(cir_bar2) +  #, aes(x=as.factor(id), y=count, fill=Amp)     # Note that id is a factor. If x is numeric, there is some space between the first bar
  # geom_bar(stat="identity", alpha=0.5) +
  geom_bar(aes(x=as.factor(id), y=count, fill=Name), stat="identity", alpha=0.5) +
  scale_fill_manual(values = manual_color)+
  ylim(-100,max(label_data$tot, na.rm=T)+10) +
  theme_minimal() +
  theme(
    # legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(5,4), "cm"),
    legend.position = c(0.5,0.5)
  ) +
  coord_polar() + #circle mode
  ####grid drawing and value marking
  geom_text(data=label_data, aes(x=id, y=tot+10, label=site, hjust=hjust),
            color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  geom_segment(data=grid_data, aes(x = end-0.25, y = 0, xend = start+0.25, yend = 0), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end-0.25, y = 50, xend = start+0.25, yend = 50), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end-0.25, y = 100, xend = start+0.25, yend = 100), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end-0.25, y = 150, xend = start+0.25, yend = 150), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  ggplot2::annotate("text", x = rep(max(cir_bar2$id)+0.5,4), y = c(0, 50, 100, 150), 
                    label = c("0", "50", "100", "150") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  
  ######base line & annotation
  geom_segment(data=base_data, aes(x = start-0.5, y = -5, xend = end+0.5, yend = -5), 
               colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  + 
  geom_text(data=base_data, aes(x=id, y= -10, label=Amp, hjust=0.5),
            color="black", fontface="bold",alpha=0.6, size=5, angle=base_data$angle, inherit.aes = FALSE ) +
  scale_color_identity(
    guide = guide_legend(),
    breaks = c("red", "green", "blue"),
    labels = c("setosa", "versicolor", "virginica")
  )
# setwd("/data/project/TRIUMPH/Figure/")
# setwd("/data/project/TRIUMPH/Figure2/new_figure_manu/")
# pdf("Cir_bar3.pdf", width = 20, height = 20)
pdf("/data/project/TRIUMPH/Figure_swkim/Fig3A.pdf", width = 20, height = 20)

p2
dev.off()




#############p-values of CNV diff between sites 
vet_test <- vet2 %>% filter(final_site =='Larynx'|final_site =='Hypopharynx')
for_test <- table(vet_test$final_site, vet_test$cnv)
fisher_re <- fisher.test(for_test)
print(paste('hypo vs larynx',fisher_re$p.value, sep = ' : '))
vet_test <- vet2 %>% filter(final_site =='Lip or oral cavity'|final_site =='Hypopharynx')
for_test <- table(vet_test$final_site, vet_test$cnv)
fisher_re <- fisher.test(for_test)
print(paste('hypo vs lip_oral',fisher_re$p.value, sep = ' : '))
vet_test <- vet2 %>% filter(final_site =='oropharynx'|final_site =='Hypopharynx')
for_test <- table(vet_test$final_site, vet_test$cnv)
fisher_re <- fisher.test(for_test)
print(paste('hypo vs oro',fisher_re$p.value, sep = ' : '))
vet_test <- vet2 %>% filter(final_site =='Lip or oral cavity'|final_site =='oropharynx')
for_test <- table(vet_test$final_site, vet_test$cnv)
fisher_re <- fisher.test(for_test)
print(paste('lip_oral vs oro',fisher_re$p.value, sep = ' : '))
vet_test <- vet2 %>% filter(final_site =='Lip or oral cavity'|final_site =='Larynx')
for_test <- table(vet_test$final_site, vet_test$cnv)
fisher_re <- fisher.test(for_test)
print(paste('lip_oral vs larynx',fisher_re$p.value, sep = ' : '))
vet_test <- vet2 %>% filter(final_site =='Larynx'|final_site =='oropharynx')
for_test <- table(vet_test$final_site, vet_test$cnv)
fisher_re <- fisher.test(for_test)
print(paste('larynx vs oro',fisher_re$p.value, sep = ' : '))



############ratio
vet_test <- vet2 %>% filter(final_site =='oropharynx'|final_site =='Hypopharynx'|final_site =='Lip or oral cavity'|final_site =='Larynx')
vet_test <- vet_test %>% 
  mutate(Smoking_status = ifelse((final_smoking == 'Current smoker')|(final_smoking == 'Former smoker'),'Smoking experience','Never smoking'))

vet_test$Smoking_status
vet_test$cnv
vet_test$cnv_amp
CNV_FreqbySite <- vet_test %>% 
  group_by(final_site) %>% 
  count(cnv) %>%
  mutate(Sum = sum(n)) %>% 
  mutate(Freq = n/Sum) %>%
  filter(cnv == "Yes")
AMP_FreqbySite <- vet_test %>% 
  group_by(final_site) %>% 
  count(cnv_amp) %>%
  mutate(Sum = sum(n)) %>% 
  mutate(Freq = n/Sum) %>%
  filter(cnv_amp == "Yes")
DEL_FreqbySite <- vet_test %>% 
  group_by(final_site) %>% 
  count(cnv_del) %>%
  mutate(Sum = sum(n)) %>% 
  mutate(Freq = n/Sum) %>%
  filter(cnv_del == "Yes")

HPV_FreqbySite <- vet_test %>% 
  group_by(final_site) %>% 
  count(final_HPV) %>%
  na.omit %>%
  mutate(Sum = sum(n)) %>% 
  mutate(Freq = n/Sum) %>%
  filter(final_HPV == "Negative")

TP53_FreqbySite <- vet_test %>% 
  group_by(final_site) %>% 
  count(TP53) %>%
  mutate(Sum = sum(n)) %>% 
  mutate(Freq = n/Sum) %>%
  filter(TP53 == TRUE)


Smoking_FreqbySite <- vet_test %>% 
  group_by(final_site) %>% 
  count(Smoking_status) %>%
  na.omit %>%
  mutate(Sum = sum(n)) %>% 
  mutate(Freq = n/Sum) %>%
  filter(Smoking_status == "Smoking experience")


Age_FreqbySite <- vet_test %>% 
  group_by(final_site) %>%
  summarise(mn = mean(final_age,na.rm = TRUE)) %>%
  mutate(Freq = mn/100)

####original
Total_FreqbySite <- rbind(CNV_FreqbySite$Freq,HPV_FreqbySite$Freq,TP53_FreqbySite$Freq)
colnames(Total_FreqbySite) <- TP53_FreqbySite$final_site
rownames(Total_FreqbySite) <- c("CNV","HPV","TP53")

Total_FreqbySite_df <- as.data.frame(melt(Total_FreqbySite))

Total_FreqbySite <- rbind(CNV_FreqbySite$Freq,HPV_FreqbySite$Freq,TP53_FreqbySite$Freq, Smoking_FreqbySite$Freq, Age_FreqbySite$Freq)
colnames(Total_FreqbySite) <- TP53_FreqbySite$final_site
rownames(Total_FreqbySite) <- c("CNV","HPV","TP53", "Smoking status", "Age")

Total_FreqbySite_df <- as.data.frame(melt(Total_FreqbySite))

####
Total_FreqbySite <- rbind(AMP_FreqbySite$Freq, DEL_FreqbySite$Freq, HPV_FreqbySite$Freq,TP53_FreqbySite$Freq, Smoking_FreqbySite$Freq, Age_FreqbySite$Freq)
colnames(Total_FreqbySite) <- TP53_FreqbySite$final_site
rownames(Total_FreqbySite) <- c("AMP", "DEL", "HPV","TP53", "Smoking status", "Age")

Total_FreqbySite_df <- as.data.frame(melt(Total_FreqbySite))


p3 <- ggplot(Total_FreqbySite_df, aes(x=Var2, y=value, group=Var1, color=Var1)) +
  geom_line(size=0.5)+
  geom_point(size=1)+
  scale_y_continuous(limits= c(0,0.95))+
  theme_classic(base_size = 3)+
  # geom_vline(aes(xintercept = Var2),
  #            size = 30, colour = alpha("grey", 0.3))+
  labs(title = "Ratio of CNV,HPV, and TP53 by Site", x = "Site", y= "Ratio") +
  theme(plot.title = element_text(hjust = 0.5,size=6,face="bold"),
        axis.text=element_text(size=5),
        axis.title=element_text(size=6))

# 
# pdf("Line_plot.pdf",width = 15, height = 10)
p3
# dev.off()

# 
# setwd("/data/project/TRIUMPH/Figure2")
# setwd("/data/project/TRIUMPH/Figure2/new_figure_manu/")
# pdf("Line_plot2.pdf",width = 3.7, height = 2.6)
# pdf("/data/project/TRIUMPH/Figure_swkim/Fig3C.pdf",width = 3.7, height = 2.6)
pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Fig3C_updateampdel.pdf",width = 3.7, height = 2.6)
p3
dev.off()


library(writexl)
# Load required libraries
library(dplyr)
# install.packages("charlatan")
library(charlatan)
##24.7.6 for rebutal letter
write_xlsx(total_clinical, path = "/Users/hwangso_macmini/Library/CloudStorage/GoogleDrive-hswbulls@gmail.com/My Drive/Work/project/HNSCC/data_rebutal.xlsx")
# Combine into a single data frame (if desired)
combined_df <- bind_rows(total_clinical)

# Identify columns that are lists
list_cols <- sapply(combined_df, is.list)

# Unlist those columns
combined_df[list_cols] <- lapply(combined_df[list_cols], unlist)

# # Function to generate a list of unique fake names
# generate_unique_names <- function(num_names) {
#   unique_names <- character(0)
#   while (length(unique_names) < num_names) {
#     new_names <- ch_name(num_names - length(unique_names))
#     unique_names <- unique(c(unique_names, new_names))
#   }
#   return(unique_names)
# }
# 
# # Generate a list of unique fake names
# num_rows <- nrow(combined_df)
# fake_names <- generate_unique_names(num_rows)
# 
# # Shuffle the fake names to ensure randomness
# set.seed(123)  # for reproducibility
# fake_names <- sample(fake_names)
# 
# # Add the fake_name column to the dataframe
# combined_df <- combined_df %>%
#   mutate(fake_name = fake_names)



# combined_df$fake_name
# Save as TSV
write.table(combined_df, 
            file = "/Users/hwangso_macmini/Library/CloudStorage/GoogleDrive-hswbulls@gmail.com/My Drive/Work/project/HNSCC/data_rebutal.tsv", 
            sep = "\t", 
            row.names = FALSE)
