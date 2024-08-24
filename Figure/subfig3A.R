source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")

source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/nanostring_import_data.R")
#################################################################################
#################################################################################
#############################Chi test on FC data#################################
#################################################################################
#################################################################################


amp_and_FC <- cnvkit_0.4[cnvkit_0.4[,1]%in%rownames(fil_FC),]
amp_and_FC <- amp_and_FC[(substr(amp_and_FC$sample,1,7)%in%colnames(fil_FC)),]
amp_and_FC <- unique(amp_and_FC)
amp_and_FC$amp_data <- apply(amp_and_FC, 1, function(x) fil_FC[x[1],substr(x[2],1,7)])

file_list <- sort(reboot_coln)
formed <- paste0(file_list,"-FP")
gene_list <- reboot_rown
cnvkit_fil <- cnvkit_0.4_onlyhnscc[cnvkit_0.4_onlyhnscc[,2]%in%formed,]
cnvkit_fil <- cnvkit_fil[cnvkit_fil[,1]%in%gene_list,]
cnvkit_fil <- unique(cnvkit_fil)
cnvkit_fil$amp_data <- apply(cnvkit_fil, 1, function(x) fil_FC[x[1], substr(x[2], 1, 7)])

for_chi <- function(x){
  total <- fil_FC[x,]
  total <- rbind(total, rep(x), paste0(names(fil_FC[x,]), "-FP"))
  amp <- cnvkit_fil[cnvkit_fil$gene == x,]
  amp <- amp[!(amp$amp_del == "Deletion"),]
  cnvkit_data <- cnvkit_0.4[cnvkit_0.4[,1]==x,]
  total_no_amp <- total[,!total[3,]%in%cnvkit_data[,2]]
  total_no_amp <- rbind(total_no_amp, rep("NA"))
  total_no_amp <- as.data.frame(t(total_no_amp))[,c(2,3,4,1)]
  colnames(total_no_amp) <- c("gene", "sample", "amp_del", "amp_data")
  colnames(amp) <- c("gene", "sample", "amp_del", "amp_data")
  
  ampvsno_amp <- rbind(amp, total_no_amp)
  ampvsno_amp$amp_del <- as.factor(ampvsno_amp$amp_del)
  ampvsno_amp$amp_data <- as.numeric(ampvsno_amp$amp_data)
  ampvsno_amp$real <- ifelse(ampvsno_amp$amp_data >1, 1, 0)
  ##ampvsno_amp <- ampvsno_amp[!(ampvsno_amp$amp_del =="Deletion"),]
  
  table <- table(ampvsno_amp[,c(3,5)])
  results <- list(table, ampvsno_amp)
  return(results)
}
intersect(row.names(fil_FC), cnvkit_0.4[,1])


for_chi_EGFR <- for_chi("EGFR")[[1]]
chisq.test(for_chi_EGFR, correct = T)  ### 2.2e-16
p1 <- fisher.test(for_chi_EGFR) ## <2.2e-16 lowest value
t1 <- for_chi("EGFR")[[2]]

for_chi_ERBB2 <- for_chi("ERBB2")[[1]]
chisq.test(for_chi_ERBB2, correct = T) ###  < 2.2e-16 
p2 <- fisher.test(for_chi_ERBB2) ##2.632e-09
t2 <- for_chi("ERBB2")[[2]]

for_chi_FGFR1 <- for_chi("FGFR1")[[1]]
chisq.test(for_chi_FGFR1, correct = T) ###no Amplification patient over logFC  > 1  ###0.7922
p3 <- fisher.test(for_chi_FGFR1) ##1
t3 <- for_chi("FGFR1")[[2]]

for_chi_FGFR2 <- for_chi("FGFR2")[[1]]
chisq.test(for_chi_FGFR2, correct = T) ###no NA patient over logFC  > 1  ###< 2.2e-16
p4 <- fisher.test(for_chi_FGFR2) ##0.0001583
t4 <- for_chi("FGFR2")[[2]]

for_chi_FGFR3 <- for_chi("FGFR3")[[1]]
chisq.test(for_chi_FGFR3, correct = T) ###only two Amplification patient  ###2.47e-06
p5 <- fisher.test(for_chi_FGFR3) ##0.01116
t5 <- for_chi("FGFR3")[[2]]

for_chi_PIK3CA <- for_chi("PIK3CA")[[1]]
chisq.test(for_chi_PIK3CA, correct = T)  ###0.0002393
p6 <- fisher.test(for_chi_PIK3CA) ##0.00109
t6 <- for_chi("PIK3CA")[[2]]

p.val <- cbind(p1$p.value,p2$p.value,p3$p.value,p6$p.value)
#test <- rbind(t1,t2,t3,t4,t5,t6)
test <- rbind(t1,t2,t3,t6)
a.adj <- p.adjust(p.val, method = "BH")
p <- function(x){
  a = wilcox.test(amp_data~amp_del, data = x)
  b = a$p.value
  return(b)
}
symbol <- function(x){
  if(x<=0.001){
    c = "***"
  }else if(b<=0.01){
    c = "**"
  }else if(b<=0.05){
    c = "*"
  }else{
    c = "NS"
  }
  return(c)
}

data = data.frame(xmin=c(0.81,1.81,2.81,3.81), xmax=c(1.19,2.19,3.19,4.19),
                  y_position=rep(6.35), p.v = p.adjust(c(p(t1),p(t2),p(t3),p(t6))))
data$q.v = p.adjust(data$p.v, method = "BH")
data$t.t = apply(data,1,function(x) if(as.numeric(x[5])<=0.001){c = "***"}
                 else if(as.numeric(x[5])<=0.01){c = "**"}
                 else if(as.numeric(x[5])<=0.05){c = "*"}else{c = "NS"})

fig8 <- ggplot(test, aes(x=gene, y=amp_data, fill = amp_del))+
  geom_boxplot(aes())+theme_bw() +   
  geom_signif(mapping=aes(x=xmin,xend=xmax, y=y_position, yend=y_position, annotation = t.t, group=c(1,2,3,4)), data = data,
              manual = T, tip_length = 0.01, inherit.aes = F, stat = "identity")+
  geom_text(aes(x= gene,y=amp_data,label=..count.., colour = amp_del), stat='count', size=4, y=6.0,
            position = position_dodge2(width = 0.75))+
  # scale_fill_grey(start = grey1, end = grey2)+
  # scale_color_grey(start = grey3, end = grey4)+
  labs(fill = "Amp or Not")+
  theme_classic()



pdf("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/Subfig3A.pdf", height = 10, width = 22)   ##########wilcox & Benjamini and Hochberg
fig8
dev.off()

  