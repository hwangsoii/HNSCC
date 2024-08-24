source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim/Figure/import_data_for_Oncoprint_final.R")

library(data.table)
library(curl)
library(tableone)
# Load the packages
# install.packages("officer")
# install.packages("magrittr")
# install.packages("flextable")
library(survival) # only needed for the dataset in this example
library(dplyr) # to modify the needed dataframe
library(tibble) # for rownames_to_column() function
library(stringr) # for str_squish()
library(flextable)
library(officer)
library(forcats)
library(tableone)
library(kableExtra)
library(tableone)
library(flextable)
library(officer)
source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/customtab.R")
All_cnvkit
total_clinical
clinical_hpv <- total_clinical %>% filter(!is.na(final_HPV))
clinical_HPV_positive <-subset(clinical_hpv, clinical_hpv$final_HPV == 'Positive')
clinical_HPV_negative <- subset(clinical_hpv, clinical_hpv$final_HPV == 'Negative')
clinical_HPV_positive$final_HPV
clinical_HPV_negative$final_HPV

HPV_positive_barcode <- clinical_HPV_positive$Tumor_Sample_Barcode
HPV_negative_barcode <- clinical_HPV_negative$Tumor_Sample_Barcode

maf_HPV = subsetMaf(maf = All_cnvkit,tsb = HPV_positive_barcode)
somaticInteractions(maf = maf_HPV, top = 25, pvalue = c(0.05, 0.1))
t <- somaticInteractions(maf = maf_HPV, top = 25, pvalue = c(0.05, 0.1))
t <-  as.data.frame(t)
t
# Set Table header
header <- str_squish(str_remove("Supplementary Table 1. Mutually exclusive or co-occurring set of genes on patients with p16(+) detected by pair-wise Fisher’s Exact test.", "\n"))

# Set Table footer
footer <- str_squish(str_remove("Numbers are No. (%) unless otherwise noted", "\n"))

# Set custom_tab() defaults
customtab_defaults()
# Create the flextable object
flextable_2 <- custom_tab(t, header,footer)
# tab2_df %>%
#   kbl(caption = "Patient chracteristics") %>%
#   kable_classic(full_width = F, html_font = "Cambria")
save_as_docx(flextable_2, path = "/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/supplementary table 1_p16_podi.docx", 
             pr_section = 
               prop_section(page_size = page_size(orient = "portrait"), 
                            type = "continuous"))

# Set Table header
header <- str_squish(str_remove("Supplementary Table 2. Mutually exclusive or co-occurring set of genes on patients with p16(-) detected by pair-wise Fisher’s Exact test.", "\n"))

# Set Table footer
footer <- str_squish(str_remove("Numbers are No. (%) unless otherwise noted", "\n"))


maf_HPV_nega = subsetMaf(maf = All_cnvkit,tsb = HPV_negative_barcode)
somaticInteractions(maf = maf_HPV_nega, top = 25, pvalue = c(0.05, 0.1))
t <- somaticInteractions(maf = maf_HPV_nega, top = 25, pvalue = c(0.05, 0.1))
t <-  as.data.frame(t)
flextable_2 <- custom_tab(t, header,footer)
# tab2_df %>%
#   kbl(caption = "Patient chracteristics") %>%
#   kable_classic(full_width = F, html_font = "Cambria")
save_as_docx(flextable_2, path = "/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/supplementary table 1_p16_nega.docx", 
             pr_section = 
               prop_section(page_size = page_size(orient = "portrait"), 
                            type = "continuous"))
å
