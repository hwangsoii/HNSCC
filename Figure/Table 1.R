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
source("D:/황신원/프로젝트/HNSCC/R/customtab.R")
getwd()
# install.packages(c("officer", "ReporteRs"))

source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")
clinical_practice <- total_clinical

library(kableExtra)
library(tableone)
library(flextable)
library(officer)
clinical_practice$Stage  [clinical_practice$Stage   == '0'] <- "I"
clinical_practice$Stage  [clinical_practice$Stage   == '1'] <- "I"
clinical_practice$Stage  [clinical_practice$Stage   == '2'] <- "II"
clinical_practice$Stage  [clinical_practice$Stage   == '3'] <- "III"
clinical_practice$Stage  [clinical_practice$Stage   == '4'] <- "IVA"
clinical_practice$Stage  [clinical_practice$Stage   == '4a'] <- "IVA"
clinical_practice$Stage  [clinical_practice$Stage   == '4A'] <- "IVA"
clinical_practice$Stage  [clinical_practice$Stage   == '4B'] <- "IVB"
clinical_practice$Stage  [clinical_practice$Stage   == '4C'] <- "IVC"
clinical_practice[is.na(clinical_practice$final_smoking),]$final_smoking <- "Not available"
clinical_practice[is.na(clinical_practice$final_HPV),]$final_HPV <- "Not available"
dt <- clinical_practice
dt$Sex <- as.factor(dt$Sex)
dt$Stage
# dt[is.na(dt$Stage),]$Stage<- 'IVA'
dt$Stage
listVars <- c("final_age", "Sex", "Stage", "final_smoking", "final_HPV")
catVars <- c("Sex", "Stage", "final_smoking",  "final_HPV")
tab1 <- dt %>% CreateTableOne(vars = listVars, 
                              data = . , 
                              factorVars = catVars, 
                              # strata = "Hormone treatment", 
                              addOverall = T,
                              test = T)
tab1
tab1_word <- print(tab1, 
                   # nonnormal = c("Progesterone receptors, fmol/L (median [IQR])", 
                   #               "Estrogen receptors, fmol/L (median [IQR])"),
                   quote = F, 
                   noSpaces = T, 
                   # smd = T, 
                   # missing = T, 
                   test = F, 
                   contDigits = 1, 
                   printToggle = F,
                   dropEqual = T, 
                   explain = F
                   , showAllLevels = T)
tab1_df <- as.data.frame(tab1_word) %>% rownames_to_column(var = "Variable")
# Rename first variable from n to No.
tab1_df$Variable[1] <- "No."

# Set Table header
header <- str_squish(str_remove("Table 1. Baseline characteristics of 419 patients enrolled in the TRIUMPH project", "\n"))

# Set Table footer
footer <- str_squish(str_remove("Numbers are No. (%) unless otherwise noted", "\n"))

# Set custom_tab() defaults
customtab_defaults()

# Create the flextable object
flextable_1 <- custom_tab(tab1_df, header, footer)
tab1_df %>%
  kbl(caption = "Patient chracteristics") %>%
  kable_classic(full_width = F, html_font = "Cambria")
save_as_docx(flextable_1, path = "/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/table_1_p16.docx", 
             pr_section = 
               prop_section(page_size = page_size(orient = "portrait"), 
                            type = "continuous"))
#####starata
listVars <- c("final_age", "Sex", "Stage", "final_smoking", "final_HPV")
catVars <- c("Sex", "Stage", "final_smoking",  "final_HPV")
dt <- clinical_practice
# dt[is.na(dt$Stage),]$Stage<- 'IVA'
tab2 <- dt %>% CreateTableOne(vars = listVars, 
                              data = . , 
                              factorVars = catVars, 
                              strata = "final_site",
                              addOverall = T,
                              test = T)
tab2
tab2_word <- print(tab2, 
                   # nonnormal = c("Progesterone receptors, fmol/L (median [IQR])", 
                   #               "Estrogen receptors, fmol/L (median [IQR])"),
                   quote = F, 
                   noSpaces = T, 
                   # smd = T, 
                   # missing = T, 
                   test = F, 
                   contDigits = 1, 
                   printToggle = F,
                   dropEqual = T, 
                   explain = F
                   , showAllLevels = T)
tab2_df <- as.data.frame(tab2_word) %>% rownames_to_column(var = "Variable")
# Rename first variable from n to No.
tab2_df$Variable[1] <- "No."

# Set Table header
header <- str_squish(str_remove("Table 1. Baseline characteristics of 419 patients enrolled in the TRIUMPH project", "\n"))

# Set Table footer
footer <- str_squish(str_remove("Numbers are No. (%) unless otherwise noted", "\n"))

# Set custom_tab() defaults
customtab_defaults()

# Create the flextable object
flextable_2 <- custom_tab(tab2_df, header, footer)
tab2_df %>%
  kbl(caption = "Patient chracteristics") %>%
  kable_classic(full_width = F, html_font = "Cambria")
save_as_docx(flextable_2, path = "/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC/수정2/추가수정4/Figure/table 1.docx", 
             pr_section = 
               prop_section(page_size = page_size(orient = "portrait"), 
                            type = "continuous"))

