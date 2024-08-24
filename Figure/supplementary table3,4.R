source("/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/R/final_R/afterswkim_p16/Figure/import_data_for_Oncoprint_final.R")
library(gtsummary) 
# install.packages("tidyr")
# install.packages("dplyr")
library(tidyr)
library(dplyr)
library(purrr)
library(tidyverse)
####################logistic regression
cnv_sample_barcode <- cnvkit_0.4_onlyhnscc$sample
total_clinical$cnv <-ifelse(total_clinical$Tumor_Sample_Barcode %in% cnv_sample_barcode,"Yes", "No")
clinical_practice <- total_clinical
clinical_practice2 <- total_clinical
clinical_practice2$cnv <- ifelse(clinical_practice2$cnv=="Yes", 1, 0)

####check smoking prevalance
clinical_practice2 <- clinical_practice2 %>% 
  mutate(Smoking_status = ifelse((final_smoking == 'Current smoker')|(final_smoking == 'Former smoker'),'Smoking experience','Never smoking'))


summary(as.factor(clinical_practice2[clinical_practice2$final_HPV == 'Positive',]$final_site))
# clinical_practice$final_site

# clinical_practice3 <- clinical_practice %>% filter(final_site !='oropharynx')
# clinical_practice3 <- clinical_practice3 %>% filter(final_site !='Lip or oral cavity')
# clinical_practice3 <- clinical_practice3 %>% 
#   mutate(Smoking_status = ifelse((final_smoking == 'Current smoker')|(final_smoking == 'Former smoker'),'Smoking experience','Never smoking'))
# clinical_practice3$cnv <- ifelse(clinical_practice3$cnv=="Yes", 1, 0)
# 
result <- glm(cnv ~ Smoking_status, family = binomial, data = clinical_practice2)
result <- glm(cnv ~ final_age + TP53  + Smoking_status,family = binomial, data = clinical_practice2)
summary(result)
summary(clinical_practice2$TP53)
summary(as.factor(clinical_practice2$Smoking_status))
# result <- glm(TP53 ~ final_age  + final_smoking, family = binomial, data = clinical_practice2)
# result <- glm(TP53 ~ final_age, family = binomial, data = clinical_practice2)
# result <- glm(cnv ~ final_smoking, family = binomial, data = clinical_practice2)
# result <- glm(cnv ~ Smoking_status, family = binomial, data = clinical_practice2)
# result <- glm(cnv ~ final_HPV, family = binomial, data = clinical_practice2)
# result <- glm(cnv ~ Recur, family = binomial, data = clinical_practice2)
# result <- glm(NOTCH1 ~ final_age, family = binomial, data = clinical_practice2)
# result <- glm(NOTCH1 ~ final_HPV, family = binomial, data = clinical_practice2)

summary(result)
# result$coefficients
# exp(cbind(coef(result), confint(result)))  
# exp(result$coefficients)
# exp(confint(result))
# 
# result <- glm(cnv ~ final_age + TP53  + Smoking_status, family = binomial, data = clinical_practice3)
# result <- glm(TP53 ~ final_age  + final_smoking, family = binomial, data = clinical_practice3)
# result <- glm(TP53 ~ final_age, family = binomial, data = clinical_practice3)
# result <- glm(cnv ~ final_smoking, family = binomial, data = clinical_practice3)
# result <- glm(cnv ~ Smoking_status, family = binomial, data = clinical_practice3)
# result <- glm(cnv ~ final_HPV, family = binomial, data = clinical_practice3)
# result <- glm(cnv ~ Recur, family = binomial, data = clinical_practice3)
# summary(result)
# install.packages("gtsummary")
# library(gtsummary)
# library(gtsummary) 
###########################Univariate
####################logistic regression
clinical_practice2 <- total_clinical
clinical_practice2$cnv <- ifelse(clinical_practice2$cnv=="Yes", 1, 0)
####check smoking prevalance
clinical_practice2 <- clinical_practice2 %>% 
  mutate(Smoking_status = ifelse((final_smoking == 'Current smoker')|(final_smoking == 'Former smoker'),'Smoking experience','Never smoking'))
##p16(+)
clinical_practice2 <- clinical_practice2[(clinical_practice2$final_HPV=="Positive") & !(is.na(clinical_practice2$final_HPV)),]
clinical_practice2$final_HPV

##p16(-)
clinical_practice2 <- total_clinical
clinical_practice2$cnv <- ifelse(clinical_practice2$cnv=="Yes", 1, 0)
####check smoking prevalance
clinical_practice2 <- clinical_practice2 %>% 
  mutate(Smoking_status = ifelse((final_smoking == 'Current smoker')|(final_smoking == 'Former smoker'),'Smoking experience','Never smoking'))
clinical_practice2 <- clinical_practice2[(clinical_practice2$final_HPV=="Negative") & !(is.na(clinical_practice2$final_HPV)),]
clinical_practice2$final_HPV
clinical_practice2$TP53_total <- as.factor(clinical_practice2$TP53_total)
clinical_practice2$Smoking_status <- as.factor(clinical_practice2$Smoking_status)
summary(clinical_practice2$Smoking_status)
summary(clinical_practice2$TP53_total)


trial_subset <-
  clinical_practice2 %>%
  select(cnv, Smoking_status, final_age, TP53) 
# trial_subset <-
#   clinical_practice2 %>%
#   select(cnv, Smoking_status, final_age, TP53) 

df_uv <-
  trial_subset %>%
  # group by trt, and nest data within group
  # group_by(cnv) %>%
  tidyr::nest() %>%
  # build univariate logistic regression models separately within grouping variable
  mutate(
    tbl_uv = map(
      data,
      ~tbl_uvregression(
        data = .x, 
        y = cnv,
        method = glm, 
        method.args = list(family = binomial),
        exponentiate = TRUE
      )
    )
  )


tbl_merge(
  tbls = df_uv$tbl_uv, # list of the univariate logistic regression tables
  # tab_spanner = paste0("**", df_uv$trt, "**") # adding stars to bold the header
)

# summary(result)
#multivariate
result <- glm(cnv ~ final_age + TP53  + Smoking_status, family = binomial,data = clinical_practice2)
tbl_regression(result,exponentiate = TRUE)

# 
# 
# tbl1 <- 
#   c("Smoking_status+final_age + TP53") %>%            # vector of covariates
#   map(
#     ~ paste("cnv", .x, sep = " ~ ") %>% # build character formula
#       as.formula() %>%                  # convert to proper formula
#       glm(data = clinical_practice2) %>%             # build linear regression model
#       tbl_regression()                  # display table with gtsummary
#   ) 
# 
# 
# tbl1
# # %>%
# # # merge tables into single table
# # tbl_merge(
# #   tab_spanner = c("**Univariate**", "**Multivariable**")
# # )
# 
# tbl2 <- tbl_merge(
#   tbls = tbl1 # list of the univariate logistic regression tables
#   # tab_spanner = paste0("**", df_uv$trt, "**") # adding stars to bold the header
# )
# tbl2
# 
# tbl3 <- 
#   c("Smoking_status", "final_age","TP53","CDKN2A","CDKN2A_cnv", "final_HPV", "Alcohol", "final_gender") %>%            # vector of covariates
#   map(
#     ~ paste("cnv", .x, sep = " ~ ") %>% # build character formula
#       as.formula() %>%                  # convert to proper formula
#       lm(data = clinical_practice2) %>%             # build linear regression model
#       tbl_regression()                  # display table with gtsummary
#   ) 
# tbl3 <- 
#   c("Smoking_status", "final_age","TP53","CDKN2A") %>%            # vector of covariates
#   map(
#     ~ paste("cnv", .x, sep = " ~ ") %>% # build character formula
#       as.formula() %>%                  # convert to proper formula
#       lm(data = clinical_practice2) %>%             # build linear regression model
#       tbl_regression()                  # display table with gtsummary
#   ) 
# tbl3
# # merge tables into single table
# tbl_merge(tbls = c(tbl3),
#           tab_spanner = c("**Univariable**", "**Multivariable**")
# )
# 
# tbl
# 
# 
# fisher.test(table(clinical_practice2$final_gender, clinical_practice2$Smoking_status))