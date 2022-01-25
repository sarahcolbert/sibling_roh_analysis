###################################
############## PREP ###############
###################################

## load packages
library(tidyverse)
library(lmerTest)

## load froh data 
froh_data <- read.table(paste(Sys.getenv("processed_dir"),"within_sibs_froh_data.txt", sep=""), header = TRUE)

## scale froh 
scaled_froh <- froh_data %>% mutate_at(froh, scale)

## load phenotype data
within_phenotype_data <- read.table(paste(Sys.getenv("processed_dir"),"within_sibs_pheno_data.txt", header = TRUE)
btwn_phenotype_data <- read.table(paste(Sys.getenv("processed_dir"),"btwn_sibs_pheno_data.txt", sep=""), header = TRUE)

## load covariate data
message(paste("loading ",(Sys.getenv("covar_file")), sep=""))
covar_data <- read.table(paste(Sys.getenv("covar_file")), header = TRUE)
message(paste("done loading ",(Sys.getenv("covar_file")), sep=""))

## scale covariate data
scaled_covar <- covar_data %>% mutate_at(c(2,4:13), scale)

## combine covariate and froh data
froh_covar <- merge(scaled_froh, scaled_covar, by=="IID")

## make df for between analysis
btwn_data <- merge(btwn_phenotype_data, froh_covar, by=="IID")

## make df for within analysis
within_data <- merge(within_phenotype_data, froh_covar, by=="IID")

###################################
######### BTWN ANALYSIS ###########
###################################

## copy df
btwn_data1 <- btwn_data

## calculate number of phenotypes
num_phenos_btwn <- ncol(btwn_data1)

for(k in X:num_phenos_btwn){
## copy df
btwn_data1 <- btwn_data
####### Step 1: regress phenotypes on covariates to get residuals
  pheno_model <- lmer(paste(colnames(btwn_data1)[k]) ~ age + sex + PC1 + PC2 + Pc3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = btwn_data1)
####### Step 2: pull residuals from the model
  pheno_resids <- add_residuals(btwn_data1, pheno_model, var = "resids")
####### Step 3: regress residuals on froh
    resids_model <- lmer(resids ~ froh, data = pheno_resids)
####### Step 4: save results
summary(resids_model)
}
