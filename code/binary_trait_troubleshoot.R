###################################
############## PREP ###############
###################################

## load packages
library(tidyverse)
library(predictmeans)
library(modelr)
library(lmerTest)

## load froh data
froh_data <- read.table(paste(Sys.getenv("processed_dir"),"within_sibs_froh_data.txt", sep=""), header = TRUE)

## load phenotype data
within_phenotype_data <- read.table(paste(Sys.getenv("processed_dir"),"within_sibs_pheno_data.txt", sep=""), header = TRUE)
btwn_phenotype_data <- read.table(paste(Sys.getenv("processed_dir"),"btwn_sibs_pheno_data.txt", sep=""), header = TRUE) %>% select(-age,-FID)

## load covariate data
message(paste("loading ",(Sys.getenv("covar_file")), sep=""))
covar_data <- read.table(paste(Sys.getenv("covar_file")), header = TRUE)
message(paste("done loading ",(Sys.getenv("covar_file")), sep=""))


## combine covariate and froh data
froh_covar <- merge(froh_data, covar_data, by="IID")

## make df for between analysis
btwn_data <- merge(froh_covar, btwn_phenotype_data, by="IID")

## make df for within analysis
within_data <- merge(froh_covar, within_phenotype_data, by="IID")

###################################
######### BTWN ANALYSIS ###########
###################################

## copy df
btwn_data1 <- btwn_data

## calculate number of phenotypes
num_phenos_btwn <- ncol(btwn_data1)

## create empty df to hold descriptive stats
all_results <- data.frame(matrix(ncol = 5, nrow = 0))

## loop through phenotypes
for(k in 19:num_phenos_btwn){
  ## copy df and remove NAs for the phenotype
  btwn_data3 <- btwn_data %>% drop_na(paste(colnames(btwn_data[k])))
  ## scale froh and covariates
  btwn_data1 <- btwn_data3 %>% mutate_at(c("froh", "age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), scale)
  ## determine if binary or continuous phenotype
    if((btwn_data1[k]==0 | btwn_data1[k]==1)){
      ## binary phenotype calculations
      ## run linear regression regressing froh on the phenotype
      pheno_model <- lmer(formula(paste('froh ~ ', colnames(btwn_data1)[k],'+ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + (1 | FID)')), data = btwn_data1)
      ## create base model to use for permutations
      base_model <- lmer(froh ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + (1 | FID), data = btwn_data1)
      ## run permutations
      perm_results<- permlmer(base_model, pheno_model, nperm=10000, plot=FALSE, seed=58)
      ## save coefficients
      phenotype <- colnames(btwn_data1)[k]
      beta <- summary(pheno_model)$coefficients[2,1]
      se <- summary(pheno_model)$coefficients[2,2]
      p <- summary(pheno_model)$coefficients[2,4]
      perm_p <- perm_results[2,9]
      type <- "linear_flipped"
            }else{
              ## continuous phenotypes calculations
              ####### Step 1: regress phenotypes on covariates to get residuals (Clark et al. equation 11)
              pheno_model <- lmer(formula(paste(colnames(btwn_data1)[k],'~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + (1 | FID)')), data = btwn_data1)
              ####### Step 2: pull residuals from the model
              pheno_resids <- add_residuals(btwn_data1, pheno_model, var = "resids")
              ####### Step 3: regress residuals on froh (Clark et al. equation 13a)
              resids_model <- lm(resids ~ froh, data = pheno_resids)
              ####### Step 4: save coefficients
              phenotype <- colnames(btwn_data1)[k]
              beta <- summary(resids_model)$coefficients[2,1]
              se <- summary(resids_model)$coefficients[2,2]
              p <- summary(resids_model)$coefficients[2,4]
              perm_p <- NA
              type <- "linear"
               }
    ## save results
    results <- as.data.frame(cbind(phenotype, beta, se, p, perm_p, type))
    all_results <- rbind(all_results, results)
          }

write.table(all_results, "all_results_test2.txt", row.names = FALSE)


###################################
######## WITHIN ANALYSIS ##########
###################################

## MAKE SURE TO REMOVE PHENOTYPES WITH < 250 FAMILIES
