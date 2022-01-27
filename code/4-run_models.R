###################################
############## PREP ###############
###################################

message("Loading packages")
library(tidyverse)
library(modelr)
library(lmerTest)
message("Done loading packages")

message(paste("Loading ",Sys.getenv("processed_dir"),"within_sibs_froh_data.txt", sep=""))
froh_data <- read.table(paste(Sys.getenv("processed_dir"),"within_sibs_froh_data.txt", sep=""), header = TRUE)
message(paste("Done loading ",Sys.getenv("processed_dir"),"within_sibs_froh_data.txt", sep=""))


message("Loading within and between phenotype data...")
within_phenotype_data <- read.table(paste(Sys.getenv("processed_dir"),"within_sibs_pheno_data.txt", sep=""), header = TRUE) %>% select(-age,-FID)
btwn_phenotype_data <- read.table(paste(Sys.getenv("processed_dir"),"btwn_sibs_pheno_data.txt", sep=""), header = TRUE) %>% select(-age,-FID)
message("Done loading within and between phenotype data")

## load covariate data
message(paste("Loading ",(Sys.getenv("covar_file")), sep=""))
covar_data <- read.table(paste(Sys.getenv("covar_file")), header = TRUE)
message(paste("Done loading ",(Sys.getenv("covar_file")), sep=""))

message("Cleaning data...")

## combine covariate and froh data
froh_covar <- merge(froh_data, covar_data, by="IID")

## make df for between analysis
btwn_data <- merge(froh_covar, btwn_phenotype_data, by="IID")

## make df for within analysis
within_data <- merge(froh_covar, within_phenotype_data, by="IID")

message("Done cleaning data")

###################################
######### BTWN ANALYSIS ###########
###################################

message("Starting between family analysis")

## copy df
btwn_data1 <- btwn_data

## calculate number of phenotypes
num_phenos_btwn <- ncol(btwn_data1)

## create empty df to hold regression estimates
all_results_btwn <- data.frame(matrix(ncol = 5, nrow = 0))

## loop through phenotypes
for(k in 19:num_phenos_btwn){
message(paste("Calculating betafroh in between family models for",colnames(btwn_data)[k], sep=" "))
  ## copy df and remove NAs for the phenotype
  btwn_data3 <- btwn_data %>% drop_na(paste(colnames(btwn_data[k])))
  ## scale froh and covariates
  btwn_data1 <- btwn_data3 %>% mutate_at(c("froh", "age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), scale)
  ## determine if binary or continuous phenotype
    if(all(btwn_data1[k]==0 | btwn_data1[k]==1)){
    message(paste("Reading",colnames(btwn_data)[k], "as a binary trait and running logistic regression", sep=" "))
      ## binary phenotype calculations
      ## Can just run logistic regression (Clark et al. equation 16)
      pheno_model <- glmer(formula(paste(colnames(btwn_data1)[k],'~ froh + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + (1 | FID)')), data = btwn_data1, family = binomial(link = 'logit'), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
      ## save coefficients
      phenotype <- colnames(btwn_data1)[k]
      beta <- summary(pheno_model)$coefficients[2,1]
      se <- summary(pheno_model)$coefficients[2,2]
      p <- summary(pheno_model)$coefficients[2,4]
      type <- "btwn_logistic"
            }else{
            message(paste("Reading",colnames(btwn_data)[k], "as a continuous trait and running linear regression", sep=" "))
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
              type <- "btwn_linear"
               }
    ## save results
    results <- as.data.frame(cbind(phenotype, beta, se, p, type))
    all_results_btwn <- rbind(all_results_btwn, results)
          }
message(paste("Regressions complete and writing results to ", Sys.getenv("output_dir"),Sys.getenv("output_name"),"_btwn_fam_analysis_results.csv", sep = ""))

## write results to csv file to be returned
write.csv(all_results_btwn, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_btwn_fam_analysis_results.csv", sep=""), row.names = FALSE)

###################################
######## WITHIN ANALYSIS ##########
###################################

message("Starting between family analysis")

## copy df
within_data1 <- within_data

## calculate number of phenotypes
num_phenos_within <- ncol(within_data1)

## create empty df to hold regression estimates
all_results_within <- data.frame(matrix(ncol = 4, nrow = 0))

## loop through phenotypes
for(k in 19:num_phenos_within){
message(paste("Calculating betafroh in within sibling models for",colnames(within_data1)[k], sep=" "))
  ## copy df and remove NAs for the phenotype
  within_data2 <- within_data1 %>% drop_na(paste(colnames(within_data1[k])))
  ## only run analysis if there's more than 250 families
  if(length(unique(within_data2$FID))>250){
  ## scale froh and covariates
  within_data3 <- within_data2 %>% mutate_at(c("froh_sibs", "age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), scale)
  ## already calculated within sibling values (Clark et al. equations 17 and 18) so just run regression
  pheno_model <- lm(formula(paste(colnames(within_data3)[k],'~ froh_sibs + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10')), data = within_data3)
      ## save coefficients
      phenotype <- colnames(within_data3)[k]
      beta <- summary(pheno_model)$coefficients[2,1]
      se <- summary(pheno_model)$coefficients[2,2]
      p <- summary(pheno_model)$coefficients[2,4]
      type <- "within"
    ## save results
    results <- as.data.frame(cbind(phenotype, beta, se, p, type))
    all_results_within <- rbind(all_results_within, results)
    }else{
    message(paste("WARNING:", colnames(within_data1)[k], " was not analyzed because phenotype data was not available for > 250 families", sep = ""))
    ## save coefficients
    phenotype <- colnames(within_data3)[k]
    beta <- NA
    se <- NA
    p <- NA
    type <- "error: results not available because Nfam < 250"
    ## save results
    results <- as.data.frame(cbind(phenotype, beta, se, p, type))
    all_results_within <- rbind(all_results_within, results)
    }
          }

message(paste("Regressions complete and writing results to ", Sys.getenv("output_dir"),Sys.getenv("output_name"),"_within_sibs_analysis_results.csv", sep = ""))


## write results to csv file to be returned
write.csv(all_results_within, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_within_sibs_analysis_results.csv", sep=""), row.names = FALSE)
