###################################
############## PREP ###############
###################################

## load packages
library(tidyverse)
library(lmerTest)

## load phenotype data
message(paste("Loading ",(Sys.getenv("pheno_file")), sep=""))
pheno1 <- read.table(paste(Sys.getenv("pheno_file")), header = TRUE)
message(paste("Done loading ",(Sys.getenv("pheno_file")), sep=""))

## load roh file data to merge with phenotype data
message(paste("Loading ",Sys.getenv("processed_dir"),"all_indivs_roh_covs.txt", sep=""))
roh_file <- read.table(paste(Sys.getenv("processed_dir"),"all_indivs_roh_covs.txt", sep=""), header = TRUE)
message(paste("Done loading ",Sys.getenv("processed_dir"),"all_indivs_roh_covs.txt", sep=""))

## clean phenotype data
## remove duplicate individuals
message("Cleaning data")
pheno2 <- pheno1 %>% distinct(IID, .keep_all = TRUE)

## merge froh and phenotype data for all individuals (between analysis)
all_data_btwn <- merge(roh_file, pheno2, by ="IID")
## also make dataframe that only has phenotype data
id_fid_age_btwn <- all_data_btwn %>% select(IID, FID, age)
pheno_data_btwn <- merge(id_fid_age_btwn, pheno2, by = "IID")

message("Done cleaning data")

###################################
### CALC STATS FOR BTWN SAMPLE ####
###################################

## create descriptive tables for each phenotype that give sample size and then other stats depending on if it is a binary or quantitative trait

## create empty df to hold descriptive stats
message("Initializing between analysis description dataframe")
btwn_pheno_descrip <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(btwn_pheno_descrip) <- c("pheno_name", "sample_size", "ncase", "ncontrols", "mean_phen", "median_phen", "sd_phen", "mean_age", "sd_age")

## for each phenotype remove the NAs to get the number of IIDs and FIDs, as well as other stats
num_phenos_btwn <- ncol(pheno_data_btwn) ## will be used to calculate number of phenotypes
message("Recording between analysis descriptive statistics for...")
for(k in 4:num_phenos_btwn){
message(paste(colnames(pheno_data_btwn)[k]))
test1 <- pheno_data_btwn %>% select(FID, IID, age, colnames(pheno_data_btwn)[k])
test2 <- test1 %>% drop_na()
sample_size <- length(unique(test2$IID))
mean_age <- mean(test2$age)
sd_age <- sd(test2$age)

## assesses if trait is binary or quantitative and pulls the correct info depending
if(all(test2[,4] %in% c(0,1))){
  ncase <- length(which(test2[,4]==1))
  ncontrols <- length(which(test2[,4]==0))
  mean_phen <- NA
  median_phen <- NA
  sd_phen <- NA
}else{
  ncase <- NA
  ncontrols <- NA
  mean_phen <- mean(test2[,4])
  median_phen <- median(test2[,4])
  sd_phen <- sd(test2[,4])}

pheno_name <- paste(colnames(test2)[4])
pheno_descrip <- cbind(pheno_name, sample_size, ncase, ncontrols, mean_phen, median_phen, sd_phen, mean_age, sd_age)
btwn_pheno_descrip <- rbind(btwn_pheno_descrip, pheno_descrip)
}

## write to file
message(paste("done recording between analysis descriptive statistics and writing to ",Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_btwn.csv", sep=""))
write.csv(btwn_pheno_descrip, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_btwn.csv", sep=""), row.names = FALSE)

###################################
####### RUN BTWN FAM MODELS #######
###################################

message("Starting between family analysis")

## create empty df to hold regression estimates
message("Initializing between analysis results dataframe")
all_results_btwn <- data.frame(matrix(ncol = 8, nrow = 0))
## make covariates df to attach to cleaned pheno df in for loop
covars <- roh_file %>% select(IID, age, sex, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

## run the between analysis
for(k in 4:num_phenos_btwn){
  message(paste("Running between family models for",colnames(pheno_data_btwn)[k], sep=" "))
  test2 <- pheno_data_btwn %>% select(FID, IID, colnames(pheno_data_btwn)[k]) %>% drop_na()

  ####### Step 1: assesses if trait is binary or continuous
  ## make pheno column
  test2$pheno <- test2[,3]
  ## name type based on if binary or continuous
  type <- ifelse(all(test2[,3] %in% c(0,1)), "binary", "continuous")

  ####### Step 2: attach covariates to df
  test3 <- merge(test2, covars, by = "IID")
  ## read sex as factor
  test3$sex <- as.factor(test3$sex)

  ####### Step 3: regress scaled phenotypes on covariates to get residuals
  cov_model <- lm(pheno ~ scale(age) + (scale(age))^2 + sex + scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7) + scale(PC8) + scale(PC9) + scale(PC10), data = test3)

  ####### Step 4: pull residuals from the model
  pheno_resids <- test3
  pheno_resids$resids <- resid(summary(cov_model))

  ####### Step 5: calculate family pheno mean using residuals
  test4 <- pheno_resids %>%
  group_by(FID) %>%
  summarise_at(vars(resids), list(fam_pheno = mean))

  ####### Step 6: calculate family froh mean values
  iid_fid_froh <- roh_file %>% select(IID, FID, froh, Fhat3) %>% filter(IID %in% pheno_resids$IID)
  test5 <- iid_fid_froh %>%
  group_by(FID) %>%
  summarise_at(vars(froh, Fhat3), list(fam = mean))

  ####### Step 8: combine mean phenotype and froh values into df
  means_df <- merge(test4, test5, by = "FID")

  ####### Step 9: run models
  froh_btwn_model <- lm(fam_pheno ~ froh_fam, data = means_df)
  fhat3_btwn_model <- lm(fam_pheno ~ Fhat3_fam, data = means_df)
  multi_btwn_model <- lm(fam_pheno ~ froh_fam + Fhat3_fam, data = means_df)

  ####### Step 10: save values
  phenotype <- colnames(pheno_data_btwn)[k]
  beta_froh <- summary(froh_btwn_model)$coefficients[2,1]
  se_froh <- summary(froh_btwn_model)$coefficients[2,2]
  p_froh <- summary(froh_btwn_model)$coefficients[2,4]
  beta_fhat3 <- summary(fhat3_btwn_model)$coefficients[2,1]
  se_fhat3 <- summary(fhat3_btwn_model)$coefficients[2,2]
  p_fhat3 <- summary(fhat3_btwn_model)$coefficients[2,4]
  beta_multi_froh <- summary(multi_btwn_model)$coefficients[2,1]
  se_multi_froh <- summary(multi_btwn_model)$coefficients[2,2]
  p_multi_froh <- summary(multi_btwn_model)$coefficients[2,4]
  beta_multi_fhat3 <- summary(multi_btwn_model)$coefficients[3,1]
  se_multi_fhat3 <- summary(multi_btwn_model)$coefficients[3,2]
  p_multi_fhat3 <- summary(multi_btwn_model)$coefficients[3,4]

  ####### Step 10: save results
  results <- as.data.frame(cbind(phenotype, beta_froh, se_froh, p_froh, beta_fhat3, se_fhat3, p_fhat3, beta_multi_froh, se_multi_froh, p_multi_froh, beta_multi_fhat3, se_multi_fhat3, p_multi_fhat3, type))
  all_results_btwn <- rbind(all_results_btwn, results)
}

message(paste("Between family analysis completed and writing results to ", Sys.getenv("output_dir"),Sys.getenv("output_name"),"_btwn_fam_analysis_results.csv", sep = ""))

## write results to csv file to be returned
write.csv(all_results_btwn, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_btwn_fam_analysis_results.csv", sep=""), row.names = FALSE)
