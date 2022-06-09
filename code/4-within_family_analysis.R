###################################
############## PREP ###############
###################################

## load packages
library(tidyverse)
library(lmerTest)
library(lfe)

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

## remove individuals that don't have within sibling froh values to create df to use in within analysis
all_data_within <- all_data_btwn %>% filter(wsib==1)
## also make dataframe that only has phenotype data
id_fid_age_within <- all_data_within %>% select(IID, FID, age)
pheno_data_wsibs <- merge(id_fid_age_within, pheno2, by = "IID")

## if trait is binary, check for N case = 500, and remove any phenos below that
## make vector with phenotype names
pheno_names <- colnames(pheno_data_wsibs)[4:ncol(pheno_data_wsibs)]
message("Checking for phenotypes with N case < 500")
for(s in 1:length(pheno_names)){
  ## pull phenotype values into a vector
  pheno_vals1 <- pheno_data_wsibs[,pheno_names[s]]
  ## remove NAs
  pheno_vals1<-pheno_vals1[!is.na(pheno_vals1)]
  ## test if binary, and if binary check how many cases (=1)
  res <- ifelse(all(pheno_vals1 %in% c(0,1)), sum(pheno_vals1), NA)
  ## if its binary and res is number, check that number isn't < 500
  if(is.na(res)){ ## do nothing
    }else if (res < 500){
  message(paste("WARNING:", pheno_names[s], "has less than 500 cases and is being removed from analysis", sep=" "))
  pheno_data_wsibs[,pheno_names[s]] <- NULL}
}

message("Done cleaning data")

###################################
### CALC STATS FOR WSIBS SAMPLE ###
###################################

## create descriptive tables for each phenotype that give sample size and then other stats depending on if it is a binary or quantitative trait

## create empty df to hold descriptive stats
message("Initializing within sibs description dataframe")
wsibs_pheno_descrip <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(wsibs_pheno_descrip) <- c("pheno_name", "sample_size", "nfams", "ncase", "ncontrols", "mean_phen", "median_phen", "sd_phen", "mean_age", "sd_age")

## initialize empty df to hold cleaned values
sibs_pheno_data <- as.data.frame(pheno_data_wsibs %>% select(IID,FID,age))

## for each phenotype remove the NAs to get the number of IIDs and FIDs, as well as other stats
num_phenos_wsibs <- ncol(pheno_data_wsibs) ## will be used to calculate number of phenotypes
message("Recording within sibs analysis descriptive statistics for...")
for(k in 4:num_phenos_wsibs){
message(paste(colnames(pheno_data_wsibs)[k]))
test1 <- pheno_data_wsibs

## make empty column to hold results
test1$pheno_cleaned <- NA
## use for loop to get value for each individual
for(i in 1:length(test1$IID)){
spec_FID <- test1$FID[i]
fam_vals <- test1 %>% filter(FID==spec_FID)
## add value to cleaned phenotype column
## this will change an individual's value to NA if they have no siblings that survived to this point or they are the only sibling with data (i.e. all other sibs = NA)
test1$pheno_cleaned[i] <- ifelse(length(fam_vals[,k])<=1 | sum(!is.na(fam_vals[,k])) <= 1, NA, test1[i,k])
}
test2 <- test1 %>% select(IID, FID, age, pheno_cleaned) %>% drop_na()
## calculate values now that df is cleaned
sample_size <- length(unique(test2$IID))
nfams <- length(unique(test2$FID))
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

## gather descriptive info into df
pheno_name <- paste(colnames(test2)[4])
pheno_descrip <- cbind(pheno_name, sample_size, nfams, ncase, ncontrols, mean_phen, median_phen, sd_phen, mean_age, sd_age)
wsibs_pheno_descrip <- rbind(wsibs_pheno_descrip, pheno_descrip)

## create df with cleaned phenotypes
sibs_within <- test1 %>% select(IID, pheno_cleaned)
sibs_pheno_data <- merge(sibs_pheno_data, sibs_within, by = "IID", all = TRUE)
newcolumn <- which(colnames(sibs_pheno_data)=="pheno_cleaned")
colnames(sibs_pheno_data)[newcolumn] <- paste(colnames(pheno_data_wsibs)[k],"_sibs",sep="")
}

## write to file
message(paste("done recording within sibs analysis descriptive statistics and writing to ",Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_wsibs.csv", sep=""))
write.csv(wsibs_pheno_descrip, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_wsibs.csv", sep=""), row.names = FALSE)

###################################
######## RUN WSIBS MODELS #########
###################################
message("Starting within sibling analysis")

## create empty df to hold regression estimates
message("Initializing within sibs analysis results dataframe")
all_results_wsibs <- data.frame(matrix(ncol = 8, nrow = 0))

## make covariates df to attach to cleaned pheno df in for loop
covars <- roh_file %>% select(IID, age, sex, froh, Fhat3)

## run models for each phenotype
for(j in 4:ncol(sibs_pheno_data)){
    message(paste("Running within sibs models for",colnames(sibs_pheno_data)[j], sep=" "))
  ####### Step 1: create dataframe with just cleaned phenotype
  test1 <- sibs_pheno_data %>% select(FID, IID, colnames(sibs_pheno_data)[j]) %>% drop_na()
  ####### Step 2: combine with covariates and froh data
  test2 <- merge(test1, covars, by = "IID")
  ####### Step 3: read FID and sex as a factor
  test2$FID <- as.factor(test2$FID)
  test2$sex <- as.factor(test2$sex)
  ####### Step 4: run regression using lfe package
  froh_wsibs_model <- felm(formula(paste(colnames(sibs_pheno_data)[j],'~ froh + scale(age) + (scale(age))^2 + sex | FID')), data = test2)
  fhat3_wsibs_model <- felm(formula(paste(colnames(sibs_pheno_data)[j],'~ Fhat3 + scale(age) + (scale(age))^2 + sex | FID')), data = test2)
  multi_wsibs_model <- felm(formula(paste(colnames(sibs_pheno_data)[j],'~ froh + Fhat3 + scale(age) + (scale(age))^2 + sex | FID')), data = test2)
  ####### Step 5: save values
  phenotype <- colnames(sibs_pheno_data)[j]
  beta_froh <- summary(froh_wsibs_model)$coefficients[2,1]
  se_froh <- summary(froh_wsibs_model)$coefficients[2,2]
  p_froh <- summary(froh_wsibs_model)$coefficients[2,4]
  beta_fhat3 <- summary(fhat3_wsibs_model)$coefficients[2,1]
  se_fhat3 <- summary(fhat3_wsibs_model)$coefficients[2,2]
  p_fhat3 <- summary(fhat3_wsibs_model)$coefficients[2,4]
  beta_multi_froh <- summary(multi_wsibs_model)$coefficients[2,1]
  se_multi_froh <- summary(multi_wsibs_model)$coefficients[2,2]
  p_multi_froh <- summary(multi_wsibs_model)$coefficients[2,4]
  beta_multi_fhat3 <- summary(multi_wsibs_model)$coefficients[3,1]
  se_multi_fhat3 <- summary(multi_wsibs_model)$coefficients[3,2]
  p_multi_fhat3 <- summary(multi_wsibs_model)$coefficients[3,4]
  ####### Step 6: save results
  results <- as.data.frame(cbind(phenotype, beta_froh, se_froh, p_froh, beta_fhat3, se_fhat3, p_fhat3, beta_multi_froh, se_multi_froh, p_multi_froh, beta_multi_fhat3, se_multi_fhat3, p_multi_fhat3))
  all_results_wsibs <- rbind(all_results_wsibs, results)
          }

message(paste("Within sibling analysis completed and writing results to ", Sys.getenv("output_dir"),Sys.getenv("output_name"),"_wsibs_analysis_results.csv", sep = ""))

## write results to csv file to be returned
write.csv(all_results_wsibs, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_wsibs_analysis_results.csv", sep=""), row.names = FALSE)
