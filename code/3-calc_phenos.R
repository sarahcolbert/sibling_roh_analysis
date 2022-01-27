###################################
############## PREP ###############
###################################

## load packages
library(tidyverse)

## load phenotype data
message(paste("loading ",(Sys.getenv("pheno_file")), sep=""))
pheno1 <- read.table(paste(Sys.getenv("pheno_file")), header = TRUE)
message(paste("done loading ",(Sys.getenv("pheno_file")), sep=""))

## load covariate data
message(paste("loading ",(Sys.getenv("covar_file")), sep=""))
covar_data <- read.table(paste(Sys.getenv("covar_file")), header = TRUE)
message(paste("done loading ",(Sys.getenv("covar_file")), sep=""))

## load roh file data to merge with phenotype data
message(paste("loading ",Sys.getenv("processed_dir"),"all_froh_data.txt", sep=""))
fam_file <- read.table(paste(Sys.getenv("processed_dir"),"all_froh_data.txt", sep=""), header = TRUE)
message(paste("done loading ",Sys.getenv("processed_dir"),"all_froh_data.txt", sep=""))

## clean phenotype data
## remove duplicate individuals
message("cleaning data")
pheno2 <- pheno1 %>% distinct(IID, .keep_all = TRUE)

## merge froh and phenotype data for all individuals (between analysis)
all_data_btwn <- merge(fam_file, pheno2, by ="IID")
## also make dataframe that only has phenotype data
id_fid_btwn <- all_data_btwn %>% select(IID, FID, age)
pheno_data_btwn <- merge(id_fid_btwn, pheno2, by = "IID")

## remove individuals that don't have within sibling froh values to create df to use in within analysis
all_data_within <- all_data_btwn %>% drop_na(froh_sibs)
## also make dataframe that only has phenotype data
id_fid_within <- all_data_within %>% select(IID, FID, age)
pheno_data_within <- merge(id_fid_within, pheno2, by = "IID")

message("done cleaning data")

###################################
### CALC STATS FOR BTWN SAMPLE ####
###################################

## create descriptive tables for each phenotype that give sample size and then other stats depending on if it is a binary or quantitative trait

## create empty df to hold descriptive stats
message("initializing between analysis description dataframe")
btwn_pheno_descrip <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(btwn_pheno_descrip) <- c("pheno_name", "sample_size", "ncase", "ncontrols", "mean_phen", "median_phen", "sd_phen", "mean_age", "sd_age")

## for each phenotype remove the NAs to get the number of IIDs and FIDs, as well as other stats
num_phenos_btwn <- ncol(pheno_data_btwn) ## will be used to calculate number of phenotypes
message("recording between analysis descriptive statistics for...")
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

message(paste("done recording between analysis descriptive statistics and writing to ",Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_btwn.csv", sep=""))

## write to file
write.csv(btwn_pheno_descrip, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_btwn.csv", sep=""), row.names = FALSE)

## write results to table so that easy to import for regression analysis
all_pheno_results <- pheno_data_btwn
write.table(all_pheno_results, paste(Sys.getenv("processed_dir"),"btwn_sibs_pheno_data.txt", sep=""), quote = FALSE, row.names=FALSE)

######################################
#### CALC WITHIN PHENOTYPE VALUES ####
######################################

message("calculating within siblings phenotype values for...")
## copy dataframe
within_sibs_phenos <- pheno_data_within
num_phenos_within <- ncol(pheno_data_within) ## will be used to calculate number of phenotypes

## calculate value of phenotype to family mean
## use for loop to go through each of the phenotype columns
for(j in 4:num_phenos_within){
message(paste(colnames(within_sibs_phenos)[j]))
  ## make empty column to hold results
  within_sibs_phenos$newcol <- NA
  newcolumn <- which(colnames(within_sibs_phenos)=="newcol")
  ## use for loop to get value for each individual
for(i in 1:length(within_sibs_phenos$IID)){
  spec_FID <- within_sibs_phenos$FID[i]
  fam_vals <- within_sibs_phenos %>% filter(FID==spec_FID)
  ## calculate residual phenotype value
  within_sibs_phenos$newcol[i] <- within_sibs_phenos[i,j]-mean(fam_vals[,j])
  ## this will change an individuals value to NA if any of their sibs are missing data
  within_sibs_phenos[i,j] <- ifelse(any(is.na(fam_vals[,j])), NA, within_sibs_phenos[i,j])
}
colnames(within_sibs_phenos)[newcolumn] <- paste(colnames(within_sibs_phenos)[j],"_sibs",sep="")
}

message(paste("done calculating within siblings phenotype values and writing to ",Sys.getenv("processed_dir"),"within_sibs_pheno_data.txt", sep=""))

## write results to table so that easy to import for regression analysis
sib_pheno_results <- within_sibs_phenos
sib_pheno_results[,4:num_phenos_within] <- NULL
write.table(sib_pheno_results, paste(Sys.getenv("processed_dir"),"within_sibs_pheno_data.txt", sep=""), quote = FALSE, row.names=FALSE)

## also save results to df with updated NAs for phenotypes not in relation to siblings
within_sample_phenos <- within_sibs_phenos
within_sample_phenos[,num_phenos_within+1:ncol(within_sample_phenos)] <- NULL

######################################
#### CALC STATS FOR WITHIN SAMPLE ####
######################################

message("initializing within analysis description dataframe")

## create descriptive tables for each phenotype that give sample size and then other stats depending on if it is a binary or quantitative trait

## create empty df to hold descriptive stats
within_pheno_descrip <- data.frame(matrix(ncol = 11, nrow = 0))

message("recording within analysis descriptive statistics for...")

## for each phenotype remove the NAs to get the number of IIDs and FIDs, as well as other stats
for(k in 4:num_phenos_within){
message(paste(colnames(within_sample_phenos)[k]))
test1 <- within_sample_phenos %>% select(FID, IID, age, colnames(within_sample_phenos)[k])
test2 <- test1 %>% drop_na()
sample_size <- length(unique(test2$IID))
mean_age <- mean(test2$age)
sd_age <- sd(test2$age)
nfams <- length(unique(test2$FID))

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
pheno_descrip <- cbind(pheno_name, sample_size, nfams, ncase, ncontrols, mean_phen, median_phen, sd_phen, mean_age, sd_age)
within_pheno_descrip <- rbind(within_pheno_descrip, pheno_descrip)
}

message(paste("done recording within analysis descriptive statistics and writing to ",Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_within.csv", sep=""))

## write to file
write.csv(within_pheno_descrip, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_within.csv", sep=""), row.names = FALSE)
