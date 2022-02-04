###################################
############## PREP ###############
###################################

## load packages
library(tidyverse)
library(lmerTest)
library(modelr)

## load phenotype data
message(paste("Loading ",(Sys.getenv("pheno_file")), sep=""))
pheno1 <- read.table(paste(Sys.getenv("pheno_file")), header = TRUE)
message(paste("Done loading ",(Sys.getenv("pheno_file")), sep=""))

## load covariate data
message(paste("Loading ",(Sys.getenv("covar_file")), sep=""))
covar_data <- read.table(paste(Sys.getenv("covar_file")), header = TRUE)
message(paste("Done loading ",(Sys.getenv("covar_file")), sep=""))

## load roh file data to merge with phenotype data
message(paste("Loading ",Sys.getenv("processed_dir"),"all_froh_data.txt", sep=""))
fam_file <- read.table(paste(Sys.getenv("processed_dir"),"all_froh_data.txt", sep=""), header = TRUE)
message(paste("Done loading ",Sys.getenv("processed_dir"),"all_froh_data.txt", sep=""))

## clean phenotype data
## remove duplicate individuals
message("Cleaning data")
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

message(paste("done recording between analysis descriptive statistics and writing to ",Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_btwn.csv", sep=""))

## write to file
write.csv(btwn_pheno_descrip, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_btwn.csv", sep=""), row.names = FALSE)

## write results to table so that easy to import for regression analysis
all_pheno_results <- pheno_data_btwn
write.table(all_pheno_results, paste(Sys.getenv("processed_dir"),"btwn_sibs_pheno_data.txt", sep=""), quote = FALSE, row.names=FALSE)

######################################
#### CALC WITHIN PHENOTYPE VALUES ####
######################################

message("Calculating within phenotype values")

## merge in covariate data
pheno_data_within1 <- pheno_data_within %>% select(-age)
test_df <- merge(covar_data,pheno_data_within1, by = "IID")

## copy dataframe
within_data1 <- test_df
num_phenos_within <- ncol(within_data1) ## will be used to calculate number of phenotypes

## STEP 1: CALCULATE TRAIT RESIDUALS

## create empty df to hold residuals
trait_residuals_within <- as.data.frame(pheno_data_within %>% select(FID, IID, age))

message("Step 1: Calculating trait residuals...")

for(m in 15:ncol(within_data1)){
  ## copy df and remove NAs for the phenotype
  within_data2 <- within_data1 %>% drop_na(paste(colnames(within_data1[m])))
  ## scale covariates
  within_data3 <- within_data2 %>% mutate_at(c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), scale)
  ## determine if binary or continuous phenotype
    if(all(within_data3[m]==0 | within_data3[m]==1)){
    message(paste("Reading",colnames(within_data3)[m], "as a binary trait so will not calculate trait residuals", sep=" "))
      ## binary phenotype stays the same
      within_data3$newcol <- within_data3[,m]
      binary_resids <- within_data3 %>% select(IID, newcol)
      ## add binary phenotype value to "residuals" column
      trait_residuals_within <- merge(trait_residuals_within, binary_resids, by = "IID", all = TRUE)
            }else{
            message(paste("Reading",colnames(within_data3)[m], "as a continuous trait and calculating trait residuals", sep=" "))
              ## continuous phenotypes calculations
              ####### Step 1: regress phenotypes on covariates to get residuals (Clark et al. equation 11)
              resids_model <- lm(formula(paste(colnames(within_data3)[m],'~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10')), data = within_data3)
              ####### Step 2: pull residuals from the model
              pheno_resids <- add_residuals(within_data3, resids_model, var = "newcol")
              cont_resids <- pheno_resids %>% select(IID, newcol)
              ## add residual to residuals dataframe
              trait_residuals_within <- merge(trait_residuals_within, cont_resids, by = "IID", all = TRUE)
               }
    ## rename residuals column to trait it corresponds to
    newcolumn <- which(colnames(trait_residuals_within)=="newcol")
    colnames(trait_residuals_within)[newcolumn] <- paste(colnames(within_data3)[m],"_resids",sep="")
          }

message("Done with Step 1")

message("Step 2: Calculating trait residuals relative to family mean for...")

## create empty df to hold residuals
sib_resids_within <- as.data.frame(trait_residuals_within %>% select(IID,FID, age))

for(j in 4:ncol(trait_residuals_within)){
  ## copy df and remove NAs for the phenotype
  sib_resids_within2 <- trait_residuals_within
  ## print phenotype
  message(paste(colnames(sib_resids_within2)[j]))
  ## make empty column to hold results
  sib_resids_within2$newcol <- NA
  ## use for loop to get value for each individual
for(i in 1:length(sib_resids_within2$IID)){
  spec_FID <- sib_resids_within2$FID[i]
  fam_vals <- sib_resids_within2 %>% filter(FID==spec_FID)
  ## calculate residual phenotype value
  sib_resids_within2$newcol[i] <- sib_resids_within2[i,j]-mean(fam_vals[,j])
  ## this will change an individual's value to NA if any of their sibs are missing data
  sib_resids_within2[i,j] <- ifelse(any(is.na(fam_vals[,j])), NA, sib_resids_within2[i,j])
}
sib_resids_within3 <- sib_resids_within2 %>% select(IID, newcol)
sib_resids_within <- merge(sib_resids_within, sib_resids_within3, by = "IID", all = TRUE)
newcolumn <- which(colnames(sib_resids_within)=="newcol")
colnames(sib_resids_within)[newcolumn] <- paste(colnames(trait_residuals_within)[j],"_sibs",sep="")
          }

message(paste("Done with step 2 and writing sibling phenotype values to ",Sys.getenv("processed_dir"),"within_sibs_pheno_data.txt", sep=""))

## write results to table so that easy to import for regression analysis
## remove any single individuals that made it through because sib had genetic data but not pheno data
sib_pheno_results <- sib_resids_within[sib_resids_within$FID %in% sib_resids_within$FID[duplicated(sib_resids_within$FID)],]
write.table(sib_pheno_results, paste(Sys.getenv("processed_dir"),"within_sibs_pheno_data.txt", sep=""), quote = FALSE, row.names=FALSE)

#####################
## need to make df with raw phenotype values but that considers if there's an NA in the family so that our descriptive stats are accurate
## combine the residual values and raw phenotype values
within_complete_df <- merge(pheno_data_within, sib_pheno_results, by = c("IID", "FID", "age"), all = TRUE)
## find number of phenotypes
num_phenos <- length(15:ncol(within_data1))
for(s in 4:(3+num_phenos)){
  for(i in 1:length(within_complete_df$IID)){
    ## if trait_sibs is NA then make trait NA
    within_complete_df[i,s] <- ifelse(is.na(within_complete_df[i,s+num_phenos]), NA, within_complete_df[i,s])
  }
}
within_sample_phenos <- within_complete_df
within_sample_phenos[,num_phenos+3:ncol(within_sample_phenos)] <- NULL

######################################
#### CALC STATS FOR WITHIN SAMPLE ####
######################################

message("Initializing within analysis description dataframe")

## create descriptive tables for each phenotype that give sample size and then other stats depending on if it is a binary or quantitative trait

## create empty df to hold descriptive stats
within_pheno_descrip <- data.frame(matrix(ncol = 11, nrow = 0))

message("Recording within analysis descriptive statistics for...")

## for each phenotype remove the NAs to get the number of IIDs and FIDs, as well as other stats
for(k in 4:ncol(within_sample_phenos)){
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

message(paste("Done recording within analysis descriptive statistics and writing to ",Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_within.csv", sep=""))

## write to file
write.csv(within_pheno_descrip, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_within.csv", sep=""), row.names = FALSE)
