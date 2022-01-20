###################################
############## PREP ###############
###################################

## load packages
library(tidyverse)

## load phenotype data
pheno1 <- read.table(paste(Sys.getenv("pheno_file")), header = TRUE)
## clean phenotype data
## remove duplicate individuals
pheno2 <- pheno1 %>% distinct(IID, .keep_all = TRUE)

## load covariate data
message(paste("loading ",(Sys.getenv("covar_file")), sep=""))
covar_data <- read.table(paste(Sys.getenv("covar_file")), header = TRUE)
message(paste("done loading ",(Sys.getenv("covar_file")), sep=""))

## load roh file data to merge with phenotype data
fam_file <- read.table(paste(Sys.getenv("processed_dir"),"all_froh_data.txt", sep=""), header = TRUE)
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

###################################
### CALC STATS FOR BTWN SAMPLE ####
###################################

## create descriptive tables for each phenotype that give sample size and then other stats depending on if it is a binary or quantitative trait

## create empty df to hold descriptive stats
btwn_pheno_descrip <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(btwn_pheno_descrip) <- c("pheno_name", "sample_size", "ncase", "ncontrols", "mean_phen", "median_phen", "sd_phen", "mean_age", "sd_age")

## for each phenotype remove the NAs to get the number of IIDs and FIDs, as well as other stats
num_phenos_btwn <- ncol(pheno_data_btwn) ## will be used to calculate number of phenotypes
for(k in 4:num_phenos_btwn){
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
write.csv(btwn_pheno_descrip, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats_btwn.csv", sep=""), row.names = FALSE)




################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################

##############################################################################################################################
################################################################################################
## still in development


################################
## old code
## load FID and IID data from the sibling roh file to merge with phenotype data
fam_file <- read.table(paste(Sys.getenv("processed_dir"),"all_froh_data.txt", sep=""), header = TRUE)
fam_data <- fam_file %>% drop_na(froh_sibs) %>% select(FID, IID)
## merge
phenotype_data_sibs <- merge(fam_data, pheno2, by ="IID")
all_data <- merge(phenotype_data_sibs, covar_data, by ="IID")

## calculate value of phenotype to family mean
## use for loop to go through each of the phenotype columns
num_phenos <- ncol(phenotype_data_sibs)
for(j in 3:num_phenos){
  ## make empty column to hold results
  phenotype_data_sibs$newcol <- NA
  newcolnum <- which(colnames(phenotype_data)=="newcol")
  ## use for loop to get value for each individual
for(i in 1:length(phenotype_data_sibs$IID)){
  spec_FID <- phenotype_data_sibs$FID[i]
  fam_vals <- phenotype_data_sibs %>% filter(FID==spec_FID)
  phenotype_data_sibs$newcol[i] <- phenotype_data_sibs[i,j]-mean(fam_vals[,j])
}
colnames(phenotype_data_sibs)[newcolnum] <- paste(colnames(phenotype_data_sibs)[j],"_sibs",sep="")
}
 
## write results to table so that easy to import for regression analysis
sib_pheno_results <- phenotype_data_sibs
sib_pheno_results[,3:num_phenos] <- NULL
write.table(sib_pheno_results, paste(Sys.getenv("processed_dir"),"within_sibs_pheno_data.txt", sep=""), quote = FALSE, row.names=FALSE)

#######################################
## create descriptive tables for each phenotype that give number of individuals, number of families, and then other stats depending on if it is a binary or quantitative trait

## create empty df to hold descriptive stats
all_pheno_descrip <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(all_pheno_descrip) <- c("pheno_name", "within_num_inds", "within_num_fams", "ncase", "ncontrols", "mean_phen", "median_phen", "sd_phen", "mean_age", "sd_age")

## add age to phenotype data
id_age <- all_data %>% select(IID, age)
phenotype_data_sibs2 <- merge(phenotype_data_sibs, id_age, by = "IID")

## for each phenotype remove the NAs to get the number of IIDs and FIDs, as well as other stats
for(k in 3:num_phenos){
test1 <- phenotype_data_sibs2 %>% select(FID, IID, age, colnames(phenotype_data_sibs2)[k], colnames(phenotype_data_sibs2)[(num_phenos-2)+k])
test2 <- test1 %>% drop_na()
within_num_inds <- length(unique(test2$IID))
within_num_fams <- length(unique(test2$FID))
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
pheno_descrip <- cbind(pheno_name, within_num_inds, within_num_fams, ncase, ncontrols, mean_phen, median_phen, sd_phen, mean_age, sd_age)
all_pheno_descrip <- rbind(all_pheno_descrip, pheno_descrip)
}

## write to file
write.csv(all_pheno_descrip, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_pheno_stats.csv", sep=""), row.names = FALSE)
