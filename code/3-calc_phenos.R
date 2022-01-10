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

## load FID and IID data from the sibling roh file to merge with phenotype data
fam_file <- read.table(paste(Sys.getenv("processed_dir"),"within_sibs_froh_data.txt", sep=""), header = TRUE)
fam_data <- fam_file %>% select(FID, IID)
## merge
phenotype_data <- merge(fam_data, pheno2, by ="IID")
all_data <- merge(phenotype_data, covar_data, by ="IID")

## calculate value of phenotype to family mean
## use for loop to go through each of the phenotype columns
num_phenos <- ncol(phenotype_data)
for(j in 3:num_phenos){
  ## make empty column to hold results
  phenotype_data$newcol <- NA
  newcolnum <- which(colnames(phenotype_data)=="newcol")
  ## use for loop to get value for each individual
for(i in 1:length(phenotype_data$IID)){
  spec_FID <- phenotype_data$FID[i]
  fam_vals <- phenotype_data %>% filter(FID==spec_FID)
  phenotype_data$newcol[i] <- phenotype_data[i,j]-mean(fam_vals[,j])
}
colnames(phenotype_data)[newcolnum] <- paste(colnames(phenotype_data)[j],"_sibs",sep="")
}
 
## write results to table so that easy to import for regression analysis
sib_pheno_results <- phenotype_data
sib_pheno_results[,3:num_phenos] <- NULL
write.table(sib_pheno_results, paste(Sys.getenv("processed_dir"),"within_sibs_pheno_data.txt", sep=""), quote = FALSE, row.names=FALSE)

#######################################
## create descriptive tables for each phenotype that give number of individuals, number of families, and then other stats depending on if it is a binary or quantitative trait


## create empty df to hold descriptive stats
all_pheno_descrip <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(all_pheno_descrip) <- c("pheno_name", "num_inds", "num_fams", "ncase", "ncontrols", "mean_phen", "median_phen", "sd_phen", "mean_age", "sd_age")

## for each phenotype remove the NAs to get the number of IIDs and FIDs, as well as other stats
for(k in 3:num_phenos){
test1 <- phenotype_data %>% select(FID, IID,colnames(phenotype_data)[k], colnames(phenotype_data)[(num_phenos-2)+k])
test2 <- test1 %>% drop_na()
num_inds <- length(unique(test2$IID))
num_fams <- length(unique(test2$FID))
mean_age <- mean(test2$age)
sd_age <- sd(test2$age)  
  
## assesses if trait is binary or quantitative and pulls the correct info depending
if(all(test2[,3] %in% c(0,1))){
  ncase <- length(which(test2[,3]==1))
  ncontrols <- length(which(test2[,3]==0)) 
  mean_phen <- NA
  median_phen <- NA
  sd_phen <- NA
}else{
  ncase <- NA
  ncontrols <- NA
  mean_phen <- mean(test2[,3])
  median_phen <- median(test2[,3])
  sd_phen <- sd(test2[,3])}

pheno_name <- paste(colnames(test2)[3])
pheno_descrip <- cbind(pheno_name, num_inds, num_fams, ncase, ncontrols, mean_phen, median_phen, sd_phen, mean_age, sd_age)
all_pheno_descrip <- rbind(all_pheno_descrip, pheno_descrip)
}

## write to file
write.csv(all_pheno_descrip, paste(Sys.getenv("output_dir"),Sys.getenv("input_prefix"),"_descriptive_pheno_stats.csv", sep=""), row.names = FALSE)
