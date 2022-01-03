## load packages
library(tidyverse)

## load phenotype data
pheno1 <- read.table(paste(Sys.getenv("pheno_file")), header = TRUE)
## clean phenotype data
## remove duplicates
pheno2 <- pheno1 %>% distinct(IID, .keep_all = TRUE)


## load FID and IID data from the sibglig roh file to merge with phenotype data
fam_file <- read.table(paste(Sys.getenv("processed_dir"),"within_sibs_froh_data.txt", sep=""), header = TRUE)
fam_data <- fam_file %>% select(FID, IID)
## merge
phenotype_data <- merge(fam_data, pheno2, by ="IID")

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

