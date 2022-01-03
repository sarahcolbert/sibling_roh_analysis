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
## make empty column to hold results
phenotype_data$col2 <- NA
## use for loop to get value for each individual
for(i in 1:length(phenotype_data$IID){
  spec_FID <- phenotype_data$FID[i]
  fam_vals <- phenotype_data %>% filter(FID==spec_FID)
  phenotype_data$col2[i] <- phenotype_data$school_2_y[i]-mean(fam_vals$school_2_y)
}
