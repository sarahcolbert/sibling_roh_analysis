## load packages
library(tidyverse)

## load data
phenotype_data <- read.table("${phenotype_file.txt}", header = TRUE)

## calculate value of phenotype (height as example) to family mean
## make empty column to hold results
phenotype_data$height_sibs <- NA
## use for loop to get value for each individual
for(i in 1:length(phenotype_data$IID){
  spec_FID <- phenotype_data$FID[i]
  fam_vals <- phenotype_data %>% filter(FID=spec_FID)
  phenotype_data$height_sibs[i] <- phenotype_data$height[i]-mean(fam_vals$height)
}

write.table(phenotype_data, "within_sibs_phenotype_data.txt", row.names=FALSE, quote = FALSE)
