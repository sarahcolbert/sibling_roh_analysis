## load packages and data
library(tidyverse)
roh_data1 <- read.table(paste(Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_roh.hom.indiv", sep=""), header = TRUE)
covar_data <- read.table(paste(Sys.getenv("covar_file")), header = TRUE)

## merge
roh_data <- merge(roh_data1, covar_data, by = "IID")

## calc froh
roh_data$froh <- roh_data$KB/(2.77*10^6)

## filter to necessary columns
froh_data2 <- roh_data %>% select(FID, IID, NSEG, KB, froh)

## remove any individuals that aren't in a sibling pair
froh_data <- froh_data2[!grepl('mis_abc', froh_data2$FID),]

## calculate value of froh relative to family mean
## make empty column to hold results
froh_data$froh_sibs <- NA
## use for loop to get value for each individual
for(i in 1:length(froh_data$IID)){
  spec_FID <- froh_data$FID[i]
  fam_vals <- froh_data %>% filter(FID==spec_FID)
  froh_data$froh_sibs[i] <- froh_data$froh[i]-mean(fam_vals$froh)
}

## save table with froh values just incase
write.table(froh_data, paste(Sys.getenv("processed_dir"),"within_sibs_froh_data.txt", sep=""), row.names=FALSE, quote = FALSE)


## calculate samples minimum, maximum, mean and median for NSEG, KB and froh
min_vals <- as.data.frame(apply(select(froh_data, c("NSEG","KB","froh")), 2, FUN = min, na.rm = TRUE))
colnames(min_vals)[1] <- "min"
mean_vals <- as.data.frame(apply(select(froh_data, c("NSEG","KB","froh")), 2, FUN = mean, na.rm = TRUE))
colnames(mean_vals)[1] <- "mean"
max_vals <- as.data.frame(apply(select(froh_data, c("NSEG","KB","froh")), 2, FUN = max, na.rm = TRUE))
colnames(max_vals)[1] <- "max" 
med_vals <- as.data.frame(apply(select(froh_data, c("NSEG","KB","froh")), 2, FUN = median, na.rm = TRUE))
colnames(med_vals)[1] <- "median" 

## put descriptive roh stats in data table and save as file (this file will be returned to us)
descript_roh_stats <- cbind(min_vals, max_vals, mean_vals, med_vals)
descript_roh_stats$n_indivs <- rep(length(froh_data$IID), length(descript_roh_stats$min))
descript_roh_stats$n_sibgroups <- rep(length(unique(froh_data$FID)), length(descript_roh_stats$min))

write.csv(descript_roh_stats, paste(Sys.getenv("output_dir"),Sys.getenv("input_prefix"),"_descriptive_roh_stats.csv", sep=""), row.names = TRUE)
