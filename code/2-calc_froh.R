## load packages and data
library(tidyverse)
message(paste("loading ",Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_roh.hom.indiv", sep=""))
roh_data1 <- read.table(paste(Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_roh.hom.indiv", sep=""), header = TRUE)
message(paste("done loading ",Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_roh.hom.indiv", sep=""))
message(paste("loading ",(Sys.getenv("covar_file")), sep=""))
covar_data <- read.table(paste(Sys.getenv("covar_file")), header = TRUE)
message(paste("done loading ",(Sys.getenv("covar_file")), sep=""))

## merge
message("merging roh and covariate data")
roh_data <- merge(roh_data1, covar_data, by = "IID")
message("done merging roh and covariate data")

######################################
##### CALCULATE FROH WITHIN SIBS #####
######################################

## calc froh
message("calculating Froh")
roh_data$froh <- roh_data$KB/(2.77*10^6)

## filter to necessary columns
froh_data2 <- roh_data %>% select(FID, IID, NSEG, KB, froh)

## remove any individuals that aren't in a sibling pair
froh_data3 <- froh_data2[froh_data2$FID %in% froh_data2$FID[duplicated(froh_data2$FID)],]

## remove any duplicate individuals
froh_data <- froh_data3 %>% distinct(IID, .keep_all = TRUE)

## calculate value of froh relative to family mean
## make empty column to hold results
froh_data$froh_sibs <- NA
## use for loop to get value for each individual
for(i in 1:length(froh_data$IID)){
  spec_FID <- froh_data$FID[i]
  fam_vals <- froh_data %>% filter(FID==spec_FID)
  froh_data$froh_sibs[i] <- froh_data$froh[i]-mean(fam_vals$froh)
}

message("done calculating Froh")

## combine info for all individuals
id_sibs_froh <- froh_data %>% select(IID, froh_sibs)
all_froh_data <- merge(id_sibs_froh, roh_data, by = "IID", all = TRUE)

## save table with froh values just incase
write.table(all_froh_data, paste(Sys.getenv("processed_dir"),"all_froh_data.txt", sep=""), row.names=FALSE, quote = FALSE)
message(paste("wrote all Froh estimates to ",Sys.getenv("processed_dir"),"all_froh_data.txt", sep=""))

## save table with froh values just incase
write.table(froh_data, paste(Sys.getenv("processed_dir"),"within_sibs_froh_data.txt", sep=""), row.names=FALSE, quote = FALSE)
message(paste("wrote within sibling Froh estimates to ",Sys.getenv("processed_dir"),"within_sibs_froh_data.txt", sep=""))

################################
##### GET BASIC FROH STATS #####
################################

## calculate samples minimum, maximum, mean and median for NSEG and froh
message("calculating sample's minimum, maximum, mean and median for NSEG and Froh")
min_vals <- as.data.frame(apply(select(all_froh_data, c("NSEG","froh")), 2, FUN = min, na.rm = TRUE))
colnames(min_vals)[1] <- "min"
mean_vals <- as.data.frame(apply(select(all_froh_data, c("NSEG","froh")), 2, FUN = mean, na.rm = TRUE))
colnames(mean_vals)[1] <- "mean"
max_vals <- as.data.frame(apply(select(all_froh_data, c("NSEG","froh")), 2, FUN = max, na.rm = TRUE))
colnames(max_vals)[1] <- "max" 
sd_vals <- as.data.frame(apply(select(all_froh_data, c("NSEG","froh")), 2, FUN = sd, na.rm = TRUE))
colnames(sd_vals)[1] <- "median" 
message("done calculating sample's minimum, maximum, mean and standard deviation for NSEG and Froh")

## put descriptive roh stats in data table and save as file (this file will be returned to us)
message ("writing table with Froh stats")
descript_roh_stats <- cbind(min_vals, max_vals, mean_vals, sd_vals)
descript_roh_stats$n_indivs <- rep(length(all_froh_data$IID), length(descript_roh_stats$min))
descript_roh_stats$n_sibgroups <- rep(length(unique(froh_data$FID)), length(descript_roh_stats$min))
## check how many individuals have froh = 0
zero_inds <- sum(all_froh_data$froh==0)
descript_roh_stats$n_zero_inds <- rep(zero_inds, length(descript_roh_stats$min))
## check how many individuals are offspring of various relative relationships
half_sib_inds <- sum(all_froh_data$froh > 0.125)
descript_roh_stats$n_half_sib_inds <- rep(half_sib_inds, length(descript_roh_stats$min))
first_cous_inds <- sum(all_froh_data$froh > 0.0625)
descript_roh_stats$n_first_cous_inds <- rep(first_cous_inds, length(descript_roh_stats$min))
half_cous_inds <- sum(all_froh_data$froh > 0.03125)
descript_roh_stats$n_half_cous_inds <- rep(half_cous_inds, length(descript_roh_stats$min))
sec_cous_inds <- sum(all_froh_data$froh > 0.0156)
descript_roh_stats$n_sec_cous_inds <- rep(sec_cous_inds, length(descript_roh_stats$min))

write.csv(descript_roh_stats, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_roh_stats.csv", sep=""), row.names = TRUE)
message(paste("wrote table with Froh stats to ",Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_roh_stats.csv", sep=""))
