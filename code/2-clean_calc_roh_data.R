##------------------------------------------------
## Set-up
##------------------------------------------------

## load packages
library(tidyverse)

## load roh data
message(paste("Loading ",Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_roh.hom.indiv", sep=""))
roh_data1 <- read.table(paste(Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_roh.hom.indiv", sep=""), header = TRUE)
message(paste("Done loading ",Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_roh.hom.indiv", sep=""))

## load Fhat3 data
message(paste("Loading ",Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_ibc.ibc", sep=""))
fhat3_data1 <- read.table(paste(Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_ibc.ibc", sep=""), header = TRUE)
message(paste("Done loading ",Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_ibc.ibc", sep=""))

## load covariates
message(paste("Loading ",(Sys.getenv("covar_file")), sep=""))
covar_data1 <- read.table(paste(Sys.getenv("covar_file")), header = TRUE) %>% distinct()
covar_data <- covar_data1[!duplicated(covar_data1$IID), ]
message(paste("Done loading ",(Sys.getenv("covar_file")), sep=""))

##------------------------------------------------
## Make dataset with all data for all individuals
##------------------------------------------------

## add covariate data to roh and Fhat3 dfs
message("Merging roh, fhat3 and covariate data")
roh_fhat3 <- merge(roh_data1, fhat3_data1, by = c("IID","FID"))
roh_df1 <- merge(roh_fhat3, covar_data, by = "IID")
message("Done merging data")

## since individuals need all the data in this dataframe to be included, remove anyone who is missing these roh estimates or covariates
roh_df <- roh_df1 %>% drop_na()

##------------------------------------------------
## Calculations
##------------------------------------------------

## calculate froh
message("Calculating Froh")
roh_df$froh <- roh_df$KB/(2.77*10^6)
message("Done calculating Froh")

##------------------------------------------------
## Determine if individual is in sibling analysis
##------------------------------------------------

## make list of people who have sibling in dataset with genetic data
## remove any individuals that aren't in a sibling pair
## and remove duplicates
sib_list <- roh_df[roh_df$FID %in% roh_df$FID[duplicated(roh_df$FID)],] %>% distinct(IID, .keep_all = TRUE) %>% select(IID, FID)

## add column that marks if individual is included in sib analysis
roh_df$wsib <- ifelse(roh_df$IID %in% sib_list$IID, 1, 0)

## save dataset
write.table(roh_df, paste(Sys.getenv("processed_dir"),"all_indivs_roh_covs.txt", sep=""), row.names=FALSE, quote = FALSE)

##------------------------------------------------
## Get ROH descriptive statistics
##------------------------------------------------

## calculate samples minimum, maximum, mean and standard deviation for NSEG and froh
message("Calculating sample's minimum, maximum, mean and standard deviation for NSEG, Froh and Fhat3")
min_vals <- as.data.frame(apply(select(roh_df, c("NSEG","froh","Fhat3")), 2, FUN = min, na.rm = TRUE))
colnames(min_vals)[1] <- "min"
mean_vals <- as.data.frame(apply(select(roh_df, c("NSEG","froh","Fhat3")), 2, FUN = mean, na.rm = TRUE))
colnames(mean_vals)[1] <- "mean"
max_vals <- as.data.frame(apply(select(roh_df, c("NSEG","froh","Fhat3")), 2, FUN = max, na.rm = TRUE))
colnames(max_vals)[1] <- "max"
sd_vals <- as.data.frame(apply(select(roh_df, c("NSEG","froh","Fhat3")), 2, FUN = sd, na.rm = TRUE))
colnames(sd_vals)[1] <- "sd"
message("Done calculating sample's minimum, maximum, mean and standard deviation for NSEG, Froh and Fhat3")

## put descriptive roh stats in data table and save as file (this file will be returned to us)
message ("Writing table with Froh stats")
descript_roh_stats <- cbind(min_vals, max_vals, mean_vals, sd_vals)
descript_roh_stats$n_indivs <- rep(length(roh_df$IID), length(descript_roh_stats$min))
## number of sib groups
descript_roh_stats$n_sibgroups <- rep(length(unique(sib_list$FID)), length(descript_roh_stats$min))
## check how many individuals have froh = 0
zero_inds <- sum(roh_df$froh==0)
descript_roh_stats$n_zero_inds <- rep(zero_inds, length(descript_roh_stats$min))
## check how many individuals are offspring of various relative relationships
half_sib_inds <- sum(roh_df$froh > 0.125)
descript_roh_stats$n_half_sib_inds <- rep(half_sib_inds, length(descript_roh_stats$min))
first_cous_inds <- sum(roh_df$froh > 0.0625)
descript_roh_stats$n_first_cous_inds <- rep(first_cous_inds, length(descript_roh_stats$min))
half_cous_inds <- sum(roh_df$froh > 0.03125)
descript_roh_stats$n_half_cous_inds <- rep(half_cous_inds, length(descript_roh_stats$min))
sec_cous_inds <- sum(roh_df$froh > 0.0156)
descript_roh_stats$n_sec_cous_inds <- rep(sec_cous_inds, length(descript_roh_stats$min))

write.csv(descript_roh_stats, paste(Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_roh_stats.csv", sep=""), row.names = TRUE)
message(paste("Wrote table with Froh stats to ",Sys.getenv("output_dir"),Sys.getenv("output_name"),"_descriptive_roh_stats.csv", sep=""))
