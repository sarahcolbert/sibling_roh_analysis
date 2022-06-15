## create objects to store errors and warnings
  errorlist <- list()
  warninglist <- list()

## load packages
library(data.table)
library(tidyverse)

## load bim file
args <- (commandArgs(TRUE));
bim_file <- as.character(args[1])
bim <- as.data.frame(fread(bim_file, h=F))
bim_file2 <- as.character(args[3])
bim2 <- as.data.frame(fread(bim_file2, h=F))

## load fam file
fam_file <- as.character(args[2])
fam <- read.table(fam_file,header=F,stringsAsFactors=F)
fam_file2 <- as.character(args[4])
fam2 <- read.table(fam_file2,header=F,stringsAsFactors=F)

## check fam files have matching FIDs
message("Checking that IID and FID match in imputed and non-imputed data")
# merge fam files
fam_nimp <- fam %>% select(V1,V2) %>% rename(FID1=V1)
fam_imp <- fam2 %>% select(V1,V2) %>% rename(FID2=V1)
fam3 <- merge(fam_nimp, fam_imp, by = "V2")
if(all(fam3$FID1 != fam3$FID2))
  {
  	msg <- paste0("Family IDs do not match in imputed and non-imputed fam files. Please fix this before going on.")
  	warninglist <- c(errorlist, msg)
  	warning("ERROR: ", msg)
  }

#### non-imputed data
## check SNPs
## check for duplicate SNPs
message("Checking for duplicate SNPs in non-imputed data")
if(any(duplicated(bim[,2])))
  {
  	msg <- "duplicate SNPs in bim file. Please remove."
  	errorlist <- c(errorlist, msg)
  	warning("ERROR: ", msg)
  }

## check # of SNPs
message("Checking non-imputed bim file: ", bim_file)
bim <- as.data.frame(fread(bim_file, h=F))
message("Number of SNPs: ", nrow(bim))
## give warning if bim file has < 250k SNPs which is a requirement
if(nrow(bim)<250000)
{
  msg <- "Less than 250k SNPs found in non-imputed bim file. Your data may not meet the study eligibility requirements."
  errorlist <- c(warninglist, msg)
  warning("Warning: ", msg)
}
## give warning if bim file has > 1.5 mil SNPs and may indicate that analyst has provided imputed data
if(nrow(bim)>1500000)
{
  msg <- "Over 1.5 million SNPs found in non-imputed bim file. Please confirm that you correctly specified non-imputed and imputed data in config file."
  warninglist <- c(warninglist, msg)
  warning("Warning: ", msg)
}

## check fam file for duplicate individuals
message("Checking fam file: ", fam_file)
if(any(duplicated(fam[,2])))
  {
  	msg <- paste0("Individual identifier is not unique. Please fix this before going on.")
  	warninglist <- c(errorlist, msg)
  	warning("ERROR: ", msg)
  }

##### imputed data
## check SNPs
## check for duplicate SNPs
message("Checking for duplicate SNPs in imputed data")
if(any(duplicated(bim2[,2])))
  {
  	msg <- "duplicate SNPs in bim file. Please remove."
  	errorlist <- c(errorlist, msg)
  	warning("ERROR: ", msg)
  }
## check # of SNPs
message("Checking non-imputed bim file: ", bim_file2)
  bim <- as.data.frame(fread(bim_file2, h=F))
  message("Number of SNPs: ", nrow(bim2))
  ## give warning if bim file has < 1.5 mil SNPs which might suggest mix up
  if(nrow(bim2)<1500000)
  {
    msg <- "Less than 1.5 million SNPs found in imputed bim file. Please confirm that you correctly specified non-imputed and imputed data in config file."
    errorlist <- c(warninglist, msg)
    warning("Warning: ", msg)
  }

## print warnings
if(length(warninglist) > 0)
  {
  	message("\n\nPlease take note of the following warnings, and fix and re-run the data check if you see fit:")
  	null <- sapply(warninglist, function(x)
  	{
  		message("- ", x)
  	})
  }

## pring errors
if(length(errorlist) > 0)
  {
  	message("\n\nThe following errors were encountered, and must be addressed before continuing:")
  	null <- sapply(errorlist, function(x)
  	{
  		message("- ", x)
  	})
  	q(status=1)
  }
  message("\n\n")
