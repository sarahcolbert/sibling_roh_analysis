## create objects to store errors and warnings
  errorlist <- list()
  warninglist <- list()

## load packages
library(data.table)

## load bim file
args <- (commandArgs(TRUE));
bim_file <- as.character(args[1])
bim <- as.data.frame(fread(bim_file, h=F))

## load fam file
fam_file <- as.character(args[2])
fam <- read.table(fam_file,header=F,stringsAsFactors=F)

## check SNPs
## check for duplicate SNPs
message("Checking for duplicate SNPs")
if(any(duplicated(bim[,2])))
  {
  	msg <- "duplicate SNPs in bim file. Please remove."
  	errorlist <- c(errorlist, msg)
  	warning("ERROR: ", msg)
  }

## check # of SNPs
message("Checking bim file: ", bim_file)
bim <- as.data.frame(fread(bim_file, h=F))
message("Number of SNPs: ", nrow(bim))
## give warning if bim file has < 250k SNPs which is a requirement
if(nrow(bim)<250000)
{
  msg <- "Less than 250k SNPs found in bim file. Your data may not meet the study eligibility requirements."
  errorlist <- c(warninglist, msg)
  warning("Warning: ", msg)
}
## give warning if bim file has > 1.5 mil SNPs and may indicate that analyst has provided imputed data
if(nrow(bim)>1500000)
{
  msg <- "Over 1.5 million SNPs found in bim file. Please confirm that you are not using imputed data."
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
