errorlist <- list()
warninglist <- list()

library(data.table)
suppressMessages(library(matrixStats))

args <- (commandArgs(TRUE));
cov_file <- as.character(args[1]);
fam_file <- as.character(args[2])

message("Checking covariates: ", cov_file)

if(cov_file == "NULL")
{
	msg <- paste0("No covariate file has been provided.")
	warninglist <- c(warninglist, msg)
	message("WARNING: ", msg)
	q()
}

cov<-fread(cov_file, h=T)
cov1<-dim(cov)[1]
cov2<-dim(cov)[2]
cov3<-dim(cov)[3]
cov4<-dim(cov)[4]
cov5<-dim(cov)[5]
cov6<-dim(cov)[6]


fam<-read.table(fam_file, h=F, stringsAsFactors=F)

## check number of covariates
if(ncol(cov)>13)
	{
	msg <- paste0("There appears to be additional columns/covariates in your file. Please only include IID, age, sex, and PCs 1-10.")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(ncol(cov)<13)
  {
	msg <- paste0("There appears to be missing columns/covariates in your file. Please make sure to include IID, age, sex, and PCs 1-10.")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}



## check IDs, column orders and name
if(names(cov)[1] !="IID")
	{
	msg <- paste0("First column in covariate file should be the sample identifier with the name IID")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(names(cov)[2] !="age")
	{
	msg <- paste0("Second column in covariate file should be age.")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}



## check age
  age <- names(cov)[-1][names(cov)[-1] %in% c("age")]
  if(length(age)<1)
  	{
  	msg <-paste0("Age is not present in the covariate file. Please check that the columns are labelled correctly")
  	errorlist <-c(errorlist,msg)
  	warning("ERROR: ", msg)
  	}

    if(any(is.na(cov$age)))
    {
    msg<-paste0("Missing values for age. Please make sure all individuals have data for this covariate.")
    errorlist <-c(errorlist, msg)
    warning("ERROR: ", msg)
    }

    if(any(cov$age<0))
    {
     msg<-paste0("Negative values in the age column.")
    errorlist <-c(errorlist, msg)
    warning("ERROR: ", msg)
    }

    if(mean(cov$age, na.rm=T)>100)
    {
     msg<-paste0("Average age is above 100, please make sure age is provided in years.")
    errorlist <-c(errorlist, msg)
    warning("ERROR: ", msg)
    }



## check sex
if(names(cov)[3] !="sex")
  {
  msg <- paste0("Third column in covariate file should be sex")
  errorlist <-c(errorlist, msg)
  warning("ERROR: ", msg)
  }

  sex <- names(cov)[-1][names(cov)[-1] %in% c("sex")]
  if(length(sex)<1)
  	{
  	msg <-paste0("Sex is not present in the covariate file. Please check that the columns are labelled correctly")
  	errorlist <-c(errorlist,msg)
  	warning("ERROR: ", msg)
  	}

  if(any(is.na(cov$sex)))
    {
    msg<-paste0("Missing values for sex. Please make sure all individuals have data for this covariate.")
    errorlist <-c(errorlist, msg)
    warning("ERROR: ", msg)
    }

  index<-cov$sex %in% c("1", "0")
  if(any(!index))
    {
    msg<-paste0("There are some values in the sex column that are neither 0 (F) nor 1 (M). Please categorise Males as 1 and Females as 0")
    errorlist<-c(errorlist, msg)
    warning("ERROR: ", msg)
    }



## check principal components
if(names(cov)[4] !="PC1")
	{
	msg <- paste0("Fourth column in covariate file should be PC1")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(names(cov)[5] !="PC2")
	{
	msg <- paste0("Fifth column in covariate file should be PC2")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(names(cov)[6] !="PC3")
	{
	msg <- paste0("Sixth column in covariate file should be PC3")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(names(cov)[7] !="PC4")
	{
	msg <- paste0("Seventh column in covariate file should be PC4")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(names(cov)[8] !="PC5")
	{
	msg <- paste0("Eigth column in covariate file should be PC5")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(names(cov)[9] !="PC6")
	{
	msg <- paste0("Ninth column in covariate file should be PC6")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(names(cov)[10] !="PC7")
	{
	msg <- paste0("Tenth column in covariate file should be PC7")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(names(cov)[11] !="PC8")
	{
	msg <- paste0("Eleventh column in covariate file should be PC8")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(names(cov)[12] !="PC9")
	{
	msg <- paste0("Twelfth column in covariate file should be PC9")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

if(names(cov)[13] !="PC10")
	{
	msg <- paste0("Thirteenth column in covariate file should be PC10")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}

pccheck <- grep("PC", names(cov))
  if(length(pccheck)<10)
    {
    msg <-paste0("The first 10 principal components are not present in the covariate file. Please check that the columns are labelled correctly and you have included PCs 1-10.")
    errorlist <-c(errorlist,msg)
    warning("ERROR: ", msg)
    }

  if(length(pccheck)>10)
    {
    msg <-paste0("You have included more than the first 10 principal components in the covariate file. Please check that the columns are labelled correctly and only inlcude PCs 1-10. Please contact us if you would like to include more.")
    errorlist <-c(errorlist,msg)
    warning("ERROR: ", msg)
  	}




## check that at least 1,000 individuals have covariate and genetic data
## (not going to check pheno data right now because it will differ cohort to cohort which phenotypes they have)
commonids_cpg <- Reduce(intersect, list(cov$IID, fam[,2]))
message("Number of samples with covariate and genetic data: ", length(commonids_cpg))

if(length(commonids_cpg) < 1000){
	msg <- paste0("must have at least 1000 individuals with covariate and genetic data.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}


## done with checks, print warnings and errors
message("\n\nCompleted checks\n")

if(length(warninglist) > 0)
{
	message("\n\nPlease take note of the following warnings, and fix and re-run the data check if you see fit:")
	null <- sapply(warninglist, function(x)
	{
		message("- ", x)
	})
}

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
