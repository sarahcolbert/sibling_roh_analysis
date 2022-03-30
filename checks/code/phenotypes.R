errorlist <- list()
warninglist <- list()

require(data.table)
require(plyr)

args <- (commandArgs(TRUE));
phenotype_file <- as.character(args[1]);
cov_file <- as.character(args[2]);

message("Checking phenotypes: ", phenotype_file)

if(phenotype_file == "NULL")
{
	msg <- paste0("No phenotype file has been provided.")
	warninglist <- c(warninglist, msg)
	message("WARNING: ", msg)
	q()
}


#Read in files

ph <- fread(phenotype_file, h=T)
p1 <- dim(ph)[1]
p2 <- dim(ph)[2]
cov <-fread(cov_file, h=T)

#Check phenotype file

if(names(ph)[1] !="IID")
	{
	msg <- paste0("First column in phenotype file should be the sample identified with the name IID. Please check that your phenotype file does not include an
  FID column")
	errorlist <-c(errorlist, msg)
	warning("ERROR: ", msg)
	}


nom <- names(ph)[-1][names(ph)[-1] %in% c("BMI", "Height")]


#Phenotype checks


if("Height" %in% nom)
{
	message("Checking Height")
	m1 <- mean(ph$Height,na.rm=T)
	age.mean<-mean(cov$age,na.rm=T)
	if((m1<100|m1>250)&age.mean>10)
	{
	msg <- paste0("please convert Height units to centimetres")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
	}
}


if("BMI" %in% nom)
{
	message("Checking BMI")
	m1<-mean(ph$BMI,na.rm=T)
	age.mean<-mean(cov$Age,na.rm=T)
	if((m1<10|m1>40)&age.mean>2)
	{
	msg <- paste0("please convert BMI units to kg/m2")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
	}
}


message("\n\nCompleted checks\n")


if(length(warninglist) > 0)
{
	message("\n\nPlease take note of the following warnings, and fix and re-run the phenotypic data check if you see fit:")
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
