errorlist <- list()
warninglist <- list()

require(data.table)

# import arguments
arguments <- commandArgs(trailingOnly = T)
relfile <- arguments[1]

## read in the relatedness output from PLINK
rel <- fread(relfile, sep=" ")


## check for monozyotic twins and give warning if family only has MZ twins
relMZ<-rel[which(rel$PI_HAT>0.98),]
MZ<-nrow(relMZ)
message("Number of monozygotic twin-pairs in sample: ", MZ)

if (MZ > 0) {
for (i in 1:nrow(relMZ)) {
  data<-rel[which(rel$FID1==relMZ$FID1[i] | rel$FID2==relMZ$FID1[i]),]
if(nrow(data)<2)
	{
	msg <-paste0("Identified family with only Monozygotic twins: Please check or remove the following families: ID ", data[1,1])
	errorlist<-c(errorlist, msg)

	warning("ERROR: ", msg)
	}
}
}

## check that first degree relative relationship is sibs and NOT parents (give warning if identified)
relPARENTS<-rel[which(rel$Z1>0.98),]

if(nrow(relPARENTS)>0)
	{
	for (i in 1:nrow(relPARENTS)) {
	msg <-paste0("Identified Parent-offspring pairs: Please check or remove the following families: IDs ", relPARENTS[i,1])
	errorlist<-c(errorlist, msg)
	warning("ERROR: ", msg)
	}
	}

## check for low IBD and give error if identified
relLOW<-rel[which(rel$PI_HAT<0.3),]
if(nrow(relLOW)>0)
	{
	for (i in 1:nrow(relPARENTS)) {
	msg <-paste0("Identified pair in same family with IBD estimate <0.3: Please check or remove the following families: IDs ", relLOW[i,1])
	errorlist<-c(errorlist, msg)
	warning("ERROR: ", msg)
	}
	}


## print errors if any are identified
message("\n\nCompleted checks\n")

if(length(errorlist) > 0)
{
	message("\n\nThe following errors were encountered, and must be addressed before continuing:")
	null <- sapply(errorlist, function(x)
	{
		message("- ", x)
	})
	q(status=1)
}
