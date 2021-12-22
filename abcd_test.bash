#########################
##### CONFIGURATION #####
#########################

## navigate to your projects directory

## making project directory for analysis
export WORKING_DIR="$PWD"
export project_dir="${WORKING_DIR}/sibling_roh_analysis/" ## will make this $1
mkdir ${project_dir} ${project_dir}output ${project_dir}log
echo "created ${project_dir}"

## get locations for data
export data_dir="/scratch/aalab/suri/data/raw_data/" ## will make this $2

############################
##### CONFIRM SIBLINGS #####
##### NOT PART OF OUR CODE #####
############ (I'm just doing this so we have ABCD data ready to test, but won't provide in our code since we link to the sibling GWAS github how to)
############################


export king="/home/sarah.c/software/king"
export genotyped_data="${data_dir}abc_abc1_white_no_hisp_sc-qc2.bed"
export out="${project_dir}output"

sbatch 1-king.sbatch

### separate bash script called 1-king.sbatch
$king -b $genotyped_data \
	--related \
        --degree 2 \
	--prefix $out/sibs \
	--rplot \
       	|& tee $out/sibs.txt
	
#Format sibling file from long to wide
##################### R CODE

# Load dependencies
require(data.table)

#Read in file and restrict to siblings
kin <- read.delim('./sibling_roh_analysis/output/sibs.kin0', header=T, as.is=T)

# Restrict to full siblings
kin2 <- kin[kin$InfType=='FS',]

#Converting from pairwise format to individual level format
#sibs<-input2[,c(1:3)]
sibs <- kin2[,c("FID1", "ID1", "ID2")]
names(sibs) <- c("FID", "IID", "Sibling")
sibs$FID <- seq.int(nrow(sibs))

#sibs2<-input2[,c(1:3)]
sibs2 <- kin2[,c("FID1", "ID1", "ID2")]
names(sibs2) <- c("FID", "Sibling", "IID")
sibs2$FID <- seq.int(nrow(sibs2))

sibs3<-rbind(sibs,sibs2)
sibs4<-sibs3[order(sibs3$FID),]

#Generate duplicated and non-duplicated sets of siblings
nodup<-sibs4[!duplicated(sibs4$IID) & !duplicated(sibs4$Sibling),]
dup<-sibs4[duplicated(sibs4$IID) | duplicated(sibs4$Sibling),]
dup$FID<-NULL
names(dup)<-c("IID", "Sibling2")

#Merge to get same family IDs for siblings
merge<-merge(nodup,dup,by="IID")
merge2<-merge[order(merge$FID),]
merge3<-merge2[,c(2,3)]
merge4<-merge2[,c(2,4)]
names(merge4)<-c("FID", "Sibling")
merge5<-rbind(merge3,merge4)
names(merge5)<-c("FID", "IID")
sibtrios<-merge5[!duplicated(merge5),]

#Merge everything together
nodup$Sibling<-NULL
final<-rbind(nodup,sibtrios)
final2<-final[order(final$FID),]
final3<-final2[!duplicated(final2),]

#Output
output1<-final3[,c(1,2)]
names(output1)<-c("FID", "IID")
output2<-final3[,c(2,2)]
names(output2)<-c("FID", "IID")
write.table(output1, "./sibling_roh_analysis/output/Siblings-FID.fam", quote=F, row.names=F, col.names=F, sep=" ")	

# R script to update FID, keeping the order of individuals the same as the original PLINK file.
library(dplyr)

# Input
sibFID <- read.delim2("./sibling_roh_analysis/output/Siblings-FID.fam", header=F, sep=" ", as.is=T)
plinkFAM <- read.delim2("/scratch/aalab/suri/data/raw_data/abc_abc1_white_no_hisp_sc-qc2.fam", header=F, sep=" ", as.is=T)

# There are some duplicate IID with different FIDS in "./sibling_roh_analysis/output/Siblings-FID.fam"

# Delete the duplicate FIDs
sibFID2 <- sibFID[!duplicated(sibFID$V2),]

# Set order - must keep the same order as in the PLINK file
plinkFAM$id <- 1:nrow(plinkFAM)

# Merge
merged <- merge(x = sibFID2, y = plinkFAM, by = "V2", sort = FALSE, all=T)

# Order
ordered <- merged[order(merged$id), ]

# Check dim
dim(plinkFAM)
dim(ordered)

# rename columns
colnames(ordered)[1:3] <- c("IID", "FID", "kingID")

# Add IID as FID for missing
ordered$FID[is.na(ordered$FID)] <- ordered$kingID[is.na(ordered$FID)]

# Format the file
ordered$kingID <- NULL
ordered$id <- NULL

# reorder
ordered2 <- ordered %>% select(FID,IID,V3,V4,V5,V6)

# Output
write.table(ordered2, "./sibling_roh_analysis/output/update.fam", quote=F, row.names=F, col.names=F, sep=" ")


############################
#####  #####
############################
