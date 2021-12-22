############################
##### CONFIRM SIBLINGS #####
##### NOT PART OF OUR CODE #####
############ (I'm just doing this so we have ABCD data ready to test, but won't provide in our code since we link to the sibling GWAS github how to)
############################

## see siblings_extra.bash


#########################
##### CONFIGURATION #####
#########################

## navigate to your projects directory

## making directories for analysis
export WORKING_DIR="$PWD"
export project_dir="${WORKING_DIR}/sibling_roh_analysis/"
export processed_dir="${project_dir}processed_files/"
export output_dir="${project_dir}output/"
mkdir ${project_dir} ${processed_dir} ${output_dir}
echo "created ${project_dir}"

## get locations for data
export data_dir="/scratch/aalab/suri/data/cleaned_data/" ## will make this $1
echo "reading data from ${data_dir}"
## get name of plink files
export input_prefix="abcd_eur_sibs" ## will make this $2

############################
##### QC AND CALL ROHS #####
############################

## use plink to exclude SNPs with >3% missingess or MAF < 5% AND exclude individuals with >3% missing data
plink --bfile ${data_dir}${input_prefix} --maf 0.05 --geno 0.03 --mind 0.03 --make-bed --out ${processed_dir}${input_prefix}_filtered

## use plink to identify continuous ROH SNPs
plink --bfile ${processed_dir}/${input_prefix}_filtered --homozyg-window-snp 50 --homozyg-snp 50  --homozyg-kb 1500  --homozyg-gap 1000  --homozyg-density 50 --homozyg-window-missing 5 --homozyg-window-het 1 --out ${processed_dir}/${input_prefix}_roh

############################
## use R to calculate froh

## load packages and data
library(tidyverse)
roh_data <- read.table(paste(Sys.getenv("processed_dir"),Sys.getenv("input_prefix"),"_roh.hom.indiv", sep=""), header = TRUE)

## calc froh
roh_data$froh <- roh_data$KB/(2.77*10^6)

## filter to necessary columns
froh_data <- roh_data %>% select(IID, NSEG, KB, froh)

## calculate minimum, maximum and mean for NSEG, KB and froh
min_vals <- as.data.frame(apply(froh_data[,2:4], 2, FUN = min, na.rm = TRUE))
colnames(min_vals)[1] <- "min"
mean_vals <- as.data.frame(apply(froh_data[,2:4], 2, FUN = mean, na.rm = TRUE))
colnames(mean_vals)[1] <- "mean"
max_vals <- as.data.frame(apply(froh_data[,2:4], 2, FUN = max, na.rm = TRUE))
colnames(max_vals)[1] <- "max" 
##########
########## ADD MEDIAN
##########

