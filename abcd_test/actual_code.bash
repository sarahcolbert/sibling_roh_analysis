############################
##### CONFIRM SIBLINGS #####
##### NOT PART OF OUR CODE #####
############ (I'm just doing this so we have ABCD data ready to test, but won't provide in our code since we link to the sibling GWAS github how to)
############################

## see 0-siblings_extra.bash


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

## give location for directory with plink files
export data_dir="/scratch/aalab/suri/data/cleaned_data/" ## will make this $1
echo "reading data from ${data_dir}"
## give name of plink files
export input_prefix="abcd_eur_sibs" ## will make this $2
## give location for covariate file
export covar_file="/scratch/aalab/suri/data/raw_data/abcd_covariates_eur.txt" ## will make this $3
## give location for phenotype file 
export pheno_file="/scratch/aalab/suri/data/raw_data/abcd_phenotypes_eur.txt" ## will make this $4


############################
##### QC AND CALL ROHS #####
############################

## use plink to perform qc and roh calling
bash 1-qc_and_call.bash

######################################
##### CALCULATE FROH WITHIN SIBS #####
##### GET DESCRIPTIVE STATISTICS FOR SAMPLE OF SIBLINGS
######################################

## use R to calculate froh
Rscript 2-calc_froh.R

############################################
##### CALCULATE PHENOTYPES WITHIN SIBS #####
##### GET DESCRIPTIVE STATISTICS FOR PHENOTYPES
############################################

## use R to calculate froh
Rscript 3-calc_phenos.R
