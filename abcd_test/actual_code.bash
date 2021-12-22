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

## see siblings_extra.bash

############################
#####  #####
############################
