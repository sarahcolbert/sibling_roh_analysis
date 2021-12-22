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
############################

# (I'm just doing this so we have ABCD data ready, but won't provide in our code since we link to the sibling GWAS github how to)

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
