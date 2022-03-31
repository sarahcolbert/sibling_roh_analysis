set -e
source ./config

mkdir -p ${project_dir}checks/results ${project_dir}checks/results/logs

section_01_logfile="${project_dir}checks/results/logs/log.txt"

exec &> >(tee ${section_01_logfile})

containsElement () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo "There is no method for ${1}."
	echo "Please run:"
	echo "${code_dir}0-checks.bash [arg]"
	echo "where arg is an optional argument that can be one of:"
	printf '%s\n' ${@:2}
	return 1
}

arg="all"
declare -a sections=('all' 'config' 'requirements' 'genetics' 'rel' 'siblings' 'covariates' 'phenotypes' 'skipsib')


if [ -n "${1}" ]; then
	arg="${1}"
	containsElement ${1} ${sections[@]}
fi

section_message () {
	echo "-----------------------------------------------"
	echo ""
	echo "$1 section"
	echo ""
	echo "to run this part on its own type:"
	echo "$ ${code_dir}0-checks.bash $1"
	echo ""
	echo "-----------------------------------------------"
	echo ""
	echo ""
}

#Check study name

if [ "$arg" = "config" ] || [ "$arg" = "all" ] || [ "$arg" = "skipsib" ]
then
	section_message "config"
	if ! [[ "$study_name" =~ [^a-zA-Z0-9_\ ] ]] && ! [ "$study_name" = "" ]
	then
		echo "The study name '${study_name}' will be used for this analysis. Change this in the config file if necessary."
		echo ""
	else
		echo "The study name '${study_name}' is invalid. Please use only alphanumeric or underscore characters, with no spaces or special characters etc."
		exit
	fi
fi

#Check PLINK version, R version and packages

if [ "$arg" = "requirements" ] || [ "$arg" = "all" ] || [ "$arg" = "skipsib" ]
then
	section_message "requirements"
    Rscript ${project_dir}checks/code/packages.R
    echo "Checking Plink version"
    plinkv=`${plink} | awk 'NR==1{print $1,$2}'`
    if echo ${plinkv} | grep -q '1.9' ; then
    echo "You are running ${plinkv}"
else
    echo "You are running ${plinkv}. Please switch to PLINK version 1.9 before continuing."
fi
fi

## check genotype file

if [ "$arg" = "genetics" ] || [ "$arg" = "all" ] || [ "$arg" = "skipsib" ]
then
	section_message "genetics"
	Rscript ${project_dir}checks/code/genetic_data.R \
		${data_dir}${input_prefix}.bim \
		${data_dir}${input_prefix}.fam
fi

## generate relatedness checks

if [ "$arg" = "rel" ] || [ "$arg" = "all" ]
then
	section_message "rel"
	${plink} \
	--bfile ${data_dir}${input_prefix} \
	--genome rel-check \
	--out ${project_dir}checks/results/relatedness
fi

## perform sibling checks

if [ "$arg" = "siblings" ] || [ "$arg" = "all" ]
then
	section_message "siblings"
	Rscript ${project_dir}checks/code/siblings.R \
		${project_dir}checks/results/relatedness.genome
fi

## check covariate file

if [ "$arg" = "covariates" ] || [ "$arg" = "all" ] || [ "$arg" = "skipsib" ]
then
	section_message "covariates"
	Rscript ${project_dir}checks/code/covariates.R \
		${covar_file} \
		${data_dir}${input_prefix}.fam
fi

## check phenotype file

if [ "$arg" = "phenotypes" ] || [ "$arg" = "all" ] || [ "$arg" = "skipsib" ]
then
	section_message "phenotypes"
	Rscript ${project_dir}checks/code/phenotypes.R \
		${pheno_file} \
		${covar_file}
fi

## finish checks

if [ "$arg" = "all" ]
then
        echo ""
        echo ""
        echo "You successfully performed all data checks!"
fi
