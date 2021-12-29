# NOTE TO SELF/TO DO:
* pre-step 4: phenotype data - ideas: make data frame with empty columns for all phenotypes we want and then will have individuals move their data into it?
* steps 4 and 5

# Runs of homozygosity in the Within Family Consortium Analyses
This repository details a Standard Operating Procedure for the data analysts that will be performing the WFC runs of homozygosity analyses in individual cohorts. We provide a fully automated analysis pipeline that will analyze and format the results to be returned. We thank blabla for allowing us to adapt code and methods from their [sibling GWAS github](https://github.com/LaurenceHowe/SiblingGWAS).


## Pre-check: Pipeline Prerequisites and Requirements
The data requirements for the pipeline are as follows:

1) Sibling data (see Pre-Step 2 for how to define siblings).

2) Called genotypes in plink binary format. Imputed data is not permitted. (see Pre-Step 4 for information on input file requirements and scripts for file conversion).

3) Complete covariate data (see Pre-Step 4 for information on what covariates should be included and the file format).

4) Phenotype data for various traits (see Pre-Step 4 for more information on phenotypes requested, coding phenotypes, and phenotype file format).

The software requirements for the pipeline are as follows:

* Plink 1.9
* R (tidyverse installed, )
* KING (if siblings are not already defined)

## Pre-Step 1: Downloading and running the pipeline

Navigate to the directory where you want to download the repository. The repository can then be downloaded using git: <br>
> git clone https://github.com/sarahcolbert/sibling_roh_analysis <br>


## Pre-Step 2: Defining siblings
The analysis pipeline requires data on siblings. We follow the suggestion of Howe et al., which is to include "all siblings from families with one or more pairs of genotyped dizygotic siblings. For example, in a family with a pair of monozygotic twins and an additional sibling, include both MZ twins and the sibling. The inclusion of both MZ twins should (very) modestly improve power by accounting for variation in the phenotypic outcome. If siblings have not been previously identified in the dataset, we suggest using [KING](https://www.kingrelatedness.com/) to infer siblings."

Instructions for defining siblings are provided by Howe et al. [here](https://github.com/LaurenceHowe/SiblingGWAS/wiki/0.1_Siblings).

## Pre-Step 3: Editing the config file
To edit the config file navigate to the sibling_roh_analysis directory and add the required information to the file. You can then run the config file using:  <br>
> source ./config <br>

## Pre-Step 4: Input Files
### Genotype data
You will need genotype data in PLINK binary format. The pipeline requires the input (.bed, .bim, .fam) files to satisfy the following requirements:

a) PLINK binary format (.bed .bim .fam) files. The first two columns must contain family IDs (FID) and individual IDs (IIDs). FIDs should be common between siblings (but unique between sets of siblings) and IIDs should be unique for each individual.

b) add requirements as needed.

### Covariate data
A covariate file should be provided containing the following columns, if available:
FID, IID, age, sex (male = 1, female = 0), batch, first 10 genomic principal components, anything else? 

### Phenotype data 
how closely do we want to follow c and d here? https://github.com/LaurenceHowe/SiblingGWAS/wiki/0.4_Phenotypes

## Steps 1 + 2: QC and ROH Calling
Genotype files must first be filtered to meet the following requirements:
1) exclude SNPs with >3% missingess
2) exclude SNPs with MAF < 5% 
3) exclude individuals with >3% missing data

Continuous ROH SNPs are identified using PLINK with the following parameters:
1) homozyg-window-snp 50
2) homozyg-snp 50
3) homozyg-kb 1500
4) homozyg-gap 1000
5) homozyg-density 50
6) homozyg-window-missing 5
7) homozyg-window-het 1

Running the code below will perform QC and then use the QCed files to call ROHs.
Before running you will need to make sure that you have plink 1.9 installed and that it can be called with the command "plink" as is done in the [1-qc_and_call.bash script](https://github.com/sarahcolbert/sibling_roh_analysis/blob/main/code/1-qc_and_call.bash).

```
bash ${code_dir}1-qc_and_call.bash
```

## Step 3: Calculate Froh (within siblings) and give descriptive statistics

Running the code below will calculate Froh for each individual, then calculate Froh within siblings and finally, create tables that describe the sample (which will be included in the return of results).

Before runnning you will need to make sure that you have R installed and your version of R should include the tidyverse package.

```
Rscript ${code_dir}2-calc_froh.R
```

## Step 4: Calculate phenotypes (within siblings)
Follow method from Clark et al. 

NTS: edit code to work for all phenotypes...probably going to have to use for loop and paste
```
## load packages
library(tidyverse)

## load data
phenotype_data <- read.table("${phenotype_file.txt}", header = TRUE)

## calculate value of phenotype (height as example) to family mean
## make empty column to hold results
phenotype_data$height_sibs <- NA
## use for loop to get value for each individual
for(i in 1:length(phenotype_data$IID){
  spec_FID <- phenotype_data$FID[i]
  fam_vals <- phenotype_data %>% filter(FID=spec_FID)
  phenotype_data$height_sibs[i] <- phenotype_data$height[i]-mean(fam_vals$height)
}

write.table(phenotype_data, "within_sibs_phenotype_data.txt", row.names=FALSE, quote = FALSE)
```

## Step 5: Run within sibling models
What will we need as covariates?? 

Results to return: beta, se, p-val, standardized beta...anything else? 

```
## load packages
library(tidyverse)
library(lmerTest)

## load data
froh_data <- read.table("within_sibs_froh_data.txt", header = TRUE)
phenotype_data <- read.table("within_sibs_phenotype_data.txt", header = TRUE)
## merge
froh_phenotype_wsibs <- merge(froh_data, phenotype_data, by = "IID")

## calculate the effect of FROH on height (as an example) within-full-siblings
height_wsibs <- lmer(height_sibs ~ froh_sibs + covars, data = froh_phenotype_wsibs)
```
