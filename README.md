# Runs of homozygosity in the Within Family Consortium Analyses
This repository details a Standard Operating Procedure for the data analysts that will be performing the WFC runs of homozygosity analyses in individual cohorts. We provide a fully automated analysis pipeline that will analyze and format the results to be returned. We thank blabla for allowing us to adapt code and methods from their [sibling GWAS github](https://github.com/LaurenceHowe/SiblingGWAS).


## Pre-check: Pipeline Prerequisites and Requirements
The data requirements for the pipeline are as follows:

1) Sibling data (see Pre-Step 3 for how to define siblings).

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


## Pre-Step 2: Editing the config file
To edit the config file navigate to the sibling_roh_analysis directory and add the required information to the file. You can then run the config file using:  <br>
> chmod u+x config
> source ./config <br>


## Pre-Step 3: Defining siblings
The analysis pipeline requires data on siblings. We follow the suggestion of Howe et al., which is to include "all siblings from families with one or more pairs of genotyped dizygotic siblings. For example, in a family with a pair of monozygotic twins and an additional sibling, include both MZ twins and the sibling. The inclusion of both MZ twins should (very) modestly improve power by accounting for variation in the phenotypic outcome. If siblings have not been previously identified in the dataset, we suggest using [KING](https://www.kingrelatedness.com/) to infer siblings."

Instructions for defining siblings are provided by Howe et al. [here](https://github.com/LaurenceHowe/SiblingGWAS/wiki/0.1_Siblings).

*** provide example of what this file should look like

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

## OTHER PRE-STEPS- GOING TO ADD

## Step 1: Prep plink files
Genotype files must first be filtered to meet the following requirements:
1) exclude SNPs with >3% missingess
2) exclude SNPs with MAF < 5% 
3) exclude individuals with >3% missing data

### For this and all the other code chunks I will make actual scripts, just putting it down like this now so I have it all in one place
```
plink --bfile ${input_prefix} --maf 0.05 --geno 0.03 --mind 0.03 --make-bed --out ${cleaned_dir}/${input_prefix}_filtered
```

## Step 2: ROH Calling
Continuous ROH SNPs can be identified using PLINK with the following parameters:
1) homozyg-window-snp 50
2) homozyg-snp 50
3) homozyg-kb 1500
4) homozyg-gap 1000
5) homozyg-density 50
6) homozyg-window-missing 5
7) homozyg-window-het 1

```
plink --bfile ${cleaned_dir}/${input_prefix}_filtered --homozyg-window-snp 50 --homozyg-snp 50  --homozyg-kb 1500  --homozyg-gap 1000  --homozyg-density 50 --homozyg-window-missing 5 --homozyg-window-het 1 --out ${out_dir}/${input_prefix}_roh
```

## Step 3: Calculate Froh and give descriptive statistics
This will use R code 

This will create files to be returned to us: 
1) ${return_dir}/${input_prefix}_descriptive_roh_stats.txt
2) ${return_dir}/${input_prefix}_descriptive_sample_stats.txt

NTS: make script that does this
```
## load packages
library(tidyverse)

## load data
roh_data <- read.table("${out_dir}/${input_prefix}_roh.hom.indiv", header = TRUE)

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


## put descriptive roh stats in data table and save as file (this file will be returned to us)
descript_roh_stats <- cbind(min_vals, max_vals, mean_vals)
write.table(descript_roh_stats, "${return_dir}/${input_prefix}_descriptive_roh_stats.txt", row.names = TRUE, quote = FALSE)

## check how many individuals have froh = 0
zero_inds <- sum(froh_data$froh==0)
## check how many total individuals roh calls were performed for
total_inds <- length(froh_data$IID)
## check how many individuals are offspring of various relative relationships
half_sib_inds <- sum(froh_data$froh > 0.125)
first_cous_inds <- sum(froh_data$froh > 0.0625)
half_cous_inds <- sum(froh_data$froh > 0.03125)
sample <- "${sample_name}"

## put descriptive sample stats in data table and save as file (this file will be returned to us)
descript_sample_stats <- as.data.frame(cbind(sample,total_inds, zero_inds, half_sib_inds, first_cous_inds, half_cous_inds))
write.table(descript_sample_stats, "${return_dir}/${input_prefix}_descriptive_sample_stats.txt",row.names = FALSE, quote = FALSE)
```

## Step 4: Calculate Froh within siblings
Follow method from Clark et al. 

NTS: make script that does this
```
## load packages
library(tidyverse)

## load data
roh_data1 <- read.table("${out_dir}/${input_prefix}_roh.hom.indiv", header = TRUE)
covariate_data <- read.table("${covariate_file.txt}", header = TRUE)
## merge
roh_data <- merge(roh_data1, covariate_data, by = "IID")

## calc froh
roh_data$froh <- roh_data$KB/(2.77*10^6)

## filter to necessary columns
froh_data <- roh_data %>% select(FID, IID, NSEG, KB, froh)

## calculate value of froh relative to family mean
## make empty column to hold results
froh_data$froh_sibs <- NA
## use for loop to get value for each individual
for(i in 1:length(froh_data$IID){
  spec_FID <- froh_data$FID[i]
  fam_vals <- froh_data %>% filter(FID=spec_FID)
  froh_data$froh_sibs[i] <- froh_data$froh[i]-mean(fam_vals$froh)
}

write.table(froh_data, "within_sibs_froh_data.txt", row.names=FALSE, quote = FALSE)
```

## Step 5: Calculate phenotype within siblings
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

## Step 6: Run within sibling models
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

