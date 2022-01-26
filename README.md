** fix intro to give proper credit


# Runs of homozygosity in the Within Family Consortium Analyses
This repository details a Standard Operating Procedure for the data analysts that will be performing the WFC runs of homozygosity analyses in individual cohorts. We provide a fully automated analysis pipeline that will analyze and format the results to be returned. We thank blabla for allowing us to adapt code and methods from their [sibling GWAS github](https://github.com/LaurenceHowe/SiblingGWAS).

The following analyses should be repeated for each subsample of genetic ancestry available in your sample. Please be sure to indicate the genetic ancestry of your sample in the config file (see Pre-Step 3).


## Pre-check: Pipeline Prerequisites and Requirements
The data requirements for the pipeline are as follows:

1) Sibling data (see Pre-Step 2 for how to define siblings).

2) Called genotypes in plink binary format. Imputed data is not permitted. (see Pre-Step 4 for information on input file requirements and scripts for file conversion).

3) Complete covariate data (see Pre-Step 4 for information on what covariates should be included and the file format).

4) Phenotype data for various traits (see Pre-Step 4 for more information on phenotypes requested, coding phenotypes, and phenotype file format).

The software requirements for the pipeline are as follows:

* Plink 1.9
* R (tidyverse, lmerTest and modelr installed)
* KING (if siblings are not already defined)

## Pre-Step 1: Downloading and running the pipeline

Navigate to the directory where you want to download the repository. The repository can then be downloaded using git: <br>
> git clone https://github.com/sarahcolbert/sibling_roh_analysis <br>


## Pre-Step 2: Defining siblings
The within-sibling analysis requires data on siblings. If siblings have not been previously identified in the dataset, we suggest using [KING](https://www.kingrelatedness.com/) to infer siblings.

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
A covariate file should be provided containing the following columns:
IID, age (defined as 2022 (or 2021?) minus birth year), sex (male = 1, female = 0), first 10 genomic principal components.


Column names should exactly match the example below:

```
IID age sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
1001 31 0 0.0067 0.0042 0.0019 -8e-04 -0.0043 0.0324 -2e-04 -0.0033 -0.0061 0.0118
1002 26 1 -0.0289 -0.0035 0.0034 0.0288 0.0021 -0.0187 -0.0013 -0.0054 0.005 0.0211
1011 12 1 0.0044 0.0022 -0.0141 -0.011 -0.0064 0.0248 -0.0157 0.008 0.0089 0.0028
1012 18 0 0.0122 0.0019 0.0034 -0.017 0.001 0.003 -0.0217 -0.0016 -0.0126 -3e-04
1013 14 0 0.007 5e-04 0.0157 0.0017 -0.0115 -0.0079 -0.0083 -0.017 0.0147 0.0227
```

### Phenotype data
We will analyse a wide range of medical, social and behavioral phenotypes based on previous findings in the literature and hypothesized relationships with autozygosity. The [analysis plan](https://docs.google.com/document/d/1weNXniAY8X03ZYm1k-TmZqH4j4h7vloHl_meiSpceII/edit#bookmark=id.nu9tiucjoq87) outlines the phenotypes that we propose to include.  

The first column in the phenotype file should be IID, then followed by the phenotypes available in your dataset. Column names for the phenotypes do not matter. Missing data should be coded as NA. For example, your phenotype file should look something like this:

```
IID   Pheno1  Pheno2   Pheno3
1001  7       21.6     0
1002  3       24.3     1
1011  9       17.7     1
1012  12      19.8     1
1013  12      19.3     0
1041  8       NA       0
1042  5       22.4     1
```

The code will ensure that analyses do not include families with missing phenotype data for at least one sibling (e.g., family 104 would be removed from the analyses which use pheno3). A csv file named pheno_descriptions_STUDYNAME.csv (replace STUDYNAME) should be returned that includes a description of the phenotypes in your dataset. Please format this csv file with quotes. We provide an [example](https://github.com/sarahcolbert/sibling_roh_analysis/blob/main/pheno_descriptions_STUDYNAME.csv) of what this might include.


## Steps 1: QC and ROH Calling
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

## Step 2: Calculate Froh (+ within siblings) and give descriptive statistics

Running the code below will calculate Froh for each individual, then calculate Froh within siblings and finally, create tables that describe the sample (which will be included in the return of results).

Before runnning you will need to make sure that you have R installed and your version of R should include the tidyverse package.

```
Rscript ${code_dir}2-calc_froh.R
```

## Step 3: Calculate phenotypes (within siblings)
Follow method from Clark et al., this code is used to calculate phenotypes within siblings. This script will also create tables that describe the distribution of the phenotypes in the sample (which will be included in the return of results).

```
Rscript ${code_dir}3-calc_phenos.R
```

## Step 4: Run within sibling models and between family models

will need to scale some variables

<strike>NTS: may need to prep phenotype files more before this (so raw phenotypes work with btwn fam analysis)</strike> wait to do this until phenotype appendix done

```
Rscript ${code_dir}4-run_models.R

```
