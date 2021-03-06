# Runs of homozygosity in the Within Family Consortium Analyses
This repository details a Standard Operating Procedure for the data analysts that will be performing the WFC runs of homozygosity analyses in individual cohorts. We provide a fully automated analysis pipeline that will analyze and format the results to be returned. We thank Howe et al. for allowing us to adapt code and methods from their [sibling GWAS github](https://github.com/LaurenceHowe/SiblingGWAS).

The following analyses should be repeated for each subsample of genetic ancestry available in your sample. Please be sure to indicate the genetic ancestry of your sample in the config file (see Pre-Step 3).


## Pre-check: Pipeline Prerequisites and Requirements
The data requirements for the pipeline are as follows:

1) Sibling data (see Pre-Step 2 for how to define siblings).

2) Imputed and non-imputed genotype data in plink binary format. (see Pre-Step 4 for information on input file requirements).

3) Complete covariate data (see Pre-Step 4 for information on what covariates should be included and the file format).

4) Phenotype data for various traits (see Pre-Step 4 for more information on phenotypes requested, coding phenotypes, and phenotype file format).

The software requirements for the pipeline are as follows:

* Plink 1.9 (if you are running into bugs make sure you are using the [newest stable version available](https://www.cog-genomics.org/plink2/))
* R >= 3.3
  * [tidyverse](https://github.com/tidyverse/tidyverse)
  * [lmerTest](https://cran.r-project.org/web/packages/lmerTest/index.html)
  * [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
  * [lfe](https://cran.r-project.org/web/packages/lfe/index.html)
* KING (if siblings are not already defined)

## Pre-Step 1: Downloading and running the pipeline

Navigate to the directory where you want to download the repository. The repository can then be downloaded using git:

```
git clone https://github.com/sarahcolbert/sibling_roh_analysis
```


## Pre-Step 2: Defining siblings
The within-sibling analysis requires data on siblings. If siblings have not been previously identified in the dataset, we suggest using [KING](https://www.kingrelatedness.com/) to infer siblings.

Instructions for defining siblings are provided by Howe et al. [here](https://github.com/LaurenceHowe/SiblingGWAS/wiki/0.1_Siblings).

We also provide an [example script](https://github.com/sarahcolbert/sibling_roh_analysis/blob/main/code/siblings_method_2.R) from David Clark with an alternate method for defining siblings.

## Pre-Step 3: Editing the config file
Navigate to the sibling_roh_analysis directory and edit the config file by adding the required information. If you would like to see an example of what a completed config file might look like please refer to the "myconfig" file (but be sure to edit the "config" file and not this one). You can then run the config file using:

```
source ./config
```

## Pre-Step 4: Input Files

In  Pre-Step 5 we provide instructions for running a script that will check to make sure these and other files meet our requirements, but please first use this information to check your data to the best of your ability.

### Non-imputed genotype data
You will need non-imputed genotype data in PLINK binary format for all F<sub>ROH</sub> analyses. **Analyses can only be run using SNPs on chromosomes 1-22. Any SNPs not on chromosomes 1-22 will be removed in step 1.** The pipeline requires the input files to satisfy the following requirements:

a) PLINK binary format (.bed .bim .fam) files.

b) The first two columns must contain family IDs (FID) and individual IDs (IIDs).

c) FIDs should be common between siblings (but unique between sets of siblings) and IIDs should be unique for each individual.

### Imputed genotype data
Imputed data will be needed to calculate Fhat3. Imputed data should meet the same requirements as non-imputed data, with a few additional requirements:

a) must be filtered to INFO > 0.6

b) must be filtered to MAF > 0.01

c) FID and IID values must be identical to those in the non-imputed genotype data. If they do not match, you will receive an error during the 'genetics' check. 

If your imputed data is in VCF or BGEN format, Howe et. al. [provide information](https://github.com/LaurenceHowe/SiblingGWAS/wiki/0.2_ImputedGenotypeData) for how to convert these files to the plink binary format.

### Covariate data
A covariate file should be provided containing the following columns:
IID, age (defined as 2022 minus birth year), sex (male = 1, female = 0), first 10 genomic principal components. Please only include one row per individual. Duplicate rows should not exist in the covariate file, but will be removed in step 2 if they do exist.


Column names should exactly match the example below:

```
IID age sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
1001 31 0 0.0067 0.0042 0.0019 -8e-04 -0.0043 0.0324 -2e-04 -0.0033 -0.0061 0.0118
1002 26 1 -0.0289 -0.0035 0.0034 0.0288 0.0021 -0.0187 -0.0013 -0.0054 0.005 0.0211
1011 12 1 0.0044 0.0022 -0.0141 -0.011 -0.0064 0.0248 -0.0157 0.008 0.0089 0.0028
1012 18 0 0.0122 0.0019 0.0034 -0.017 0.001 0.003 -0.0217 -0.0016 -0.0126 -3e-04
1013 14 0 0.007 5e-04 0.0157 0.0017 -0.0115 -0.0079 -0.0083 -0.017 0.0147 0.0227
```
Please reach out to us at sarah.colbert@wustl.edu or emma.c.johnson@wustl.edu if you have additional covariates you would like to include for your sample.


### Phenotype data
We will analyze a wide range of medical, social and behavioral phenotypes based on previous findings in the literature and hypothesized relationships with autozygosity. The [analysis plan](https://docs.google.com/document/d/1weNXniAY8X03ZYm1k-TmZqH4j4h7vloHl_meiSpceII/edit#bookmark=id.nu9tiucjoq87) outlines the phenotypes that we propose to include and the desired coding formats for each phenotype.  

The first column in the phenotype file should be IID, then followed by the phenotypes available in your dataset. Column names for the phenotypes do not matter, **except for height, which should be labelled "Height" and BMI, which should be labelled "BMI"**. Missing data should be coded as NA. For binary phenotypes, please use a binary coding of control = 0 and case = 1. For example, your phenotype file should look something like this:

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

A csv file named pheno_descriptions_STUDYNAME.csv (replace STUDYNAME) should be returned that includes a description of the phenotypes in your dataset. This should include the column name of the phenotype, which phenotype in the analysis plan it corresponds to and any derivations from the preferred coding we outline in the analysis plan google doc. Please format this csv file with quotes.

## Pre-Step 5: Checks

The set-up script runs checks to ensure that the input files are in the correct format and checks the software requirements. This script has 7 steps:

1. Checks basic information such as the study name.
2. Checks installation of PLINK, R and necessary R packages.
3. Checks the genotype file formatting for compatibility with the pipeline.
4. Estimates pairwise IBD between individuals with the same FID.
5. Evaluates the sibling pairs by running checks for parent-offspring pairs and monozygotic twins.
6. Checks the covariate file formatting for compatibility with the pipeline.
7. Checks the phenotype file formatting for compatibility with the pipeline.

Each step of the script can be run separately using the following as arguments:
'config', 'requirements', 'genetics', 'rel', 'siblings', 'covariates', 'phenotypes'

The steps should be run in order. If you are running multiple steps and encounter a warning or error at a specific step, the script will fail and subsequent steps will not be run. You should address this warning or error, then re-run the step it occurred at. Once you pass this step you may continue onto the next steps.

To run through all steps in order:

```
bash ${code_dir}0-checks.bash
```

For example, to run only the phenotypes section:

```
bash ${code_dir}0-checks.bash phenotypes
```

Steps 4 and 5 are included as sanity checks for derived siblings. Some cohorts which are confident in their derived pedigrees may wish to run through skipping these steps as they can take a longer amount of time compared to the other checks (for example, in a sample of ~5500 individuals, these steps took ~10 CPU mins). To skip these steps you can use the code:

```
bash ${code_dir}0-checks.bash skipsib
```

## Step 1: QC and ROH Calling
Genotype files must first be filtered to meet the following requirements:
1) exclude SNPs with >3% missingess
2) exclude SNPs with MAF < 5%
3) exclude SNPs not on chromosomes 1-22
4) exclude individuals with >3% missing data

Continuous ROH SNPs are identified using PLINK with the following parameters:
1) homozyg-window-snp 50
2) homozyg-snp 50
3) homozyg-kb 1500
4) homozyg-gap 1000
5) homozyg-density 50
6) homozyg-window-missing 5
7) homozyg-window-het 1

Fhat3 is calculated using the --ibc flag in PLINK. 

Running the code below will perform QC on the _non-imputed data_ and then use the QCed _non-imputed_ files to call ROHs and the provided _imputed_ files to calculate Fhat3.
Before running you will need to make sure that you have installed PLINK 1.9 and that you specified the location of the executable file during configuration [1-qc_and_call.bash script](https://github.com/sarahcolbert/sibling_roh_analysis/blob/main/code/1-qc_and_call.bash).

```
bash ${code_dir}1-qc_and_call.bash
```

## Step 2: Calculate Froh and gather descriptive statistics

Running the code below will calculate Froh for each individual and create tables that describe the sample (which will be included in the return of results).

```
Rscript ${code_dir}2-clean_calc_roh_data.R
```

## Step 3: Between family models
This code is used to run between family models for Froh and Fhat3 both separately and together. This script will also clean the phenotype data for the between family analyses and create tables that describe the distribution of the phenotypes in the sample (which will be included in the return of results). Please check the log output after running this code so that any warnings or errors can be reported.

```
Rscript ${code_dir}3-btwn_family_analysis.R
```

## Step 4: Run within sibling models
This code is used to run within sibling models for Froh and Fhat3 both separately and together. This script will also clean the phenotype data for the within sibling analyses and create tables that describe the distribution of the phenotypes in the sample (which will be included in the return of results). Please check the log output after running this code so that any warnings or errors can be reported.

```
Rscript ${code_dir}4-within_family_analysis.R
```

## Step 5: Return your results
To prepare a folder with the files to be uploaded to the box link that has been sent to you, please use the following code inside the sibling_roh_analysis directory:
```
tar -zcvf ${output_name}_results.tar.gz ${output_dir}
```
This will create a compressed file, following the naming scheme "STUDY_NAME_ANCESTRY_ANALYST_INITIALS_DATE_results.tar.gz", containing a directory with the results to return.
