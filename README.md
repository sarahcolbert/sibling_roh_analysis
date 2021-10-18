# Runs of homozygosity in the Within Family Consortium Analyses
This repository details a Standard Operating Procedure for the data analysts that will be performing the WFC runs of homozygosity analyses in individual cohorts. We provide a fully automated analysis pipeline that will analyze and format the results to be returned. We thank blabla for allowing us to adapt code and methods from their [sibling GWAS github](https://github.com/LaurenceHowe/SiblingGWAS).


## Pre-check: Pipeline Prerequisites and Requirements
The data requirements for the pipeline are as follows:

1) Sibling data (see Pre-Step 3 for how to define siblings).

2) Genotype data in plink binary format (see Pre-Step 4 for information on input file requirements and scripts for file conversion).

3) Complete covariate data (see Pre-Step 4 for information on what covariates should be included and the file format).

4) Phenotype data for various traits (see Pre-Step 4 for more information on phenotypes requested, coding phenotypes, and phenotype file format).

The software requirements for the pipeline are as follows:

* Plink 1.9
* R
* KING (if siblings are not already defined)

## Pre-Step 1: Downloading and running the pipeline

Navigate to the directory where you want to download the repository. The repository can then be downloaded using git: <br>
> git clone https://github.com/sarahcolbert/siblingRohAnalyses <br>
<br>
Once the repository is downloaded, run the following command to check that files have downloaded properly: <br>

> head ./siblingRohAnalyses/resources/parameters <br>

<br>


## Pre-Step 2: Editing the config file

working on this...not sure what will need to be in it


## Pre-Step 3: Defining siblings
The analysis pipeline requires data on siblings. We follow the suggestion of Howe et al., which is to include "all siblings from families with one or more pairs of genotyped dizygotic siblings. For example, in a family with a pair of monozygotic twins and an additional sibling, include both MZ twins and the sibling. The inclusion of both MZ twins should (very) modestly improve power by accounting for variation in the phenotypic outcome. If siblings have not been previously identified in the dataset, we suggest using [KING](https://www.kingrelatedness.com/) to infer siblings."

Instructions for defining siblings are provided by Howe et al. [here](https://github.com/LaurenceHowe/SiblingGWAS/wiki/0.1_Siblings).

## Pre-Step 4: Input Files
### Genotype data
You will need genotype data in PLINK binary format. The pipeline requires the input (.bed, .bim, .fam) files to satisfy the following requirements:

a) PLINK binary format (.bed .bim .fam) files. The first two columns must contain family IDs (FID) and individual IDs (IIDs). FIDs should be common between siblings (but unique between sets of siblings) and IIDs should be unique for each individual.

b) add stuff as we need it

### Covariate data
A covariate file should be provided containing the following columns, if available:

IID, age, sex, batch, first 10 genomic principal components, anything else? 

# NOTES/OLD FRAMEWORK

## Pre-Step 0.1: Clone and set up working directory
make some kind of configuration file...what kind of variables will we need?
variables:
* working directory
* phenotype file(s?) (can offer code for gathering phenotypes if we need to..or some code for renaming phenos to be uniform? or is that too much to ask? I am fine changing pheno names in house)
* covariates file: WHAT DO WE WANT AS COVS? (will we need to add code for scaling covariates?)
* genotype files (hopefully plink format, but need to offer code for converting or merging? Do we need plink files not split by chr? That's what I was using for our ABCD analyses but can it work split up by chr?)
* directories to be made

## Pre-Step 0.2: Summary
adapt from: https://github.com/LaurenceHowe/SiblingGWAS/wiki/2.0_summary

## Step 1: Sibling identification
I see three options: follow sibling identification method from sibling gwas, follow method from clark et al., or let individual cohorts submit however they already have sibs coded? Or we can come up with our own method?

## Step 2: QC and file prep (could make master script to run this and step 3 all at once)
NTS: should be able to use most code from htcf /suri/projects/froh_abcd/white/code/4-prep_v2.bash
exclude SNPs with >3% missingess or MAF < 5% AND exclude individuals with >3% missing data

## Step 3: Calculate ROH, Froh, Num segments, etc. using plink
NTS: should be able to use most code from htcf /suri/projects/froh_abcd/white/code/5-calc_froh_clark_v2.R

## Step 4: Get descriptive statistics
NTS: should be able to use most code from local ./Desktop/projects/froh_abcd/white/scripts/1-autozygosity_variables_distributions.Rmd
We don't want any figures from the cohorts, correct? I'm assuming cohort specific figs won't be helpful and we'll just make figs after we meta-analyze. We just won't be able to make any figs like the violin plot type since that requires individual data points, right?

## Step 5: Run models for phenotypes
Will need to look at clark method a little more but can also use code from local ./Desktop/projects/froh_abcd/white/scripts/5-child_cog_clark.Rmd
