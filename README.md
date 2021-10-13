# Runs of homozygosity in the Within Family Consortium Analyses
This repository details a Standard Operating Procedure for the data analysts that will be performing the WFC runs of homozygosity analyses in individual cohorts. We provide a fully automated analysis pipeline that will analyze and format the results to be returned.

## Pre-Step 0.1: Clone and set up working directory
make some kind of configuration file...what kind of variables will we need?
variables:
* working directory
* phenotype files?
* genotype files (hopefully plink format)
* directories to be made

## Step 1: Sibling identification
I see three options: follow sibling identification method from sibling gwas, follow method from clark et al., or let individual cohorts submit however they already have sibs coded? Or we can come up with our own method?

## Step 2: QC and file prep
NTS: should be able to use most code from /froh_abcd/eur/code/1-prep.bash

## Step 3: Calculate ROH and Froh, Num segments, etc. using plink
NTS: should be able to use most code from /froh_abcd/eur/code/2-cacl_froh_v1.R
