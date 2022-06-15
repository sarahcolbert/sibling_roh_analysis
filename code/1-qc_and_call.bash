source ./config

## step 1: variant QC for non-imputed data
## use plink to exclude SNPs with >3% missingess, with MAF < 5% and not on chr1-22
${plink} --bfile ${data_dir}${input_prefix} --chr 1-22 --maf 0.05 --geno 0.03 --make-bed --out ${processed_dir}${input_prefix}_variant_qc

## step 2: sample QC for non-imputed data
## exclude individuals with >3% missing data
${plink} --bfile ${processed_dir}${input_prefix}_variant_qc --mind 0.03 --make-bed --out ${processed_dir}${input_prefix}_filtered

## use plink to identify continuous ROH SNPs using non-imputed data
${plink} --bfile ${processed_dir}${input_prefix}_filtered --homozyg-window-snp 50 --homozyg-snp 50  --homozyg-kb 1500  --homozyg-gap 1000  --homozyg-density 50 --homozyg-window-missing 5 --homozyg-window-het 1 --out ${processed_dir}${input_prefix}_roh

## use plink to calculate Fhat3 from imputed data
${plink} --bfile ${processed_dir}${imputed_prefix} --autosome --ibc --out ${processed_dir}${input_prefix}_ibc
