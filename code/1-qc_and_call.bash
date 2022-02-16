## use plink to exclude SNPs with >3% missingess or MAF < 5% AND exclude individuals with >3% missing data
plink --bfile ${data_dir}${input_prefix} --maf 0.05 --geno 0.03 --mind 0.03 --make-bed --out ${processed_dir}${input_prefix}_filtered

## use plink to identify continuous ROH SNPs
plink --bfile ${processed_dir}${input_prefix}_filtered --homozyg-window-snp 50 --homozyg-snp 50  --homozyg-kb 1500  --homozyg-gap 1000  --homozyg-density 50 --homozyg-window-missing 5 --homozyg-window-het 1 --out ${processed_dir}${input_prefix}_roh

## use plink to calculate F_GRM
${gcta} --bfile ${processed_dir}${input_prefix}_filtered --autosome --ibc --out ${processed_dir}${input_prefix}_ibc
