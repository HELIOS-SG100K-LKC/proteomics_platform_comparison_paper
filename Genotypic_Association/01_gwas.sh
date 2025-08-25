#!/bin/bash

# Move to working directory
cd /your/directory/path || exit 1

# === Run fastGWA using GCTA ===
gcta64 --bfile your_data_folder/genotype_data_prefix \
       --grm-sparse your_data_folder/genotype_data_prefix \
       --fastGWA-mlm \
       --pheno your_pheno_folder/your_pheno.txt \ ## One Phenotype at a time (FID,IID, Pheno-Col) | Rank Inverse Normal transformed
       --qcovar your_covariate_folder/quantitative_covariates.txt \
       --covar your_covariate_folder/categorical_covariates.txt \
       --thread-num 5\
       --out your_output_folder/output_prefix

# === Clump significant associations using PLINK ===
plink --bfile your_data_folder/genotype_data_prefix \
      --clump your_output_folder/output_prefix.fastGWA \
      --clump-p1 5e-08 \
      --clump-r2 0.01 \
      --clump-kb 500 \
      --clump-field P \
      --clump-snp-field SNP \
      --out your_clumped_folder/output_prefix \
      --threads 5
