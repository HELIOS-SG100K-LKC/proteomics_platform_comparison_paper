#!/bin/bash

# Navigate to your working directory
cd /your/directory/path || exit 1

fine_map() {
  i=$1


## File containing [S.No Protein_ID SNP_ID Chromosome Position]

  protein=$(awk '{if($1=="'${i}'") print $2}' your_filename.txt)
  snp=$(awk '{if($1=="'${i}'") print $5}' your_filename.txt)
  chr=$(awk '{if($1=="'${i}'") print $3}' your_filename.txt)
  lower=$(awk '{if($1=="'${i}'") print $4-1500000}' your_filename.txt)
  upper=$(awk '{if($1=="'${i}'") print $4+1500000}' your_filename.txt)

  echo "extracting ${protein}-${snp} : chr${chr}:${lower}-${upper}"

  awk -v chr="$chr" -v lower="$lower" -v upper="$upper" 'BEGIN{OFS=FS="\t"} NR==1 || ($1 == chr && $3 >= lower && $3 <= upper)' \
    your/sumstat_directory/your_summary_file_prefix_${protein}.fastGWA > tmp-sumstats/${protein}-${snp}.sumstats

  SuSiEx \
    --sst_file=tmp-sumstats/${protein}-${snp}.sumstats \
    --n_gwas=197 \
    --ref_file=your/reference/path/genomic_data_all \
    --ld_file=your/ld_directory/${protein}-${snp} \
    --out_dir=Results \
    --out_name=olink_${protein}-${snp} \
    --chr=${chr} \
    --bp=${lower},${upper} \
    --chr_col=1 \
    --snp_col=2 \
    --bp_col=3 \
    --a1_col=4 \
    --a2_col=5 \
    --eff_col=8 \
    --se_col=9 \
    --pval_col=10 \
    --plink=~/bin/plink \
    --threads=4 \
    --level=0.90
}

export -f fine_map

# Run in parallel for all entries (adjust range if needed)
parallel -j 10 fine_map ::: {1..305}
