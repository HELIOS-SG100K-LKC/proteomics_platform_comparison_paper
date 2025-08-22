# remove environment variables
rm(list = ls())

# load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)


# read eqtl data
eqtlgen <- fread("eqtlgen_grch38_all.sumstats", header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)


# read annotated pqtl summary data
df <- fread("pQTLS_gene-annot_all_p5e08_annotated.txt", 
               header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)


# read list of unique uniprot-analyte pairs
trim_ids <- fread("cis-pQTL_unique-uniprot-analyte-pair.txt", 
            header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)


# extract cis pqtls
df %>% 
  filter(pqtl == "cis") -> cis_pqtl


# lookup cis-eQTLs
cis_pqtl %>%
  filter(analyte_ID %in% trim_ids$analyte_ID) %>% 
  dplyr::rename("Gene" = "ensembl_gene_id") %>% 
  inner_join(eqtlgen, by = c("Chr", "BP", "Gene")) -> cis_pqtl_eqtlgen

nrow(cis_pqtl_eqtlgen)

cis_pqtl_eqtlgen %>% 
  filter(p<0.05) -> cis_pqtl_eqtlgen

