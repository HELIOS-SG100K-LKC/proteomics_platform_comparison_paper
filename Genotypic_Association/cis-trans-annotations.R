# remove environment variables
rm(list = ls())

# load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(biomaRt)


# read pqtl summary statistics
df <- fread("pQTLS_gene-annot_all_p5e08.txt", 
            header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

colnames(df)[1] <- "ID"
length(unique(df$ID))

# read platform specific annotation file
annot <- fread("Annotation_All.csv", 
                 header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

colnames(annot)[1] <- "analyte_ID"
annot %>% 
  dplyr::select(analyte_ID, UniProt) -> annot


# annotate pqtl data with uniprot id in annotation file
df %>% 
  mutate(analyte_ID = sub(".*_", "", ID)) %>% 
  inner_join(annot, by = "analyte_ID") -> df_annot


# get ensembl id gene positions
ensembl <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", mirror = "useast")

ids <- unique(df_annot$UniProt)
ensembl_df <- getBM(attributes=c("uniprot_gn_id", "ensembl_gene_id", "chromosome_name", "start_position","end_position", "external_gene_name"),
                    filters="uniprot_gn_id",values = ids, mart=ensembl)
ensembl_df %>% 
  filter(chromosome_name %in% as.character(c(1:21))) -> ensembl_df


# merge ensembl data with pqtl data 
ensembl_df %>% 
  dplyr::rename(UniProt = uniprot_gn_id) -> ensembl_df

df_annot %>% 
  inner_join(ensembl_df, by = "UniProt") -> df_annot_ensembl


# annotate cis and trans
df_annot_ensembl %>% 
  mutate(Chr = as.character(Chr)) %>% 
  mutate(pqtl = ifelse(Chr == chromosome_name & BP > (start_position - 1000000) & BP < (end_position + 1000000), "cis", "trans")) -> df_annot_ensembl_pqtl


# save file
# write.table(df_annot_ensembl_pqtl, "pQTLS_gene-annot_all_p5e08_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# number of cis associations
df_annot_ensembl_pqtl %>% 
  filter(P < 5e-08  & pqtl == "cis") -> pqtl_cis

nrow(pqtl_cis) #total cis pqtls
length(unique(pqtl_cis$UniProt)) #total uniprots with a cis pqtl

# number of trans association
n <- length(unique(annot$UniProt))
pthresh <- 5e-08/n

df_annot_ensembl_pqtl %>% 
  filter(P < pthresh  & pqtl == "trans") -> pqtl_trans

nrow(pqtl_trans) #total trans pqtls
length(unique(pqtl_trans$UniProt)) #total uniprots with a trans pqtl

