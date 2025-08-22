# remove environment variables
rm(list = ls())

# load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)


# read list of unique uniprot-analyte pairs
trim_ids <- fread("cis-pQTL_unique-uniprot-analyte-pair.txt", 
                  header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)


# read finemapping results and filter for uniprots
df_finemap <- fread("finemap/All_cs_olink.txt", 
               header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

df_finemap %>% 
  separate(Protein_SNP_Pair, into = c("analyte_ID", "lead_SNP"), sep = "-") %>% 
  mutate(analyte_ID = sub(".*_", "", analyte_ID)) %>% 
  filter(analyte_ID %in% trim_ids$analyte_ID) -> df_finemap_trim
  
length(unique(df_finemap_trim$lead_SNP))
length(unique(df_finemap_trim$analyte_ID))


# get causal variants
df_finemap_trim %>% 
  filter(CS_PIP>0.5) -> finemap_causal

length(unique(finemap_causal$lead_SNP))

