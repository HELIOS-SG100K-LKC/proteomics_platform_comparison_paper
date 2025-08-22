#Generate Unique Files


#SomaLogic pre-ANML

#SomaLogic Data
Dat <- read.csv("SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label_SomaLogic <- Dat[,c(1:20)]
ProtMatrix_SomaLogic <- Dat[,c(21:ncol(Dat))]
rm(Dat)

ProtMatrix_SomaLogic <- ProtMatrix_SomaLogic[which(Label_SomaLogic$SampleType == "Sample"),]
Label_SomaLogic <- Label_SomaLogic[which(Label_SomaLogic$SampleType == "Sample"),]
ProtMatrix_SomaLogic <- log2(ProtMatrix_SomaLogic)

#Analyte Annotation
Analyte_Annotation_SomaLogic <- read.csv("SomaLogic_pre-anml_Preprocessing/Somalogic_Analyte_Annotation_All.csv")
rownames(Analyte_Annotation_SomaLogic) <- Analyte_Annotation_SomaLogic$SeqId

colnames(ProtMatrix_SomaLogic) <- Analyte_Annotation_SomaLogic$SeqId
rm(Analyte_Annotation_SomaLogic)

#LoD
LoD <- read.csv("SomaLogic_pre-ANML_LoD/Somalogic_LOD.csv")

library(dplyr)
df_unique <- LoD %>%
  group_by(UniProt) %>%
  slice_min(PropBelowLOD, with_ties = FALSE) %>%
  ungroup()
rm(LoD)

write.csv(df_unique, "SomaLogic_Analyte_Annotation_All_UniProt_pre-ANML.csv", row.names = FALSE)

ProtMatrix_SomaLogic <- ProtMatrix_SomaLogic[,df_unique$SeqID]
colnames(ProtMatrix_SomaLogic) <- df_unique$UniProt
rm(df_unique)

A <- cbind(Label_SomaLogic, ProtMatrix_SomaLogic)

write.csv(A, "SomaLogic_Merged_All_UniProt_pre-ANML.csv", row.names = FALSE)

rm(A)
rm(Label_SomaLogic)
rm(ProtMatrix_SomaLogic)


#SomaLogic ANML

#SomaLogic Data
Dat <- read.csv("Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label_SomaLogic <- Dat[,c(1:20)]
ProtMatrix_SomaLogic <- Dat[,c(21:ncol(Dat))]
rm(Dat)

ProtMatrix_SomaLogic <- ProtMatrix_SomaLogic[which(Label_SomaLogic$SampleType == "Sample"),]
Label_SomaLogic <- Label_SomaLogic[which(Label_SomaLogic$SampleType == "Sample"),]
ProtMatrix_SomaLogic <- log2(ProtMatrix_SomaLogic)

#Analyte Annotation
Analyte_Annotation_SomaLogic <- read.csv("Somalogic_Analyte_Annotation_All.csv")
rownames(Analyte_Annotation_SomaLogic) <- Analyte_Annotation_SomaLogic$SeqId

colnames(ProtMatrix_SomaLogic) <- Analyte_Annotation_SomaLogic$SeqId
rm(Analyte_Annotation_SomaLogic)

#LoD
LoD <- read.csv("Somalogic_LOD.csv")

library(dplyr)
df_unique <- LoD %>%
  group_by(UniProt) %>%
  slice_min(PropBelowLOD, with_ties = FALSE) %>%
  ungroup()
rm(LoD)

write.csv(df_unique, "SomaLogic_Analyte_Annotation_All_UniProt.csv", row.names = FALSE)

ProtMatrix_SomaLogic <- ProtMatrix_SomaLogic[,df_unique$SeqID]
colnames(ProtMatrix_SomaLogic) <- df_unique$UniProt
rm(df_unique)

A <- cbind(Label_SomaLogic, ProtMatrix_SomaLogic)

write.csv(A, "SomaLogic_Merged_All_UniProt.csv", row.names = FALSE)

rm(A)
rm(Label_SomaLogic)
rm(ProtMatrix_SomaLogic)


#Olink

#Olink Data
Dat <- read.csv("Olink_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label_Olink <- Dat[,c(1:19)]
ProtMatrix_Olink <- Dat[,c(20:ncol(Dat))]
rm(Dat)

ProtMatrix_Olink <- ProtMatrix_Olink[which(Label_Olink$SampleType == "SAMPLE"),]
Label_Olink <- Label_Olink[which(Label_Olink$SampleType == "SAMPLE"),]

#LoD
LoD <- read.csv("Olink_LOD.csv")

library(dplyr)
df_unique <- LoD %>%
  group_by(UniProt) %>%
  slice_min(PropBelowLOD, with_ties = FALSE) %>%
  ungroup()
rm(LoD)

write.csv(df_unique, "Olink_Assay_Annotation_All_UniProt.csv", row.names = FALSE)

ProtMatrix_Olink <- ProtMatrix_Olink[,df_unique$OlinkID]
colnames(ProtMatrix_Olink) <- df_unique$UniProt
rm(df_unique)

A <- cbind(Label_Olink, ProtMatrix_Olink)

write.csv(A, "Olink_Merged_All_UniProt.csv", row.names = FALSE)

rm(A)
rm(Label_Olink)
rm(ProtMatrix_Olink)


#Thermo Fisher

#Thermo Fisher Data
library(readr)
Dat <- read_csv("Thermo_Fisher_Merged_All.csv")
Label_Thermo_Fisher <- Dat[,c(1:10)]
ProtMatrix_Thermo_Fisher <- Dat[,c(11:ncol(Dat))]
rm(Dat)

Dat <- read_csv("Thermo_Fisher_Merged_All_Imputed.csv")
Label_Thermo_Fisher2 <- Dat[,c(1:10)]
ProtMatrix_Thermo_Fisher2 <- Dat[,c(11:ncol(Dat))]
rm(Dat)

#LoD
LoD <- read.csv("Thermo_Fisher_ProteinMissingness.csv")

library(dplyr)
df_unique <- LoD %>%
  group_by(UniProt) %>%
  slice_min(Missing_Rate_Protein, with_ties = FALSE) %>%
  ungroup()
rm(LoD)

write.csv(df_unique, "Thermo_Fisher_Protein_Annotation_All_UniProt.csv", row.names = FALSE)

ProtMatrix_Thermo_Fisher <- ProtMatrix_Thermo_Fisher[,df_unique$Protein]
colnames(ProtMatrix_Thermo_Fisher) <- df_unique$UniProt
ProtMatrix_Thermo_Fisher2 <- ProtMatrix_Thermo_Fisher2[,df_unique$Protein]
colnames(ProtMatrix_Thermo_Fisher2) <- df_unique$UniProt
rm(df_unique)

A <- cbind(Label_Thermo_Fisher, ProtMatrix_Thermo_Fisher)

write.csv(A, "Thermo_Fisher_Merged_All_UniProt.csv", row.names = FALSE)

B <- cbind(Label_Thermo_Fisher2, ProtMatrix_Thermo_Fisher2)

write.csv(B, "Thermo_Fisher_Merged_All_UniProt_Imputed.csv", row.names = FALSE)

rm(A)
rm(B)
rm(Label_Thermo_Fisher)
rm(ProtMatrix_Thermo_Fisher)
rm(Label_Thermo_Fisher2)
rm(ProtMatrix_Thermo_Fisher2)



