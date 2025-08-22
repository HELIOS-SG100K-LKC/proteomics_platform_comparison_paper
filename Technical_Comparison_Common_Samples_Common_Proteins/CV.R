#CV


#SomaLogic pre-ANML

#SomaLogic Data
Dat <- read.csv("Somalogic_Merged_All_UniProt_pre-ANML.csv")
rownames(Dat) <- Dat$UniqueID
Label_SomaLogic <- Dat[,c(1:20)]
ProtMatrix_SomaLogic <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation_SomaLogic <- read.csv("Somalogic_Analyte_Annotation_All_UniProt_pre-ANML.csv")
rownames(Analyte_Annotation_SomaLogic) <- Analyte_Annotation_SomaLogic$UniProt

#SomaLogic ANML

#SomaLogic Data
Dat <- read.csv("Somalogic_Merged_All_UniProt.csv")
rownames(Dat) <- Dat$UniqueID
Label_SomaLogic2 <- Dat[,c(1:20)]
ProtMatrix_SomaLogic2 <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation_SomaLogic2 <- read.csv("Somalogic_Analyte_Annotation_All_UniProt.csv")
rownames(Analyte_Annotation_SomaLogic2) <- Analyte_Annotation_SomaLogic2$UniProt

#Olink

#Olink Data
Dat <- read.csv("Olink_Merged_All_UniProt.csv")
rownames(Dat) <- Dat$UniqueID
Label_Olink <- Dat[,c(1:19)]
ProtMatrix_Olink <- Dat[,c(20:ncol(Dat))]
rm(Dat)

#Assay Annotation
Assay_Annotation_Olink <- read.csv("Olink_Assay_Annotation_All_UniProt.csv")
rownames(Assay_Annotation_Olink) <- Assay_Annotation_Olink$UniProt

#Thermo Fisher

#Thermo Fisher Data
Dat <- read.csv("Thermo_Fisher_Merged_All_UniProt_Imputed.csv")
rownames(Dat) <- Dat$biosample_id
Label_TF <- Dat[,c(1:10)]
ProtMatrix_TF <- Dat[,c(11:ncol(Dat))]
rm(Dat)

#Protein Annotation
Protein_Annotation_TF <- read.csv("Thermo_Fisher_Protein_Annotation_All_UniProt.csv")
rownames(Protein_Annotation_TF) <- Protein_Annotation_TF$UniProt

#Common Samples
SampleID_SomaLogic <- paste(Label_SomaLogic$FREG0_PID, Label_SomaLogic$tech_rep, sep = "_")
rownames(Label_SomaLogic) <- SampleID_SomaLogic
rownames(ProtMatrix_SomaLogic) <- SampleID_SomaLogic

SampleID_SomaLogic2 <- paste(Label_SomaLogic2$FREG0_PID, Label_SomaLogic2$tech_rep, sep = "_")
rownames(Label_SomaLogic2) <- SampleID_SomaLogic2
rownames(ProtMatrix_SomaLogic2) <- SampleID_SomaLogic2

SampleID_Olink <- paste(Label_Olink$FREG0_PID, Label_Olink$tech_rep, sep = "_")
rownames(Label_Olink) <- SampleID_Olink
rownames(ProtMatrix_Olink) <- SampleID_Olink

SampleID_TF <- paste(Label_TF$FREG0_PID, Label_TF$tech_rep, sep = "_")
rownames(Label_TF) <- SampleID_TF
rownames(ProtMatrix_TF) <- SampleID_TF

Intersection <- intersect(SampleID_SomaLogic, intersect(SampleID_SomaLogic2, intersect(SampleID_Olink, SampleID_TF)))

Label_SomaLogic <- Label_SomaLogic[Intersection,]
ProtMatrix_SomaLogic <- ProtMatrix_SomaLogic[Intersection,]

Label_SomaLogic2 <- Label_SomaLogic2[Intersection,]
ProtMatrix_SomaLogic2 <- ProtMatrix_SomaLogic2[Intersection,]

Label_Olink <- Label_Olink[Intersection,]
ProtMatrix_Olink <- ProtMatrix_Olink[Intersection,]

Label_TF <- Label_TF[Intersection,]
ProtMatrix_TF <- ProtMatrix_TF[Intersection,]

rm(SampleID_SomaLogic)
rm(SampleID_SomaLogic2)
rm(SampleID_Olink)
rm(SampleID_TF)
rm(Intersection)

#Common Proteins
Intersection <- intersect(colnames(ProtMatrix_SomaLogic), intersect(colnames(ProtMatrix_SomaLogic2), intersect(colnames(ProtMatrix_Olink), colnames(ProtMatrix_TF))))

Analyte_Annotation_SomaLogic <- Analyte_Annotation_SomaLogic[Intersection,]
ProtMatrix_SomaLogic <- ProtMatrix_SomaLogic[,Intersection]

Analyte_Annotation_SomaLogic2 <- Analyte_Annotation_SomaLogic2[Intersection,]
ProtMatrix_SomaLogic2 <- ProtMatrix_SomaLogic2[,Intersection]

Assay_Annotation_Olink <- Assay_Annotation_Olink[Intersection,]
ProtMatrix_Olink <- ProtMatrix_Olink[,Intersection]

Protein_Annotation_TF <- Protein_Annotation_TF[Intersection,]
ProtMatrix_TF <- ProtMatrix_TF[,Intersection]

rm(Intersection)


#CV
CV <- function(x) {
  return(100 * sqrt(exp((log(2) * sd(x))^2) - 1))
}

#ID
ID <- as.data.frame(table(Label_TF$FREG0_PID))
ID <- as.character(ID$Var1[which(ID$Freq == 2)])

#SomaLogic pre-ANML

#Overall CV between Technical Replicates
CV_SomaLogic_Tech <- matrix(NA, nrow = ncol(ProtMatrix_SomaLogic), ncol = 45)
for (i in 1:45) {
  CV_SomaLogic_Tech[,i] <- apply(ProtMatrix_SomaLogic[which(Label_SomaLogic$FREG0_PID == ID[i]),], 2, CV)
}
rm(i)

summary(apply(CV_SomaLogic_Tech, 1, mean))
rm(CV_SomaLogic_Tech)

rm(ProtMatrix_SomaLogic)
rm(Label_SomaLogic)
rm(Analyte_Annotation_SomaLogic)

#SomaLogic ANML

#Overall CV between Technical Replicates
CV_SomaLogic_Tech2 <- matrix(NA, nrow = ncol(ProtMatrix_SomaLogic2), ncol = 45)
for (i in 1:45) {
  CV_SomaLogic_Tech2[,i] <- apply(ProtMatrix_SomaLogic2[which(Label_SomaLogic2$FREG0_PID == ID[i]),], 2, CV)
}
rm(i)

summary(apply(CV_SomaLogic_Tech2, 1, mean))
rm(CV_SomaLogic_Tech2)

rm(ProtMatrix_SomaLogic2)
rm(Label_SomaLogic2)
rm(Analyte_Annotation_SomaLogic2)

#Olink

#Overall CV between Technical Replicates
CV_Olink_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Olink), ncol = 45)
for (i in 1:45) {
  CV_Olink_Tech[,i] <- apply(ProtMatrix_Olink[which(Label_Olink$FREG0_PID == ID[i]),], 2, CV)
}
rm(i)

summary(apply(CV_Olink_Tech, 1, mean))
rm(CV_Olink_Tech)

rm(ProtMatrix_Olink)
rm(Label_Olink)
rm(Assay_Annotation_Olink)

#Thermo

#Overall CV between Technical Replicates
CV_TF_Tech <- matrix(NA, nrow = ncol(ProtMatrix_TF), ncol = 45)
for (i in 1:45) {
  CV_TF_Tech[,i] <- apply(ProtMatrix_TF[which(Label_TF$FREG0_PID == ID[i]),], 2, CV)
}
rm(i)

summary(apply(CV_TF_Tech, 1, mean))
rm(CV_TF_Tech)

rm(ProtMatrix_TF)
rm(Label_TF)
rm(Protein_Annotation_TF)

rm(ID)
rm(CV)
gc()


