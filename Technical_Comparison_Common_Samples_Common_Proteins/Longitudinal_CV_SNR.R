#Longitudinal CV SNR


#SomaLogic pre-ANML

#SomaLogic Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/5_Common_Samples_Common_Proteins/1_Preprocessing/Somalogic_Merged_All_UniProt_pre-ANML.csv")
rownames(Dat) <- Dat$UniqueID
Label_SomaLogic <- Dat[,c(1:20)]
ProtMatrix_SomaLogic <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation_SomaLogic <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/5_Common_Samples_Common_Proteins/1_Preprocessing/Somalogic_Analyte_Annotation_All_UniProt_pre-ANML.csv")
rownames(Analyte_Annotation_SomaLogic) <- Analyte_Annotation_SomaLogic$UniProt

#SomaLogic ANML

#SomaLogic Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/5_Common_Samples_Common_Proteins/1_Preprocessing/Somalogic_Merged_All_UniProt.csv")
rownames(Dat) <- Dat$UniqueID
Label_SomaLogic2 <- Dat[,c(1:20)]
ProtMatrix_SomaLogic2 <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation_SomaLogic2 <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/5_Common_Samples_Common_Proteins/1_Preprocessing/Somalogic_Analyte_Annotation_All_UniProt.csv")
rownames(Analyte_Annotation_SomaLogic2) <- Analyte_Annotation_SomaLogic2$UniProt

#Olink

#Olink Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/5_Common_Samples_Common_Proteins/1_Preprocessing/Olink_Merged_All_UniProt.csv")
rownames(Dat) <- Dat$UniqueID
Label_Olink <- Dat[,c(1:19)]
ProtMatrix_Olink <- Dat[,c(20:ncol(Dat))]
rm(Dat)

#Assay Annotation
Assay_Annotation_Olink <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/5_Common_Samples_Common_Proteins/1_Preprocessing/Olink_Assay_Annotation_All_UniProt.csv")
rownames(Assay_Annotation_Olink) <- Assay_Annotation_Olink$UniProt

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

Intersection <- intersect(SampleID_SomaLogic, intersect(SampleID_SomaLogic2, SampleID_Olink))

Label_SomaLogic <- Label_SomaLogic[Intersection,]
ProtMatrix_SomaLogic <- ProtMatrix_SomaLogic[Intersection,]

Label_SomaLogic2 <- Label_SomaLogic2[Intersection,]
ProtMatrix_SomaLogic2 <- ProtMatrix_SomaLogic2[Intersection,]

Label_Olink <- Label_Olink[Intersection,]
ProtMatrix_Olink <- ProtMatrix_Olink[Intersection,]

rm(SampleID_SomaLogic)
rm(SampleID_SomaLogic2)
rm(SampleID_Olink)
rm(Intersection)

#Common Proteins
Intersection <- intersect(colnames(ProtMatrix_SomaLogic), intersect(colnames(ProtMatrix_SomaLogic2), colnames(ProtMatrix_Olink)))

Analyte_Annotation_SomaLogic <- Analyte_Annotation_SomaLogic[Intersection,]
ProtMatrix_SomaLogic <- ProtMatrix_SomaLogic[,Intersection]

Analyte_Annotation_SomaLogic2 <- Analyte_Annotation_SomaLogic2[Intersection,]
ProtMatrix_SomaLogic2 <- ProtMatrix_SomaLogic2[,Intersection]

Assay_Annotation_Olink <- Assay_Annotation_Olink[Intersection,]
ProtMatrix_Olink <- ProtMatrix_Olink[,Intersection]

rm(Intersection)


#CV
CV <- function(x) {
  return(100 * sqrt(exp((log(2) * sd(x))^2) - 1))
}


#SomaLogic pre-ANML

#Intra-Plate Technical CV
ProtMatrix_Tech <- ProtMatrix_SomaLogic[which(Label_SomaLogic$tech_rep %in% c(1:2)),]
Label_Tech <- Label_SomaLogic[which(Label_SomaLogic$tech_rep %in% c(1:2)),]
ID <- unique(Label_Tech$FREG0_PID)

CV_SomaLogic_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Tech), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Tech[,i] <- apply(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Tech <- as.data.frame(CV_SomaLogic_Tech)
rownames(CV_SomaLogic_Tech) <- colnames(ProtMatrix_Tech)

CV_Intra_ANML <- apply(CV_SomaLogic_Tech, 1, mean)
round(quantile(CV_Intra_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Intra-plate Longitudinal CV
ProtMatrix_Bio <- ProtMatrix_SomaLogic[union(which(Label_SomaLogic$tech_rep == 3), which(Label_SomaLogic$bio_rep_id == "Y")),]
Label_Bio <- Label_SomaLogic[union(which(Label_SomaLogic$tech_rep == 3), which(Label_SomaLogic$bio_rep_id == "Y")),]
ID <- unique(Label_Bio$FREG0_PID)

CV_SomaLogic_Bio <- matrix(NA, nrow = ncol(ProtMatrix_Bio), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Bio[,i] <- apply(t(ProtMatrix_Bio[which(Label_Bio$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Bio <- as.data.frame(CV_SomaLogic_Bio)
rownames(CV_SomaLogic_Bio) <- colnames(ProtMatrix_Bio)

CV_Intra_Longi_ANML <- apply(CV_SomaLogic_Bio, 1, mean, na.rm = TRUE)
round(quantile(CV_Intra_Longi_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Bio)
rm(ID)
rm(Label_Bio)
rm(ProtMatrix_Bio)

#Intra-plate SNR
SNR_Intra_ANML <- CV_Intra_Longi_ANML / CV_Intra_ANML
round(quantile(SNR_Intra_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Intra_ANML)
rm(CV_Intra_Longi_ANML)
rm(SNR_Intra_ANML)

rm(ProtMatrix_SomaLogic)
rm(Label_SomaLogic)
rm(Analyte_Annotation_SomaLogic)


#SomaLogic ANML

#Intra-Plate Technical CV
ProtMatrix_Tech <- ProtMatrix_SomaLogic2[which(Label_SomaLogic2$tech_rep %in% c(1:2)),]
Label_Tech <- Label_SomaLogic2[which(Label_SomaLogic2$tech_rep %in% c(1:2)),]
ID <- unique(Label_Tech$FREG0_PID)

CV_SomaLogic_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Tech), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Tech[,i] <- apply(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Tech <- as.data.frame(CV_SomaLogic_Tech)
rownames(CV_SomaLogic_Tech) <- colnames(ProtMatrix_Tech)

CV_Intra_ANML <- apply(CV_SomaLogic_Tech, 1, mean)
round(quantile(CV_Intra_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Intra-plate Longitudinal CV
ProtMatrix_Bio <- ProtMatrix_SomaLogic2[union(which(Label_SomaLogic2$tech_rep == 3), which(Label_SomaLogic2$bio_rep_id == "Y")),]
Label_Bio <- Label_SomaLogic2[union(which(Label_SomaLogic2$tech_rep == 3), which(Label_SomaLogic2$bio_rep_id == "Y")),]
ID <- unique(Label_Bio$FREG0_PID)

CV_SomaLogic_Bio <- matrix(NA, nrow = ncol(ProtMatrix_Bio), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Bio[,i] <- apply(t(ProtMatrix_Bio[which(Label_Bio$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Bio <- as.data.frame(CV_SomaLogic_Bio)
rownames(CV_SomaLogic_Bio) <- colnames(ProtMatrix_Bio)

CV_Intra_Longi_ANML <- apply(CV_SomaLogic_Bio, 1, mean, na.rm = TRUE)
round(quantile(CV_Intra_Longi_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Bio)
rm(ID)
rm(Label_Bio)
rm(ProtMatrix_Bio)

#Intra-plate SNR
SNR_Intra_ANML <- CV_Intra_Longi_ANML / CV_Intra_ANML
round(quantile(SNR_Intra_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Intra_ANML)
rm(CV_Intra_Longi_ANML)
rm(SNR_Intra_ANML)

rm(ProtMatrix_SomaLogic2)
rm(Label_SomaLogic2)
rm(Analyte_Annotation_SomaLogic2)


#Olink

#Intra-Plate Technical CV
ProtMatrix_Tech <- ProtMatrix_Olink[which(Label_Olink$tech_rep %in% c(1:2)),]
Label_Tech <- Label_Olink[which(Label_Olink$tech_rep %in% c(1:2)),]
ID <- unique(Label_Tech$FREG0_PID)

CV_Olink_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Tech), ncol = length(ID))
for (i in 1:100) {
  CV_Olink_Tech[,i] <- apply(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_Olink_Tech <- as.data.frame(CV_Olink_Tech)
rownames(CV_Olink_Tech) <- colnames(ProtMatrix_Tech)

CV_Intra_ANML <- apply(CV_Olink_Tech, 1, mean)
round(quantile(CV_Intra_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Olink_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Intra-plate Longitudinal CV
ProtMatrix_Bio <- ProtMatrix_Olink[union(which(Label_Olink$tech_rep == 3), which(Label_Olink$bio_rep_id == "Y")),]
Label_Bio <- Label_Olink[union(which(Label_Olink$tech_rep == 3), which(Label_Olink$bio_rep_id == "Y")),]
ID <- unique(Label_Bio$FREG0_PID)

CV_Olink_Bio <- matrix(NA, nrow = ncol(ProtMatrix_Bio), ncol = length(ID))
for (i in 1:100) {
  CV_Olink_Bio[,i] <- apply(t(ProtMatrix_Bio[which(Label_Bio$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_Olink_Bio <- as.data.frame(CV_Olink_Bio)
rownames(CV_Olink_Bio) <- colnames(ProtMatrix_Bio)

CV_Intra_Longi_ANML <- apply(CV_Olink_Bio, 1, mean, na.rm = TRUE)
round(quantile(CV_Intra_Longi_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Olink_Bio)
rm(ID)
rm(Label_Bio)
rm(ProtMatrix_Bio)

#Intra-plate SNR
SNR_Intra_ANML <- CV_Intra_Longi_ANML / CV_Intra_ANML
round(quantile(SNR_Intra_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Intra_ANML)
rm(CV_Intra_Longi_ANML)
rm(SNR_Intra_ANML)

rm(ProtMatrix_Olink)
rm(Label_Olink)
rm(Assay_Annotation_Olink)


rm(CV)
gc()

