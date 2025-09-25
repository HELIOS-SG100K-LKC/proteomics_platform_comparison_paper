#Longitudinal CV SNR


#CV
CV <- function(x) {
  return(100 * sqrt(exp((log(2) * sd(x))^2) - 1))
}


#SomaLogic ANML

#SomaLogic Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/1_SomaLogic_anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:20)]
ProtMatrix <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)


#Intra-Plate Technical CV
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep %in% c(1:2)),]
Label_Tech <- Label[which(Label$tech_rep %in% c(1:2)),]
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
ProtMatrix_Bio <- ProtMatrix[union(which(Label$tech_rep == 3), which(Label$bio_rep_id == "Y")),]
Label_Bio <- Label[union(which(Label$tech_rep == 3), which(Label$bio_rep_id == "Y")),]
ID <- unique(Label_Bio$FREG0_PID)

CV_SomaLogic_Bio <- matrix(NA, nrow = ncol(ProtMatrix_Bio), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Bio[,i] <- apply(t(ProtMatrix_Bio[which(Label_Bio$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Bio <- as.data.frame(CV_SomaLogic_Bio)
rownames(CV_SomaLogic_Bio) <- colnames(ProtMatrix_Bio)

CV_Intra_Longi_ANML <- apply(CV_SomaLogic_Bio, 1, mean)
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


#Inter-Plate Technical CV
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep %in% c(2:3)),]
Label_Tech <- Label[which(Label$tech_rep %in% c(2:3)),]
ID <- unique(Label_Tech$FREG0_PID)

CV_SomaLogic_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Tech), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Tech[,i] <- apply(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Tech <- as.data.frame(CV_SomaLogic_Tech)
rownames(CV_SomaLogic_Tech) <- colnames(ProtMatrix_Tech)

CV_Inter_ANML <- apply(CV_SomaLogic_Tech, 1, mean)
round(quantile(CV_Inter_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Inter-plate Longitudinal CV
ProtMatrix_Bio <- ProtMatrix[union(which(Label$tech_rep == 2), which(Label$bio_rep_id == "Y")),]
Label_Bio <- Label[union(which(Label$tech_rep == 2), which(Label$bio_rep_id == "Y")),]
ID <- unique(Label_Bio$FREG0_PID)

CV_SomaLogic_Bio <- matrix(NA, nrow = ncol(ProtMatrix_Bio), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Bio[,i] <- apply(t(ProtMatrix_Bio[which(Label_Bio$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Bio <- as.data.frame(CV_SomaLogic_Bio)
rownames(CV_SomaLogic_Bio) <- colnames(ProtMatrix_Bio)

CV_Inter_Longi_ANML <- apply(CV_SomaLogic_Bio, 1, mean)
round(quantile(CV_Inter_Longi_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Bio)
rm(ID)
rm(Label_Bio)
rm(ProtMatrix_Bio)

#Inter-plate SNR
SNR_Inter_ANML <- CV_Inter_Longi_ANML / CV_Inter_ANML
round(quantile(SNR_Inter_ANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Inter_ANML)
rm(CV_Inter_Longi_ANML)
rm(SNR_Inter_ANML)


rm(ProtMatrix)
rm(Label)
gc()


#SomaLogic pre-ANML

#SomaLogic Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:20)]
ProtMatrix <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)


#Intra-Plate Technical CV
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep %in% c(1:2)),]
Label_Tech <- Label[which(Label$tech_rep %in% c(1:2)),]
ID <- unique(Label_Tech$FREG0_PID)

CV_SomaLogic_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Tech), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Tech[,i] <- apply(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Tech <- as.data.frame(CV_SomaLogic_Tech)
rownames(CV_SomaLogic_Tech) <- colnames(ProtMatrix_Tech)

CV_Intra_preANML <- apply(CV_SomaLogic_Tech, 1, mean)
round(quantile(CV_Intra_preANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Intra-plate Longitudinal CV
ProtMatrix_Bio <- ProtMatrix[union(which(Label$tech_rep == 3), which(Label$bio_rep_id == "Y")),]
Label_Bio <- Label[union(which(Label$tech_rep == 3), which(Label$bio_rep_id == "Y")),]
ID <- unique(Label_Bio$FREG0_PID)

CV_SomaLogic_Bio <- matrix(NA, nrow = ncol(ProtMatrix_Bio), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Bio[,i] <- apply(t(ProtMatrix_Bio[which(Label_Bio$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Bio <- as.data.frame(CV_SomaLogic_Bio)
rownames(CV_SomaLogic_Bio) <- colnames(ProtMatrix_Bio)

CV_Intra_Longi_preANML <- apply(CV_SomaLogic_Bio, 1, mean)
round(quantile(CV_Intra_Longi_preANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Bio)
rm(ID)
rm(Label_Bio)
rm(ProtMatrix_Bio)

#Intra-plate SNR
SNR_Intra_preANML <- CV_Intra_Longi_preANML / CV_Intra_preANML
round(quantile(SNR_Intra_preANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Intra_preANML)
rm(CV_Intra_Longi_preANML)
rm(SNR_Intra_preANML)


#Inter-Plate Technical CV
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep %in% c(2:3)),]
Label_Tech <- Label[which(Label$tech_rep %in% c(2:3)),]
ID <- unique(Label_Tech$FREG0_PID)

CV_SomaLogic_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Tech), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Tech[,i] <- apply(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Tech <- as.data.frame(CV_SomaLogic_Tech)
rownames(CV_SomaLogic_Tech) <- colnames(ProtMatrix_Tech)

CV_Inter_preANML <- apply(CV_SomaLogic_Tech, 1, mean)
round(quantile(CV_Inter_preANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Inter-plate Longitudinal CV
ProtMatrix_Bio <- ProtMatrix[union(which(Label$tech_rep == 2), which(Label$bio_rep_id == "Y")),]
Label_Bio <- Label[union(which(Label$tech_rep == 2), which(Label$bio_rep_id == "Y")),]
ID <- unique(Label_Bio$FREG0_PID)

CV_SomaLogic_Bio <- matrix(NA, nrow = ncol(ProtMatrix_Bio), ncol = length(ID))
for (i in 1:100) {
  CV_SomaLogic_Bio[,i] <- apply(t(ProtMatrix_Bio[which(Label_Bio$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_SomaLogic_Bio <- as.data.frame(CV_SomaLogic_Bio)
rownames(CV_SomaLogic_Bio) <- colnames(ProtMatrix_Bio)

CV_Inter_Longi_preANML <- apply(CV_SomaLogic_Bio, 1, mean)
round(quantile(CV_Inter_Longi_preANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Bio)
rm(ID)
rm(Label_Bio)
rm(ProtMatrix_Bio)

#Inter-plate SNR
SNR_Inter_preANML <- CV_Inter_Longi_preANML / CV_Inter_preANML
round(quantile(SNR_Inter_preANML, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Inter_preANML)
rm(CV_Inter_Longi_preANML)
rm(SNR_Inter_preANML)


rm(ProtMatrix)
rm(Label)
gc()


#Olink

#Olink Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/3_Olink_Intensity_Preprocessing/Olink_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:19)]
ProtMatrix <- Dat[,c(20:ncol(Dat))]
rm(Dat)

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "SAMPLE"),]
Label <- Label[which(Label$SampleType == "SAMPLE"),]


#Intra-Plate Technical CV
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep %in% c(1:2)),]
Label_Tech <- Label[which(Label$tech_rep %in% c(1:2)),]
ID <- unique(Label_Tech$FREG0_PID)

CV_Olink_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Tech), ncol = length(ID))
for (i in 1:100) {
  CV_Olink_Tech[,i] <- apply(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_Olink_Tech <- as.data.frame(CV_Olink_Tech)
rownames(CV_Olink_Tech) <- colnames(ProtMatrix_Tech)

CV_Intra_Olink <- apply(CV_Olink_Tech, 1, mean)
round(quantile(CV_Intra_Olink, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Olink_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Intra-Plate Longitudinal CV
ProtMatrix_Bio <- ProtMatrix[union(which(Label$tech_rep == 3), which(Label$bio_rep_id == "Y")),]
Label_Bio <- Label[union(which(Label$tech_rep == 3), which(Label$bio_rep_id == "Y")),]
ID <- unique(Label_Bio$FREG0_PID)

CV_Olink_Bio <- matrix(NA, nrow = ncol(ProtMatrix_Bio), ncol = length(ID))
for (i in 1:100) {
  CV_Olink_Bio[,i] <- apply(t(ProtMatrix_Bio[which(Label_Bio$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_Olink_Bio <- as.data.frame(CV_Olink_Bio)
rownames(CV_Olink_Bio) <- colnames(ProtMatrix_Bio)

CV_Intra_Longi_Olink <- apply(CV_Olink_Bio, 1, mean, na.rm = T)
round(quantile(CV_Intra_Longi_Olink, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Olink_Bio)
rm(ID)
rm(Label_Bio)
rm(ProtMatrix_Bio)

#Intra-plate SNR
SNR_Intra_Olink <- CV_Intra_Longi_Olink / CV_Intra_Olink
round(quantile(SNR_Intra_Olink, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Intra_Olink)
rm(CV_Intra_Longi_Olink)
rm(SNR_Intra_Olink)


#Inter-Plate Technical CV
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep %in% c(2:3)),]
Label_Tech <- Label[which(Label$tech_rep %in% c(2:3)),]
ID <- unique(Label_Tech$FREG0_PID)

CV_Olink_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Tech), ncol = length(ID))
for (i in 1:100) {
  CV_Olink_Tech[,i] <- apply(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_Olink_Tech <- as.data.frame(CV_Olink_Tech)
rownames(CV_Olink_Tech) <- colnames(ProtMatrix_Tech)

CV_Inter_Olink <- apply(CV_Olink_Tech, 1, mean, na.rm = T)
round(quantile(CV_Inter_Olink, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Olink_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Inter-Plate Longitudinal CV
ProtMatrix_Bio <- ProtMatrix[union(which(Label$tech_rep == 2), which(Label$bio_rep_id == "Y")),]
Label_Bio <- Label[union(which(Label$tech_rep == 2), which(Label$bio_rep_id == "Y")),]
ID <- unique(Label_Bio$FREG0_PID)

CV_Olink_Bio <- matrix(NA, nrow = ncol(ProtMatrix_Bio), ncol = length(ID))
for (i in 1:100) {
  CV_Olink_Bio[,i] <- apply(t(ProtMatrix_Bio[which(Label_Bio$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_Olink_Bio <- as.data.frame(CV_Olink_Bio)
rownames(CV_Olink_Bio) <- colnames(ProtMatrix_Bio)

CV_Inter_Longi_Olink <- apply(CV_Olink_Bio, 1, mean, na.rm = T)
round(quantile(CV_Inter_Longi_Olink, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Olink_Bio)
rm(ID)
rm(Label_Bio)
rm(ProtMatrix_Bio)

#Inter-plate SNR
SNR_Inter_Olink <- CV_Inter_Longi_Olink / CV_Inter_Olink
round(quantile(SNR_Inter_Olink, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Inter_Olink)
rm(CV_Inter_Longi_Olink)
rm(SNR_Inter_Olink)


rm(ProtMatrix)
rm(Label)


rm(CV)
gc()

