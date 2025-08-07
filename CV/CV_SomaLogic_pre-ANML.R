#CV SomaLogic pre-ANML


#Data Loading

#SomaLogic Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:20)]
ProtMatrix <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Analyte_Annotation_All.csv")
rownames(Analyte_Annotation) <- Analyte_Annotation$AptName

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)

#CV
CV <- function(x) {
  return(100 * sqrt(exp((log(2) * sd(x))^2) - 1))
}


#Overall Intra-Plate CV
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

CV_Intra <- apply(CV_SomaLogic_Tech, 1, mean)
summary(CV_Intra)
round(quantile(CV_Intra, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Overall Inter-Plate CV
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

CV_Inter <- apply(CV_SomaLogic_Tech, 1, mean)
summary(CV_Inter)
round(quantile(CV_Inter, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_SomaLogic_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)


#Stratify by LoD
LOD_SomaLogic <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/3_Limit_of_Detection/4_SomaLogic_pre-ANML_LoD/SomaLogic_LOD.csv")
identical(LOD_SomaLogic$SeqID, Analyte_Annotation$SeqId)

#Intra

#Above
median(CV_Intra[which(LOD_SomaLogic$PropBelowLOD <= 0.2)])

#Below
median(CV_Intra[which(LOD_SomaLogic$PropBelowLOD > 0.2)])

#Inter

#Above
median(CV_Inter[which(LOD_SomaLogic$PropBelowLOD <= 0.2)])

#Below
median(CV_Inter[which(LOD_SomaLogic$PropBelowLOD > 0.2)])

rm(LOD_SomaLogic)


#Different Sets

#Set 1: SomaScan 5K
library(readxl)
ProteinList5K <- read_excel("SomaScan_Menu_1.3K_5K_7K.xlsx")
ProteinList5K <- ProteinList5K[-c(1:2),10]
ProteinList5K <- unique(ProteinList5K$...10)

#Set2 : SomaScan 7K - SomaScan 5K
ProteinList7K <- read_excel("SomaScan_Menu_1.3K_5K_7K.xlsx")
ProteinList7K <- ProteinList7K[-c(1:2),1]
ProteinList7K <- unique(ProteinList7K$`Supplementary Data 1. List of all SOMAmers in the plasma 7k, 5k, and 1.3k SomaScan assays.`)
ProteinList2 <- setdiff(ProteinList7K, ProteinList5K)

#Set 3: SomaScan 11K - SomaScan 7K
ProteinList3 <- setdiff(Analyte_Annotation$SeqId, ProteinList7K)
rm(ProteinList7K)

#Intra

#Set 1
median(CV_Intra[which(Analyte_Annotation$SeqId %in% ProteinList5K)])

#Set 2
median(CV_Intra[which(Analyte_Annotation$SeqId %in% ProteinList2)])

#Set 3
median(CV_Intra[which(Analyte_Annotation$SeqId %in% ProteinList3)])

#Inter

#Set 1
median(CV_Inter[which(Analyte_Annotation$SeqId %in% ProteinList5K)])

#Set 2
median(CV_Inter[which(Analyte_Annotation$SeqId %in% ProteinList2)])

#Set 3
median(CV_Inter[which(Analyte_Annotation$SeqId %in% ProteinList3)])

rm(ProteinList5K)
rm(ProteinList2)
rm(ProteinList3)


#Different Dilution Groups

#Intra

#1:5
median(CV_Intra[which(Analyte_Annotation$Dilution2 == 1/5)])

#1:200
median(CV_Intra[which(Analyte_Annotation$Dilution2 == 1/200)])

#1:20000
median(CV_Intra[which(Analyte_Annotation$Dilution2 == 1/20000)])

#Inter

#1:5
median(CV_Inter[which(Analyte_Annotation$Dilution2 == 1/5)])

#1:200
median(CV_Inter[which(Analyte_Annotation$Dilution2 == 1/200)])

#1:20000
median(CV_Inter[which(Analyte_Annotation$Dilution2 == 1/20000)])

rm(CV_Intra, CV_Inter)


rm(CV)
rm(ProtMatrix)
rm(Label)
rm(Analyte_Annotation)
gc()

