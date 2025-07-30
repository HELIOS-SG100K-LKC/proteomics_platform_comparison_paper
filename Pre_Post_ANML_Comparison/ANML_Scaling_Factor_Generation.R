#ANML Scaling Factor Generation

#Pre-ANML

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

#ANML

#SomaLogic Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/1_SomaLogic_anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label2 <- Dat[,c(1:20)]
ProtMatrix2 <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation2 <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/1_SomaLogic_anml_Preprocessing/Somalogic_Analyte_Annotation_All.csv")
rownames(Analyte_Annotation2) <- Analyte_Annotation2$AptName

#Remove Control Samples
ProtMatrix2 <- ProtMatrix2[which(Label2$SampleType == "Sample"),]
Label2 <- Label2[which(Label2$SampleType == "Sample"),]

#Mapping
identical(rownames(Label), rownames(Label2))
identical(colnames(Label), colnames(Label2))

identical(rownames(Analyte_Annotation), rownames(Analyte_Annotation2))
identical(colnames(Analyte_Annotation), colnames(Analyte_Annotation2))

identical(rownames(ProtMatrix), rownames(ProtMatrix2))
identical(colnames(ProtMatrix), colnames(ProtMatrix2))

#Scaling Factor Calculation
identical(which(Analyte_Annotation$Dilution2 == 0.2), which(Analyte_Annotation2$Dilution2 == 0.2))
Scaling_Factor1 <- ProtMatrix2[,which(Analyte_Annotation2$Dilution2 == 0.2)]/ProtMatrix[,which(Analyte_Annotation$Dilution2 == 0.2)]
Scaling_Factor_Median1 <- apply(Scaling_Factor1, 1, median)

identical(which(Analyte_Annotation$Dilution2 == 0.005), which(Analyte_Annotation2$Dilution2 == 0.005))
Scaling_Factor2 <- ProtMatrix2[,which(Analyte_Annotation2$Dilution2 == 0.005)]/ProtMatrix[,which(Analyte_Annotation$Dilution2 == 0.005)]
Scaling_Factor_Median2 <- apply(Scaling_Factor2, 1, median)

identical(which(Analyte_Annotation$Dilution2 == 0.00005), which(Analyte_Annotation2$Dilution2 == 0.00005))
Scaling_Factor3 <- ProtMatrix2[,which(Analyte_Annotation2$Dilution2 == 0.00005)]/ProtMatrix[,which(Analyte_Annotation$Dilution2 == 0.00005)]
Scaling_Factor_Median3 <- apply(Scaling_Factor3, 1, median)

#Median
Scaling_Factor <- cbind(Scaling_Factor_Median1, Scaling_Factor_Median2, Scaling_Factor_Median3)
colnames(Scaling_Factor) <- c("ANML_1_5", "ANML_1_200", "ANML_1_20000")
Scaling_Factor <- cbind(Label[,c("SOMA_ID","FREG0_PID","FREG14_Visit_number","tech_rep_id","bio_rep_id","tech_rep","UniqueID","SampleId")], Scaling_Factor)

write.csv(Scaling_Factor, "SomaLogic_ANML_Scaling_Factor.csv", row.names = FALSE)

rm(Scaling_Factor)
rm(Scaling_Factor1)
rm(Scaling_Factor2)
rm(Scaling_Factor3)
rm(Scaling_Factor_Median1)
rm(Scaling_Factor_Median2)
rm(Scaling_Factor_Median3)

rm(ProtMatrix)
rm(ProtMatrix2)
rm(Label)
rm(Label2)
rm(Analyte_Annotation)
rm(Analyte_Annotation2)
gc()
