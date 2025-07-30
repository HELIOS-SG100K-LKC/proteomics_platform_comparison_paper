#Limit of Detection


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

#Thermo Fisher

#Thermo Fisher Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/5_Common_Samples_Common_Proteins/1_Preprocessing/Thermo_Fisher_Merged_All_UniProt.csv")
rownames(Dat) <- Dat$biosample_id
Label_TF <- Dat[,c(1:10)]
ProtMatrix_TF <- Dat[,c(11:ncol(Dat))]
rm(Dat)

#Protein Annotation
Protein_Annotation_TF <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/5_Common_Samples_Common_Proteins/1_Preprocessing/Thermo_Fisher_Protein_Annotation_All_UniProt.csv")
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


#LoD

#SomaLogic pre-ANML

#SomaLogic Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:20)]
ProtMatrix <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Analyte_Annotation_All.csv")
rownames(Analyte_Annotation) <- Analyte_Annotation$AptName

#Remove Control Samples and Retain Buffer Samples
ProtMatrix_Buffer <- ProtMatrix[which(Label$SampleType == "Buffer"),]
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)

#LOD Calculation

#LOD for Each Protein
LOD <- vector()
for (i in 1:10675) {
  LOD[i] <- median(ProtMatrix_Buffer[,i]) + 4.9 * mad(ProtMatrix_Buffer[,i])
}
rm(i)

#When comparing with LoD, the RFUs were transformed back to original scale.
LOD_TF_SomaLogic <- 2^(ProtMatrix) < matrix(rep(LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF_SomaLogic <- as.data.frame(LOD_TF_SomaLogic)
colnames(LOD_TF_SomaLogic) <- Analyte_Annotation$SeqId

rm(ProtMatrix_Buffer)
rm(ProtMatrix)
rm(Analyte_Annotation)
rm(Label)
rm(LOD)

#SomaLogic ANML

#SomaLogic Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/1_SomaLogic_anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:20)]
ProtMatrix <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/1_SomaLogic_anml_Preprocessing/Somalogic_Analyte_Annotation_All.csv")
rownames(Analyte_Annotation) <- Analyte_Annotation$AptName

#Remove Control Samples and Retain Buffer Samples
ProtMatrix_Buffer <- ProtMatrix[which(Label$SampleType == "Buffer"),]
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)

#LOD Calculation

#LOD for Each Protein
LOD <- vector()
for (i in 1:10675) {
  LOD[i] <- median(ProtMatrix_Buffer[,i]) + 4.9 * mad(ProtMatrix_Buffer[,i])
}
rm(i)

#When comparing with LoD, the RFUs were transformed back to original scale.
LOD_TF_SomaLogic2 <- 2^(ProtMatrix) < matrix(rep(LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF_SomaLogic2 <- as.data.frame(LOD_TF_SomaLogic2)
colnames(LOD_TF_SomaLogic2) <- Analyte_Annotation$SeqId

rm(ProtMatrix_Buffer)
rm(ProtMatrix)
rm(Analyte_Annotation)
rm(Label)
rm(LOD)

#Olink

#Olink Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/3_Olink_Intensity_Preprocessing/Olink_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:19)]
ProtMatrix <- Dat[,c(20:ncol(Dat))]
rm(Dat)

#Assay Annotation
Assay_Annotation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/3_Olink_Intensity_Preprocessing/Olink_Assay_Annotation_All.csv")
rownames(Assay_Annotation) <- Assay_Annotation$OlinkID

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "SAMPLE"),]
Label <- Label[which(Label$SampleType == "SAMPLE"),]

#Limit of Detection (LoD) Analysis

#Package
library(OlinkAnalyze)
library(arrow)

#Data Loading
my_NPX_data <- as.data.frame(read_NPX("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/3_Olink_Intensity_Preprocessing/Q-14739_Chambers_NPX_2024-08-05.parquet"))

#LOD
my_NPX_data <- olink_lod(my_NPX_data, lod_method = "NCLOD")

#Matching and Reshape
my_NPX_data <- my_NPX_data[which(my_NPX_data$OlinkID %in% colnames(ProtMatrix)),]
my_NPX_data <- my_NPX_data[which(my_NPX_data$SampleID %in% Label$SampleID.1),]
library(reshape2)
ProtMatrix_LOD <- dcast(my_NPX_data, SampleID ~ OlinkID, value.var = "LOD")
rownames(ProtMatrix_LOD) <- ProtMatrix_LOD$SampleID
ProtMatrix_LOD <- ProtMatrix_LOD[,-1]
sum(is.na(ProtMatrix_LOD))
ProtMatrix_LOD <- ProtMatrix_LOD[Label$SampleID.1, colnames(ProtMatrix)]
rm(my_NPX_data)
gc()

#Proportion below LOD
LOD_TF_Olink <- ProtMatrix < ProtMatrix_LOD
LOD_TF_Olink <- as.data.frame(LOD_TF_Olink)

rm(ProtMatrix)
rm(Assay_Annotation)
rm(Label)
rm(ProtMatrix_LOD)


#SomaLogic pre-ANML
LOD_TF_SomaLogic <- LOD_TF_SomaLogic[Label_SomaLogic$UniqueID, Analyte_Annotation_SomaLogic$SeqID]
sum(LOD_TF_SomaLogic)/(nrow(LOD_TF_SomaLogic) * ncol(LOD_TF_SomaLogic))

#Sample
Rate_LoD_Sample_SomaLogic <- vector()
for (i in 1:nrow(LOD_TF_SomaLogic)) {
  Rate_LoD_Sample_SomaLogic[i] <- sum(LOD_TF_SomaLogic[i,])/ncol(LOD_TF_SomaLogic)
}
rm(i)
summary(Rate_LoD_Sample_SomaLogic)

#Protein
Rate_LoD_Protein_SomaLogic <- vector()
for (i in 1:ncol(LOD_TF_SomaLogic)) {
  Rate_LoD_Protein_SomaLogic[i] <- sum(LOD_TF_SomaLogic[,i])/nrow(LOD_TF_SomaLogic)
}
rm(i)

rm(LOD_TF_SomaLogic)
rm(Analyte_Annotation_SomaLogic)
rm(Label_SomaLogic)
rm(ProtMatrix_SomaLogic)


#SomaLogic ANML
LOD_TF_SomaLogic2 <- LOD_TF_SomaLogic2[Label_SomaLogic2$UniqueID, Analyte_Annotation_SomaLogic2$SeqID]
sum(LOD_TF_SomaLogic2)/(nrow(LOD_TF_SomaLogic2) * ncol(LOD_TF_SomaLogic2))

#Sample
Rate_LoD_Sample_SomaLogic2 <- vector()
for (i in 1:nrow(LOD_TF_SomaLogic2)) {
  Rate_LoD_Sample_SomaLogic2[i] <- sum(LOD_TF_SomaLogic2[i,])/ncol(LOD_TF_SomaLogic2)
}
rm(i)
summary(Rate_LoD_Sample_SomaLogic2)

#Protein
Rate_LoD_Protein_SomaLogic2 <- vector()
for (i in 1:ncol(LOD_TF_SomaLogic2)) {
  Rate_LoD_Protein_SomaLogic2[i] <- sum(LOD_TF_SomaLogic2[,i])/nrow(LOD_TF_SomaLogic2)
}
rm(i)

rm(LOD_TF_SomaLogic2)
rm(Analyte_Annotation_SomaLogic2)
rm(Label_SomaLogic2)
rm(ProtMatrix_SomaLogic2)


#Olink
LOD_TF_Olink <- LOD_TF_Olink[Label_Olink$UniqueID, Assay_Annotation_Olink$OlinkID]
sum(LOD_TF_Olink)/(nrow(LOD_TF_Olink) * ncol(LOD_TF_Olink))

#Sample
Rate_LoD_Sample_Olink <- vector()
for (i in 1:nrow(LOD_TF_Olink)) {
  Rate_LoD_Sample_Olink[i] <- sum(LOD_TF_Olink[i,])/ncol(LOD_TF_Olink)
}
rm(i)
summary(Rate_LoD_Sample_Olink)

#Protein
Rate_LoD_Protein_Olink <- vector()
for (i in 1:ncol(LOD_TF_Olink)) {
  Rate_LoD_Protein_Olink[i] <- sum(LOD_TF_Olink[,i])/nrow(LOD_TF_Olink)
}
rm(i)

rm(LOD_TF_Olink)
rm(Assay_Annotation_Olink)
rm(Label_Olink)
rm(ProtMatrix_Olink)


#Thermo Fisher
Rate_LoD_Sample_TF <- vector()
for (i in 1:nrow(ProtMatrix_TF)) {
  Rate_LoD_Sample_TF[i] <- sum(is.na(ProtMatrix_TF[i,]))/ncol(ProtMatrix_TF)
}
rm(i)
summary(Rate_LoD_Sample_TF)

#Protein
Rate_LoD_Protein_TF <- vector()
for (i in 1:ncol(ProtMatrix_TF)) {
  Rate_LoD_Protein_TF[i] <- sum(is.na(ProtMatrix_TF[,i]))/nrow(ProtMatrix_TF)
}
rm(i)

rm(Protein_Annotation_TF)
rm(Label_TF)
rm(ProtMatrix_TF)

#Visualization

#Sample
Class <- 1:length(c(Rate_LoD_Sample_SomaLogic, Rate_LoD_Sample_SomaLogic2, Rate_LoD_Sample_Olink, Rate_LoD_Sample_TF))
Class[1:length(Rate_LoD_Sample_SomaLogic)] <- "SomaLogic pre-ANML"
Class[(length(Rate_LoD_Sample_SomaLogic) + 1):(length(Rate_LoD_Sample_SomaLogic) + length(Rate_LoD_Sample_SomaLogic2))] <- "SomaLogic ANML"
Class[(length(Rate_LoD_Sample_SomaLogic) + length(Rate_LoD_Sample_SomaLogic2) + 1):(length(Rate_LoD_Sample_SomaLogic) + length(Rate_LoD_Sample_SomaLogic2) + length(Rate_LoD_Sample_Olink))] <- "Olink"
Class[(length(Rate_LoD_Sample_SomaLogic) + length(Rate_LoD_Sample_SomaLogic2) + length(Rate_LoD_Sample_Olink) + 1):(length(Rate_LoD_Sample_SomaLogic) + length(Rate_LoD_Sample_SomaLogic2) + length(Rate_LoD_Sample_Olink) + length(Rate_LoD_Sample_TF))] <- "Thermo"
Data <- data.frame(Rate_LoD_Sample = c(Rate_LoD_Sample_SomaLogic, Rate_LoD_Sample_SomaLogic2, Rate_LoD_Sample_Olink, Rate_LoD_Sample_TF), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Rate_LoD_Sample, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("SomaLogic pre-ANML" = "#0072B5FF", "SomaLogic ANML" = "#E18727FF", "Olink" = "#20854EFF", "Thermo" = "#BC3C29FF")) +
  scale_color_manual(values = c("SomaLogic pre-ANML" = "#0072B5FF", "SomaLogic ANML" = "#E18727FF", "Olink" = "#20854EFF", "Thermo" = "#BC3C29FF")) +
  scale_x_discrete(limits = c("SomaLogic pre-ANML","SomaLogic ANML","Olink", "Thermo")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  labs(y = "Sample missingness", x = "Platform") +
  ggtitle("Sample missingness")
rm(Class)
rm(Data)
rm(Rate_LoD_Sample_SomaLogic)
rm(Rate_LoD_Sample_SomaLogic2)
rm(Rate_LoD_Sample_Olink)
rm(Rate_LoD_Sample_TF)

#Protein
Rate_SomaLogic <- vector()
Rate_SomaLogic[1] <- length(which(Rate_LoD_Protein_SomaLogic == 0))
Rate_SomaLogic[2] <- length(which(Rate_LoD_Protein_SomaLogic > 0.2))
Rate_SomaLogic[3] <- length(which(Rate_LoD_Protein_SomaLogic > 0.5))
Rate_SomaLogic[4] <- length(which(Rate_LoD_Protein_SomaLogic == 1))

Rate_SomaLogic2 <- vector()
Rate_SomaLogic2[1] <- length(which(Rate_LoD_Protein_SomaLogic2 == 0))
Rate_SomaLogic2[2] <- length(which(Rate_LoD_Protein_SomaLogic2 > 0.2))
Rate_SomaLogic2[3] <- length(which(Rate_LoD_Protein_SomaLogic2 > 0.5))
Rate_SomaLogic2[4] <- length(which(Rate_LoD_Protein_SomaLogic2 == 1))

Rate_Olink <- vector()
Rate_Olink[1] <- length(which(Rate_LoD_Protein_Olink == 0))
Rate_Olink[2] <- length(which(Rate_LoD_Protein_Olink > 0.2))
Rate_Olink[3] <- length(which(Rate_LoD_Protein_Olink > 0.5))
Rate_Olink[4] <- length(which(Rate_LoD_Protein_Olink == 1))

Rate_Thermo <- vector()
Rate_Thermo[1] <- length(which(Rate_LoD_Protein_TF == 0))
Rate_Thermo[2] <- length(which(Rate_LoD_Protein_TF > 0.2))
Rate_Thermo[3] <- length(which(Rate_LoD_Protein_TF > 0.5))
Rate_Thermo[4] <- length(which(Rate_LoD_Protein_TF == 1))

Data <- data.frame(Rate = c(Rate_SomaLogic, Rate_SomaLogic2, Rate_Olink, Rate_Thermo), Threshold = rep(c("1.=0%", "2.>20%", "3.>50%", "4.=100%"), 4), Platform = rep(c("SomaLogic pre-ANML", "SomaLogic ANML", "Olink", "Thermo Fisher"), each = 4))
ggplot(Data, aes(x = Platform, y = Rate, fill = Threshold, color = Threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(limits = c("SomaLogic pre-ANML", "SomaLogic ANML", "Olink", "Thermo Fisher")) +
  scale_fill_manual(values = c("1.=0%" = "#7876B1FF", "2.>20%" = "#6F99ADFF", "3.>50%" = "#FFDC91FF", "4.=100%" = "#EE4C97FF")) +
  scale_color_manual(values = c("1.=0%" = "#7876B1FF", "2.>20%" = "#6F99ADFF", "3.>50%" = "#FFDC91FF", "4.=100%" = "#EE4C97FF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  labs(x = "Platform", y = "Number of proteins") +
  geom_text(aes(label = Rate), color = "black", position = position_dodge(width = 1), vjust = -0.5, size = 5.5) +
  ggtitle("Protein missingness")
rm(Data)
rm(Rate_SomaLogic)
rm(Rate_SomaLogic2)
rm(Rate_Olink)
rm(Rate_Thermo)

rm(Rate_LoD_Protein_SomaLogic)
rm(Rate_LoD_Protein_SomaLogic2)
rm(Rate_LoD_Protein_Olink)
rm(Rate_LoD_Protein_TF)
gc()

