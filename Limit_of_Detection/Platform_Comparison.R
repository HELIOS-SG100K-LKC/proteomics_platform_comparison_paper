#Limit of Detection Comparison among SomaLogic pre-ANML, SomaLogic ANML, Olink and Thermo Fisher


#Results
Proportion_Unified_Formula <- vector()
Overall_Proportion <- vector()


#SomaLogic pre-ANML

#SomaLogic Data
Dat <- read.csv("SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:20)]
ProtMatrix <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation <- read.csv("SomaLogic_pre-anml_Preprocessing/Somalogic_Analyte_Annotation_All.csv")
rownames(Analyte_Annotation) <- Analyte_Annotation$AptName

#Remove Control Samples and Retain Buffer Samples
ProtMatrix_Buffer <- ProtMatrix[which(Label$SampleType == "Buffer"),]
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)

#500
Index <- union(which(Label$tech_rep %in% c(2:3)), which(Label$bio_rep_id == "Y"))
ProtMatrix <- ProtMatrix[-Index,]
Label <- Label[-Index,]
rm(Index)

#LOD Calculation

#LOD for Each Protein
LOD <- vector()
for (i in 1:10675) {
  LOD[i] <- median(ProtMatrix_Buffer[,i]) + 4.9 * mad(ProtMatrix_Buffer[,i])
}
rm(i)

#LoD
LOD <- data.frame(SeqID = Analyte_Annotation$SeqId, LOD)
rownames(LOD) <- LOD$SeqID

#Proportion below LOD
rownames(Analyte_Annotation) <- Analyte_Annotation$SeqId
colnames(ProtMatrix) <- Analyte_Annotation$SeqId
colnames(ProtMatrix_Buffer) <- Analyte_Annotation$SeqId

#When comparing with LoD, the RFUs were transformed back to original scale.
LOD_TF <- 2^(ProtMatrix) < matrix(rep(LOD$LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
(Overall_Proportion[1] <- sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF)))

#Proportion below LoD for each sample
Rate_LoD_Sample_SomaLogic <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample_SomaLogic[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample_SomaLogic)

#Proportion below LoD for each protein
Rate_LoD_Protein_SomaLogic <- vector()
for (i in 1:ncol(LOD_TF)) {
  Rate_LoD_Protein_SomaLogic[i] <- sum(LOD_TF[,i])/nrow(LOD_TF)
}
rm(i)
rm(LOD_TF)
rm(LOD)

#Formula 1

#LOD for Each Protein
LOD <- vector()
for (i in 1:10675) {
  LOD[i] <- mean(ProtMatrix_Buffer[,i]) + 3.3 * sd(ProtMatrix_Buffer[,i])
}
rm(i)

#When comparing with LoD, the RFUs were transformed back to original scale.
LOD_TF <- 2^(ProtMatrix) < matrix(rep(LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
(Proportion_Unified_Formula[1] <- sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF)))

rm(LOD_TF)
rm(LOD)

#Formula 2

#LOD for Each Protein
ProtMatrix <- 2^ProtMatrix
Quantile <- vector()
for (i in 1:10675) {
  Quantile[i] <- quantile(ProtMatrix[,i], 0.05)
}
rm(i)

LOD <- vector()
for (i in 1:10675) {
  Index <- which(ProtMatrix[,i] < Quantile[i])
  LOD[i] <- mean(ProtMatrix_Buffer[,i]) + 1.65 * sd(ProtMatrix_Buffer[,i]) + 1.65 * sd(ProtMatrix[Index,i])
}
rm(i)
rm(Index)
rm(Quantile)

#LOD
LOD_TF <- ProtMatrix < matrix(rep(LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
(Proportion_Unified_Formula[2] <- sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF)))
rm(LOD_TF)
rm(LOD)

rm(ProtMatrix)
rm(Analyte_Annotation)
rm(Label)
rm(ProtMatrix_Buffer)
gc()


#SomaLogic ANML

#SomaLogic Data
Dat <- read.csv("Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:20)]
ProtMatrix <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation <- read.csv("Somalogic_Analyte_Annotation_All.csv")
rownames(Analyte_Annotation) <- Analyte_Annotation$AptName

#Remove Control Samples and Retain Buffer Samples
ProtMatrix_Buffer <- ProtMatrix[which(Label$SampleType == "Buffer"),]
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)

#500
Index <- union(which(Label$tech_rep %in% c(2:3)), which(Label$bio_rep_id == "Y"))
ProtMatrix <- ProtMatrix[-Index,]
Label <- Label[-Index,]
rm(Index)

#LOD Calculation

#LOD for Each Protein
LOD <- vector()
for (i in 1:10675) {
  LOD[i] <- median(ProtMatrix_Buffer[,i]) + 4.9 * mad(ProtMatrix_Buffer[,i])
}
rm(i)

#LoD
LOD <- data.frame(SeqID = Analyte_Annotation$SeqId, LOD)
rownames(LOD) <- LOD$SeqID

#Proportion below LOD
rownames(Analyte_Annotation) <- Analyte_Annotation$SeqId
colnames(ProtMatrix) <- Analyte_Annotation$SeqId
colnames(ProtMatrix_Buffer) <- Analyte_Annotation$SeqId

#When comparing with LoD, the RFUs were transformed back to original scale.
LOD_TF <- 2^(ProtMatrix) < matrix(rep(LOD$LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
(Overall_Proportion[2] <- sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF)))

#Proportion below LoD for each sample
Rate_LoD_Sample_SomaLogic2 <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample_SomaLogic2[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample_SomaLogic2)

#Proportion below LoD for each protein
Rate_LoD_Protein_SomaLogic2 <- vector()
for (i in 1:ncol(LOD_TF)) {
  Rate_LoD_Protein_SomaLogic2[i] <- sum(LOD_TF[,i])/nrow(LOD_TF)
}
rm(i)
rm(LOD_TF)
rm(LOD)

#Formula 1

#LOD for Each Protein
LOD <- vector()
for (i in 1:10675) {
  LOD[i] <- mean(ProtMatrix_Buffer[,i]) + 3.3 * sd(ProtMatrix_Buffer[,i])
}
rm(i)

#When comparing with LoD, the RFUs were transformed back to original scale.
LOD_TF <- 2^(ProtMatrix) < matrix(rep(LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
(Proportion_Unified_Formula[3] <- sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF)))

rm(LOD_TF)
rm(LOD)

#Formula 2

#LOD for Each Protein
ProtMatrix <- 2^ProtMatrix
Quantile <- vector()
for (i in 1:10675) {
  Quantile[i] <- quantile(ProtMatrix[,i], 0.05)
}
rm(i)

LOD <- vector()
for (i in 1:10675) {
  Index <- which(ProtMatrix[,i] < Quantile[i])
  LOD[i] <- mean(ProtMatrix_Buffer[,i]) + 1.65 * sd(ProtMatrix_Buffer[,i]) + 1.65 * sd(ProtMatrix[Index,i])
}
rm(i)
rm(Index)
rm(Quantile)

#LOD
LOD_TF <- ProtMatrix < matrix(rep(LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
(Proportion_Unified_Formula[4] <- sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF)))
rm(LOD_TF)
rm(LOD)

rm(ProtMatrix)
rm(Analyte_Annotation)
rm(Label)
rm(ProtMatrix_Buffer)
gc()


#Olink

#Olink Data
Dat <- read.csv("Olink_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:19)]
ProtMatrix <- Dat[,c(20:ncol(Dat))]
rm(Dat)

#Assay Annotation
Assay_Annotation <- read.csv("Olink_Assay_Annotation_All.csv")
rownames(Assay_Annotation) <- Assay_Annotation$OlinkID

#Remove Control Samples
ProtMatrix_Negative_Control <- ProtMatrix[which(Label$SampleType == "NEGATIVE_CONTROL"),]
ProtMatrix <- ProtMatrix[which(Label$SampleType == "SAMPLE"),]
Label <- Label[which(Label$SampleType == "SAMPLE"),]

#216
Index <- union(which(Label$tech_rep %in% c(2:3)), which(Label$bio_rep_id == "Y"))
ProtMatrix <- ProtMatrix[-Index,]
Label <- Label[-Index,]
rm(Index)

#Limit of Detection (LoD) Analysis

#Package
library(OlinkAnalyze)
library(arrow)

#Data Loading
my_NPX_data <- as.data.frame(read_NPX("Q-14739_Chambers_NPX_2024-08-05.parquet"))

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
LOD_TF <- ProtMatrix < ProtMatrix_LOD
LOD_TF <- as.data.frame(LOD_TF)

#Overall
(Overall_Proportion[3] <- sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF)))

#Proportion below LoD for each sample
Rate_LoD_Sample_Olink <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample_Olink[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample_Olink)

#Proportion below LoD for each protein
Rate_LoD_Protein_Olink <- vector()
for (i in 1:ncol(LOD_TF)) {
  Rate_LoD_Protein_Olink[i] <- sum(LOD_TF[,i])/nrow(LOD_TF)
}
rm(i)

rm(LOD_TF)
rm(ProtMatrix_LOD)

#Formula 1

#LOD for Each Protein
LOD <- vector()
for (i in 1:5405) {
  LOD[i] <- mean(ProtMatrix_Negative_Control[,i]) + 3.3 * sd(ProtMatrix_Negative_Control[,i])
}
rm(i)

LOD_TF <- ProtMatrix < matrix(rep(LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
(Proportion_Unified_Formula[5] <- sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF)))
rm(LOD_TF)
rm(LOD)

#Formula 2

#LOD for Each Protein
Quantile <- vector()
for (i in 1:5405) {
  Quantile[i] <- quantile(ProtMatrix[,i], 0.05)
}
rm(i)

LOD <- vector()
for (i in 1:5405) {
  Index <- which(ProtMatrix[,i] < Quantile[i])
  LOD[i] <- mean(ProtMatrix_Negative_Control[,i]) + 1.65 * sd(ProtMatrix_Negative_Control[,i]) + 1.65 * sd(ProtMatrix[Index,i])
}
rm(i)

rm(Index)
rm(Quantile)

LOD_TF <- ProtMatrix < matrix(rep(LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
(Proportion_Unified_Formula[6] <- sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF)))
rm(LOD_TF)
rm(LOD)

rm(ProtMatrix)
rm(Label)
rm(Assay_Annotation)
rm(ProtMatrix_Negative_Control)
gc()


#Thermo Fisher

#Thermo Fisher Data
Dat <- read.csv("Thermo_Fisher_Merged_All.csv")
rownames(Dat) <- Dat$biosample_id
Label <- Dat[,c(1:10)]
ProtMatrix <- Dat[,c(11:ncol(Dat))]
min(ProtMatrix, na.rm = T)
rm(Dat)

#Protein Annotation
Protein_Annotation <- read.csv("Thermo_Fisher_Protein_Annotation_All.csv")
rownames(Protein_Annotation) <- Protein_Annotation$ProteinID

#46
Index <- which(Label$tech_rep %in% 2)
ProtMatrix <- ProtMatrix[-Index,]
Label <- Label[-Index,]
rm(Index)

#Missingness Calculation

#Overall
(Overall_Proportion[4] <- sum(is.na(ProtMatrix))/(nrow(ProtMatrix) * ncol(ProtMatrix)))

#Sample Missingness
Missing_Rate_Sample <- vector()
for (i in 1:nrow(ProtMatrix)) {
  Missing_Rate_Sample[i] <- sum(is.na(ProtMatrix[i,]))/ncol(ProtMatrix)
}
rm(i)
summary(Missing_Rate_Sample)

#Protein Missingness
Missing_Rate_Protein <- vector()
for (i in 1:ncol(ProtMatrix)) {
  Missing_Rate_Protein[i] <- sum(is.na(ProtMatrix[,i]))/nrow(ProtMatrix)
}
rm(i)

rm(ProtMatrix)
rm(Protein_Annotation)
rm(Label)


#Comparison of overall LoD proportion using unified and recommended formulas
library(ggplot2)
Overall <- data.frame(Platform = c(rep("1.SomaLogic pre-ANML", 2), rep("2.SomaLogic ANML", 2), rep("3.Olink", 2), c("1.SomaLogic pre-ANML","2.SomaLogic ANML","3.Olink","4.Thermo Fisher")), Overall = c(Proportion_Unified_Formula, Overall_Proportion), Formula = c(rep(c("1","2"), 3), rep("Recommended",4)))
ggplot(data = Overall, aes(x = Formula, y = Overall, fill = Platform, color = Platform)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("1.SomaLogic pre-ANML" = "#0072B5FF", "2.SomaLogic ANML" = "#E18727FF", "3.Olink" = "#20854EFF", "4.Thermo Fisher" = "#BC3C29FF")) +
  scale_color_manual(values = c("1.SomaLogic pre-ANML" = "#0072B5FF", "2.SomaLogic ANML" = "#E18727FF", "3.Olink" = "#20854EFF", "4.Thermo Fisher" = "#BC3C29FF")) +
  theme_classic() +
  scale_x_discrete(limits = c("Recommended","1","2")) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  labs(y = "Overall missingness", x = "Formula") +
  geom_text(aes(label = round(Overall, 2)), vjust = -0.3, hjust = 0.5, size = 4) +
  ggtitle("Overall missingness")
rm(Overall)
rm(Proportion_Unified_Formula)
rm(Overall_Proportion)


#Sample
Class <- 1:length(c(Rate_LoD_Sample_SomaLogic, Rate_LoD_Sample_SomaLogic2, Rate_LoD_Sample_Olink, Missing_Rate_Sample))
Class[1:length(Rate_LoD_Sample_SomaLogic)] <- "SomaLogic pre-ANML"
Class[(length(Rate_LoD_Sample_SomaLogic) + 1):(length(Rate_LoD_Sample_SomaLogic) + length(Rate_LoD_Sample_SomaLogic2))] <- "SomaLogic ANML"
Class[(length(Rate_LoD_Sample_SomaLogic) + length(Rate_LoD_Sample_SomaLogic2) + 1):(length(Rate_LoD_Sample_SomaLogic) + length(Rate_LoD_Sample_SomaLogic2) + length(Rate_LoD_Sample_Olink))] <- "Olink"
Class[(length(Rate_LoD_Sample_SomaLogic) + length(Rate_LoD_Sample_SomaLogic2) + length(Rate_LoD_Sample_Olink) + 1):(length(Rate_LoD_Sample_SomaLogic) + length(Rate_LoD_Sample_SomaLogic2) + length(Rate_LoD_Sample_Olink) + length(Missing_Rate_Sample))] <- "Thermo Fisher"
Data <- data.frame(Rate_LoD_Sample = c(Rate_LoD_Sample_SomaLogic, Rate_LoD_Sample_SomaLogic2, Rate_LoD_Sample_Olink, Missing_Rate_Sample), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Rate_LoD_Sample, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("SomaLogic pre-ANML" = "#0072B5FF", "SomaLogic ANML" = "#E18727FF", "Olink" = "#20854EFF", "Thermo Fisher" = "#BC3C29FF")) +
  scale_color_manual(values = c("SomaLogic pre-ANML" = "#0072B5FF", "SomaLogic ANML" = "#E18727FF", "Olink" = "#20854EFF", "Thermo Fisher" = "#BC3C29FF")) +
  scale_x_discrete(limits = c("SomaLogic pre-ANML","SomaLogic ANML","Olink","Thermo Fisher")) +
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
rm(Missing_Rate_Sample)


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
Rate_Thermo[1] <- length(which(Missing_Rate_Protein == 0))
Rate_Thermo[2] <- length(which(Missing_Rate_Protein > 0.2))
Rate_Thermo[3] <- length(which(Missing_Rate_Protein > 0.5))
Rate_Thermo[4] <- length(which(Missing_Rate_Protein == 1))

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
  geom_hline(yintercept = 10675, color = "#0072B5FF") +
  geom_hline(yintercept = 5405, color = "#20854EFF") +
  geom_hline(yintercept = 5945, color = "#BC3C29FF") +
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
rm(Missing_Rate_Protein)
gc()


#pre-ANML & ANMl Stratification by Sets
pre <- read.csv("Sample_Missingness_pre-ANML_Sets.csv")
ANML <- read.csv("Sample_Missingness_ANML_Sets.csv")
Data <- rbind(pre, ANML)

library(ggplot2)
ggplot(data = Data, aes(x = Normalization, y = Rate_LoD_Sample, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_color_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_x_discrete(limits = c("pre-ANML","ANML")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ylim(0,0.9) +
  facet_wrap(~Class) +
  labs(y = "Sample missingness", x = "Normalization")
rm(pre)
rm(ANML)
rm(Data)

#pre-ANML & ANMl Stratification by Sets
pre <- read.csv("Protein_Missingness_pre-ANML_Sets.csv")
ANML <- read.csv("Protein_Missingness_ANML_Sets.csv")
Data <- rbind(pre, ANML)

library(ggplot2)
ggplot(data = Data, aes(x = Normalization, y = Rate, fill = Set, color = Set)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_color_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_x_discrete(limits = c("pre-ANML","ANML")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ylim(0,0.9) +
  facet_wrap(~Set) +
  labs(y = "Protein missingness", x = "Normalization") +
  geom_text(aes(label = paste0(round(Rate * 100, 1), "%")), color = "black", position = position_dodge(width = 1), vjust = -0.5, size = 4)

rm(pre)
rm(ANML)
rm(Data)


#pre-ANML & ANMl Stratification by Dilutions
pre <- read.csv("Sample_Missingness_pre-ANML_Dilution.csv")
ANML <- read.csv("Sample_Missingness_ANML_Dilution.csv")
Data <- rbind(pre, ANML)
Data$Class[which(Data$Class == "1:5")] <- "1.1:5"
Data$Class[which(Data$Class == "1:200")] <- "2.1:200"
Data$Class[which(Data$Class == "1:20000")] <- "3.1:20000"

library(ggplot2)
ggplot(data = Data, aes(x = Normalization, y = Rate_LoD_Sample, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("1.1:5" = "#0072B5FF", "2.1:200" = "#E18727FF", "3.1:20000" = "#20854EFF")) +
  scale_color_manual(values = c("1.1:5" = "#0072B5FF", "2.1:200" = "#E18727FF", "3.1:20000" = "#20854EFF")) +
  scale_x_discrete(limits = c("pre-ANML","ANML")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ylim(0,0.9) +
  facet_wrap(~Class) +
  labs(y = "Sample missingness", x = "Normalization")
rm(pre)
rm(ANML)
rm(Data)

#pre-ANML & ANMl Stratification by Dilutions
pre <- read.csv("Protein_Missingness_pre-ANML_Dilution.csv")
ANML <- read.csv("Protein_Missingness_ANML_Dilution.csv")
Data <- rbind(pre, ANML)
Data$Set[which(Data$Set == "1:5")] <- "1.1:5"
Data$Set[which(Data$Set == "1:200")] <- "2.1:200"
Data$Set[which(Data$Set == "1:20000")] <- "3.1:20000"

library(ggplot2)
ggplot(data = Data, aes(x = Normalization, y = Rate, fill = Set, color = Set)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("1.1:5" = "#0072B5FF", "2.1:200" = "#E18727FF", "3.1:20000" = "#20854EFF")) +
  scale_color_manual(values = c("1.1:5" = "#0072B5FF", "2.1:200" = "#E18727FF", "3.1:20000" = "#20854EFF")) +
  scale_x_discrete(limits = c("pre-ANML","ANML")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ylim(0,0.9) +
  facet_wrap(~Set) +
  labs(y = "Protein missingness", x = "Normalization") +
  geom_text(aes(label = paste0(round(Rate * 100, 1), "%")), color = "black", position = position_dodge(width = 1), vjust = -0.5, size = 4)

rm(pre)
rm(ANML)
rm(Data)


