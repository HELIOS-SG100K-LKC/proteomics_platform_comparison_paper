#SomaLogic ANML LOD


#Data Loading

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

#Official LOD for Each Protein
library(readxl)
LOD_Official <- read_excel("SomaScan_11K_Annotated_Content.xlsx")
LOD_Official <- as.data.frame(LOD_Official[-c(1:4),c(1,10)])
colnames(LOD_Official) <- c("SeqID","LOD")
rownames(LOD_Official) <- LOD_Official$SeqID

#Comparison
LOD <- data.frame(SeqID = Analyte_Annotation$SeqId, LOD)
rownames(LOD) <- LOD$SeqID
LOD_Official <- LOD_Official[rownames(LOD),]
LOD_Official$LOD <- as.numeric(LOD_Official$LOD)
cor(LOD$LOD, LOD_Official$LOD, method = "pearson")
rm(LOD_Official)

#Proportion below LOD
rownames(Analyte_Annotation) <- Analyte_Annotation$SeqId
colnames(ProtMatrix) <- Analyte_Annotation$SeqId
colnames(ProtMatrix_Buffer) <- Analyte_Annotation$SeqId

#When comparing with LoD, the RFUs were transformed back to original scale
LOD_TF <- 2^(ProtMatrix) < matrix(rep(LOD$LOD, each = nrow(ProtMatrix)), nrow = nrow(ProtMatrix), ncol = ncol(ProtMatrix), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF))

#Proportion below LoD for each sample
Rate_LoD_Sample <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample)

#Proportion below LoD for each protein
Rate_LoD_Protein <- vector()
for (i in 1:ncol(LOD_TF)) {
  Rate_LoD_Protein[i] <- sum(LOD_TF[,i])/nrow(LOD_TF)
}
rm(i)
rm(LOD_TF)

#Sample
Class <- 1:800
Class[1:800] <- ""
Data <- data.frame(Rate_LoD_Sample = Rate_LoD_Sample, Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Rate_LoD_Sample)) +
  geom_violin(trim = TRUE, fill = "#E18727FF", linewidth = 0.8, colour = "#E18727FF") +
  geom_boxplot(fill = "white", width = 0.3, linewidth = 1, outlier.size = 1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  ggtitle("SomaLogic ANML") +
  labs(y = "Sample missingness", x = "Violin plot")
rm(Class)
rm(Data)

#Protein
Rate <- vector()
Rate[1] <- length(Rate_LoD_Protein)
Rate[2] <- length(which(Rate_LoD_Protein == 0))
Rate[3] <- length(which(Rate_LoD_Protein > 0.2))
Rate[4] <- length(which(Rate_LoD_Protein > 0.5))
Rate[5] <- length(which(Rate_LoD_Protein == 1))
Rate_Percentage <- Rate / length(Rate_LoD_Protein) * 100
Data <- data.frame(Rate, Rate_Percentage, Threshold = c("All", "=0%", ">20%", ">50%", "=100%"))

ggplot(Data, aes(x = Threshold, y = Rate)) +
  geom_bar(stat = "identity", position = "dodge", color = "#E18727FF", fill = "#E18727FF") +
  geom_text(aes(label = paste0(Rate, " (", sprintf("%.0f", Rate_Percentage), "%)")),
            hjust = 0.4, vjust = -0.4, size = 3.5, color = "black") +
  scale_x_discrete(limits = Data$Threshold) +
  scale_y_continuous(
    name = "Number of proteins",
    sec.axis = sec_axis(~ . / length(Rate_LoD_Protein) * 100, name = "Proportion of all proteins (%)")
  ) +
  theme_classic() +
  ggtitle(paste0("SomaLogic ANML")) +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  labs(x = "Protein missingness", y = "Number of proteins")

rm(Rate)
rm(Data)
rm(Rate_Percentage)


#Different Sets

#First set (SomaScan 5K)
library(readxl)
ProteinList5K <- read_excel("SomaScan_Menu_1.3K_5K_7K.xlsx")
ProteinList5K <- ProteinList5K[-c(1:2),10]
ProteinList5K <- unique(ProteinList5K$...10)

ProtMatrix1 <- ProtMatrix[,which(Analyte_Annotation$SeqId %in% ProteinList5K)]
Analyte_Annotation1 <- Analyte_Annotation[which(Analyte_Annotation$SeqId %in% ProteinList5K),]

#Proportion below LOD
LOD_TF <- 2^(ProtMatrix1) < matrix(rep(LOD$LOD[which(Analyte_Annotation$SeqId %in% ProteinList5K)], each = nrow(ProtMatrix1)), nrow = nrow(ProtMatrix1), ncol = ncol(ProtMatrix1), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF))

#Proportion below LoD for each sample
Rate_LoD_Sample1 <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample1[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample1)
rm(LOD_TF)

#Proportion below LoD for each protein
Rate_LoD_Protein1 <- Rate_LoD_Protein[which(Analyte_Annotation$SeqId %in% ProteinList5K)]

nrow(Analyte_Annotation1)
length(which(Rate_LoD_Protein1 > 0.2))
length(which(Rate_LoD_Protein1 > 0.2))/nrow(Analyte_Annotation1)

rm(Analyte_Annotation1)
rm(ProtMatrix1)
rm(Rate_LoD_Protein1)

#Second set (SomaScan 7K - SomaScan 5K)
ProteinList7K <- read_excel("SomaScan_Menu_1.3K_5K_7K.xlsx")
ProteinList7K <- ProteinList7K[-c(1:2),1]
ProteinList7K <- unique(ProteinList7K$`Supplementary Data 1. List of all SOMAmers in the plasma 7k, 5k, and 1.3k SomaScan assays.`)

ProteinList2 <- setdiff(ProteinList7K, ProteinList5K)
ProtMatrix2 <- ProtMatrix[,which(Analyte_Annotation$SeqId %in% ProteinList2)]
Analyte_Annotation2 <- Analyte_Annotation[which(Analyte_Annotation$SeqId %in% ProteinList2),]

#Proportion below LOD
LOD_TF <- 2^(ProtMatrix2) < matrix(rep(LOD$LOD[which(Analyte_Annotation$SeqId %in% ProteinList2)], each = nrow(ProtMatrix2)), nrow = nrow(ProtMatrix2), ncol = ncol(ProtMatrix2), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF))

#Proportion below LoD for each sample
Rate_LoD_Sample2 <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample2[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample2)
rm(LOD_TF)

#Proportion below LoD for each protein
Rate_LoD_Protein2 <- Rate_LoD_Protein[which(Analyte_Annotation$SeqId %in% ProteinList2)]

nrow(Analyte_Annotation2)
length(which(Rate_LoD_Protein2 > 0.2))
length(which(Rate_LoD_Protein2 > 0.2))/nrow(Analyte_Annotation2)

rm(Analyte_Annotation2)
rm(ProtMatrix2)
rm(Rate_LoD_Protein2)

#Third set (SomaScan 11K - SomaScan 7K)
ProteinList3 <- setdiff(Analyte_Annotation$SeqId, ProteinList7K)
ProtMatrix3 <- ProtMatrix[,which(Analyte_Annotation$SeqId %in% ProteinList3)]
Analyte_Annotation3 <- Analyte_Annotation[which(Analyte_Annotation$SeqId %in% ProteinList3),]

#Proportion below LOD
LOD_TF <- 2^(ProtMatrix3) < matrix(rep(LOD$LOD[which(Analyte_Annotation$SeqId %in% ProteinList3)], each = nrow(ProtMatrix3)), nrow = nrow(ProtMatrix3), ncol = ncol(ProtMatrix3), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF))

#Proportion below LoD for each sample
Rate_LoD_Sample3 <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample3[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample3)
rm(LOD_TF)

#Proportion below LoD for each protein
Rate_LoD_Protein3 <- Rate_LoD_Protein[which(Analyte_Annotation$SeqId %in% ProteinList3)]

nrow(Analyte_Annotation3)
length(which(Rate_LoD_Protein3 > 0.2))
length(which(Rate_LoD_Protein3 > 0.2))/nrow(Analyte_Annotation3)

rm(Analyte_Annotation3)
rm(ProtMatrix3)
rm(Rate_LoD_Protein3)

#Comparison

#Sample
Class <- 1:length(c(Rate_LoD_Sample1, Rate_LoD_Sample2, Rate_LoD_Sample3))
Class[1:length(Rate_LoD_Sample1)] <- "Set 1"
Class[(length(Rate_LoD_Sample1) + 1):(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2))] <- "Set 2"
Class[(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2) + 1):(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2) + length(Rate_LoD_Sample3))] <- "Set 3"
Data <- data.frame(Rate_LoD_Sample = c(Rate_LoD_Sample1, Rate_LoD_Sample2, Rate_LoD_Sample3), Class = Class)
Data$Class <- as.factor(Data$Class)
ggplot(data = Data, aes(x = Class, y = Rate_LoD_Sample, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_color_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_x_discrete(limits = c("Set 1","Set 2","Set 3")) +
  theme_classic() +
  ggtitle("SomaLogic ANML") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ylim(0,0.9) +
  labs(y = "Sample missingness", x = "Different sets")

Data$Normalization <- rep("ANML", 2400)
write.csv(Data, "Sample_Missingness_ANML_Sets.csv", row.names = F)

rm(Class)
rm(Data)

#Protein
Data <- data.frame(Rate = c(0.0417004,0.03042328,0.04066917), Set = c("Set 1","Set 2","Set 3"))
ggplot(Data, aes(x = Set, y = Rate, fill = Set, color = Set)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_color_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_x_discrete(limits = c("Set 1","Set 2","Set 3")) +
  scale_y_continuous(limits = c(0, 0.9)) +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(legend.position = "none") +
  ggtitle("SomaLogic ANML") +
  labs(y = "Proportion of proteins below LoD", x = "Different sets") +
  geom_text(aes(label = paste0(round(Rate * 100, 1), "%")), position = position_dodge(width = 1), vjust = -0.5, size = 4)

Data$Normalization <- rep("ANML", 3)
write.csv(Data, "Protein_Missingness_ANML_Sets.csv", row.names = F)

rm(Data)

rm(ProteinList2)
rm(ProteinList3)
rm(ProteinList5K)
rm(ProteinList7K)
rm(Rate_LoD_Sample1)
rm(Rate_LoD_Sample2)
rm(Rate_LoD_Sample3)


#Dilution Series

#1:5
ProtMatrix1 <- ProtMatrix[,which(Analyte_Annotation$Dilution2 == 1/5)]
Analyte_Annotation1 <- Analyte_Annotation[which(Analyte_Annotation$Dilution2 == 1/5),]

#Proportion below LOD
LOD_TF <- 2^(ProtMatrix1) < matrix(rep(LOD$LOD[which(Analyte_Annotation$Dilution2 == 1/5)], each = nrow(ProtMatrix1)), nrow = nrow(ProtMatrix1), ncol = ncol(ProtMatrix1), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF))

#Proportion below LoD for each sample
Rate_LoD_Sample1 <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample1[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample1)
rm(LOD_TF)

#Proportion below LoD for each protein
Rate_LoD_Protein1 <- Rate_LoD_Protein[which(Analyte_Annotation$Dilution2 == 1/5)]

nrow(Analyte_Annotation1)
length(which(Rate_LoD_Protein1 > 0.2))
length(which(Rate_LoD_Protein1 > 0.2))/nrow(Analyte_Annotation1)

rm(Analyte_Annotation1)
rm(ProtMatrix1)
rm(Rate_LoD_Protein1)

#1:200
ProtMatrix2 <- ProtMatrix[,which(Analyte_Annotation$Dilution2 == 1/200)]
Analyte_Annotation2 <- Analyte_Annotation[which(Analyte_Annotation$Dilution2 == 1/200),]

#Proportion below LOD
LOD_TF <- 2^(ProtMatrix2) < matrix(rep(LOD$LOD[which(Analyte_Annotation$Dilution2 == 1/200)], each = nrow(ProtMatrix2)), nrow = nrow(ProtMatrix2), ncol = ncol(ProtMatrix2), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF))

#Proportion below LoD for each sample
Rate_LoD_Sample2 <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample2[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample2)
rm(LOD_TF)

#Proportion below LoD for each protein
Rate_LoD_Protein2 <- Rate_LoD_Protein[which(Analyte_Annotation$Dilution2 == 1/200)]

nrow(Analyte_Annotation2)
length(which(Rate_LoD_Protein2 > 0.2))
length(which(Rate_LoD_Protein2 > 0.2))/nrow(Analyte_Annotation2)

rm(Analyte_Annotation2)
rm(ProtMatrix2)
rm(Rate_LoD_Protein2)

#1:20000
ProtMatrix3 <- ProtMatrix[,which(Analyte_Annotation$Dilution2 == 1/20000)]
Analyte_Annotation3 <- Analyte_Annotation[which(Analyte_Annotation$Dilution2 == 1/20000),]

#Proportion below LOD
LOD_TF <- 2^(ProtMatrix3) < matrix(rep(LOD$LOD[which(Analyte_Annotation$Dilution2 == 1/20000)], each = nrow(ProtMatrix3)), nrow = nrow(ProtMatrix3), ncol = ncol(ProtMatrix3), byrow = F)
LOD_TF <- as.data.frame(LOD_TF)

#Overall
sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF))

#Proportion below LoD for each sample
Rate_LoD_Sample3 <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample3[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample3)
rm(LOD_TF)

#Proportion below LoD for each protein
Rate_LoD_Protein3 <- Rate_LoD_Protein[which(Analyte_Annotation$Dilution2 == 1/20000)]

nrow(Analyte_Annotation3)
length(which(Rate_LoD_Protein3 > 0.2))
length(which(Rate_LoD_Protein3 > 0.2))/nrow(Analyte_Annotation3)

rm(Analyte_Annotation3)
rm(ProtMatrix3)
rm(Rate_LoD_Protein3)

#Comparison

#Sample
Class <- 1:length(c(Rate_LoD_Sample1, Rate_LoD_Sample2, Rate_LoD_Sample3))
Class[1:length(Rate_LoD_Sample1)] <- "1:5"
Class[(length(Rate_LoD_Sample1) + 1):(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2))] <- "1:200"
Class[(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2) + 1):(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2) + length(Rate_LoD_Sample3))] <- "1:20000"
Data <- data.frame(Rate_LoD_Sample = c(Rate_LoD_Sample1, Rate_LoD_Sample2, Rate_LoD_Sample3), Class = Class)
Data$Class <- as.factor(Data$Class)
ggplot(data = Data, aes(x = Class, y = Rate_LoD_Sample, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("1:5" = "#0072B5FF", "1:200" = "#E18727FF", "1:20000" = "#20854EFF")) +
  scale_color_manual(values = c("1:5" = "#0072B5FF", "1:200" = "#E18727FF", "1:20000" = "#20854EFF")) +
  scale_x_discrete(limits = c("1:5","1:200","1:20000")) +
  theme_classic() +
  ggtitle("SomaLogic ANML") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ylim(0,0.75) +
  labs(y = "Sample missingness", x = "Dilution series")

Data$Normalization <- rep("ANML", 2400)
write.csv(Data, "Sample_Missingness_ANML_Dilution.csv", row.names = F)

rm(Class)
rm(Data)

rm(Rate_LoD_Sample1)
rm(Rate_LoD_Sample2)
rm(Rate_LoD_Sample3)

#Protein
Data <- data.frame(Rate = c(0.04688593,0.006349206,0.009478673), Set = c("1:5","1:200","1:20000"))
ggplot(Data, aes(x = Set, y = Rate, fill = Set, color = Set)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c("1:5" = "#0072B5FF", "1:200" = "#E18727FF", "1:20000" = "#20854EFF")) +
  scale_color_manual(values = c("1:5" = "#0072B5FF", "1:200" = "#E18727FF", "1:20000" = "#20854EFF")) +
  scale_x_discrete(limits = c("1:5","1:200","1:20000")) +
  scale_y_continuous(limits = c(0, 0.75)) +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(legend.position = "none") +
  ggtitle("SomaLogic ANML") +
  labs(y = "Proportion of proteins below LoD", x = "Dilution") +
  geom_text(aes(label = paste0(round(Rate * 100, 1), "%")), position = position_dodge(width = 1), vjust = -0.5, size = 4)

Data$Normalization <- rep("ANML", 3)
write.csv(Data, "Protein_Missingness_ANML_Dilution.csv", row.names = F)

rm(Data)


#Proteins with proportion below LoD > 20%
Dat <- data.frame(Protein = colnames(ProtMatrix), Rate_LoD = Rate_LoD_Protein)
Dat <- Dat[Dat$Rate_LoD > 0.2,]
Dat <- cbind(Dat, Analyte_Annotation[Dat$Protein,])
write.csv(Dat, "SomaLogic_HighProportionProteins_0.2.csv", row.names = F)
rm(Dat)

rm(Rate_LoD_Sample)
gc()

#Output
write.csv(data.frame(SeqID = Analyte_Annotation$SeqId, UniProt = Analyte_Annotation$UniProt, log2LOD = log2(LOD$LOD), PropBelowLOD = Rate_LoD_Protein), "SomaLogic_LOD.csv", row.names = F)
rm(Rate_LoD_Protein)


#Remove All
rm(list = ls())
gc()

