#Olink LoD


#Data Loading

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
ProtMatrix_Negative_Control <- ProtMatrix[which(Label$SampleType == "NEGATIVE_CONTROL"),]
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

#Official LOD for Each Protein
library(readxl)
LOD_Official <- read_excel("olink-explore-ht-validation-data-results.xlsx", sheet = 2)
LOD_Official <- as.data.frame(LOD_Official[-1,c(1,6)])
colnames(LOD_Official) <- c("UNIPROT","LOD")
LOD_Official <- LOD_Official[match(unique(LOD_Official$UNIPROT), LOD_Official$UNIPROT),]
rownames(LOD_Official) <- LOD_Official$UNIPROT

#Comparison
LOD <- apply(ProtMatrix_LOD, 2, median)
LOD <- data.frame(OlinkID = Assay_Annotation$OlinkID, LOD)

LOD2 <- data.frame(UniProt = Assay_Annotation$UniProt, LOD$LOD)
LOD2 <- LOD2[match(unique(LOD2$UniProt), LOD2$UniProt),]
rownames(LOD2) <- LOD2$UniProt

LOD_Official <- LOD_Official[intersect(rownames(LOD2), LOD_Official$UNIPROT),]
LOD_Official$LOD <- as.numeric(LOD_Official$LOD)
LOD2 <- LOD2[rownames(LOD_Official),]
cor(LOD2$LOD, LOD_Official$LOD, method = "pearson", use = "pairwise.complete.obs")
rm(LOD_Official)
rm(LOD2)

#Proportion below LOD
LOD_TF <- ProtMatrix < ProtMatrix_LOD
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

LOD <- data.frame(OlinkID = LOD$OlinkID, UniProt = Assay_Annotation$UniProt, MedianLOD = LOD$LOD, PropBelowLOD = Rate_LoD_Protein)
write.csv(LOD, "Olink_LoD.csv", row.names = F)
rm(LOD)

#Sample
Class <- 1:514
Class[1:514] <- ""
Data <- data.frame(Rate_LoD_Sample = Rate_LoD_Sample, Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Rate_LoD_Sample)) +
  geom_violin(trim = TRUE, fill = "#20854EFF", linewidth = 0.8, colour = "#20854EFF") +
  geom_boxplot(fill = "white", width = 0.3, linewidth = 1, outlier.size = 1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  ggtitle("Olink") +
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
  geom_bar(stat = "identity", position = "dodge", color = "#20854EFF", fill = "#20854EFF") +
  geom_text(aes(label = paste0(Rate, " (", sprintf("%.0f", Rate_Percentage), "%)")),
            hjust = 0.4, vjust = -0.4, size = 3.5, color = "black") +
  scale_x_discrete(limits = Data$Threshold) +
  scale_y_continuous(
    name = "Number of proteins",
    sec.axis = sec_axis(~ . / length(Rate_LoD_Protein) * 100, name = "Proportion of all proteins (%)")
  ) +
  theme_classic() +
  ggtitle(paste0("Olink")) +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  labs(x = "Protein missingness", y = "Number of proteins")

rm(Rate)
rm(Data)
rm(Rate_Percentage)


#Different Sets

#First set: Olink Explore 1536
ProteinList1536 <- read_excel("explore-1536-assay-list-20210227-web-1.xlsx")
colnames(ProteinList1536) <- ProteinList1536[1,]
ProteinList1536 <- ProteinList1536[-1,]

ProtMatrix1 <- ProtMatrix[,which(Assay_Annotation$UniProt %in% ProteinList1536$`Uniprot ID`)]
ProtMatrix_LOD1 <- ProtMatrix_LOD[,which(Assay_Annotation$UniProt %in% ProteinList1536$`Uniprot ID`)]
Assay_Annotation1 <- Assay_Annotation[which(Assay_Annotation$UniProt %in% ProteinList1536$`Uniprot ID`),]

#Proportion below LOD
LOD_TF <- ProtMatrix1 < ProtMatrix_LOD1
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
Rate_LoD_Protein1 <- Rate_LoD_Protein[which(Assay_Annotation$UniProt %in% ProteinList1536$`Uniprot ID`)]

length(Rate_LoD_Protein1)
length(which(Rate_LoD_Protein1 > 0.2))
length(which(Rate_LoD_Protein1 > 0.2))/length(Rate_LoD_Protein1)

rm(Rate_LoD_Protein1)
rm(ProtMatrix1)
rm(ProtMatrix_LOD1)
rm(Assay_Annotation1)

#Second set: Olink Explore 3072 - Olink Explore 1536
ProteinList3072 <- read_excel("olink-explore-3072-validation-data-results.xlsx")
colnames(ProteinList3072) <- ProteinList3072[1,]
ProteinList3072 <- ProteinList3072[-1,]
ProteinList3072 <- ProteinList3072[-which(ProteinList3072$UniProt %in% ProteinList1536$`Uniprot ID`),]

ProtMatrix2 <- ProtMatrix[,which(Assay_Annotation$UniProt %in% ProteinList3072$UniProt)]
ProtMatrix_LOD2 <- ProtMatrix_LOD[,which(Assay_Annotation$UniProt %in% ProteinList3072$UniProt)]
Assay_Annotation2 <- Assay_Annotation[which(Assay_Annotation$UniProt %in% ProteinList3072$UniProt),]

#Proportion below LOD
LOD_TF <- ProtMatrix2 < ProtMatrix_LOD2
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
Rate_LoD_Protein2 <- Rate_LoD_Protein[which(Assay_Annotation$UniProt %in% ProteinList3072$UniProt)]

length(Rate_LoD_Protein2)
length(which(Rate_LoD_Protein2 > 0.2))
length(which(Rate_LoD_Protein2 > 0.2))/length(Rate_LoD_Protein2)

rm(Rate_LoD_Protein2)
rm(ProtMatrix2)
rm(ProtMatrix_LOD2)
rm(Assay_Annotation2)

#Third set: Olink Explore HT - Olink Explore 3072
ProteinList3072 <- read_excel("olink-explore-3072-validation-data-results.xlsx")
colnames(ProteinList3072) <- ProteinList3072[1,]
ProteinList3072 <- ProteinList3072[-1,]
ProtList <- Assay_Annotation[-which(Assay_Annotation$UniProt %in% ProteinList3072$UniProt),]

ProtMatrix3 <- ProtMatrix[,which(Assay_Annotation$UniProt %in% ProtList$UniProt)]
ProtMatrix_LOD3 <- ProtMatrix_LOD[,which(Assay_Annotation$UniProt %in% ProtList$UniProt)]
Assay_Annotation3 <- Assay_Annotation[which(Assay_Annotation$UniProt %in% ProtList$UniProt),]

#Proportion below LOD
LOD_TF <- ProtMatrix3 < ProtMatrix_LOD3
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
Rate_LoD_Protein3 <- Rate_LoD_Protein[which(Assay_Annotation$UniProt %in% ProtList$UniProt)]

length(Rate_LoD_Protein3)
length(which(Rate_LoD_Protein3 > 0.2))
length(which(Rate_LoD_Protein3 > 0.2))/length(Rate_LoD_Protein3)

rm(Rate_LoD_Protein3)
rm(ProtMatrix3)
rm(ProtMatrix_LOD3)
rm(Assay_Annotation3)

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
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ggtitle("Olink") +
  ylim(0,0.9) +
  labs(y = "Sample missingness", x = "Different sets")
rm(Class)
rm(Data)

#Protein
Data <- data.frame(Rate = c(0.1839161,0.5236749,0.8710938), Set = c("Set 1","Set 2","Set 3"))
ggplot(Data, aes(x = Set, y = Rate, fill = Set, color = Set)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_color_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_x_discrete(limits = c("Set 1","Set 2","Set 3")) +
  scale_y_continuous(limits = c(0,0.9)) +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(legend.position = "none") +
  ggtitle("Olink") +
  labs(y = "Proportion of proteins below LoD", x = "Different sets") +
  geom_text(aes(label = paste0(round(Rate * 100, 1), "%")), position = position_dodge(width = 1), vjust = -0.5, size = 4)
rm(Data)

rm(Rate_LoD_Sample1)
rm(Rate_LoD_Sample2)
rm(Rate_LoD_Sample3)
rm(ProteinList1536)
rm(ProteinList3072)
rm(ProtList)


#Dilution series
Assay_Annotation$Dilution <- Assay_Annotation$Block
Assay_Annotation$Dilution[which(Assay_Annotation$Block %in% c(1,2,3,4))] <- "1:1"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 5)] <- "1:10"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 6)] <- "1:100"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 7)] <- "1:1000"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 8)] <- "1:100000"

#1:1
ProtMatrix1 <- ProtMatrix[,which(Assay_Annotation$Dilution == "1:1")]
ProtMatrix_LOD1 <- ProtMatrix_LOD[,which(Assay_Annotation$Dilution == "1:1")]
Assay_Annotation1 <- Assay_Annotation[which(Assay_Annotation$Dilution == "1:1"),]

#Proportion below LOD
LOD_TF <- ProtMatrix1 < ProtMatrix_LOD1
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
Rate_LoD_Protein1 <- Rate_LoD_Protein[which(Assay_Annotation$Dilution == "1:1")]

length(Rate_LoD_Protein1)
length(which(Rate_LoD_Protein1 > 0.2))
length(which(Rate_LoD_Protein1 > 0.2))/length(Rate_LoD_Protein1)

rm(Rate_LoD_Protein1)
rm(ProtMatrix1)
rm(ProtMatrix_LOD1)
rm(Assay_Annotation1)

#1:10
ProtMatrix2 <- ProtMatrix[,which(Assay_Annotation$Dilution == "1:10")]
ProtMatrix_LOD2 <- ProtMatrix_LOD[,which(Assay_Annotation$Dilution == "1:10")]
Assay_Annotation2 <- Assay_Annotation[which(Assay_Annotation$Dilution == "1:10"),]

#Proportion below LOD
LOD_TF <- ProtMatrix2 < ProtMatrix_LOD2
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
Rate_LoD_Protein2 <- Rate_LoD_Protein[which(Assay_Annotation$Dilution == "1:10")]

length(Rate_LoD_Protein2)
length(which(Rate_LoD_Protein2 > 0.2))
length(which(Rate_LoD_Protein2 > 0.2))/length(Rate_LoD_Protein2)

rm(Rate_LoD_Protein2)
rm(ProtMatrix2)
rm(ProtMatrix_LOD2)
rm(Assay_Annotation2)

#1:100
ProtMatrix3 <- ProtMatrix[,which(Assay_Annotation$Dilution == "1:100")]
ProtMatrix_LOD3 <- ProtMatrix_LOD[,which(Assay_Annotation$Dilution == "1:100")]
Assay_Annotation3 <- Assay_Annotation[which(Assay_Annotation$Dilution == "1:100"),]

#Proportion below LOD
LOD_TF <- ProtMatrix3 < ProtMatrix_LOD3
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
Rate_LoD_Protein3 <- Rate_LoD_Protein[which(Assay_Annotation$Dilution == "1:100")]

length(Rate_LoD_Protein3)
length(which(Rate_LoD_Protein3 > 0.2))
length(which(Rate_LoD_Protein3 > 0.2))/length(Rate_LoD_Protein3)

rm(Rate_LoD_Protein3)
rm(ProtMatrix3)
rm(ProtMatrix_LOD3)
rm(Assay_Annotation3)

#1:1000
ProtMatrix4 <- ProtMatrix[,which(Assay_Annotation$Dilution == "1:1000")]
ProtMatrix_LOD4 <- ProtMatrix_LOD[,which(Assay_Annotation$Dilution == "1:1000")]
Assay_Annotation4 <- Assay_Annotation[which(Assay_Annotation$Dilution == "1:1000"),]

#Proportion below LOD
LOD_TF <- ProtMatrix4 < ProtMatrix_LOD4
LOD_TF <- as.data.frame(LOD_TF)

#Overall
sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF))

#Proportion below LoD for each sample
Rate_LoD_Sample4 <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample4[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample4)
rm(LOD_TF)

#Proportion below LoD for each protein
Rate_LoD_Protein4 <- Rate_LoD_Protein[which(Assay_Annotation$Dilution == "1:1000")]

length(Rate_LoD_Protein4)
length(which(Rate_LoD_Protein4 > 0.2))
length(which(Rate_LoD_Protein4 > 0.2))/length(Rate_LoD_Protein4)

rm(Rate_LoD_Protein4)
rm(ProtMatrix4)
rm(ProtMatrix_LOD4)
rm(Assay_Annotation4)

#1:100000
ProtMatrix5 <- ProtMatrix[,which(Assay_Annotation$Dilution == "1:100000")]
ProtMatrix_LOD5 <- ProtMatrix_LOD[,which(Assay_Annotation$Dilution == "1:100000")]
Assay_Annotation5 <- Assay_Annotation[which(Assay_Annotation$Dilution == "1:100000"),]

#Proportion below LOD
LOD_TF <- ProtMatrix5 < ProtMatrix_LOD5
LOD_TF <- as.data.frame(LOD_TF)

#Overall
sum(LOD_TF)/(nrow(LOD_TF) * ncol(LOD_TF))

#Proportion below LoD for each sample
Rate_LoD_Sample5 <- vector()
for (i in 1:nrow(LOD_TF)) {
  Rate_LoD_Sample5[i] <- sum(LOD_TF[i,])/ncol(LOD_TF)
}
rm(i)
summary(Rate_LoD_Sample5)
rm(LOD_TF)

#Proportion below LoD for each protein
Rate_LoD_Protein5 <- Rate_LoD_Protein[which(Assay_Annotation$Dilution == "1:100000")]

length(Rate_LoD_Protein5)
length(which(Rate_LoD_Protein5 > 0.2))
length(which(Rate_LoD_Protein5 > 0.2))/length(Rate_LoD_Protein5)

rm(Rate_LoD_Protein5)
rm(ProtMatrix5)
rm(ProtMatrix_LOD5)
rm(Assay_Annotation5)

#Sample
Class <- 1:length(c(Rate_LoD_Sample1, Rate_LoD_Sample2, Rate_LoD_Sample3, Rate_LoD_Sample4, Rate_LoD_Sample5))
Class[1:length(Rate_LoD_Sample1)] <- "1:1"
Class[(length(Rate_LoD_Sample1) + 1):(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2))] <- "1:10"
Class[(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2) + 1):(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2) + length(Rate_LoD_Sample3))] <- "1:100"
Class[(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2) + length(Rate_LoD_Sample3) + 1):(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2) + length(Rate_LoD_Sample3) + length(Rate_LoD_Sample4))] <- "1:1000"
Class[(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2) + length(Rate_LoD_Sample3) + length(Rate_LoD_Sample4) + 1):(length(Rate_LoD_Sample1) + length(Rate_LoD_Sample2) + length(Rate_LoD_Sample3) + length(Rate_LoD_Sample4) + length(Rate_LoD_Sample5))] <- "1:100000"
Data <- data.frame(Rate_LoD_Sample = c(Rate_LoD_Sample1, Rate_LoD_Sample2, Rate_LoD_Sample3, Rate_LoD_Sample4, Rate_LoD_Sample5), Class = Class)
Data$Class <- as.factor(Data$Class)
ggplot(data = Data, aes(x = Class, y = Rate_LoD_Sample, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("1:1" = "#0072B5FF", "1:10" = "#E18727FF", "1:100" = "#20854EFF", "1:1000" = "#BC3C29FF", "1:100000" = "#7876B1FF")) +
  scale_color_manual(values = c("1:1" = "#0072B5FF", "1:10" = "#E18727FF", "1:100" = "#20854EFF", "1:1000" = "#BC3C29FF", "1:100000" = "#7876B1FF")) +
  scale_x_discrete(limits = c("1:1","1:10","1:100","1:1000","1:100000")) +
  theme_classic() +
  ggtitle("Olink") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ylim(0,0.75) +
  labs(y = "Sample missingness", x = "Dilution series")
rm(Class)
rm(Data)

rm(Rate_LoD_Sample1)
rm(Rate_LoD_Sample2)
rm(Rate_LoD_Sample3)
rm(Rate_LoD_Sample4)
rm(Rate_LoD_Sample5)

#Protein
Data <- data.frame(Rate = c(0.723932,0.08605852,0.06666667,0.0530303,0.1029412), Set = c("1:1","1:10","1:100","1:1000","1:100000"))
ggplot(Data, aes(x = Set, y = Rate, fill = Set, color = Set)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c("1:1" = "#0072B5FF", "1:10" = "#E18727FF", "1:100" = "#20854EFF", "1:1000" = "#BC3C29FF", "1:100000" = "#7876B1FF")) +
  scale_color_manual(values = c("1:1" = "#0072B5FF", "1:10" = "#E18727FF", "1:100" = "#20854EFF", "1:1000" = "#BC3C29FF", "1:100000" = "#7876B1FF")) +
  scale_x_discrete(limits = c("1:1","1:10","1:100","1:1000","1:100000")) +
  scale_y_continuous(limits = c(0,0.75)) +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(legend.position = "none") +
  ggtitle("Olink") +
  labs(y = "Proportion of proteins below LoD", x = "Dilution") +
  geom_text(aes(label = paste0(round(Rate * 100, 1), "%")), position = position_dodge(width = 1), vjust = -0.5, size = 4)
rm(Data)


#Proteins with proportion below LoD > 20%
Dat <- data.frame(Protein = colnames(ProtMatrix), Rate_LoD = Rate_LoD_Protein)
Dat <- Dat[Dat$Rate_LoD > 0.2,]
Dat <- cbind(Dat, Assay_Annotation[Dat$Protein,])
write.csv(Dat, "Olink_HighProportionProteins_0.2.csv", row.names = F)
rm(Dat)

rm(Rate_LoD_Protein)
rm(Rate_LoD_Sample)
rm(ProtMatrix_LOD)
gc()


#Remove All
rm(list = ls())
gc()

