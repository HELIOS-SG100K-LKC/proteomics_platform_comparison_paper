#Correlation Olink


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
ProtMatrix <- ProtMatrix[which(Label$SampleType == "SAMPLE"),]
Label <- Label[which(Label$SampleType == "SAMPLE"),]


#Overall Correlation between Technical Replicates
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"),]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Visualization
Class <- 1:length(Mean_Cor_Olink_Tech)
Class[1:length(Mean_Cor_Olink_Tech)] <- ""
Data <- data.frame(Mean_Cor_Tech = Mean_Cor_Olink_Tech, Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor_Tech, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8, color = "#20854EFF", fill = "#20854EFF") +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Violin plot") +
  ggtitle("Olink Sample-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_Olink_Tech)


#Stratify by LoD
LOD_Olink <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/3_Limit_of_Detection/2_Olink_LoD/Olink_LOD.csv")

#Tech 1
G1 <- which(LOD_Olink$PropBelowLOD <= 0.2)
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), G1]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech_G1 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech_G1[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech_G1)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Tech 2
G2 <- which(LOD_Olink$PropBelowLOD > 0.2)
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), G2]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech_G2 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech_G2[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech_G2)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Visualization
Class <- 1:length(c(Mean_Cor_Olink_Tech_G1, Mean_Cor_Olink_Tech_G2))
Class[1:length(Mean_Cor_Olink_Tech_G1)] <- "AboveLoD"
Class[(length(Mean_Cor_Olink_Tech_G1) + 1):(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2))] <- "BelowLoD"
Data <- data.frame(Mean_Cor_Tech = c(Mean_Cor_Olink_Tech_G1, Mean_Cor_Olink_Tech_G2), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor_Tech, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("AboveLoD" = "#0072B5FF", "BelowLoD" = "#E18727FF")) +
  scale_color_manual(values = c("AboveLoD" = "#0072B5FF", "BelowLoD" = "#E18727FF")) +
  scale_x_discrete(limits = c("AboveLoD","BelowLoD")) +
  theme_classic() +
  ylim(0.15, 1) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Proteins used for calculation") +
  ggtitle("Olink Sample-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_Olink_Tech_G1)
rm(Mean_Cor_Olink_Tech_G2)
rm(LOD_Olink)


#Different Sets

#Set 1: Olink Explore 1536
library(readxl)
ProteinList1536 <- read_excel("explore-1536-assay-list-20210227-web-1.xlsx")
colnames(ProteinList1536) <- ProteinList1536[1,]
ProteinList1536 <- ProteinList1536[-1,]
V1 <- which(Assay_Annotation$UniProt %in% ProteinList1536$`Uniprot ID`)

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), V1]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech_V1 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech_V1[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech_V1)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Set 2: Olink Explore 3072 - Olink Explore 1536
ProteinList3072 <- read_excel("olink-explore-3072-validation-data-results.xlsx")
colnames(ProteinList3072) <- ProteinList3072[1,]
ProteinList3072 <- ProteinList3072[-1,]
ProteinList3072 <- ProteinList3072[-which(ProteinList3072$UniProt %in% ProteinList1536$`Uniprot ID`),]
V2 <- which(Assay_Annotation$UniProt %in% ProteinList3072$UniProt)

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), V2]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech_V2 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech_V2[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech_V2)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Set 3: Olink Explore HT - Olink Explore 3072
ProteinList3072 <- read_excel("olink-explore-3072-validation-data-results.xlsx")
colnames(ProteinList3072) <- ProteinList3072[1,]
ProteinList3072 <- ProteinList3072[-1,]
V3 <- setdiff(1:5405, which(Assay_Annotation$UniProt %in% ProteinList3072$UniProt))

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), V3]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech_V3 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech_V3[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech_V3)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Visualization
Class <- 1:length(c(Mean_Cor_Olink_Tech_V1, Mean_Cor_Olink_Tech_V2, Mean_Cor_Olink_Tech_V3))
Class[1:length(Mean_Cor_Olink_Tech_V1)] <- "Set 1"
Class[(length(Mean_Cor_Olink_Tech_V1) + 1):(length(Mean_Cor_Olink_Tech_V1) + length(Mean_Cor_Olink_Tech_V2))] <- "Set 2"
Class[(length(Mean_Cor_Olink_Tech_V1) + length(Mean_Cor_Olink_Tech_V2) + 1):(length(Mean_Cor_Olink_Tech_V1) + length(Mean_Cor_Olink_Tech_V2) + length(Mean_Cor_Olink_Tech_V3))] <- "Set 3"
Data <- data.frame(Mean_Cor_Tech = c(Mean_Cor_Olink_Tech_V1, Mean_Cor_Olink_Tech_V2, Mean_Cor_Olink_Tech_V3), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor_Tech, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_color_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_x_discrete(limits = c("Set 1","Set 2","Set 3")) +
  theme_classic() +
  ylim(0.1, 1) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Different sets") +
  ggtitle("Olink Sample-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_Olink_Tech_V1)
rm(Mean_Cor_Olink_Tech_V2)
rm(Mean_Cor_Olink_Tech_V3)

rm(ProteinList1536)
rm(ProteinList3072)


#Dilution series
Assay_Annotation$Dilution <- Assay_Annotation$Block
Assay_Annotation$Dilution[which(Assay_Annotation$Block %in% c(1,2,3,4))] <- "1:1"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 5)] <- "1:10"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 6)] <- "1:100"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 7)] <- "1:1000"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 8)] <- "1:100000"

#1:1
D1 <- which(Assay_Annotation$Dilution == "1:1")

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), D1]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech_V1 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech_V1[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech_V1)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#1:10
D2 <- which(Assay_Annotation$Dilution == "1:10")

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), D2]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech_V2 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech_V2[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech_V2)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#1:100
D3 <- which(Assay_Annotation$Dilution == "1:100")

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), D3]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech_V3 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech_V3[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech_V3)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#1:1000
D4 <- which(Assay_Annotation$Dilution == "1:1000")

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), D4]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech_V4 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech_V4[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech_V4)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#1:100000
D5 <- which(Assay_Annotation$Dilution == "1:100000")

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), D5]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_Olink_Tech_V5 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech_V5[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech_V5)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Visualization
Class <- 1:length(c(Mean_Cor_Olink_Tech_V1, Mean_Cor_Olink_Tech_V2, Mean_Cor_Olink_Tech_V3, Mean_Cor_Olink_Tech_V4, Mean_Cor_Olink_Tech_V5))
Class[1:length(Mean_Cor_Olink_Tech_V1)] <- "1:1"
Class[(length(Mean_Cor_Olink_Tech_V1) + 1):(length(Mean_Cor_Olink_Tech_V1) + length(Mean_Cor_Olink_Tech_V2))] <- "1:10"
Class[(length(Mean_Cor_Olink_Tech_V1) + length(Mean_Cor_Olink_Tech_V2) + 1):(length(Mean_Cor_Olink_Tech_V1) + length(Mean_Cor_Olink_Tech_V2) + length(Mean_Cor_Olink_Tech_V3))] <- "1:100"
Class[(length(Mean_Cor_Olink_Tech_V1) + length(Mean_Cor_Olink_Tech_V2) + length(Mean_Cor_Olink_Tech_V3) + 1):(length(Mean_Cor_Olink_Tech_V1) + length(Mean_Cor_Olink_Tech_V2) + length(Mean_Cor_Olink_Tech_V3) + length(Mean_Cor_Olink_Tech_V4))] <- "1:1000"
Class[(length(Mean_Cor_Olink_Tech_V1) + length(Mean_Cor_Olink_Tech_V2) + length(Mean_Cor_Olink_Tech_V3) + length(Mean_Cor_Olink_Tech_V4) + 1):(length(Mean_Cor_Olink_Tech_V1) + length(Mean_Cor_Olink_Tech_V2) + length(Mean_Cor_Olink_Tech_V3) + length(Mean_Cor_Olink_Tech_V4) + length(Mean_Cor_Olink_Tech_V5))] <- "1:100000"
Data <- data.frame(Mean_Cor_Tech = c(Mean_Cor_Olink_Tech_V1, Mean_Cor_Olink_Tech_V2, Mean_Cor_Olink_Tech_V3, Mean_Cor_Olink_Tech_V4, Mean_Cor_Olink_Tech_V5), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor_Tech, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("1:1" = "#0072B5FF", "1:10" = "#E18727FF", "1:100" = "#20854EFF", "1:1000" = "#BC3C29FF", "1:100000" = "#7876B1FF")) +
  scale_color_manual(values = c("1:1" = "#0072B5FF", "1:10" = "#E18727FF", "1:100" = "#20854EFF", "1:1000" = "#BC3C29FF", "1:100000" = "#7876B1FF")) +
  scale_x_discrete(limits = c("1:1","1:10","1:100","1:1000","1:100000")) +
  theme_classic() +
  ylim(0, 1) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Dilution series") +
  ggtitle("Olink Sample-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_Olink_Tech_V1)
rm(Mean_Cor_Olink_Tech_V2)
rm(Mean_Cor_Olink_Tech_V3)
rm(Mean_Cor_Olink_Tech_V4)
rm(Mean_Cor_Olink_Tech_V5)


#Protein-Level

#Olink
Mean_Cor_Olink_Tech <- vector()
for (i in 1:5405) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,i])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech)

#Visualization
Class <- 1:5405
Class[1:5405] <- ""
library(ggplot2)
Data <- data.frame(Mean_Cor = Mean_Cor_Olink_Tech, Class = Class)
Data$Class <- as.factor(Data$Class)
ggplot(data = Data, aes(x = Class, y = Mean_Cor)) +
  geom_violin(trim = TRUE, fill = "#20854EFF", linewidth = 0.8, colour = "#20854EFF") +
  geom_boxplot(fill = "white", width = 0.3, linewidth = 1, outlier.size = 1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  ggtitle("Olink Protein-Level") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Violin plot")
rm(Class)
rm(Data)
rm(Mean_Cor_Olink_Tech)


#Stratify by LoD

#Above LOD
Mean_Cor_Olink_Tech_G1 <- vector()
for (i in 1:length(G1)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,G1[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech_G1[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech_G1)

#Below LOD
Mean_Cor_Olink_Tech_G2 <- vector()
for (i in 1:length(G2)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,G2[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech_G2[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech_G2)

#Visualization
Class <- 1:length(c(Mean_Cor_Olink_Tech_G1, Mean_Cor_Olink_Tech_G2))
Class[1:length(Mean_Cor_Olink_Tech_G1)] <- "AboveLOD"
Class[(length(Mean_Cor_Olink_Tech_G1) + 1):(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2))] <- "BelowLOD"
Data <- data.frame(Mean_Cor_Tech = c(Mean_Cor_Olink_Tech_G1, Mean_Cor_Olink_Tech_G2), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor_Tech, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("AboveLOD" = "#0072B5FF", "BelowLOD" = "#E18727FF")) +
  scale_color_manual(values = c("AboveLOD" = "#0072B5FF", "BelowLOD" = "#E18727FF")) +
  scale_x_discrete(limits = c("AboveLOD","BelowLOD")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ylim(-0.4,1) +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Proteins used for calculation") +
  ggtitle("Olink Protein-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_Olink_Tech_G1)
rm(Mean_Cor_Olink_Tech_G2)
rm(G1)
rm(G2)


#Different Sets

#Set 1
Mean_Cor_Olink_Tech_G1 <- vector()
for (i in 1:length(V1)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,V1[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech_G1[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech_G1)

#Set 2
Mean_Cor_Olink_Tech_G2 <- vector()
for (i in 1:length(V2)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,V2[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech_G2[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech_G2)

#Set 3
Mean_Cor_Olink_Tech_G3 <- vector()
for (i in 1:length(V3)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,V3[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech_G3[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech_G3)

#Visualization
Class <- 1:length(c(Mean_Cor_Olink_Tech_G1, Mean_Cor_Olink_Tech_G2, Mean_Cor_Olink_Tech_G3))
Class[1:length(Mean_Cor_Olink_Tech_G1)] <- "Set 1"
Class[(length(Mean_Cor_Olink_Tech_G1) + 1):(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2))] <- "Set 2"
Class[(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2) + 1):(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2) + length(Mean_Cor_Olink_Tech_G3))] <- "Set 3"
Data <- data.frame(Mean_Cor = c(Mean_Cor_Olink_Tech_G1, Mean_Cor_Olink_Tech_G2, Mean_Cor_Olink_Tech_G3), Class = Class)
Data$Class <- as.factor(Data$Class)
ggplot(data = Data, aes(x = Class, y = Mean_Cor, fill = Class, color = Class)) +
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
  ggtitle("Olink Protein-Level") +
  ylim(-0.2,1) +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Different sets")
rm(Class)
rm(Data)

rm(Mean_Cor_Olink_Tech_G1)
rm(Mean_Cor_Olink_Tech_G2)
rm(Mean_Cor_Olink_Tech_G3)
rm(V1)
rm(V2)
rm(V3)


#Dilution series

#1:1
Mean_Cor_Olink_Tech_G1 <- vector()
for (i in 1:length(D1)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,D1[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech_G1[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech_G1)

#1:10
Mean_Cor_Olink_Tech_G2 <- vector()
for (i in 1:length(D2)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,D2[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech_G2[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech_G2)

#1:100
Mean_Cor_Olink_Tech_G3 <- vector()
for (i in 1:length(D3)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,D3[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech_G3[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech_G3)

#1:1000
Mean_Cor_Olink_Tech_G4 <- vector()
for (i in 1:length(D4)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,D4[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech_G4[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech_G4)

#1:100000
Mean_Cor_Olink_Tech_G5 <- vector()
for (i in 1:length(D5)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,D5[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech_G5[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_Olink_Tech_G5)

#Visualization
Class <- 1:length(c(Mean_Cor_Olink_Tech_G1, Mean_Cor_Olink_Tech_G2, Mean_Cor_Olink_Tech_G3, Mean_Cor_Olink_Tech_G4, Mean_Cor_Olink_Tech_G5))
Class[1:length(Mean_Cor_Olink_Tech_G1)] <- "1:1"
Class[(length(Mean_Cor_Olink_Tech_G1) + 1):(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2))] <- "1:10"
Class[(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2) + 1):(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2) + length(Mean_Cor_Olink_Tech_G3))] <- "1:100"
Class[(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2) + length(Mean_Cor_Olink_Tech_G3) + 1):(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2) + length(Mean_Cor_Olink_Tech_G3) + length(Mean_Cor_Olink_Tech_G4))] <- "1:1000"
Class[(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2) + length(Mean_Cor_Olink_Tech_G3) + length(Mean_Cor_Olink_Tech_G4) + 1):(length(Mean_Cor_Olink_Tech_G1) + length(Mean_Cor_Olink_Tech_G2) + length(Mean_Cor_Olink_Tech_G3) + length(Mean_Cor_Olink_Tech_G4) + length(Mean_Cor_Olink_Tech_G5))] <- "1:100000"
Data <- data.frame(Mean_Cor = c(Mean_Cor_Olink_Tech_G1, Mean_Cor_Olink_Tech_G2, Mean_Cor_Olink_Tech_G3, Mean_Cor_Olink_Tech_G4, Mean_Cor_Olink_Tech_G5), Class = Class)
Data$Class <- as.factor(Data$Class)
ggplot(data = Data, aes(x = Class, y = Mean_Cor, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("1:1" = "#0072B5FF", "1:10" = "#E18727FF", "1:100" = "#20854EFF", "1:1000" = "#BC3C29FF", "1:100000" = "#7876B1FF")) +
  scale_color_manual(values = c("1:1" = "#0072B5FF", "1:10" = "#E18727FF", "1:100" = "#20854EFF", "1:1000" = "#BC3C29FF", "1:100000" = "#7876B1FF")) +
  scale_x_discrete(limits = c("1:1","1:10","1:100","1:1000","1:100000")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ggtitle("Olink Protein-Level") +
  ylim(-0.2,1) +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Dilution series")
rm(Class)
rm(Data)

rm(Mean_Cor_Olink_Tech_G1)
rm(Mean_Cor_Olink_Tech_G2)
rm(Mean_Cor_Olink_Tech_G3)
rm(Mean_Cor_Olink_Tech_G4)
rm(Mean_Cor_Olink_Tech_G5)
rm(D1)
rm(D2)
rm(D3)
rm(D4)
rm(D5)


rm(Assay_Annotation)
rm(Label)
rm(ProtMatrix)
gc()

