#Correlation SomaLogic ANML


#Data Loading

#SomaLogic Data
Dat <- read.csv("Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:20)]
ProtMatrix <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation <- read.csv("Somalogic_Analyte_Annotation_All.csv")
rownames(Analyte_Annotation) <- Analyte_Annotation$AptName

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)


#Overall Correlation between Technical Replicates
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"),]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_SomaLogic_Tech <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Visualization
Class <- 1:length(Mean_Cor_SomaLogic_Tech)
Class[1:length(Mean_Cor_SomaLogic_Tech)] <- ""
Data <- data.frame(Mean_Cor_Tech = Mean_Cor_SomaLogic_Tech, Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor_Tech, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8, color = "#E18727FF", fill = "#E18727FF") +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Violin plot") +
  ggtitle("SomaLogic ANML Sample-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_SomaLogic_Tech)


#Stratify by LoD
LOD_SomaLogic <- read.csv("SomaLogic_LOD.csv")

#Tech 1
G1 <- which(LOD_SomaLogic$PropBelowLOD <= 0.2)
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), G1]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_SomaLogic_Tech_G1 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech_G1[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech_G1)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Tech 2
G2 <- which(LOD_SomaLogic$PropBelowLOD > 0.2)
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), G2]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_SomaLogic_Tech_G2 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech_G2[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech_G2)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Tech Visualization 2
Class <- 1:length(c(Mean_Cor_SomaLogic_Tech_G1, Mean_Cor_SomaLogic_Tech_G2))
Class[1:length(Mean_Cor_SomaLogic_Tech_G1)] <- "AboveLoD"
Class[(length(Mean_Cor_SomaLogic_Tech_G1) + 1):(length(Mean_Cor_SomaLogic_Tech_G1) + length(Mean_Cor_SomaLogic_Tech_G2))] <- "BelowLoD"
Data <- data.frame(Mean_Cor_Tech = c(Mean_Cor_SomaLogic_Tech_G1, Mean_Cor_SomaLogic_Tech_G2), Class = Class)
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
  ggtitle("SomaLogic ANML Sample-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_SomaLogic_Tech_G1)
rm(Mean_Cor_SomaLogic_Tech_G2)
rm(LOD_SomaLogic)


#Different Sets

#Set 1: SomaScan 5K
library(readxl)
ProteinList5K <- read_excel("SomaScan_Menu_1.3K_5K_7K.xlsx")
ProteinList5K <- ProteinList5K[-c(1:2),10]
ProteinList5K <- unique(ProteinList5K$...10)

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), which(Analyte_Annotation$SeqId %in% ProteinList5K)]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_SomaLogic_Tech1 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech1[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech1)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Set2 : SomaScan 7K - SomaScan 5K
ProteinList7K <- read_excel("SomaScan_Menu_1.3K_5K_7K.xlsx")
ProteinList7K <- ProteinList7K[-c(1:2),1]
ProteinList7K <- unique(ProteinList7K$`Supplementary Data 1. List of all SOMAmers in the plasma 7k, 5k, and 1.3k SomaScan assays.`)
ProteinList2 <- setdiff(ProteinList7K, ProteinList5K)

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), which(Analyte_Annotation$SeqId %in% ProteinList2)]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_SomaLogic_Tech2 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech2[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech2)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Set 3: SomaScan 11K - SomaScan 7K
ProteinList3 <- setdiff(Analyte_Annotation$SeqId, ProteinList7K)

ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), which(Analyte_Annotation$SeqId %in% ProteinList3)]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_SomaLogic_Tech3 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech3[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech3)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Visualization
Class <- 1:length(c(Mean_Cor_SomaLogic_Tech1, Mean_Cor_SomaLogic_Tech2, Mean_Cor_SomaLogic_Tech3))
Class[1:length(Mean_Cor_SomaLogic_Tech1)] <- "Set 1"
Class[(length(Mean_Cor_SomaLogic_Tech1) + 1):(length(Mean_Cor_SomaLogic_Tech1) + length(Mean_Cor_SomaLogic_Tech2))] <- "Set 2"
Class[(length(Mean_Cor_SomaLogic_Tech1) + length(Mean_Cor_SomaLogic_Tech2) + 1):(length(Mean_Cor_SomaLogic_Tech1) + length(Mean_Cor_SomaLogic_Tech2) + length(Mean_Cor_SomaLogic_Tech3))] <- "Set 3"
Data <- data.frame(Mean_Cor_SomaLogic_Tech = c(Mean_Cor_SomaLogic_Tech1, Mean_Cor_SomaLogic_Tech2, Mean_Cor_SomaLogic_Tech3), Class = Class)
Data$Class <- as.factor(Data$Class)

ggplot(data = Data, aes(x = Class, y = Mean_Cor_SomaLogic_Tech, fill = Class, color = Class)) +
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
  ggtitle("SomaLogic ANML Sample-Level") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Different sets")

Data$Normalization <- rep("ANML", 300)
write.csv(Data, "Correlation_Set_ANML_Sample.csv", row.names = F)

rm(Class)
rm(Data)

rm(Mean_Cor_SomaLogic_Tech1)
rm(Mean_Cor_SomaLogic_Tech2)
rm(Mean_Cor_SomaLogic_Tech3)
rm(ProteinList7K)


#Different Dilution Groups

#1:5
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), which(Analyte_Annotation$Dilution2 == 1/5)]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_SomaLogic_Tech1 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech1[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech1)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#1:200
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), which(Analyte_Annotation$Dilution2 == 1/200)]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_SomaLogic_Tech2 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech2[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech2)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#1:20000
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep_id == "Y"), which(Analyte_Annotation$Dilution2 == 1/20000)]
Label_Tech <- Label[which(Label$tech_rep_id == "Y"),]
ID <- unique(Label_Tech$FREG0_PID)

Mean_Cor_SomaLogic_Tech3 <- vector()
for (i in 1:100) {
  COR <- cor(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech3[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech3)

rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Visualization
Class <- 1:length(c(Mean_Cor_SomaLogic_Tech1, Mean_Cor_SomaLogic_Tech2, Mean_Cor_SomaLogic_Tech3))
Class[1:length(Mean_Cor_SomaLogic_Tech1)] <- "1:5"
Class[(length(Mean_Cor_SomaLogic_Tech1) + 1):(length(Mean_Cor_SomaLogic_Tech1) + length(Mean_Cor_SomaLogic_Tech2))] <- "1:200"
Class[(length(Mean_Cor_SomaLogic_Tech1) + length(Mean_Cor_SomaLogic_Tech2) + 1):(length(Mean_Cor_SomaLogic_Tech1) + length(Mean_Cor_SomaLogic_Tech2) + length(Mean_Cor_SomaLogic_Tech3))] <- "1:20000"
Data <- data.frame(Mean_Cor_SomaLogic_Tech = c(Mean_Cor_SomaLogic_Tech1, Mean_Cor_SomaLogic_Tech2, Mean_Cor_SomaLogic_Tech3), Class = Class)
Data$Class <- as.factor(Data$Class)

ggplot(data = Data, aes(x = Class, y = Mean_Cor_SomaLogic_Tech, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("1:5" = "#0072B5FF", "1:200" = "#E18727FF", "1:20000" = "#20854EFF")) +
  scale_color_manual(values = c("1:5" = "#0072B5FF", "1:200" = "#E18727FF", "1:20000" = "#20854EFF")) +
  scale_x_discrete(limits = c("1:5","1:200","1:20000")) +
  theme_classic() +
  ylim(0, 1) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ggtitle("SomaLogic ANML Sample-Level") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Dilution series")

Data$Normalization <- rep("ANML", 300)
write.csv(Data, "Correlation_Dilution_ANML_Sample.csv", row.names = F)

rm(Class)
rm(Data)

rm(Mean_Cor_SomaLogic_Tech1)
rm(Mean_Cor_SomaLogic_Tech2)
rm(Mean_Cor_SomaLogic_Tech3)


#Protein-Level Correlation between Technical Replicates

#SomaLogic
Mean_Cor_SomaLogic_Tech <- vector()
for (i in 1:10675) {
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
  
  Mean_Cor_SomaLogic_Tech[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_SomaLogic_Tech)

#Visualization
Class <- 1:10675
Class[1:10675] <- ""
library(ggplot2)
Data <- data.frame(Mean_Cor = Mean_Cor_SomaLogic_Tech, Class = Class)
Data$Class <- as.factor(Data$Class)
ggplot(data = Data, aes(x = Class, y = Mean_Cor)) +
  geom_violin(trim = TRUE, fill = "#E18727FF", linewidth = 0.8, colour = "#E18727FF") +
  geom_boxplot(fill = "white", width = 0.3, linewidth = 1, outlier.size = 1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  ggtitle("SomaLogic ANML Protein-Level") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Violin plot")
rm(Class)
rm(Data)
rm(Mean_Cor_SomaLogic_Tech)


#Stratify by LoD

#Above LOD
Mean_Cor_SomaLogic_Tech_G1 <- vector()
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
  
  Mean_Cor_SomaLogic_Tech_G1[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_SomaLogic_Tech_G1)

#Below LOD
Mean_Cor_SomaLogic_Tech_G2 <- vector()
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
  
  Mean_Cor_SomaLogic_Tech_G2[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_SomaLogic_Tech_G2)

#Visualization
Class <- 1:length(c(Mean_Cor_SomaLogic_Tech_G1, Mean_Cor_SomaLogic_Tech_G2))
Class[1:length(Mean_Cor_SomaLogic_Tech_G1)] <- "AboveLOD"
Class[(length(Mean_Cor_SomaLogic_Tech_G1) + 1):(length(Mean_Cor_SomaLogic_Tech_G1) + length(Mean_Cor_SomaLogic_Tech_G2))] <- "BelowLOD"
Data <- data.frame(Mean_Cor_Tech = c(Mean_Cor_SomaLogic_Tech_G1, Mean_Cor_SomaLogic_Tech_G2), Class = Class)
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
  ggtitle("SomaLogic ANML Protein-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_SomaLogic_Tech_G1)
rm(Mean_Cor_SomaLogic_Tech_G2)
rm(G1)
rm(G2)


#Different Sets

#Set 1
Mean_Cor_SomaLogic_Tech_G1 <- vector()
G1 <- which(Analyte_Annotation$SeqId %in% ProteinList5K)
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
  
  Mean_Cor_SomaLogic_Tech_G1[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(G1)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_SomaLogic_Tech_G1)

#Set 2
Mean_Cor_SomaLogic_Tech_G2 <- vector()
G2 <- which(Analyte_Annotation$SeqId %in% ProteinList2)
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
  
  Mean_Cor_SomaLogic_Tech_G2[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(G2)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_SomaLogic_Tech_G2)

#Set 3
Mean_Cor_SomaLogic_Tech_G3 <- vector()
G3 <- which(Analyte_Annotation$SeqId %in% ProteinList3)
for (i in 1:length(G3)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,G3[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_SomaLogic_Tech_G3[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(G3)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_SomaLogic_Tech_G3)

#Visualization
Class <- 1:length(c(Mean_Cor_SomaLogic_Tech_G1, Mean_Cor_SomaLogic_Tech_G2, Mean_Cor_SomaLogic_Tech_G3))
Class[1:length(Mean_Cor_SomaLogic_Tech_G1)] <- "Set 1"
Class[(length(Mean_Cor_SomaLogic_Tech_G1) + 1):(length(Mean_Cor_SomaLogic_Tech_G1) + length(Mean_Cor_SomaLogic_Tech_G2))] <- "Set 2"
Class[(length(Mean_Cor_SomaLogic_Tech_G1) + length(Mean_Cor_SomaLogic_Tech_G2) + 1):(length(Mean_Cor_SomaLogic_Tech_G1) + length(Mean_Cor_SomaLogic_Tech_G2) + length(Mean_Cor_SomaLogic_Tech_G3))] <- "Set 3"
Data <- data.frame(Mean_Cor = c(Mean_Cor_SomaLogic_Tech_G1, Mean_Cor_SomaLogic_Tech_G2, Mean_Cor_SomaLogic_Tech_G3), Class = Class)
Data$Class <- as.factor(Data$Class)

ggplot(data = Data, aes(x = Class, y = Mean_Cor, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_color_manual(values = c("Set 1" = "#0072B5FF", "Set 2" = "#E18727FF", "Set 3" = "#20854EFF")) +
  scale_x_discrete(limits = c("Set 1","Set 2","Set 3")) +
  theme_classic() +
  ylim(-0.2,1) +
  ggtitle("SomaLogic ANML Protein-Level") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Different sets")

Data$Normalization <- rep("ANML", 10675)
write.csv(Data, "Correlation_Set_ANML_Protein.csv", row.names = F)

rm(Class)
rm(Data)

rm(Mean_Cor_SomaLogic_Tech_G1)
rm(Mean_Cor_SomaLogic_Tech_G2)
rm(Mean_Cor_SomaLogic_Tech_G3)
rm(ProteinList5K)
rm(ProteinList2)
rm(ProteinList3)


#Dilution series

#1:5
Mean_Cor_SomaLogic_Tech_G1 <- vector()
G1 <- which(Analyte_Annotation$Dilution2 == 1/5)
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
  
  Mean_Cor_SomaLogic_Tech_G1[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(G1)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_SomaLogic_Tech_G1)

#1:200
Mean_Cor_SomaLogic_Tech_G2 <- vector()
G2 <- which(Analyte_Annotation$Dilution2 == 1/200)
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
  
  Mean_Cor_SomaLogic_Tech_G2[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(G2)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_SomaLogic_Tech_G2)

#1:20000
Mean_Cor_SomaLogic_Tech_G3 <- vector()
G3 <- which(Analyte_Annotation$Dilution2 == 1/20000)
for (i in 1:length(G3)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,G3[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec3 <- Vec[which(Vec$tech_rep == "3"),]
  rownames(Vec3) <- Vec3$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, intersect(Vec2$FREG0_PID, Vec3$FREG0_PID)),]
  Vec2 <- Vec2[rownames(Vec1),]
  Vec3 <- Vec3[rownames(Vec1),]
  
  Mean_Cor_SomaLogic_Tech_G3[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(G3)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_SomaLogic_Tech_G3)

#Visualization
Class <- 1:length(c(Mean_Cor_SomaLogic_Tech_G1, Mean_Cor_SomaLogic_Tech_G2, Mean_Cor_SomaLogic_Tech_G3))
Class[1:length(Mean_Cor_SomaLogic_Tech_G1)] <- "1:5"
Class[(length(Mean_Cor_SomaLogic_Tech_G1) + 1):(length(Mean_Cor_SomaLogic_Tech_G1) + length(Mean_Cor_SomaLogic_Tech_G2))] <- "1:200"
Class[(length(Mean_Cor_SomaLogic_Tech_G1) + length(Mean_Cor_SomaLogic_Tech_G2) + 1):(length(Mean_Cor_SomaLogic_Tech_G1) + length(Mean_Cor_SomaLogic_Tech_G2) + length(Mean_Cor_SomaLogic_Tech_G3))] <- "1:20000"
Data <- data.frame(Mean_Cor = c(Mean_Cor_SomaLogic_Tech_G1, Mean_Cor_SomaLogic_Tech_G2, Mean_Cor_SomaLogic_Tech_G3), Class = Class)
Data$Class <- as.factor(Data$Class)

ggplot(data = Data, aes(x = Class, y = Mean_Cor, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("1:5" = "#0072B5FF", "1:200" = "#E18727FF", "1:20000" = "#20854EFF")) +
  scale_color_manual(values = c("1:5" = "#0072B5FF", "1:200" = "#E18727FF", "1:20000" = "#20854EFF")) +
  scale_x_discrete(limits = c("1:5","1:200","1:20000")) +
  theme_classic() +
  ylim(-0.2,1) +
  ggtitle("SomaLogic ANML Protein-Level") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Dilution series")

Data$Normalization <- rep("ANML", 10675)
write.csv(Data, "Correlation_Dilution_ANML_Protein.csv", row.names = F)

rm(Class)
rm(Data)

rm(Mean_Cor_SomaLogic_Tech_G1)
rm(Mean_Cor_SomaLogic_Tech_G2)
rm(Mean_Cor_SomaLogic_Tech_G3)


rm(ProtMatrix)
rm(Label)
rm(Analyte_Annotation)
gc()


