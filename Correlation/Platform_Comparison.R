#Platform Comparison


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

rm(Analyte_Annotation)
rm(Label)
rm(ProtMatrix)


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

rm(Analyte_Annotation)
rm(Label)
rm(ProtMatrix)


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

rm(Assay_Annotation)
rm(Label)
rm(ProtMatrix)


#Thermo Fisher

#Thermo Fisher Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/4_Thermo_Fisher_Preprocessing/Thermo_Fisher_Merged_All_Imputed.csv")
rownames(Dat) <- Dat$biosample_id
Label <- Dat[,c(1:10)]
ProtMatrix <- Dat[,c(11:ncol(Dat))]
min(ProtMatrix, na.rm = T)
rm(Dat)

#Protein Annotation
Protein_Annotation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/4_Thermo_Fisher_Preprocessing/Thermo_Fisher_Protein_Annotation_All.csv")
rownames(Protein_Annotation) <- Protein_Annotation$ProteinID

#Overall Correlation between Technical Replicates
ID <- as.data.frame(table(Label$FREG0_PID))
ID <- as.character(ID$Var1[which(ID$Freq == 2)])

Mean_Cor_Thermo_Fisher_Tech <- vector()
for (i in 1:45) {
  COR <- cor(t(ProtMatrix[which(Label$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Thermo_Fisher_Tech[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Thermo_Fisher_Tech)

rm(ID)
rm(Label)
rm(ProtMatrix)
rm(Protein_Annotation)


#Visualization
Class <- 1:length(c(Mean_Cor_SomaLogic_Tech, Mean_Cor_SomaLogic_Tech2, Mean_Cor_Olink_Tech, Mean_Cor_Thermo_Fisher_Tech))
Class[1:length(Mean_Cor_SomaLogic_Tech)] <- "SomaLogic pre-ANML"
Class[(length(Mean_Cor_SomaLogic_Tech)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2))] <- "SomaLogic ANML"
Class[(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+length(Mean_Cor_Olink_Tech))] <- "Olink"
Class[(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+length(Mean_Cor_Olink_Tech)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+length(Mean_Cor_Olink_Tech)+length(Mean_Cor_Thermo_Fisher_Tech))] <- "Thermo Fisher"
Data <- data.frame(Mean_Cor = c(Mean_Cor_SomaLogic_Tech, Mean_Cor_SomaLogic_Tech2, Mean_Cor_Olink_Tech, Mean_Cor_Thermo_Fisher_Tech), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor, fill = Class, color = Class)) +
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
  ggtitle("Sample-Level") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Platform")
rm(Class)
rm(Data)
rm(Mean_Cor_SomaLogic_Tech)
rm(Mean_Cor_SomaLogic_Tech2)
rm(Mean_Cor_Olink_Tech)
rm(Mean_Cor_Thermo_Fisher_Tech)
gc()


#Platform Comparison Protein-Level

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

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)

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

rm(Analyte_Annotation)
rm(Label)
rm(ProtMatrix)


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

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)

#SomaLogic
Mean_Cor_SomaLogic_Tech2 <- vector()
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
  
  Mean_Cor_SomaLogic_Tech2[i] <- mean(c(cor(Vec1[,3], Vec2[,3], method = "pearson"), cor(Vec1[,3], Vec3[,3], method = "pearson"), cor(Vec2[,3], Vec3[,3], method = "pearson")))
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
rm(Vec3)
summary(Mean_Cor_SomaLogic_Tech2)

rm(Analyte_Annotation)
rm(Label)
rm(ProtMatrix)


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

rm(Assay_Annotation)
rm(Label)
rm(ProtMatrix)


#Thermo Fisher

#Thermo Fisher Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/4_Thermo_Fisher_Preprocessing/Thermo_Fisher_Merged_All_Imputed.csv")
rownames(Dat) <- Dat$biosample_id
Label <- Dat[,c(1:10)]
ProtMatrix <- Dat[,c(11:ncol(Dat))]
min(ProtMatrix, na.rm = T)
rm(Dat)

#Protein Annotation
Protein_Annotation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/4_Thermo_Fisher_Preprocessing/Thermo_Fisher_Protein_Annotation_All.csv")
rownames(Protein_Annotation) <- Protein_Annotation$ProteinID

#Thermo_Fisher
Mean_Cor_Thermo_Fisher_Tech <- vector()
for (i in 1:5945) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,i])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, Vec2$FREG0_PID),]
  Vec2 <- Vec2[rownames(Vec1),]
  
  Mean_Cor_Thermo_Fisher_Tech[i] <- cor(Vec1[,3], Vec2[,3], method = "pearson")
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
summary(Mean_Cor_Thermo_Fisher_Tech)

rm(Label)
rm(ProtMatrix)
rm(Protein_Annotation)


#Visualization
Class <- 1:length(c(Mean_Cor_SomaLogic_Tech, Mean_Cor_SomaLogic_Tech2, Mean_Cor_Olink_Tech, Mean_Cor_Thermo_Fisher_Tech))
Class[1:length(Mean_Cor_SomaLogic_Tech)] <- "SomaLogic pre-ANML"
Class[(length(Mean_Cor_SomaLogic_Tech)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2))] <- "SomaLogic ANML"
Class[(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+length(Mean_Cor_Olink_Tech))] <- "Olink"
Class[(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+length(Mean_Cor_Olink_Tech)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+length(Mean_Cor_Olink_Tech)+length(Mean_Cor_Thermo_Fisher_Tech))] <- "Thermo Fisher"
Data <- data.frame(Mean_Cor = c(Mean_Cor_SomaLogic_Tech, Mean_Cor_SomaLogic_Tech2, Mean_Cor_Olink_Tech, Mean_Cor_Thermo_Fisher_Tech), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor, fill = Class, color = Class)) +
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
  ggtitle("Protein-Level") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Platform")
rm(Class)
rm(Data)
rm(Mean_Cor_SomaLogic_Tech)
rm(Mean_Cor_SomaLogic_Tech2)
rm(Mean_Cor_Olink_Tech)
rm(Mean_Cor_Thermo_Fisher_Tech)
gc()


#pre-ANML & ANMl Stratification by Sets
pre <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/4_Correlation/4_Correlation_SomaLogic_pre-ANML/Correlation_Set_pre-ANML_Sample.csv")
ANML <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/4_Correlation/1_Correlation_SomaLogic/Correlation_Set_ANML_Sample.csv")
Data <- rbind(pre, ANML)

library(ggplot2)
ggplot(data = Data, aes(x = Normalization, y = Mean_Cor_SomaLogic_Tech, fill = Class, color = Class)) +
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
  ylim(0.1,1) +
  facet_wrap(~Class) +
  ggtitle("SomaLogic Sample-Level") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Normalization")
rm(pre)
rm(ANML)
rm(Data)

#pre-ANML & ANMl Stratification by Sets
pre <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/4_Correlation/4_Correlation_SomaLogic_pre-ANML/Correlation_Set_pre-ANML_Protein.csv")
ANML <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/4_Correlation/1_Correlation_SomaLogic/Correlation_Set_ANML_Protein.csv")
Data <- rbind(pre, ANML)

library(ggplot2)
ggplot(data = Data, aes(x = Normalization, y = Mean_Cor, fill = Class, color = Class)) +
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
  ylim(-0.2,1) +
  facet_wrap(~Class) +
  ggtitle("SomaLogic Protein-Level") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Normalization")
rm(pre)
rm(ANML)
rm(Data)


#pre-ANML & ANMl Stratification by Dilutions
pre <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/4_Correlation/4_Correlation_SomaLogic_pre-ANML/Correlation_Dilution_pre-ANML_Sample.csv")
ANML <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/4_Correlation/1_Correlation_SomaLogic/Correlation_Dilution_ANML_Sample.csv")
Data <- rbind(pre, ANML)
Data$Class[which(Data$Class == "1:5")] <- "1.1:5"
Data$Class[which(Data$Class == "1:200")] <- "2.1:200"
Data$Class[which(Data$Class == "1:20000")] <- "3.1:20000"

library(ggplot2)
ggplot(data = Data, aes(x = Normalization, y = Mean_Cor_SomaLogic_Tech, fill = Class, color = Class)) +
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
  ylim(0,1) +
  facet_wrap(~Class) +
  ggtitle("SomaLogic Sample-Level") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Normalization")
rm(pre)
rm(ANML)
rm(Data)

#pre-ANML & ANMl Stratification by Dilutions
pre <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/4_Correlation/4_Correlation_SomaLogic_pre-ANML/Correlation_Dilution_pre-ANML_Protein.csv")
ANML <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/4_Correlation/1_Correlation_SomaLogic/Correlation_Dilution_ANML_Protein.csv")
Data <- rbind(pre, ANML)
Data$Class[which(Data$Class == "1:5")] <- "1.1:5"
Data$Class[which(Data$Class == "1:200")] <- "2.1:200"
Data$Class[which(Data$Class == "1:20000")] <- "3.1:20000"

library(ggplot2)
ggplot(data = Data, aes(x = Normalization, y = Mean_Cor, fill = Class, color = Class)) +
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
  ylim(-0.2,1) +
  facet_wrap(~Class) +
  ggtitle("SomaLogic Protein-Level") +
  labs(y = "Mean Pearson correlation between technical replicates", x = "Normalization")
rm(pre)
rm(ANML)
rm(Data)

