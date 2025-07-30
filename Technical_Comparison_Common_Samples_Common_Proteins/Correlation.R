#Correlation


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
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/5_Common_Samples_Common_Proteins/1_Preprocessing/Thermo_Fisher_Merged_All_UniProt_Imputed.csv")
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


#Sample-Level

#ID
ID <- as.data.frame(table(Label_TF$FREG0_PID))
ID <- as.character(ID$Var1[which(ID$Freq == 2)])

#SomaLogic pre-ANML

#Overall Correlation between Technical Replicates
Mean_Cor_SomaLogic_Tech <- vector()
for (i in 1:45) {
  COR <- cor(t(ProtMatrix_SomaLogic[which(Label_SomaLogic$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech)

#SomaLogic ANML

#Overall Correlation between Technical Replicates
Mean_Cor_SomaLogic_Tech2 <- vector()
for (i in 1:45) {
  COR <- cor(t(ProtMatrix_SomaLogic2[which(Label_SomaLogic2$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_SomaLogic_Tech2[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_SomaLogic_Tech2)

#Olink

#Overall Correlation between Technical Replicates
Mean_Cor_Olink_Tech <- vector()
for (i in 1:45) {
  COR <- cor(t(ProtMatrix_Olink[which(Label_Olink$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Olink_Tech[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Olink_Tech)

#Thermo

#Overall Correlation between Technical Replicates
Mean_Cor_Thermo_Fisher_Tech <- vector()
for (i in 1:45) {
  COR <- cor(t(ProtMatrix_TF[which(Label_TF$FREG0_PID == ID[i]),]), method = "pearson")
  Mean_Cor_Thermo_Fisher_Tech[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Thermo_Fisher_Tech)

rm(ID)

#Visualization
Class <- 1:length(c(Mean_Cor_SomaLogic_Tech, Mean_Cor_SomaLogic_Tech2, Mean_Cor_Olink_Tech, Mean_Cor_Thermo_Fisher_Tech))
Class[1:length(Mean_Cor_SomaLogic_Tech)] <- "SomaLogic pre-ANML"
Class[(length(Mean_Cor_SomaLogic_Tech)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2))] <- "SomaLogic ANML"
Class[(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+length(Mean_Cor_Olink_Tech))] <- "Olink"
Class[(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+length(Mean_Cor_Olink_Tech)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic_Tech2)+length(Mean_Cor_Olink_Tech)+length(Mean_Cor_Thermo_Fisher_Tech))] <- "Thermo"

Data <- data.frame(Mean_Cor = c(Mean_Cor_SomaLogic_Tech, Mean_Cor_SomaLogic_Tech2, Mean_Cor_Olink_Tech, Mean_Cor_Thermo_Fisher_Tech), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("SomaLogic pre-ANML" = "#0072B5FF", "SomaLogic ANML" = "#E18727FF", "Olink" = "#20854EFF", "Thermo" = "#BC3C29FF")) +
  scale_color_manual(values = c("SomaLogic pre-ANML" = "#0072B5FF", "SomaLogic ANML" = "#E18727FF", "Olink" = "#20854EFF", "Thermo" = "#BC3C29FF")) +
  scale_x_discrete(limits = c("SomaLogic pre-ANML","SomaLogic ANML","Olink","Thermo")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ggtitle("Sample-Level") +
  labs(y = "Pearson correlation between technical replicates", x = "Platform")
rm(Class)
rm(Data)
rm(Mean_Cor_SomaLogic_Tech)
rm(Mean_Cor_SomaLogic_Tech2)
rm(Mean_Cor_Olink_Tech)
rm(Mean_Cor_Thermo_Fisher_Tech)
gc()


#Protein-Level

#SomaLogic pre-ANML
Mean_Cor_SomaLogic_Tech <- vector()
for (i in 1:1740) {
  Vec <- cbind(Label_SomaLogic[,c("FREG0_PID","tech_rep")], ProtMatrix_SomaLogic[,i])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, Vec2$FREG0_PID),]
  Vec2 <- Vec2[rownames(Vec1),]
  
  Mean_Cor_SomaLogic_Tech[i] <- cor(Vec1[,3], Vec2[,3], method = "pearson")
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
summary(Mean_Cor_SomaLogic_Tech)

rm(Analyte_Annotation_SomaLogic)
rm(Label_SomaLogic)
rm(ProtMatrix_SomaLogic)

#SomaLogic ANML
Mean_Cor_SomaLogic2_Tech <- vector()
for (i in 1:1740) {
  Vec <- cbind(Label_SomaLogic2[,c("FREG0_PID","tech_rep")], ProtMatrix_SomaLogic2[,i])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, Vec2$FREG0_PID),]
  Vec2 <- Vec2[rownames(Vec1),]
  
  Mean_Cor_SomaLogic2_Tech[i] <- cor(Vec1[,3], Vec2[,3], method = "pearson")
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
summary(Mean_Cor_SomaLogic2_Tech)

rm(Analyte_Annotation_SomaLogic2)
rm(Label_SomaLogic2)
rm(ProtMatrix_SomaLogic2)

#Olink
Mean_Cor_Olink_Tech <- vector()
for (i in 1:1740) {
  Vec <- cbind(Label_Olink[,c("FREG0_PID","tech_rep")], ProtMatrix_Olink[,i])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, Vec2$FREG0_PID),]
  Vec2 <- Vec2[rownames(Vec1),]
  
  Mean_Cor_Olink_Tech[i] <- cor(Vec1[,3], Vec2[,3], method = "pearson")
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
summary(Mean_Cor_Olink_Tech)

rm(Assay_Annotation_Olink)
rm(Label_Olink)
rm(ProtMatrix_Olink)

#Thermo_Fisher
Mean_Cor_Thermo_Fisher_Tech <- vector()
for (i in 1:1740) {
  Vec <- cbind(Label_TF[,c("FREG0_PID","tech_rep")], ProtMatrix_TF[,i])
  
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

rm(Protein_Annotation_TF)
rm(Label_TF)
rm(ProtMatrix_TF)

#Visualization
Class <- 1:length(c(Mean_Cor_SomaLogic_Tech, Mean_Cor_SomaLogic2_Tech, Mean_Cor_Olink_Tech, Mean_Cor_Thermo_Fisher_Tech))
Class[1:length(Mean_Cor_SomaLogic_Tech)] <- "SomaLogic pre-ANML"
Class[(length(Mean_Cor_SomaLogic_Tech)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic2_Tech))] <- "SomaLogic ANML"
Class[(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic2_Tech)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic2_Tech)+length(Mean_Cor_Olink_Tech))] <- "Olink"
Class[(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic2_Tech)+length(Mean_Cor_Olink_Tech)+1):(length(Mean_Cor_SomaLogic_Tech)+length(Mean_Cor_SomaLogic2_Tech)+length(Mean_Cor_Olink_Tech)+length(Mean_Cor_Thermo_Fisher_Tech))] <- "Thermo"

Data <- data.frame(Mean_Cor = c(Mean_Cor_SomaLogic_Tech, Mean_Cor_SomaLogic2_Tech, Mean_Cor_Olink_Tech, Mean_Cor_Thermo_Fisher_Tech), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("SomaLogic pre-ANML" = "#0072B5FF", "SomaLogic ANML" = "#E18727FF", "Olink" = "#20854EFF", "Thermo" = "#BC3C29FF")) +
  scale_color_manual(values = c("SomaLogic pre-ANML" = "#0072B5FF", "SomaLogic ANML" = "#E18727FF", "Olink" = "#20854EFF", "Thermo" = "#BC3C29FF")) +
  scale_x_discrete(limits = c("SomaLogic pre-ANML","SomaLogic ANML","Olink","Thermo")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  ggtitle("Protein-Level") +
  labs(y = "Pearson correlation between technical replicates", x = "Platform")
rm(Class)
rm(Data)
rm(Mean_Cor_SomaLogic_Tech)
rm(Mean_Cor_SomaLogic2_Tech)
rm(Mean_Cor_Olink_Tech)
rm(Mean_Cor_Thermo_Fisher_Tech)
gc()

