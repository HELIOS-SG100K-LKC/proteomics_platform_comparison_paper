#Correlation Thermo Fisher


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

#Visualization
Class <- 1:length(Mean_Cor_Thermo_Fisher_Tech)
Class[1:length(Mean_Cor_Thermo_Fisher_Tech)] <- ""
Data <- data.frame(Mean_Cor_Tech = Mean_Cor_Thermo_Fisher_Tech, Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor_Tech, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8, color = "#BC3C29FF", fill = "#BC3C29FF") +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  labs(y = "Pearson correlation between technical replicates", x = "Violin plot") +
  ggtitle("Thermo Fisher Sample-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_Thermo_Fisher_Tech)


#Stratify by LoD
Missingness <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/3_Limit_of_Detection/3_Thermo_Fisher_Missingness/Thermo_Fisher_ProteinMissingness.csv")

#Tech 1
G1 <- which(Missingness$Missing_Rate_Protein <= 0.2)
ID <- as.data.frame(table(Label$FREG0_PID))
ID <- as.character(ID$Var1[which(ID$Freq == 2)])

Mean_Cor_Thermo_Fisher_Tech_G1 <- vector()
for (i in 1:45) {
  COR <- cor(t(ProtMatrix[which(Label$FREG0_PID == ID[i]), G1]), method = "pearson")
  Mean_Cor_Thermo_Fisher_Tech_G1[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Thermo_Fisher_Tech_G1)

rm(ID)

#Tech 2
G2 <- which(Missingness$Missing_Rate_Protein > 0.2)
ID <- as.data.frame(table(Label$FREG0_PID))
ID <- as.character(ID$Var1[which(ID$Freq == 2)])

Mean_Cor_Thermo_Fisher_Tech_G2 <- vector()
for (i in 1:45) {
  COR <- cor(t(ProtMatrix[which(Label$FREG0_PID == ID[i]), G2]), method = "pearson")
  Mean_Cor_Thermo_Fisher_Tech_G2[i] <- mean(COR[upper.tri(COR)])
}
rm(i)
rm(COR)
summary(Mean_Cor_Thermo_Fisher_Tech_G2)

rm(ID)

#Visualization
Class <- 1:length(c(Mean_Cor_Thermo_Fisher_Tech_G1, Mean_Cor_Thermo_Fisher_Tech_G2))
Class[1:length(Mean_Cor_Thermo_Fisher_Tech_G1)] <- "AboveLoD"
Class[(length(Mean_Cor_Thermo_Fisher_Tech_G1) + 1):(length(Mean_Cor_Thermo_Fisher_Tech_G1) + length(Mean_Cor_Thermo_Fisher_Tech_G2))] <- "BelowLoD"
Data <- data.frame(Mean_Cor_Tech = c(Mean_Cor_Thermo_Fisher_Tech_G1, Mean_Cor_Thermo_Fisher_Tech_G2), Class = Class)
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
  labs(y = "Pearson correlation between technical replicates", x = "Proteins used for calculation") +
  ggtitle("Thermo Fisher Sample-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_Thermo_Fisher_Tech_G1)
rm(Mean_Cor_Thermo_Fisher_Tech_G2)
rm(Missingness)


#Protein-Level

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

#Visualization
Class <- 1:5945
Class[1:5945] <- ""
library(ggplot2)
Data <- data.frame(Mean_Cor = Mean_Cor_Thermo_Fisher_Tech, Class = Class)
Data$Class <- as.factor(Data$Class)
ggplot(data = Data, aes(x = Class, y = Mean_Cor)) +
  geom_violin(trim = TRUE, fill = "#BC3C29FF", linewidth = 0.8, colour = "#BC3C29FF") +
  geom_boxplot(fill = "white", width = 0.3, linewidth = 1, outlier.size = 1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  ggtitle("Thermo Fisher Protein-Level") +
  labs(y = "Pearson correlation between technical replicates", x = "Violin plot")
rm(Class)
rm(Data)
rm(Mean_Cor_Thermo_Fisher_Tech)


#Stratify by LoD

#Above LoD
Mean_Cor_Thermo_Fisher_Tech_G1 <- vector()
for (i in 1:length(G1)) {
  
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,G1[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, Vec2$FREG0_PID),]
  Vec2 <- Vec2[rownames(Vec1),]

  Mean_Cor_Thermo_Fisher_Tech_G1[i] <- cor(Vec1[,3], Vec2[,3], method = "pearson")
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
summary(Mean_Cor_Thermo_Fisher_Tech_G1)

#Below LoD
Mean_Cor_Thermo_Fisher_Tech_G2 <- vector()
for (i in 1:length(G2)) {
  Vec <- cbind(Label[,c("FREG0_PID","tech_rep")], ProtMatrix[,G2[i]])
  
  Vec1 <- Vec[which(Vec$tech_rep == "1"),]
  rownames(Vec1) <- Vec1$FREG0_PID
  
  Vec2 <- Vec[which(Vec$tech_rep == "2"),]
  rownames(Vec2) <- Vec2$FREG0_PID
  
  Vec1 <- Vec1[intersect(Vec1$FREG0_PID, Vec2$FREG0_PID),]
  Vec2 <- Vec2[rownames(Vec1),]

  Mean_Cor_Thermo_Fisher_Tech_G2[i] <- cor(Vec1[,3], Vec2[,3], method = "pearson")
}
rm(i)
rm(Vec)
rm(Vec1)
rm(Vec2)
summary(Mean_Cor_Thermo_Fisher_Tech_G2)

#Visualization
Class <- 1:length(c(Mean_Cor_Thermo_Fisher_Tech_G1, Mean_Cor_Thermo_Fisher_Tech_G2))
Class[1:length(Mean_Cor_Thermo_Fisher_Tech_G1)] <- "AboveLOD"
Class[(length(Mean_Cor_Thermo_Fisher_Tech_G1) + 1):(length(Mean_Cor_Thermo_Fisher_Tech_G1) + length(Mean_Cor_Thermo_Fisher_Tech_G2))] <- "BelowLOD"
Data <- data.frame(Mean_Cor_Tech = c(Mean_Cor_Thermo_Fisher_Tech_G1, Mean_Cor_Thermo_Fisher_Tech_G2), Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Mean_Cor_Tech, fill = Class, color = Class)) +
  geom_violin(trim = TRUE, linewidth = 0.8) +
  geom_boxplot(fill = "white", width = 0.25, linewidth = 0.5, outlier.size = 1, color = "black") +
  scale_fill_manual(values = c("AboveLOD" = "#0072B5FF", "BelowLOD" = "#E18727FF")) +
  scale_color_manual(values = c("AboveLOD" = "#0072B5FF", "BelowLOD" = "#E18727FF")) +
  scale_x_discrete(limits = c("AboveLOD","BelowLOD")) +
  theme_classic() +
  ylim(-0.4,1) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  theme(legend.position = "none") +
  labs(y = "Pearson correlation between technical replicates", x = "Proteins used for calculation") +
  ggtitle("Thermo Fisher Protein-Level")
rm(Class)
rm(Data)
rm(Mean_Cor_Thermo_Fisher_Tech_G1)
rm(Mean_Cor_Thermo_Fisher_Tech_G2)
rm(G1)
rm(G2)


rm(ProtMatrix)
rm(Label)
rm(Protein_Annotation)
gc()

