#Overall Outlier Detection


#SomaLogic pre-ANML

#SomaLogic Data
Dat <- read.csv("SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:20)]
ProtMatrix <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)

#Outlier Detection

#PCA
set.seed(1)
PCA <- prcomp(ProtMatrix)
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PC_Mean <- apply(PCA[,c(1:2)], 2, mean)
PC_SD <- apply(PCA[,c(1:2)], 2, sd)
library(ggplot2)
p1 <- ggplot(PCA, aes(x = PC1, y = PC2)) +
  geom_point(size = 2.5, alpha = 1, color = "#0072B5FF", fill = "#0072B5FF") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  labs(x = "PC1 (24.3%)", y = "PC2 (5.5%)") +
  ggtitle("SomaLogic pre-ANML") +
  geom_vline(xintercept = PC_Mean[1] + 5 * PC_SD[1], linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = PC_Mean[1] - 5 * PC_SD[1], linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = PC_Mean[2] + 5 * PC_SD[2], linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = PC_Mean[2] - 5 * PC_SD[2], linetype = "dashed", color = "black", linewidth = 1)
p1

which(PCA$PC1 > PC_Mean[1] + 5 * PC_SD[1])
which(PCA$PC1 < PC_Mean[1] - 5 * PC_SD[1])

which(PCA$PC2 > PC_Mean[2] + 5 * PC_SD[2])
which(PCA$PC2 < PC_Mean[2] - 5 * PC_SD[2])

rm(PC_Mean)
rm(PC_SD)
rm(PCA)

#IQR
IQR <- apply(ProtMatrix, 1, IQR)
which(IQR > mean(IQR) + 5 * sd(IQR))
which(IQR < mean(IQR) - 5 * sd(IQR))
B <- names(which(IQR < mean(IQR) - 5 * sd(IQR)))
B
#"Plate 3_D4"

#Median
Median <- apply(ProtMatrix, 1, median)
which(Median > mean(Median) + 5 * sd(Median))
which(Median < mean(Median) - 5 * sd(Median))

#Visualization
IQR_Median <- data.frame(IQR, Median)
library(ggplot2)
p2 <- ggplot(IQR_Median, aes(x = Median, y = IQR)) +
  geom_point(size = 2.5, alpha = 1, color = "#0072B5FF", fill = "#0072B5FF") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  labs(x = "Median", y = "IQR") +
  ggtitle("SomaLogic pre-ANML") +
  geom_vline(xintercept = mean(Median) + 5 * sd(Median), linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = mean(Median) - 5 * sd(Median), linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = mean(IQR) + 5 * sd(IQR), linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = mean(IQR) - 5 * sd(IQR), linetype = "dashed", color = "black", linewidth = 1)
p2

rm(IQR)
rm(Median)
rm(IQR_Median)

#Correlation
Cor_Matrix <- cor(t(ProtMatrix), method = "pearson")
Mean_Cor <- vector()
for (i in 1:ncol(Cor_Matrix)) {
  Mean_Cor[i] <- mean(Cor_Matrix[i,][-i])
}
rm(i)
rm(Cor_Matrix)
summary(Mean_Cor)

#Visualization
Class <- 1:800
Class[1:800] <- ""
Data <- data.frame(Mean_Cor = Mean_Cor, Class = Class)
Data$Class <- as.factor(Data$Class)
p3 <- ggplot(data = Data, aes(x = Class, y = Mean_Cor)) +
  geom_violin(trim = TRUE, fill = "#0072B5FF", linewidth = 0.8, colour = "#0072B5FF") +
  geom_boxplot(fill = "white", width = 0.3, linewidth = 1, outlier.size = 1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  ggtitle("SomaLogic pre-ANML") +
  labs(y = "Mean Pearson correlation with other samples", x = "Violin plot") +
  geom_hline(yintercept = mean(Mean_Cor) - 5 * sd(Mean_Cor), linetype = "dashed", color = "black", linewidth = 1)
p3

rm(Class)
rm(Data)

which(Mean_Cor < mean(Mean_Cor) - 5 * sd(Mean_Cor))
D <- rownames(ProtMatrix)[which(Mean_Cor < mean(Mean_Cor) - 5 * sd(Mean_Cor))]
D
#"Plate 5_A11" "Plate 9_H1"

#Outlier Output
Out <- Label[union(B,D),]
write.csv(Out, "SomaLogic_Outlier_Samples_pre_ANML.csv", row.names = F)
rm(Mean_Cor)
rm(Out)
rm(B)
rm(D)

#Combined Figure
library(patchwork)
p1 + p2 + p3

rm(ProtMatrix)
rm(Label)
rm(p1)
rm(p2)
rm(p3)


#SomaLogic ANML

#SomaLogic Data
Dat <- read.csv("Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:20)]
ProtMatrix <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "Sample"),]
Label <- Label[which(Label$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix <- log2(ProtMatrix)
min(ProtMatrix)

#Outlier Detection

#PCA
set.seed(1)
PCA <- prcomp(ProtMatrix)
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PC_Mean <- apply(PCA[,c(1:2)], 2, mean)
PC_SD <- apply(PCA[,c(1:2)], 2, sd)
library(ggplot2)
p1 <- ggplot(PCA, aes(x = PC1, y = PC2)) +
  geom_point(size = 2.5, alpha = 1, color = "#E18727FF", fill = "#E18727FF") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  labs(x = "PC1 (9.3%)", y = "PC2 (5.3%)") +
  ggtitle("SomaLogic ANML") +
  geom_vline(xintercept = PC_Mean[1] + 5 * PC_SD[1], linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = PC_Mean[1] - 5 * PC_SD[1], linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = PC_Mean[2] + 5 * PC_SD[2], linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = PC_Mean[2] - 5 * PC_SD[2], linetype = "dashed", color = "black", linewidth = 1)
p1

which(PCA$PC1 > PC_Mean[1] + 5 * PC_SD[1])
which(PCA$PC1 < PC_Mean[1] - 5 * PC_SD[1])

which(PCA$PC2 > PC_Mean[2] + 5 * PC_SD[2])
which(PCA$PC2 < PC_Mean[2] - 5 * PC_SD[2])

rm(PC_Mean)
rm(PC_SD)
rm(PCA)

#IQR
IQR <- apply(ProtMatrix, 1, IQR)
which(IQR > mean(IQR) + 5 * sd(IQR))
which(IQR < mean(IQR) - 5 * sd(IQR))

#Median
Median <- apply(ProtMatrix, 1, median)
which(Median > mean(Median) + 5 * sd(Median))
C <- names(which(Median > mean(Median) + 5 * sd(Median)))
C
#"Plate 5_A11" "Plate 9_H1"
which(Median < mean(Median) - 5 * sd(Median))

#Visualization
IQR_Median <- data.frame(IQR, Median)
library(ggplot2)
p2 <- ggplot(IQR_Median, aes(x = Median, y = IQR)) +
  geom_point(size = 2.5, alpha = 1, color = "#E18727FF", fill = "#E18727FF") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  labs(x = "Median", y = "IQR") +
  ggtitle("SomaLogic ANML") +
  geom_vline(xintercept = mean(Median) + 5 * sd(Median), linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = mean(Median) - 5 * sd(Median), linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = mean(IQR) + 5 * sd(IQR), linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = mean(IQR) - 5 * sd(IQR), linetype = "dashed", color = "black", linewidth = 1)
p2

rm(IQR)
rm(Median)
rm(IQR_Median)

#Correlation
Cor_Matrix <- cor(t(ProtMatrix), method = "pearson")
Mean_Cor <- vector()
for (i in 1:ncol(Cor_Matrix)) {
  Mean_Cor[i] <- mean(Cor_Matrix[i,][-i])
}
rm(i)
rm(Cor_Matrix)
summary(Mean_Cor)

#Visualization
Class <- 1:800
Class[1:800] <- ""
Data <- data.frame(Mean_Cor = Mean_Cor, Class = Class)
Data$Class <- as.factor(Data$Class)
p3 <- ggplot(data = Data, aes(x = Class, y = Mean_Cor)) +
  geom_violin(trim = TRUE, fill = "#E18727FF", linewidth = 0.8, colour = "#E18727FF") +
  geom_boxplot(fill = "white", width = 0.3, linewidth = 1, outlier.size = 1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  ggtitle("SomaLogic ANML") +
  labs(y = "Mean Pearson correlation with other samples", x = "Violin plot") +
  geom_hline(yintercept = mean(Mean_Cor) - 5 * sd(Mean_Cor), linetype = "dashed", color = "black", linewidth = 1)
p3

rm(Class)
rm(Data)

which(Mean_Cor < mean(Mean_Cor) - 5 * sd(Mean_Cor))
D <- rownames(ProtMatrix)[which(Mean_Cor < mean(Mean_Cor) - 5 * sd(Mean_Cor))]
D
#"Plate 5_A11" "Plate 9_H1"

#Outlier Output
Out <- Label[union(C,D),]
write.csv(Out, "SomaLogic_Outlier_Samples_ANML.csv", row.names = F)
rm(Mean_Cor)
rm(Out)
rm(C)
rm(D)

#Combined Figure
library(patchwork)
p1 + p2 + p3

rm(ProtMatrix)
rm(Label)
rm(p1)
rm(p2)
rm(p3)


#Olink

#Olink Data
Dat <- read.csv("Olink_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label <- Dat[,c(1:19)]
ProtMatrix <- Dat[,c(20:ncol(Dat))]
rm(Dat)

#Remove Control Samples
ProtMatrix <- ProtMatrix[which(Label$SampleType == "SAMPLE"),]
Label <- Label[which(Label$SampleType == "SAMPLE"),]

#Outlier Detection

#PCA
set.seed(1)
PCA <- prcomp(ProtMatrix)
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PC_Mean <- apply(PCA[,c(1:2)], 2, mean)
PC_SD <- apply(PCA[,c(1:2)], 2, sd)
library(ggplot2)
p1 <- ggplot(PCA, aes(x = PC1, y = PC2)) +
  geom_point(size = 2.5, alpha = 1, color = "#20854EFF", fill = "#20854EFF") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  labs(x = "PC1 (3.7%)", y = "PC2 (1.6%)") +
  ggtitle("Olink") +
  geom_vline(xintercept = PC_Mean[1] + 5 * PC_SD[1], linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = PC_Mean[1] - 5 * PC_SD[1], linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = PC_Mean[2] + 5 * PC_SD[2], linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = PC_Mean[2] - 5 * PC_SD[2], linetype = "dashed", color = "black", linewidth = 1)
p1

which(PCA$PC1 > PC_Mean[1] + 5 * PC_SD[1])
which(PCA$PC1 < PC_Mean[1] - 5 * PC_SD[1])

which(PCA$PC2 > PC_Mean[2] + 5 * PC_SD[2])
which(PCA$PC2 < PC_Mean[2] - 5 * PC_SD[2])

rm(PC_Mean)
rm(PC_SD)
rm(PCA)

#IQR
IQR <- apply(ProtMatrix, 1, IQR)
which(IQR > mean(IQR) + 5 * sd(IQR))
which(IQR < mean(IQR) - 5 * sd(IQR))

#Median
Median <- apply(ProtMatrix, 1, median)
which(Median > mean(Median) + 5 * sd(Median))
which(Median < mean(Median) - 5 * sd(Median))

#Visualization
IQR_Median <- data.frame(IQR, Median)
library(ggplot2)
p2 <- ggplot(IQR_Median, aes(x = Median, y = IQR)) +
  geom_point(size = 2.5, alpha = 1, color = "#20854EFF", fill = "#20854EFF") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  labs(x = "Median", y = "IQR") +
  ggtitle("Olink") +
  geom_vline(xintercept = mean(Median) + 5 * sd(Median), linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = mean(Median) - 5 * sd(Median), linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = mean(IQR) + 5 * sd(IQR), linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = mean(IQR) - 5 * sd(IQR), linetype = "dashed", color = "black", linewidth = 1)
p2

rm(IQR)
rm(Median)
rm(IQR_Median)

#Correlation
Cor_Matrix <- cor(t(ProtMatrix), method = "pearson")
Mean_Cor <- vector()
for (i in 1:ncol(Cor_Matrix)) {
  Mean_Cor[i] <- mean(Cor_Matrix[i,][-i])
}
rm(i)
rm(Cor_Matrix)
summary(Mean_Cor)

#Visualization
Class <- 1:514
Class[1:514] <- ""
Data <- data.frame(Mean_Cor = Mean_Cor, Class = Class)
Data$Class <- as.factor(Data$Class)
p3 <- ggplot(data = Data, aes(x = Class, y = Mean_Cor)) +
  geom_violin(trim = TRUE, fill = "#20854EFF", linewidth = 0.8, colour = "#20854EFF") +
  geom_boxplot(fill = "white", width = 0.3, linewidth = 1, outlier.size = 1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  ggtitle("Olink") +
  labs(y = "Mean Pearson correlation with other samples", x = "Violin plot") +
  geom_hline(yintercept = mean(Mean_Cor) - 5 * sd(Mean_Cor), linetype = "dashed", color = "black", linewidth = 1)
p3

rm(Class)
rm(Data)

which(Mean_Cor < mean(Mean_Cor) - 5 * sd(Mean_Cor))
rm(Mean_Cor)

#Combined Figure
library(patchwork)
p1 + p2 + p3

rm(ProtMatrix)
rm(Label)
rm(p1)
rm(p2)
rm(p3)


#Thermo Fisher

#Thermo Fisher Data
Dat <- read.csv("Thermo_Fisher_Merged_All_Imputed.csv")
rownames(Dat) <- Dat$biosample_id
Label <- Dat[,c(1:10)]
ProtMatrix <- Dat[,c(11:ncol(Dat))]
min(ProtMatrix, na.rm = T)
rm(Dat)

#Outlier Detection

#PCA
set.seed(1)
PCA <- prcomp(ProtMatrix)
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PC_Mean <- apply(PCA[,c(1:2)], 2, mean)
PC_SD <- apply(PCA[,c(1:2)], 2, sd)
library(ggplot2)
p1 <- ggplot(PCA, aes(x = PC1, y = PC2)) +
  geom_point(size = 2.5, alpha = 1, color = "#BC3C29FF", fill = "#BC3C29FF") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  labs(x = "PC1 (39.4%)", y = "PC2 (7.1%)") +
  ggtitle("Thermo Fisher") +
  geom_vline(xintercept = PC_Mean[1] + 5 * PC_SD[1], linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = PC_Mean[1] - 5 * PC_SD[1], linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = PC_Mean[2] + 5 * PC_SD[2], linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = PC_Mean[2] - 5 * PC_SD[2], linetype = "dashed", color = "black", linewidth = 1)
p1

which(PCA$PC1 > PC_Mean[1] + 5 * PC_SD[1])
which(PCA$PC1 < PC_Mean[1] - 5 * PC_SD[1])

which(PCA$PC2 > PC_Mean[2] + 5 * PC_SD[2])
which(PCA$PC2 < PC_Mean[2] - 5 * PC_SD[2])

rm(PC_Mean)
rm(PC_SD)
rm(PCA)

#IQR
IQR <- apply(ProtMatrix, 1, IQR, na.rm = T)
which(IQR > mean(IQR) + 5 * sd(IQR))
which(IQR < mean(IQR) - 5 * sd(IQR))

#Median
Median <- apply(ProtMatrix, 1, median, na.rm = T)
which(Median > mean(Median) + 5 * sd(Median))
which(Median < mean(Median) - 5 * sd(Median))

#Visualization
IQR_Median <- data.frame(IQR, Median)
library(ggplot2)
p2 <- ggplot(IQR_Median, aes(x = Median, y = IQR)) +
  geom_point(size = 2.5, alpha = 1, color = "#BC3C29FF", fill = "#BC3C29FF") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  labs(x = "Median", y = "IQR") +
  ggtitle("Thermo Fisher") +
  geom_vline(xintercept = mean(Median) + 5 * sd(Median), linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = mean(Median) - 5 * sd(Median), linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = mean(IQR) + 5 * sd(IQR), linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = mean(IQR) - 5 * sd(IQR), linetype = "dashed", color = "black", linewidth = 1)
p2

rm(IQR)
rm(Median)
rm(IQR_Median)

#Correlation
Cor_Matrix <- cor(t(ProtMatrix), method = "pearson")
Mean_Cor <- vector()
for (i in 1:ncol(Cor_Matrix)) {
  Mean_Cor[i] <- mean(Cor_Matrix[i,][-i])
}
rm(i)
rm(Cor_Matrix)
summary(Mean_Cor)

#Visualization
Class <- 1:95
Class[1:95] <- ""
Data <- data.frame(Mean_Cor = Mean_Cor, Class = Class)
Data$Class <- as.factor(Data$Class)
p3 <- ggplot(data = Data, aes(x = Class, y = Mean_Cor)) +
  geom_violin(trim = TRUE, fill = "#BC3C29FF", linewidth = 0.8, colour = "#BC3C29FF") +
  geom_boxplot(fill = "white", width = 0.3, linewidth = 1, outlier.size = 1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  ggtitle("Thermo Fisher") +
  labs(y = "Mean Pearson correlation with other samples", x = "Violin plot") +
  geom_hline(yintercept = mean(Mean_Cor) - 5 * sd(Mean_Cor), linetype = "dashed", color = "black", linewidth = 1)
p3

rm(Class)
rm(Data)

which(Mean_Cor < mean(Mean_Cor) - 5 * sd(Mean_Cor))
rm(Mean_Cor)

#Combined Figure
library(patchwork)
p1 + p2 + p3

rm(ProtMatrix)
rm(Label)
rm(p1)
rm(p2)
rm(p3)
gc()



