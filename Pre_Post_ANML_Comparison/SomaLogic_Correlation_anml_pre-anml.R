#SomaLogic correlation between ANML and pre-ANML normalization strategies

#ANML

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

#500
Index <- union(which(Label$tech_rep %in% c(2:3)), which(Label$bio_rep_id == "Y"))
ProtMatrix <- ProtMatrix[-Index,]
Label <- Label[-Index,]
rm(Index)

#Pre-ANML

#SomaLogic Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Label2 <- Dat[,c(1:20)]
ProtMatrix2 <- Dat[,c(21:ncol(Dat))]
rm(Dat)

#Analyte Annotation
Analyte_Annotation2 <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Analyte_Annotation_All.csv")
rownames(Analyte_Annotation2) <- Analyte_Annotation2$AptName

#Remove Control Samples
ProtMatrix2 <- ProtMatrix2[which(Label2$SampleType == "Sample"),]
Label2 <- Label2[which(Label2$SampleType == "Sample"),]

#Log2 Transformation
ProtMatrix2 <- log2(ProtMatrix2)
min(ProtMatrix2)

#500
Index <- union(which(Label2$tech_rep %in% c(2:3)), which(Label2$bio_rep_id == "Y"))
ProtMatrix2 <- ProtMatrix2[-Index,]
Label2 <- Label2[-Index,]
rm(Index)

#Mapping
Label2 <- Label2[rownames(Label), colnames(Label)]
Analyte_Annotation2 <- Analyte_Annotation2[rownames(Analyte_Annotation), colnames(Analyte_Annotation)]
ProtMatrix2 <- ProtMatrix2[rownames(ProtMatrix), colnames(ProtMatrix)]

#Pearson Correlation

#Protein-Level
Protein_Cor <- vector()
for (i in 1:10675) {
  Protein_Cor[i] <- cor(ProtMatrix[,i], ProtMatrix2[,i])
}
rm(i)
summary(Protein_Cor)

length(which(Protein_Cor < 0.8))
length(which(Protein_Cor < 0.6))

#Visualization
data <- data.frame(correlation = Protein_Cor)
data$group <- cut(data$correlation, breaks = c(-Inf, 0.6, 0.8, Inf), labels = c("< 0.6", "0.6 - 0.8", "> 0.8"))
ggplot(data, aes(x = correlation, fill = group)) +
  geom_histogram(binwidth = 0.05, color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("< 0.6" = "#0072B5FF", "0.6 - 0.8" = "#E18727FF", "> 0.8" = "#BC3C29FF")) +
  labs(title = "SomaLogic Protein-level", x = "Pearson correlation between pre-ANML and ANML", y = "Frequency", fill = "Correlation Range") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1))
rm(data)

write.csv(data.frame(Protein = colnames(ProtMatrix), Corr = Protein_Cor), "Protein_Level_Pearson_Correlation_anml_pre-anml.csv", row.names = FALSE)

rm(ProtMatrix)
rm(ProtMatrix2)
rm(Label)
rm(Label2)
rm(Analyte_Annotation)
rm(Analyte_Annotation2)

rm(Protein_Cor)
gc()
