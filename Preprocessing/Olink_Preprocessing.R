#Olink Preprocessing


#Package
library(OlinkAnalyze)
library(arrow)


#Data Loading
my_NPX_data <- as.data.frame(read_NPX("Q-14739_Chambers_NPX_2024-08-05.parquet"))

#Panel & Explore Version
unique(my_NPX_data$Panel)
#Explore_HT
unique(my_NPX_data$ExploreVersion)
#7.0.5

#Normalization Method
unique(my_NPX_data$Normalization)
#Intensity

#Sample QC
unique(my_NPX_data$SampleQC)
table(my_NPX_data$SampleID[which(my_NPX_data$SampleQC == "WARN")])
my_NPX_data$NPX[which(my_NPX_data$SampleQC == "WARN")]
#Retain with caution

table(my_NPX_data$SampleID[which(my_NPX_data$SampleQC == "FAIL")])
my_NPX_data$NPX[which(my_NPX_data$SampleQC == "FAIL")]
sum(is.na(my_NPX_data$NPX))
#All NaN for FAIL rows
#Remove these samples ("0439600684" and "2137803170")

my_NPX_data <- my_NPX_data[-which(my_NPX_data$SampleID %in% c("0439600684","2137803170")),]
sum(is.na(my_NPX_data$NPX))

#Assay QC
unique(my_NPX_data$AssayQC)
table(my_NPX_data$AssayType[which(my_NPX_data$AssayQC == "NA")])
#Control Proteins: Retain

table(my_NPX_data$AssayType[which(my_NPX_data$AssayQC == "WARN")])
table(my_NPX_data$OlinkID[which(my_NPX_data$AssayQC == "WARN")])
my_NPX_data$NPX[which(my_NPX_data$AssayQC == "WARN")]
#Retain with caution

#Protein Matrix
library(reshape2)
ProtMatrix <- dcast(my_NPX_data, SampleID ~ OlinkID, value.var = "NPX")
rownames(ProtMatrix) <- ProtMatrix$SampleID
ProtMatrix <- ProtMatrix[,-1]
sum(is.na(ProtMatrix))

#Assay Annotation
Assay_Annotation <- my_NPX_data[,c("OlinkID","UniProt","Assay","AssayType","Block","AssayQC")]
Assay_Annotation <- Assay_Annotation[match(colnames(ProtMatrix), Assay_Annotation$OlinkID),]
rownames(Assay_Annotation) <- Assay_Annotation$OlinkID

#Sample Annotation
Sample_Annotation <- my_NPX_data[,c("SampleID","SampleType","WellID","PlateID","SampleQC")]
Sample_Annotation <- Sample_Annotation[match(rownames(ProtMatrix), Sample_Annotation$SampleID),]
rownames(Sample_Annotation) <- Sample_Annotation$SampleID

#Matching
ProtMatrix <- ProtMatrix[rownames(Sample_Annotation), rownames(Assay_Annotation)]

rm(my_NPX_data)
gc()

#Protein groups which have at least 2 UNIPROT IDs
A <- data.frame(Assay_Annotation$UniProt, length = nchar(Assay_Annotation$UniProt))
which(A$length >= 13)

#Remove
Assay_Annotation <- Assay_Annotation[-which(A$length >= 13),]
rm(A)

#Protein Matrix
ProtMatrix <- as.data.frame(ProtMatrix[,Assay_Annotation$OlinkID])

table(Assay_Annotation$AssayType)
table(Sample_Annotation$SampleType)

#Now the matrix contains
#Proteins: 5405 human proteins and 24 spike-in controls
#Samples: regular samples (514 in total) and control samples (60 in total)
#Plates: 6 plates 3 batches

#Label
Label <- read.csv("olink_sample_plating_June2024_withControlBatch.csv")

#Mapping
Sample_Annotation$PlateId2 <- rep("Plate", nrow(Sample_Annotation))

#Batch 1
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateID == "Q-14739__flowcellposition_PA_SS240047_SP240208")] <- "Plate 1"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateID == "Q-14739__flowcellposition_PB_SS240047_SP240211")] <- "Plate 4"

#Batch 2
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateID == "Q-14739__flowcellposition_PA_SS240048_SP240209")] <- "Plate 2"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateID == "Q-14739__flowcellposition_PB_SS240048_SP240212")] <- "Plate 5"

#Batch 3
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateID == "Q-14739__flowcellposition_PA_SS240049_SP240210")] <- "Plate 3"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateID == "Q-14739__flowcellposition_PB_SS240049_SP240213")] <- "Plate 6"

#Unique ID
Sample_Annotation$UniqueID <- paste(Sample_Annotation$PlateId2, Sample_Annotation$WellID, sep = "_")
rownames(Sample_Annotation) <- Sample_Annotation$UniqueID
rownames(ProtMatrix) <- Sample_Annotation$UniqueID

#Label Unique ID
Label$UniqueID <- paste(Label$plate, Label$well, sep = "_")
rownames(Label) <- Label$UniqueID
intersect(Label$UniqueID, Sample_Annotation$UniqueID)

#Remove QC Fail samples
setdiff(Label$UniqueID, Sample_Annotation$UniqueID)
Label <- Label[-which(rownames(Label) %in% c("Plate 4_A1","Plate 6_A6")),]
intersect(Label$UniqueID, Sample_Annotation$UniqueID)

#Mapping
Sample_Annotation <- Sample_Annotation[rownames(Label),]
ProtMatrix <- ProtMatrix[rownames(Label),]


#Batch Effects Evaluation

#PCA
set.seed(1)
PCA <- prcomp(ProtMatrix)
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PCA <- as.data.frame(cbind(PCA, Batch = Label$Batch, Plate = Label$plate, SampleType = Sample_Annotation$SampleType))
PCA$Batch <- as.factor(PCA$Batch)
PCA$Plate <- as.factor(PCA$Plate)
PCA$SampleType <- as.factor(PCA$SampleType)

#Batch
library(ggplot2)
ggplot(PCA, aes(x = PC1, y = PC2, color = Batch, fill = Batch, shape = SampleType)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  ggtitle("Olink") +
  labs(x = "PC1 (21.5%)", y = "PC2 (6.1%)")

#Plate
library(ggplot2)
ggplot(PCA, aes(x = PC1, y = PC2, color = Plate, fill = Plate, shape = SampleType)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  ggtitle("Olink") +
  labs(x = "PC1 (21.5%)", y = "PC2 (6.1%)")
rm(PCA)


#Extract 16 control proteins and only retain experimental samples
Assay_Annotation2 <- Assay_Annotation[which(Assay_Annotation$AssayType %in% c("inc_ctrl","amp_ctrl")),]
ProtMatrix2 <- ProtMatrix[which(Sample_Annotation$SampleType == "SAMPLE"), which(Assay_Annotation$AssayType %in% c("inc_ctrl","amp_ctrl"))]
Label2 <- Label[which(Sample_Annotation$SampleType == "SAMPLE"),]
Sample_Annotation2 <- Sample_Annotation[which(Sample_Annotation$SampleType == "SAMPLE"),]

#Merge for General Analysis
Olink_Merged <- as.data.frame(cbind(Label2, Sample_Annotation2[,c("SampleID","SampleType","PlateID","SampleQC")], ProtMatrix2))
write.csv(Olink_Merged, "Olink_Merged_Control.csv", row.names = FALSE)
write.csv(Assay_Annotation2, "Olink_Assay_Annotation_Control.csv", row.names = FALSE)
rm(Olink_Merged)

rm(Assay_Annotation2)
rm(ProtMatrix2)
rm(Label2)
rm(Sample_Annotation2)

#Remove 24 control proteins
Assay_Annotation <- Assay_Annotation[which(Assay_Annotation$AssayType == "assay"),]
ProtMatrix <- ProtMatrix[,rownames(Assay_Annotation)]

#Merge for General Analysis
Olink_Merged <- as.data.frame(cbind(Label, Sample_Annotation[,c("SampleID","SampleType","PlateID","SampleQC")], ProtMatrix))
write.csv(Olink_Merged, "Olink_Merged_All.csv", row.names = FALSE)
write.csv(Assay_Annotation, "Olink_Assay_Annotation_All.csv", row.names = FALSE)
rm(Olink_Merged)

#Remove 60 control samples
ProtMatrix <- ProtMatrix[which(Sample_Annotation$SampleType == "SAMPLE"),]
Label <- Label[which(Sample_Annotation$SampleType == "SAMPLE"),]
Sample_Annotation <- Sample_Annotation[which(Sample_Annotation$SampleType == "SAMPLE"),]

#PCA
PCA <- prcomp(ProtMatrix)
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PCA <- cbind(PCA, Plate = Label$plate, Batch = Label$Batch)
PCA$Plate <- as.factor(PCA$Plate)
PCA$Batch <- as.factor(PCA$Batch)
library(ggplot2)
ggplot(PCA, aes(x = PC1, y = PC2, fill = Plate, color = Plate)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  ggtitle("Olink") +
  labs(x = "PC1 (3.7%)", y = "PC2 (1.6%)")
ggplot(PCA, aes(x = PC1, y = PC2, fill = Batch, color = Batch)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  ggtitle("Olink") +
  labs(x = "PC1 (3.7%)", y = "PC2 (1.6%)")
rm(PCA)

rm(Assay_Annotation, Label, ProtMatrix, Sample_Annotation)
gc()

