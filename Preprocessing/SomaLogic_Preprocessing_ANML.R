#SomaLogic Preprocessing


#Basic Setting
Sys.setenv(LANGUAGE = "en")

#Package
library(SomaDataIO)


#Data Loading
my_adat <- read_adat("SS-2453461_v5.0_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.20240521.adat")
print(my_adat)
attr(my_adat, "Header.Meta")
attr(my_adat, "Col.Meta")
attr(my_adat, "row_meta")

#Analyte Annotation
Analyte_Annotation <- as.data.frame(getAnalyteInfo(my_adat))
rownames(Analyte_Annotation) <- Analyte_Annotation$AptName

#Non-human proteins
table(Analyte_Annotation$Organism)
ProtIndex_Non_Human <- which(Analyte_Annotation$Organism != "Human")
table(Analyte_Annotation$Type)
ProtIndex_Non_Regular_Protein <- which(Analyte_Annotation$Type %in% c("Non-Human","Spuriomer","Non-Biotin","Non-Cleavable"))

#Protein groups which have at least 2 UNIPROT IDs
A <- data.frame(Analyte_Annotation$UniProt, length = nchar(Analyte_Annotation$UniProt))
which(A$length >= 13)

#Analyte QC
table(Analyte_Annotation$ColCheck)
#No need to remove FLAG proteins

#Remove
Analyte_Annotation <- Analyte_Annotation[-union(ProtIndex_Non_Human, union(ProtIndex_Non_Regular_Protein, which(A$length >= 13))),]
rm(ProtIndex_Non_Human, ProtIndex_Non_Regular_Protein, A)

#Sample Annotation
my_adat <- as.data.frame(my_adat)
Sample_Annotation <- as.data.frame(my_adat[,c(1:33)])

#Remove FLAG buffer sample
table(Sample_Annotation$RowCheck)
which(Sample_Annotation$RowCheck == "FLAG")
Sample_Annotation[414,]#Buffer
my_adat <- my_adat[-which(Sample_Annotation$RowCheck == "FLAG"),]
Sample_Annotation <- Sample_Annotation[-which(Sample_Annotation$RowCheck == "FLAG"),]

#Sample Note
table(Sample_Annotation$AssayNote)
#No need to remove LEAK, but should treat it with caution

#Protein Matrix
Protein_Matrix <- as.data.frame(my_adat[,Analyte_Annotation$AptName])
sum(is.na(Protein_Matrix))
rm(my_adat)

table(Sample_Annotation$SampleType)
table(Analyte_Annotation$Type)

#Now the matrix contains
#Proteins: 10675 human proteins and 12 spike-in controls
#Samples: regular samples (800 in total) and control samples (111 in total)
#Plates: 10 plates 2 batches

#Label
Label <- read.csv("soma_800sample_plating_June2024_withbatch_Plate10update.csv")

#Mapping
Sample_Annotation$PlateId2 <- rep("Plate", nrow(Sample_Annotation))

#Batch 1
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateId == "PLT30760")] <- "Plate 1"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateId == "PLT30761")] <- "Plate 2"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateId == "PLT30762")] <- "Plate 5"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateId == "PLT30763")] <- "Plate 6"

#Batch 2
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateId == "PLT30786")] <- "Plate 3"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateId == "PLT30787")] <- "Plate 4"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateId == "PLT30788")] <- "Plate 7"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateId == "PLT30789")] <- "Plate 8"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateId == "PLT30790")] <- "Plate 9"
Sample_Annotation$PlateId2[which(Sample_Annotation$PlateId == "PLT30892")] <- "Plate 10"

#Unique ID
Sample_Annotation$UniqueID <- paste(Sample_Annotation$PlateId2, Sample_Annotation$PlatePosition, sep = "_")
rownames(Sample_Annotation) <- Sample_Annotation$UniqueID
rownames(Protein_Matrix) <- Sample_Annotation$UniqueID

#Label Unique ID
Label$UniqueID <- paste(Label$plate, Label$well, sep = "_")
rownames(Label) <- Label$UniqueID
intersect(Label$UniqueID, Sample_Annotation$UniqueID)

#Remove FLAG buffer sample
setdiff(Label$UniqueID, Sample_Annotation$UniqueID)
Label <- Label[-which(rownames(Label) == "Plate 3_C3"),]
intersect(Label$UniqueID, Sample_Annotation$UniqueID)

#Mapping
Sample_Annotation <- Sample_Annotation[rownames(Label),]
Protein_Matrix <- Protein_Matrix[rownames(Label),]


#Batch Effects Evaluation

#PCA
set.seed(1)
PCA <- prcomp(log2(Protein_Matrix))
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
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  ggtitle("SomaLogic ANML") +
  labs(x = "PC1 (70.3%)", y = "PC2 (4.2%)")

#Plate
library(ggplot2)
ggplot(PCA, aes(x = PC1, y = PC2, color = Plate, fill = Plate, shape = SampleType)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF","#374E55FF","#80796BFF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  ggtitle("SomaLogic ANML") +
  labs(x = "PC1 (70.3%)", y = "PC2 (4.2%)")
rm(PCA)


#Extract 12 control proteins and only retain experimental samples
Analyte_Annotation2 <- Analyte_Annotation[which(Analyte_Annotation$Type == "Hybridization Control Elution"),]
Protein_Matrix2 <- Protein_Matrix[which(Sample_Annotation$SampleType == "Sample"), which(Analyte_Annotation$Type == "Hybridization Control Elution")]
Label2 <- Label[which(Sample_Annotation$SampleType == "Sample"),]
Sample_Annotation2 <- Sample_Annotation[which(Sample_Annotation$SampleType == "Sample"),]

Protein_Matrix2 <- log2(Protein_Matrix2)

#Merge for General Analysis
Somalogic_Merged <- as.data.frame(cbind(Label2, Sample_Annotation2[,c("PlateId","SampleId","SampleType","SampleNotes","AliquotingNotes","AssayNotes","RowCheck")], Protein_Matrix2))
write.csv(Somalogic_Merged, "Somalogic_Merged_Control.csv", row.names = FALSE)
write.csv(Analyte_Annotation2, "Somalogic_Analyte_Annotation_Control.csv", row.names = FALSE)
rm(Somalogic_Merged)

rm(Analyte_Annotation2)
rm(Protein_Matrix2)
rm(Label2)
rm(Sample_Annotation2)

#Remove 12 control proteins
Analyte_Annotation <- Analyte_Annotation[-which(Analyte_Annotation$Type == "Hybridization Control Elution"),]
Protein_Matrix <- Protein_Matrix[,rownames(Analyte_Annotation)]

#Merge for General Analysis
Somalogic_Merged <- as.data.frame(cbind(Label, Sample_Annotation[,c("PlateId","SampleId","SampleType","SampleNotes","AliquotingNotes","AssayNotes","RowCheck")], Protein_Matrix))
write.csv(Somalogic_Merged, "Somalogic_Merged_All.csv", row.names = FALSE)
write.csv(Analyte_Annotation, "Somalogic_Analyte_Annotation_All.csv", row.names = FALSE)
rm(Somalogic_Merged)

write.csv(as.data.frame(cbind(Label[,-13], Sample_Annotation)), "Somalogic_Sample_Annotation_All.csv", row.names = FALSE)

#Remove 111 control samples
Protein_Matrix <- Protein_Matrix[which(Sample_Annotation$SampleType == "Sample"),]
Label <- Label[which(Sample_Annotation$SampleType == "Sample"),]
Sample_Annotation <- Sample_Annotation[which(Sample_Annotation$SampleType == "Sample"),]

#Log2 transformation for regular samples
Protein_Matrix <- log2(Protein_Matrix)

#PCA
PCA <- prcomp(Protein_Matrix)
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PCA <- cbind(PCA, Plate = Label$plate, Batch = Label$Batch)
PCA$Plate <- as.factor(PCA$Plate)
PCA$Batch <- as.factor(PCA$Batch)
library(ggplot2)
ggplot(PCA, aes(x = PC1, y = PC2, fill = Plate, color = Plate)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF","#374E55FF","#80796BFF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  ggtitle("SomaLogic ANML") +
  labs(x = "PC1 (9.3%)", y = "PC2 (5.3%)")
ggplot(PCA, aes(x = PC1, y = PC2, fill = Batch, color = Batch)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  ggtitle("SomaLogic ANML") +
  labs(x = "PC1 (9.3%)", y = "PC2 (5.3%)")
rm(PCA)


rm(Analyte_Annotation, Label, Protein_Matrix, Sample_Annotation)
gc()

