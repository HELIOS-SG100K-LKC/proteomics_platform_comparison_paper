#Thermo Fisher Preprocessing


#Data Loading
library(readr)
Dat <- as.data.frame(read_tsv("2024US-1201-01_ThermoSingapore_ProteinGroupData_2TechReps.tsv"))
Dat$biosample_id2 <- paste(Dat$biosample_id, Dat$rep_no, Dat$plate_id, sep = "_")
Dat$biosample_id[which(Dat$biosample_id == "PPS2")] <- Dat$biosample_id2[which(Dat$biosample_id == "PPS2")]
Dat <- Dat[,-ncol(Dat)]

#Use Technical Injection Run 1
Dat <- Dat[which(Dat$run_no == "Run1"),]

#Protein Matrix
library(reshape2)
ProtMatrix <- dcast(Dat, biosample_id ~ protein_group_id, value.var = "intensity")
rownames(ProtMatrix) <- ProtMatrix$biosample_id
ProtMatrix <- ProtMatrix[,-1]
ProtMatrix <- log2(ProtMatrix)

#Sample Annotation
Sample_Annotation <- Dat[,c("biosample_id","plate_id")]
Sample_Annotation <- Sample_Annotation[match(rownames(ProtMatrix), Sample_Annotation$biosample_id),]
rownames(Sample_Annotation) <- Sample_Annotation$biosample_id
Sample_Annotation$SampleType <- c(rep("Experimental", 95), rep("Pooled", 15))

#Matching
ProtMatrix <- ProtMatrix[rownames(Sample_Annotation),]

rm(Dat)
gc()

#Remove Contamination Proteins
colnames(ProtMatrix)[grep("CON", colnames(ProtMatrix))]
ProtMatrix <- ProtMatrix[,-grep("CON", colnames(ProtMatrix))]

#Annotation
A <- unlist(strsplit(colnames(ProtMatrix), ";"))
#All human proteins
rm(A)

#Missing Value Imputation
library(SeqKnn)
ProtMatrix_SeqKNN <- SeqKNN(t(ProtMatrix), k = 10)
ProtMatrix_SeqKNN <- as.data.frame(t(ProtMatrix_SeqKNN))

#Exclude protein groups which have at least 2 different UNIPROT IDs
A <- data.frame(colnames(ProtMatrix), length = nchar(colnames(ProtMatrix)))

#Check same UNIPROT ID
check_same_uniprot_id <- function(protein_string) {
  proteins <- unlist(strsplit(protein_string, ";"))
  base_ids <- gsub("-.*", "", proteins)
  all(base_ids == base_ids[1])
}
B <- as.data.frame(sapply(colnames(ProtMatrix)[which(A$length >= 13)], check_same_uniprot_id))
colnames(B) <- "Results"
Index <- rownames(B)[which(B$Results == FALSE)]

ProtMatrix <- ProtMatrix[,-which(colnames(ProtMatrix) %in% Index)]
ProtMatrix_SeqKNN <- ProtMatrix_SeqKNN[,-which(colnames(ProtMatrix_SeqKNN) %in% Index)]

rm(A, B, Index)
rm(check_same_uniprot_id)

#Now the matrix contains
#Proteins: 5945 human protein groups
#Samples: regular samples (95 in total) and control samples (15 in total)
#Plates: 3 plates


#Batch Effects Evaluation

#PCA
set.seed(1)
PCA <- prcomp(ProtMatrix_SeqKNN)
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PCA <- as.data.frame(cbind(PCA, Plate = Sample_Annotation$plate_id, SampleType = Sample_Annotation$SampleType))
PCA$Plate <- as.factor(PCA$Plate)
PCA$SampleType <- as.factor(PCA$SampleType)

#Plate
library(ggplot2)
ggplot(PCA, aes(x = PC1, y = PC2, color = Plate, fill = Plate, shape = SampleType)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  ggtitle("Thermo Fisher") +
  labs(x = "PC1 (36.8%)", y = "PC2 (13.5%)")
rm(PCA)


#Label
Label <- read.csv("TF_Plating_File_Aug2024.csv")
rownames(Label) <- Label$TubeID
intersect(Label$TubeID, Sample_Annotation$biosample_id)

setdiff(Label$TubeID, Sample_Annotation$biosample_id)
Label <- Label[-which(rownames(Label) %in% setdiff(Label$TubeID, Sample_Annotation$biosample_id)),]
intersect(Label$TubeID, Sample_Annotation$biosample_id)

#Mapping
Sample_Annotation <- Sample_Annotation[rownames(Label),]
ProtMatrix <- ProtMatrix[rownames(Label),]
ProtMatrix_SeqKNN <- ProtMatrix_SeqKNN[rownames(Label),]

#PCA
PCA <- prcomp(ProtMatrix_SeqKNN)
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PCA <- cbind(PCA, Plate = Sample_Annotation$plate_id)
PCA$Plate <- as.factor(PCA$Plate)
library(ggplot2)
ggplot(PCA, aes(x = PC1, y = PC2, fill = Plate, color = Plate)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  ggtitle("Thermo Fisher") +
  labs(x = "PC1 (39.4%)", y = "PC2 (7.1%)")
rm(PCA)


#Merge All
Thermo_Fisher_Merged <- as.data.frame(cbind(Label[,-ncol(Label)], Sample_Annotation, ProtMatrix))
write.csv(Thermo_Fisher_Merged, "Thermo_Fisher_Merged_All.csv", row.names = FALSE)

Thermo_Fisher_Merged_Imputed <- as.data.frame(cbind(Label[,-ncol(Label)], Sample_Annotation, ProtMatrix_SeqKNN))
write.csv(Thermo_Fisher_Merged_Imputed, "Thermo_Fisher_Merged_All_Imputed.csv", row.names = FALSE)

rm(Thermo_Fisher_Merged, Thermo_Fisher_Merged_Imputed)

Protein_Annotation <- data.frame(ProteinID = colnames(ProtMatrix), length = nchar(colnames(ProtMatrix)))
C <- vector()
for (i in 1:nrow(Protein_Annotation)) {
  if (Protein_Annotation$length[i] < 13) {
    C[i] <- gsub("-.*", "", Protein_Annotation$ProteinID[i])
  } else {
    C[i] <- substr(Protein_Annotation$ProteinID[i], 1, 6)
  }
}
Protein_Annotation$UniProt <- C
rm(C)
rm(i)
Protein_Annotation$Group <- rep("", nrow(Protein_Annotation))
Protein_Annotation$Group[which(Protein_Annotation$length >= 13)] <- "Isoform_group"
Protein_Annotation$Group[intersect(which(Protein_Annotation$length < 13), grep("-", Protein_Annotation$ProteinID))] <- "Isoform"
Protein_Annotation$Group[setdiff(which(Protein_Annotation$length < 13), grep("-", Protein_Annotation$ProteinID))] <- "Protein"

write.csv(Protein_Annotation, "Thermo_Fisher_Protein_Annotation_All.csv", row.names = FALSE)

rm(Label, ProtMatrix, ProtMatrix_SeqKNN, Sample_Annotation)
rm(Protein_Annotation)
gc()

