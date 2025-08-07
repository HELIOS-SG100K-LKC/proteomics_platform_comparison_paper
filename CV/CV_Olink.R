#CV Olink


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

#CV
CV <- function(x) {
  return(100 * sqrt(exp((log(2) * sd(x))^2) - 1))
}


#Overall Intra-Plate CV
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep %in% c(1:2)),]
Label_Tech <- Label[which(Label$tech_rep %in% c(1:2)),]
ID <- unique(Label_Tech$FREG0_PID)

CV_Olink_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Tech), ncol = length(ID))
for (i in 1:100) {
  CV_Olink_Tech[,i] <- apply(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_Olink_Tech <- as.data.frame(CV_Olink_Tech)
rownames(CV_Olink_Tech) <- colnames(ProtMatrix_Tech)

CV_Intra <- apply(CV_Olink_Tech, 1, mean)
summary(CV_Intra)
round(quantile(CV_Intra, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Olink_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)

#Overall Inter-Plate CV
ProtMatrix_Tech <- ProtMatrix[which(Label$tech_rep %in% c(2:3)),]
Label_Tech <- Label[which(Label$tech_rep %in% c(2:3)),]
ID <- unique(Label_Tech$FREG0_PID)

CV_Olink_Tech <- matrix(NA, nrow = ncol(ProtMatrix_Tech), ncol = length(ID))
for (i in 1:100) {
  CV_Olink_Tech[,i] <- apply(t(ProtMatrix_Tech[which(Label_Tech$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_Olink_Tech <- as.data.frame(CV_Olink_Tech)
rownames(CV_Olink_Tech) <- colnames(ProtMatrix_Tech)

CV_Inter <- apply(CV_Olink_Tech, 1, mean, na.rm = T)
summary(CV_Inter)
round(quantile(CV_Inter, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_Olink_Tech)
rm(ID)
rm(Label_Tech)
rm(ProtMatrix_Tech)


#Stratify by LoD
LOD_Olink <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/3_Limit_of_Detection/2_Olink_LoD/Olink_LOD.csv")
identical(LOD_Olink$OlinkID, Assay_Annotation$OlinkID)

#Intra

#Above
median(CV_Intra[which(LOD_Olink$PropBelowLOD <= 0.2)])

#Below
median(CV_Intra[which(LOD_Olink$PropBelowLOD > 0.2)])

#Inter

#Above
median(CV_Inter[which(LOD_Olink$PropBelowLOD <= 0.2)])

#Below
median(CV_Inter[which(LOD_Olink$PropBelowLOD > 0.2)])

rm(LOD_Olink)


#Different Sets

#Set 1: Olink Explore 1536
library(readxl)
ProteinList1536 <- read_excel("explore-1536-assay-list-20210227-web-1.xlsx")
colnames(ProteinList1536) <- ProteinList1536[1,]
ProteinList1536 <- ProteinList1536[-1,]
V1 <- which(Assay_Annotation$UniProt %in% ProteinList1536$`Uniprot ID`)

#Set 2: Olink Explore 3072 - Olink Explore 1536
ProteinList3072 <- read_excel("olink-explore-3072-validation-data-results.xlsx")
colnames(ProteinList3072) <- ProteinList3072[1,]
ProteinList3072 <- ProteinList3072[-1,]
ProteinList3072 <- ProteinList3072[-which(ProteinList3072$UniProt %in% ProteinList1536$`Uniprot ID`),]
V2 <- which(Assay_Annotation$UniProt %in% ProteinList3072$UniProt)
rm(ProteinList1536)
rm(ProteinList3072)

#Set 3: Olink Explore HT - Olink Explore 3072
ProteinList3072 <- read_excel("olink-explore-3072-validation-data-results.xlsx")
colnames(ProteinList3072) <- ProteinList3072[1,]
ProteinList3072 <- ProteinList3072[-1,]
V3 <- setdiff(1:5405, which(Assay_Annotation$UniProt %in% ProteinList3072$UniProt))
rm(ProteinList3072)

#Intra

#Set 1
median(CV_Intra[V1])

#Set 2
median(CV_Intra[V2])

#Set 3
median(CV_Intra[V3])

#Inter

#Set 1
median(CV_Inter[V1])

#Set 2
median(CV_Inter[V2])

#Set 3
median(CV_Inter[V3])

rm(V1)
rm(V2)
rm(V3)


#Dilution series
Assay_Annotation$Dilution <- Assay_Annotation$Block
Assay_Annotation$Dilution[which(Assay_Annotation$Block %in% c(1,2,3,4))] <- "1:1"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 5)] <- "1:10"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 6)] <- "1:100"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 7)] <- "1:1000"
Assay_Annotation$Dilution[which(Assay_Annotation$Block == 8)] <- "1:100000"

#1:1
D1 <- which(Assay_Annotation$Dilution == "1:1")

#1:10
D2 <- which(Assay_Annotation$Dilution == "1:10")

#1:100
D3 <- which(Assay_Annotation$Dilution == "1:100")

#1:1000
D4 <- which(Assay_Annotation$Dilution == "1:1000")

#1:100000
D5 <- which(Assay_Annotation$Dilution == "1:100000")

#Intra

#1:1
median(CV_Intra[D1])

#1:10
median(CV_Intra[D2])

#1:100
median(CV_Intra[D3])

#1:1000
median(CV_Intra[D4])

#1:100000
median(CV_Intra[D5])

#Inter

#1:1
median(CV_Inter[D1])

#1:10
median(CV_Inter[D2])

#1:100
median(CV_Inter[D3])

#1:1000
median(CV_Inter[D4])

#1:100000
median(CV_Inter[D5])

rm(D1)
rm(D2)
rm(D3)
rm(D4)
rm(D5)

rm(CV_Intra)
rm(CV_Inter)
rm(CV)
rm(Assay_Annotation)
rm(Label)
rm(ProtMatrix)
gc()

