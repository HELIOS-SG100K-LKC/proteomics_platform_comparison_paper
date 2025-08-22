#CV Thermo Fisher


#Thermo Fisher Data
Dat <- read.csv("Thermo_Fisher_Merged_All_Imputed.csv")
rownames(Dat) <- Dat$biosample_id
Label <- Dat[,c(1:10)]
ProtMatrix <- Dat[,c(11:ncol(Dat))]
min(ProtMatrix, na.rm = T)
rm(Dat)

#Protein Annotation
Protein_Annotation <- read.csv("Thermo_Fisher_Protein_Annotation_All.csv")
rownames(Protein_Annotation) <- Protein_Annotation$ProteinID

#CV
CV <- function(x) {
  return(100 * sqrt(exp((log(2) * sd(x))^2) - 1))
}


#Overall CV
ID <- as.data.frame(table(Label$FREG0_PID))
ID <- as.character(ID$Var1[which(ID$Freq == 2)])

CV_TF_Tech <- matrix(NA, nrow = ncol(ProtMatrix), ncol = length(ID))
for (i in 1:45) {
  CV_TF_Tech[,i] <- apply(t(ProtMatrix[which(Label$FREG0_PID == ID[i]),]), 1, CV)
}
rm(i)
CV_TF_Tech <- as.data.frame(CV_TF_Tech)
rownames(CV_TF_Tech) <- colnames(ProtMatrix)

CV_Intra <- apply(CV_TF_Tech, 1, mean)
summary(CV_Intra)
round(quantile(CV_Intra, c(0.25, 0.5, 0.75, 0.9)), 2)

rm(CV_TF_Tech)
rm(ID)


#Stratify by LoD
Missingness <- read.csv("Thermo_Fisher_ProteinMissingness.csv")
identical(Missingness$Protein, Protein_Annotation$ProteinID)

#Intra

#Above
median(CV_Intra[which(Missingness$Missing_Rate_Protein <= 0.2)])

#Below
median(CV_Intra[which(Missingness$Missing_Rate_Protein > 0.2)])

rm(Missingness)


rm(CV_Intra)
rm(CV)
rm(ProtMatrix)
rm(Label)
rm(Protein_Annotation)
gc()


