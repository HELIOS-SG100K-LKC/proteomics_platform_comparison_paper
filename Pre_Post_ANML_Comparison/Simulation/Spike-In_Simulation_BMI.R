#Spike-In Simulation


#Data Preparation

#SomaLogic pre-ANML Data
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

#500 Samples
Index <- intersect(which(Label$tech_rep_id == "N"), which(Label$bio_rep_id == "N"))
Index <- union(Index, which(Label$tech_rep == 1))
ProtMatrix <- ProtMatrix[Index,]
Label <- Label[Index,]
rm(Index)
rownames(ProtMatrix) <- Label$FREG0_PID
rownames(Label) <- Label$FREG0_PID

#Largest Dilution Bin 1:5
ProtMatrix <- ProtMatrix[,which(Analyte_Annotation$Dilution2 == 1/5)]
Analyte_Annotation <- Analyte_Annotation[which(Analyte_Annotation$Dilution2 == 1/5),]

#Phenotypes
Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/1_Demographics/HELIOS_Core_v4.csv")
Phenotype <- Phenotype[which(Phenotype$FREG0_PID %in% Label$FREG0_PID),]
Phenotype <- Phenotype[which(Phenotype$FREG14_Visit_number == "1"),]
Phenotype$BMI <- Phenotype$DBI14_Weight/((Phenotype$DBI13_Height/100)^2)
rownames(Phenotype) <- Phenotype$FREG0_PID
Phenotype <- Phenotype[rownames(Label),]

rm(Label)
rm(Analyte_Annotation)


#Define the Phenotype
Phenotype <- Phenotype$BMI


#Simulation

#INT Transformation
library(RNOmni)
Phenotype <- RankNorm(Phenotype)

#Linear Regression on original data
beta <- vector()
p <- vector()
for (i in 1:ncol(ProtMatrix)) {
  lm_model <- summary(lm(Phenotype ~ ProtMatrix[,i]))
  beta[i] <- lm_model$coefficients[2,1]
  p[i] <- lm_model$coefficients[2,4]
}
rm(i)
rm(lm_model)

#Select 45 p from 10^-0.5 to 10^-9 and 5 p from 10^-10 to 10^-20, uniformly distributed in -log10 space
Desired_neg_log10p <- c(seq(0.5, 9, length.out = 45), seq(10, 20, length.out = 5))
Desired_p <- 10^(-Desired_neg_log10p)
rm(Desired_neg_log10p)

#Positive
ProtMatrix_Pos <- ProtMatrix[,which(beta > 0)]
beta_Pos <- beta[which(beta > 0)]
p_Pos <- p[which(beta > 0)]

#Find the closest p-value to the desired p-values
Index <- vector()
for (i in 1:50) {
  Index[i] <- which.min(abs(p_Pos - Desired_p[i]))
}
rm(i)
DF_Pos <- data.frame(Direction = rep("Positive", 50), Desired_P = Desired_p, Aptamer = colnames(ProtMatrix_Pos)[Index], Original_Beta = beta_Pos[Index], Original_P = p_Pos[Index])
ProtMatrix_Pos <- ProtMatrix_Pos[,Index]

rm(Index)
rm(beta_Pos)
rm(p_Pos)

#Check
length(unique(DF_Pos$Aptamer))
sum(DF_Pos$Original_Beta > 0)
data.frame(DF_Pos$Desired_P, DF_Pos$Desired_P - DF_Pos$Original_P)

#Simulation
target_rho_Pearson <- vector()
pp <- matrix(NA, 40000, 50)
for (j in 1:40000) {
  
  #Set your desired Pearson correlation in the latent (Gaussian) space
  target_rho <- j/100000
  target_rho_Pearson[j] <- target_rho
  
  print(j)
  
  #For each measured Phenotype, sample a latent Protein value
  #from the conditional distribution: N(target_rho * Phenotype, 1 - target_rho^2)
  set.seed(1)
  z_Protein <- rnorm(length(Phenotype), mean = target_rho * Phenotype, sd = sqrt(1 - target_rho^2))
  
  #Convert the new latent Protein values to uniform scale
  U_Protein <- pnorm(z_Protein)
  
  #Generate simulated proteins and calculate the new p-values
  for (i in 1:50) {
    #Map these uniform values to the actual Protein scale using the empirical quantile function
    #One simple method is to sort the measured Protein values and assign new ones accordingly
    #For each U value, select the corresponding quantile (using ceiling to index)
    Protein_new <- sort(ProtMatrix_Pos[,i])[ceiling(U_Protein * length(Phenotype))]
    
    #Linear regression
    pp[j,i] <- summary(lm(Phenotype ~ Protein_new))$coefficients[2,4]
  }
}
rm(j)
rm(i)
rm(target_rho)
rm(z_Protein)
rm(U_Protein)
rm(Protein_new)

#Select the target rho
colnames(pp) <- colnames(ProtMatrix_Pos)
rownames(pp) <- target_rho_Pearson
apply(pp, 2, max)
apply(pp, 2, min)

Row_Index <- vector()
for (i in 1:50) {
  Row_Index[i] <- which.min(abs(pp[,i] - Desired_p[i]))
}
rm(i)
DF_Pos$Target_Pearson <- target_rho_Pearson[Row_Index]
rm(Row_Index)
rm(target_rho_Pearson)
rm(pp)

#Generate final simulated proteins
pp <- vector()
betabeta <- vector()
cor_Pearson <- vector()
cor_Spearman <- vector()
Simulated_Protein <- matrix(NA, nrow = length(Phenotype), ncol = 50)
for (i in 1:50) {
  #Set your desired Pearson correlation in the latent (Gaussian) space
  target_rho <- DF_Pos$Target_Pearson[i]

  #For each measured Phenotype, sample a latent Protein value
  #from the conditional distribution: N(target_rho * Phenotype, 1 - target_rho^2)
  set.seed(1)
  z_Protein <- rnorm(length(Phenotype), mean = target_rho * Phenotype, sd = sqrt(1 - target_rho^2))
  
  #Convert the new latent Protein values to uniform scale
  U_Protein <- pnorm(z_Protein)
  
  #Generate simulated proteins and calculate the new p-values
  #Map these uniform values to the actual Protein scale using the empirical quantile function
  #One simple method is to sort the measured Protein values and assign new ones accordingly
  #For each U value, select the corresponding quantile (using ceiling to index)
  Protein_new <- sort(ProtMatrix_Pos[,i])[ceiling(U_Protein * length(Phenotype))]
  Simulated_Protein[,i] <- Protein_new
  
  #Check the new correlation (Spearman correlation is often more robust to the mapping)
  cor_Pearson[i] <- cor(Phenotype, Protein_new, method = "pearson")
  cor_Spearman[i] <- cor(Phenotype, Protein_new, method = "spearman")
  
  #Linear regression
  betabeta[i] <- summary(lm(Phenotype ~ Protein_new))$coefficients[2,1]
  pp[i] <- summary(lm(Phenotype ~ Protein_new))$coefficients[2,4]
}
rm(i)
rm(target_rho)
rm(z_Protein)
rm(U_Protein)
rm(Protein_new)

#Results Summary
DF_Pos$Simulated_P <- pp
DF_Pos$Simulated_Beta <- betabeta
DF_Pos$Simulated_Correlation_Pearson <- cor_Pearson
DF_Pos$Simulated_Correlation_Spearman <- cor_Spearman
rm(pp)
rm(betabeta)
rm(cor_Pearson)
rm(cor_Spearman)

#Check
sum(DF_Pos$Simulated_Beta > 0)
data.frame(DF_Pos$Desired_P, DF_Pos$Desired_P - DF_Pos$Simulated_P)

rownames(Simulated_Protein) <- rownames(ProtMatrix_Pos)
colnames(Simulated_Protein) <- colnames(ProtMatrix_Pos)
Simulated_Protein <- as.data.frame(Simulated_Protein)

#Check Distribution
par(mfrow = c(2,5))
for (i in 1:50) {
  plot(density(Simulated_Protein[,i]), col = "blue", main = "Simulated Protein (Blue) vs Original Protein (Red)", xlab = "Protein", lwd = 5)
  lines(density(ProtMatrix_Pos[,i]), col = "red", lwd = 2)
}
rm(i)
rm(ProtMatrix_Pos)

#Negative
ProtMatrix_Neg <- ProtMatrix[,which(beta < 0)]
beta_Neg <- beta[which(beta < 0)]
p_Neg <- p[which(beta < 0)]

rm(beta)
rm(p)

#Find the closest p-value to the desired p-values
Index <- vector()
for (i in 1:50) {
  Index[i] <- which.min(abs(p_Neg - Desired_p[i]))
}
rm(i)
DF_Neg <- data.frame(Direction = rep("Negative", 50), Desired_P = Desired_p, Aptamer = colnames(ProtMatrix_Neg)[Index], Original_Beta = beta_Neg[Index], Original_P = p_Neg[Index])
ProtMatrix_Neg <- ProtMatrix_Neg[,Index]

rm(Index)
rm(beta_Neg)
rm(p_Neg)

#Check
length(unique(DF_Neg$Aptamer))
sum(DF_Neg$Original_Beta < 0)
data.frame(DF_Neg$Desired_P, DF_Neg$Desired_P - DF_Neg$Original_P)

#Simulation
target_rho_Pearson <- vector()
pp <- matrix(NA, 40000, 50)
for (j in 1:40000) {
  
  #Set your desired Pearson correlation in the latent (Gaussian) space
  target_rho <- (-j/100000)
  target_rho_Pearson[j] <- target_rho
  
  print(j)
  
  #For each measured Phenotype, sample a latent Protein value
  #from the conditional distribution: N(target_rho * Phenotype, 1 - target_rho^2)
  set.seed(2)
  z_Protein <- rnorm(length(Phenotype), mean = target_rho * Phenotype, sd = sqrt(1 - target_rho^2))
  
  #Convert the new latent Protein values to uniform scale
  U_Protein <- pnorm(z_Protein)
  
  #Generate simulated proteins and calculate the new p-values
  for (i in 1:50) {
    #Map these uniform values to the actual Protein scale using the empirical quantile function
    #One simple method is to sort the measured Protein values and assign new ones accordingly
    #For each U value, select the corresponding quantile (using ceiling to index)
    Protein_new <- sort(ProtMatrix_Neg[,i])[ceiling(U_Protein * length(Phenotype))]
    
    #Linear regression
    pp[j,i] <- summary(lm(Phenotype ~ Protein_new))$coefficients[2,4]
  }
}
rm(j)
rm(i)
rm(target_rho)
rm(z_Protein)
rm(U_Protein)
rm(Protein_new)

#Select the target rho
colnames(pp) <- colnames(ProtMatrix_Neg)
rownames(pp) <- target_rho_Pearson
apply(pp, 2, max)
apply(pp, 2, min)

Row_Index <- vector()
for (i in 1:50) {
  Row_Index[i] <- which.min(abs(pp[,i] - Desired_p[i]))
}
rm(i)
DF_Neg$Target_Pearson <- target_rho_Pearson[Row_Index]
rm(Row_Index)
rm(target_rho_Pearson)
rm(pp)

#Generate final simulated proteins
pp <- vector()
betabeta <- vector()
cor_Pearson <- vector()
cor_Spearman <- vector()
Simulated_Protein2 <- matrix(NA, nrow = length(Phenotype), ncol = 50)
for (i in 1:50) {
  #Set your desired Pearson correlation in the latent (Gaussian) space
  target_rho <- DF_Neg$Target_Pearson[i]
  
  #For each measured Phenotype, sample a latent Protein value
  #from the conditional distribution: N(target_rho * Phenotype, 1 - target_rho^2)
  set.seed(2)
  z_Protein <- rnorm(length(Phenotype), mean = target_rho * Phenotype, sd = sqrt(1 - target_rho^2))
  
  #Convert the new latent Protein values to uniform scale
  U_Protein <- pnorm(z_Protein)
  
  #Generate simulated proteins and calculate the new p-values
  #Map these uniform values to the actual Protein scale using the empirical quantile function
  #One simple method is to sort the measured Protein values and assign new ones accordingly
  #For each U value, select the corresponding quantile (using ceiling to index)
  Protein_new <- sort(ProtMatrix_Neg[,i])[ceiling(U_Protein * length(Phenotype))]
  Simulated_Protein2[,i] <- Protein_new
  
  #Check the new correlation (Spearman correlation is often more robust to the mapping)
  cor_Pearson[i] <- cor(Phenotype, Protein_new, method = "pearson")
  cor_Spearman[i] <- cor(Phenotype, Protein_new, method = "spearman")
  
  #Linear regression
  betabeta[i] <- summary(lm(Phenotype ~ Protein_new))$coefficients[2,1]
  pp[i] <- summary(lm(Phenotype ~ Protein_new))$coefficients[2,4]
}
rm(i)
rm(target_rho)
rm(z_Protein)
rm(U_Protein)
rm(Protein_new)

#Results Summary
DF_Neg$Simulated_P <- pp
DF_Neg$Simulated_Beta <- betabeta
DF_Neg$Simulated_Correlation_Pearson <- cor_Pearson
DF_Neg$Simulated_Correlation_Spearman <- cor_Spearman
rm(pp)
rm(betabeta)
rm(cor_Pearson)
rm(cor_Spearman)

#Check
sum(DF_Neg$Simulated_Beta < 0)
data.frame(DF_Neg$Desired_P, DF_Neg$Desired_P - DF_Neg$Simulated_P)

rownames(Simulated_Protein2) <- rownames(ProtMatrix_Neg)
colnames(Simulated_Protein2) <- colnames(ProtMatrix_Neg)
Simulated_Protein2 <- as.data.frame(Simulated_Protein2)

#Check Distribution
par(mfrow = c(2,5))
for (i in 1:50) {
  plot(density(Simulated_Protein2[,i]), col = "blue", main = "Simulated Protein (Blue) vs Original Protein (Red)", xlab = "Protein", lwd = 5)
  lines(density(ProtMatrix_Neg[,i]), col = "red", lwd = 2)
}
rm(i)
rm(ProtMatrix_Neg)


#Summary
rm(Desired_p)

Simulated_Protein <- as.data.frame(cbind(Simulated_Protein, Simulated_Protein2))
rm(Simulated_Protein2)

DF <- data.frame(rbind(DF_Pos, DF_Neg))
rm(DF_Pos)
rm(DF_Neg)

Phenotype <- data.frame(Phenotype)
rownames(Phenotype) <- rownames(ProtMatrix)
write.csv(data.frame(Phenotype), file = "Simulated_Protein_BMI_Phenotype.csv", row.names = TRUE)
rm(Phenotype)

#PCA

#Protein-Level
set.seed(1)
PCA <- prcomp(t(cbind(Simulated_Protein, ProtMatrix)))
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PCA$Protein <- c(rep("Simulated Protein", 100), rep("Original Protein", 8574))
library(ggplot2)
ggplot(PCA, aes(x = PC1, y = PC2, fill = Protein, color = Protein)) +
  geom_point(size = 2.5, alpha = 0.5) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  labs(x = "PC1 (93.9%)", y = "PC2 (0.4%)")
rm(PCA)

#Sample-Level
set.seed(1)
ProtMatrix2 <- ProtMatrix[,-which(colnames(ProtMatrix) %in% colnames(Simulated_Protein))]
ProtMatrix2 <- as.data.frame(cbind(ProtMatrix2, Simulated_Protein))
ProtMatrix2 <- ProtMatrix2[,colnames(ProtMatrix)]
PCA <- prcomp(rbind(ProtMatrix, ProtMatrix2))
(summary(PCA))$importance[2,c(1:2)]
PCA <- as.data.frame(predict(PCA)[,1:2])
PCA$Sample <- c(rep("Original", 500), rep("Simulated", 500))
library(ggplot2)
ggplot(PCA, aes(x = PC1, y = PC2, fill = Sample, color = Sample)) +
  geom_point(size = 2.5, alpha = 0.5) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  labs(x = "PC1 (23.9%)", y = "PC2 (5.7%)")
rm(PCA)
rm(ProtMatrix2)
rm(ProtMatrix)


#Output
Simulated_Protein <- 2^(Simulated_Protein)

write.csv(DF, file = "Simulated_Protein_BMI_Annotation.csv", row.names = FALSE)
write.csv(Simulated_Protein, file = "Simulated_Protein_BMI.csv", row.names = TRUE)

rm(Simulated_Protein)
rm(DF)
gc()

