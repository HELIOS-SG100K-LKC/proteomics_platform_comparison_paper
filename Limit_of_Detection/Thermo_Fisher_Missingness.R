#Thermo Fisher Missingness


#Data Loading

#Thermo Fisher Data
Dat <- read.csv("Thermo_Fisher_Merged_All.csv")
rownames(Dat) <- Dat$biosample_id
Label <- Dat[,c(1:10)]
ProtMatrix <- Dat[,c(11:ncol(Dat))]
min(ProtMatrix, na.rm = T)
rm(Dat)

#Protein Annotation
Protein_Annotation <- read.csv("Thermo_Fisher_Protein_Annotation_All.csv")
rownames(Protein_Annotation) <- Protein_Annotation$ProteinID

#46
Index <- which(Label$tech_rep == 2)
ProtMatrix <- ProtMatrix[-Index,]
Label <- Label[-Index,]
rm(Index)

#Missingness Calculation

#Overall
sum(is.na(ProtMatrix))/(nrow(ProtMatrix) * ncol(ProtMatrix))

#Sample Missingness
Missing_Rate_Sample <- vector()
for (i in 1:nrow(ProtMatrix)) {
  Missing_Rate_Sample[i] <- sum(is.na(ProtMatrix[i,]))/ncol(ProtMatrix)
}
rm(i)
summary(Missing_Rate_Sample)

#Protein Missingness
Missing_Rate_Protein <- vector()
for (i in 1:ncol(ProtMatrix)) {
  Missing_Rate_Protein[i] <- sum(is.na(ProtMatrix[,i]))/nrow(ProtMatrix)
}
rm(i)

#Sample
Class <- 1:46
Class[1:46] <- ""
Data <- data.frame(Missing_Rate_Sample = Missing_Rate_Sample, Class = Class)
Data$Class <- as.factor(Data$Class)
library(ggplot2)
ggplot(data = Data, aes(x = Class, y = Missing_Rate_Sample)) +
  geom_violin(trim = TRUE, fill = "#BC3C29FF", linewidth = 0.8, colour = "#BC3C29FF") +
  geom_boxplot(fill = "white", width = 0.3, linewidth = 1, outlier.size = 1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.1)) +
  ggtitle("Thermo Fisher") +
  labs(y = "Sample missingness", x = "Violin plot")
rm(Class)
rm(Data)

#Protein
Rate <- vector()
Rate[1] <- length(Missing_Rate_Protein)
Rate[2] <- length(which(Missing_Rate_Protein == 0))
Rate[3] <- length(which(Missing_Rate_Protein > 0.2))
Rate[4] <- length(which(Missing_Rate_Protein > 0.5))
Rate[5] <- length(which(Missing_Rate_Protein == 1))
Rate_Percentage <- Rate / length(Missing_Rate_Protein) * 100
Data <- data.frame(Rate, Rate_Percentage, Threshold = c("All", "=0%", ">20%", ">50%", "=100%"))

ggplot(Data, aes(x = Threshold, y = Rate)) +
  geom_bar(stat = "identity", position = "dodge", color = "#BC3C29FF", fill = "#BC3C29FF") +
  geom_text(aes(label = paste0(Rate, " (", sprintf("%.0f", Rate_Percentage), "%)")),
            hjust = 0.4, vjust = -0.4, size = 3.5, color = "black") +
  scale_x_discrete(limits = Data$Threshold) +
  scale_y_continuous(
    name = "Number of proteins",
    sec.axis = sec_axis(~ . / length(Missing_Rate_Protein) * 100, name = "Proportion of all proteins (%)")
  ) +
  theme_classic() +
  ggtitle(paste0("Thermo Fisher")) +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
  labs(x = "Protein missingness", y = "Number of proteins")

rm(Rate)
rm(Data)
rm(Rate_Percentage)


#Proteins with missing rate > 20%
Dat <- data.frame(Protein = Protein_Annotation$ProteinID, Missing_Rate_Protein = Missing_Rate_Protein)
Dat <- Dat[Dat$Missing_Rate_Protein > 0.2,]
Dat <- cbind(Dat, Protein_Annotation[Dat$Protein,])
write.csv(Dat, "Thermo_Fisher_HighProportionProteins_0.2.csv", row.names = F)
rm(Dat)

#Proteins with missingness annotation
Dat <- data.frame(Protein = Protein_Annotation$ProteinID, Missing_Rate_Protein = Missing_Rate_Protein)
Dat <- cbind(Dat, Protein_Annotation[Dat$Protein,])
write.csv(Dat, "Thermo_Fisher_ProteinMissingness.csv", row.names = F)
rm(Dat)

rm(Missing_Rate_Sample)
rm(Missing_Rate_Protein)


#Remove All
rm(list = ls())
gc()


