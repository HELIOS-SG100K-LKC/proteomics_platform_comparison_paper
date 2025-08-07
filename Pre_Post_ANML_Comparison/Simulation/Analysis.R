#Analysis

#Basic Setting
Sys.setenv(LANGUAGE = "en")

#Package
library(SomaDataIO)


#1 BMI

#Data Loading
my_adat <- read_adat("Simulation_Data_1.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/1_BMI/Simulated_Protein_BMI.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/1_BMI/Simulated_Protein_BMI_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/1_BMI/Simulated_Protein_BMI_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b1 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("BMI\nsrho with ANML factor = -0.35\nsrho with BMI = 1") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b1

#P
library(ggplot2)
p1 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("BMI\nsrho with ANML factor = -0.35\nsrho with BMI = 1") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p1

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#2 WHR

#Data Loading
my_adat <- read_adat("Simulation_Data_2.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/2_WHR/Simulated_Protein_WHR.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/2_WHR/Simulated_Protein_WHR_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/2_WHR/Simulated_Protein_WHR_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b2 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("WHR\nsrho with ANML factor = -0.14\nsrho with BMI = 0.38") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b2

#P
library(ggplot2)
p2 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("WHR\nsrho with ANML factor = -0.14\nsrho with BMI = 0.38") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p2

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#3 Triglyceride

#Data Loading
my_adat <- read_adat("Simulation_Data_3.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/3_Triglyceride/Simulated_Protein_Triglyceride.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/3_Triglyceride/Simulated_Protein_Triglyceride_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/3_Triglyceride/Simulated_Protein_Triglyceride_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b3 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Triglyceride\nsrho with ANML factor = -0.29\nsrho with BMI = 0.40") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b3

#P
library(ggplot2)
p3 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Triglyceride\nsrho with ANML factor = -0.29\nsrho with BMI = 0.40") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p3

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#4 Insulin

#Data Loading
my_adat <- read_adat("Simulation_Data_4.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/4_Insulin/Simulated_Protein_Insulin.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/4_Insulin/Simulated_Protein_Insulin_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/4_Insulin/Simulated_Protein_Insulin_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b4 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Insulin\nsrho with ANML factor = -0.36\nsrho with BMI = 0.59") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b4

#P
library(ggplot2)
p4 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.border = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Insulin\nsrho with ANML factor = -0.36\nsrho with BMI = 0.59") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p4

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#5 Glucose

#Data Loading
my_adat <- read_adat("Simulation_Data_5.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/5_Glucose/Simulated_Protein_Glucose.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/5_Glucose/Simulated_Protein_Glucose_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/5_Glucose/Simulated_Protein_Glucose_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b5 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Glucose\nsrho with ANML factor = -0.13\nsrho with BMI = 0.27") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b5

#P
library(ggplot2)
p5 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Glucose\nsrho with ANML factor = -0.13\nsrho with BMI = 0.27") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p5

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#6 ALT

#Data Loading
my_adat <- read_adat("Simulation_Data_6.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/6_ALT/Simulated_Protein_ALT.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/6_ALT/Simulated_Protein_ALT_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/6_ALT/Simulated_Protein_ALT_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b6 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("ALT\nsrho with ANML factor = -0.17\nsrho with BMI = 0.34") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b6

#P
library(ggplot2)
p6 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("ALT\nsrho with ANML factor = -0.17\nsrho with BMI = 0.34") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p6

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#7 GGT

#Data Loading
my_adat <- read_adat("Simulation_Data_7.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/7_GGT/Simulated_Protein_GGT.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/7_GGT/Simulated_Protein_GGT_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/7_GGT/Simulated_Protein_GGT_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b7 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.05, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("GGT\nsrho with ANML factor = -0.16\nsrho with BMI = 0.34") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b7

#P
library(ggplot2)
p7 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("GGT\nsrho with ANML factor = -0.16\nsrho with BMI = 0.34") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p7

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#8 HDL

#Data Loading
my_adat <- read_adat("Simulation_Data_8.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/8_HDL/Simulated_Protein_HDL.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/8_HDL/Simulated_Protein_HDL_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/8_HDL/Simulated_Protein_HDL_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b8 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("HDL\nsrho with ANML factor = 0.22\nsrho with BMI = -0.44") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b8

#P
library(ggplot2)
p8 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("HDL\nsrho with ANML factor = 0.22\nsrho with BMI = -0.44") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p8

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#9 CRP

#Data Loading
my_adat <- read_adat("Simulation_Data_9.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/9_CRP/Simulated_Protein_CRP.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/9_CRP/Simulated_Protein_CRP_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/9_CRP/Simulated_Protein_CRP_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b9 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("CRP\nsrho with ANML factor = -0.37\nsrho with BMI = 0.52") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b9

#P
library(ggplot2)
p9 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("CRP\nsrho with ANML factor = -0.37\nsrho with BMI = 0.52") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p9

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#10 Creatinine

#Data Loading
my_adat <- read_adat("Simulation_Data_10.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/10_Creatinine/Simulated_Protein_Creatinine.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/10_Creatinine/Simulated_Protein_Creatinine_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/10_Creatinine/Simulated_Protein_Creatinine_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b10 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Creatinine\nsrho with ANML factor = 0.06\nsrho with BMI = 0.06") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b10

#P
library(ggplot2)
p10 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Creatinine\nsrho with ANML factor = 0.06\nsrho with BMI = 0.06") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p10

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#11 Uric Acid

#Data Loading
my_adat <- read_adat("Simulation_Data_11.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/11_Uric_Acid/Simulated_Protein_Uric_Acid.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/11_Uric_Acid/Simulated_Protein_Uric_Acid_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/11_Uric_Acid/Simulated_Protein_Uric_Acid_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b11 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Uric Acid\nsrho with ANML factor = -0.1\nsrho with BMI = 0.33") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b11

#P
library(ggplot2)
p11 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Uric Acid\nsrho with ANML factor = -0.1\nsrho with BMI = 0.33") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p11

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#12 SBP

#Data Loading
my_adat <- read_adat("Simulation_Data_12.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/12_SBP/Simulated_Protein_SBP.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/12_SBP/Simulated_Protein_SBP_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/12_SBP/Simulated_Protein_SBP_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b12 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("SBP\nsrho with ANML factor = -0.14\nsrho with BMI = 0.22") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b12

#P
library(ggplot2)
p12 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("SBP\nsrho with ANML factor = -0.14\nsrho with BMI = 0.22") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p12

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#13 Total Cholesterol

#Data Loading
my_adat <- read_adat("Simulation_Data_13.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/13_Total_Cholesterol/Simulated_Protein_Total_Cholesterol.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/13_Total_Cholesterol/Simulated_Protein_Total_Cholesterol_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/13_Total_Cholesterol/Simulated_Protein_Total_Cholesterol_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b13 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Total Cholesterol\nsrho with ANML factor = -0.12\nsrho with BMI = -0.04") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b13

#P
library(ggplot2)
p13 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Total Cholesterol\nsrho with ANML factor = -0.12\nsrho with BMI = -0.04") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p13

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#14 Age

#Data Loading
my_adat <- read_adat("Simulation_Data_14.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/14_Age/Simulated_Protein_Age.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/14_Age/Simulated_Protein_Age_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/14_Age/Simulated_Protein_Age_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b14 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Age\nsrho with ANML factor = 0.01\nsrho with BMI = -0.05") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b14

#P
library(ggplot2)
p14 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Age\nsrho with ANML factor = 0.01\nsrho with BMI = -0.05") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p14

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#15 Urea

#Data Loading
my_adat <- read_adat("Simulation_Data_15.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/15_Urea/Simulated_Protein_Urea.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/15_Urea/Simulated_Protein_Urea_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/15_Urea/Simulated_Protein_Urea_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b15 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Urea\nsrho with ANML factor = 0.02\nsrho with BMI = -0.01") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b15

#P
library(ggplot2)
p15 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Urea\nsrho with ANML factor = 0.02\nsrho with BMI = -0.01") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p15

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#16 Cl

#Data Loading
my_adat <- read_adat("Simulation_Data_16.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/16_Cl/Simulated_Protein_Cl.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/16_Cl/Simulated_Protein_Cl_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/16_Cl/Simulated_Protein_Cl_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b16 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Cl\nsrho with ANML factor = 0.18\nsrho with BMI = 0.07") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b16

#P
library(ggplot2)
p16 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Cl\nsrho with ANML factor = 0.18\nsrho with BMI = 0.07") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p16

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#17 Na

#Data Loading
my_adat <- read_adat("Simulation_Data_17.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/17_Na/Simulated_Protein_Na.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/17_Na/Simulated_Protein_Na_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/17_Na/Simulated_Protein_Na_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b17 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Na\nsrho with ANML factor = 0.14\nsrho with BMI = -0.06") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b17

#P
library(ggplot2)
p17 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Na\nsrho with ANML factor = 0.14\nsrho with BMI = -0.06") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p17

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#18 Albumin

#Data Loading
my_adat <- read_adat("Simulation_Data_18.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/18_Albumin/Simulated_Protein_Albumin.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/18_Albumin/Simulated_Protein_Albumin_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/18_Albumin/Simulated_Protein_Albumin_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b18 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Albumin\nsrho with ANML factor = -0.11\nsrho with BMI = -0.13") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b18

#P
library(ggplot2)
p18 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("Albumin\nsrho with ANML factor = -0.11\nsrho with BMI = -0.13") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p18

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#19 BMD_Hip

#Data Loading
my_adat <- read_adat("Simulation_Data_19.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/19_BMD_Hip/Simulated_Protein_BMD_Hip.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/19_BMD_Hip/Simulated_Protein_BMD_Hip_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/19_BMD_Hip/Simulated_Protein_BMD_Hip_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b19 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("BMD_Hip\nsrho with ANML factor = -0.21\nsrho with BMI = 0.46") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b19

#P
library(ggplot2)
p19 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("BMD_Hip\nsrho with ANML factor = -0.21\nsrho with BMI = 0.46") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p19

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#20 BMD_Lumbar

#Data Loading
my_adat <- read_adat("Simulation_Data_20.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/20_BMD_Lumbar/Simulated_Protein_BMD_Lumbar.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/2_SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID
Dat <- Dat[which(Dat$SampleType == "Sample"),]
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID),]
my_adat <- as.data.frame(my_adat)
rownames(my_adat) <- my_adat$SubjectID
my_adat <- my_adat[Dat$SOMA_ID,]
rownames(my_adat) <- Dat$FREG0_PID
rm(Dat)

my_adat <- my_adat[rownames(Simulation), colnames(Simulation)]

#Phenotypic Association
my_adat <- log2(my_adat)
Simulation <- log2(Simulation)

Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/20_BMD_Lumbar/Simulated_Protein_BMD_Lumbar_Phenotype.csv")
rownames(Phenotype) <- Phenotype$X
Phenotype <- Phenotype[rownames(Simulation),]
Phenotype <- Phenotype$Phenotype

#pre-ANML
pre_ANML_beta <- vector()
pre_ANML_p <- vector()
for (i in 1:100) {
  pre_ANML_beta[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,1]
  pre_ANML_p[i] <- summary(lm(Phenotype ~ Simulation[,i]))$coefficients[2,4]
}
rm(i)

Check <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/4_Spike-in_Simulation/1_Simulation/20_BMD_Lumbar/Simulated_Protein_BMD_Lumbar_Annotation.csv")
identical(Check$Aptamer, colnames(Simulation))
sum(abs(Check$Simulated_Beta - pre_ANML_beta))
sum(abs(Check$Simulated_P - pre_ANML_p))
rm(Check)

#ANML
ANML_beta <- vector()
ANML_p <- vector()
for (i in 1:100) {
  ANML_beta[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,1]
  ANML_p[i] <- summary(lm(Phenotype ~ my_adat[,i]))$coefficients[2,4]
}
rm(i)

#Visualization

#Beta
library(ggplot2)
library(ggpubr)
b20 <- ggplot(data.frame(Beta_pre_ANML = pre_ANML_beta, Beta_ANML = ANML_beta, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = Beta_pre_ANML, y = Beta_ANML)) +
  geom_point(size = 3, aes(fill = Direction, color = Direction)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_regline_equation(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = 0.1, label.y = -0.5, size = 4.5) +
  stat_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(-1.5, 1.6)) +
  scale_y_continuous(limits = c(-1.5, 1.6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("BMD_Lumbar\nsrho with ANML factor = -0.06\nsrho with BMI = 0.28") +
  labs(x = "Beta (Pre-ANML)", y = "Beta (ANML)")
b20

#P
library(ggplot2)
p20 <- ggplot(data.frame(P_pre_ANML = pre_ANML_p, P_ANML = ANML_p, Direction = c(rep("Positive", 50), rep("Negative", 50))), aes(x = -log10(P_pre_ANML), y = -log10(P_ANML), fill = Direction, color = Direction)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_color_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  scale_x_continuous(limits = c(0, 32)) +
  scale_y_continuous(limits = c(0, 32)) +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 13.5, face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("BMD_Lumbar\nsrho with ANML factor = -0.06\nsrho with BMI = 0.28") +
  labs(x = "-log10p (Pre-ANML)", y = "-log10p (ANML)")
p20

rm(pre_ANML_beta, ANML_beta)
rm(pre_ANML_p, ANML_p)
rm(Phenotype, my_adat, Simulation)


#Summary
library(patchwork)
b9 + b4 + b1 + b3 + b19 + b6 + b7 + b2 + b12 + b5 + b13 + b18 + b11 + b20 + b14 + b15 + b10 + b17 + b16 + b8 + plot_layout(ncol = 4)
p9 + p4 + p1 + p3 + p19 + p6 + p7 + p2 + p12 + p5 + p13 + p18 + p11 + p20 + p14 + p15 + p10 + p17 + p16 + p8 + plot_layout(ncol = 4)

b9 + b1 + b15 + b8 + plot_layout(ncol = 2)

rm(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20)
rm(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20)

Cor_ANML <- c(-0.37, -0.36, -0.35, -0.29, -0.21, -0.17, -0.16, -0.14, -0.14, -0.13, -0.12, -0.11, -0.10, -0.06, 0.01, 0.02, 0.06, 0.14, 0.18, 0.22)
Cor_BMI <- c(0.52, 0.59, 1, 0.40, 0.46, 0.34, 0.34, 0.38, 0.22, 0.27, -0.04, -0.13, 0.33, 0.28, -0.05, -0.01, 0.06, -0.06, 0.07, -0.44)
Intercept <- c(-0.32, -0.32, -0.27, -0.29, -0.19, -0.14, -0.14, -0.13, -0.13, -0.12, -0.095, -0.098, -0.11, -0.055, -0.019, 0.0025, 0.026, 0.083, 0.11, 0.18)

cor(Cor_ANML, Intercept, method = "spearman")
cor(Cor_BMI, Intercept, method = "spearman")
cor(Cor_BMI, Cor_ANML, method = "spearman")

rm(Cor_ANML, Cor_BMI, Intercept)

gc()

