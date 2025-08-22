#Analysis

#Basic Setting
Sys.setenv(LANGUAGE = "en")

#Package
library(SomaDataIO)


#1 BMI

#Data Loading
my_adat <- read_adat("Simulation_Data_1.anmlSMP.adat")

#Simulated Data
Simulation <- read.csv("Simulated_Protein_BMI.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data
Dat <- read.csv("SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
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

Phenotype <- read.csv("Simulated_Protein_BMI_Phenotype.csv")
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

Check <- read.csv("Simulated_Protein_BMI_Annotation.csv")
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

rm(b1)
rm(p1)

gc()


