#Correlation between ANML Scaling Factors and Phenotypes


#ANML Fraction
ANML_Fraction <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/1_Preprocessing/1_SomaLogic_anml_Preprocessing/Somalogic_Sample_Annotation_All.csv")
Index <- intersect(which(ANML_Fraction$tech_rep_id == "N"), which(ANML_Fraction$bio_rep_id == "N"))
Index <- union(Index, which(ANML_Fraction$tech_rep == 1))
ANML_Fraction <- ANML_Fraction[Index,]
rm(Index)

#Fraction
summary(ANML_Fraction$ANMLFractionUsed_20)
summary(ANML_Fraction$ANMLFractionUsed_0_5)
summary(ANML_Fraction$ANMLFractionUsed_0_005)
rm(ANML_Fraction)


#ANML Scaling Factors
ANML_Scaling_Factors <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/2_SomaLogic_Comparison_anml_pre-anml/2_ANML_Scaling_Factor_Generation/Somalogic_ANML_Scaling_Factor.csv")
Index <- intersect(which(ANML_Scaling_Factors$tech_rep_id == "N"), which(ANML_Scaling_Factors$bio_rep_id == "N"))
Index <- union(Index, which(ANML_Scaling_Factors$tech_rep == 1))
ANML_Scaling_Factors <- ANML_Scaling_Factors[Index,]
rm(Index)
rownames(ANML_Scaling_Factors) <- ANML_Scaling_Factors$FREG0_PID

#Phenotypes
Phenotype <- read.csv("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/1_Demographics/HELIOS_Core_v4.csv")
Phenotype <- Phenotype[which(Phenotype$FREG0_PID %in% ANML_Scaling_Factors$FREG0_PID),]
Phenotype <- Phenotype[which(Phenotype$FREG14_Visit_number == "1"),]
rownames(Phenotype) <- Phenotype$FREG0_PID
Phenotype <- Phenotype[rownames(ANML_Scaling_Factors),]

Phenotype$BMI <- Phenotype$DBI14_Weight / ((Phenotype$DBI13_Height / 100) ^ 2)
Phenotype$WHR <- c(Phenotype$FWH16_Waist1 + Phenotype$FWH17_Waist2 + Phenotype$FWH18_Waist3) / c(Phenotype$FWH19_Hip1 + Phenotype$FWH20_Hip2 + Phenotype$FWH21_Hip3) 
Phenotype$SBP <- (Phenotype$DBP9_Sbp_1 + Phenotype$DBP10_Sbp_2 + Phenotype$DBP11_Sbp_3) / 3
Phenotype$DBP <- (Phenotype$DBP12_Dbp_1 + Phenotype$DBP13_Dbp_2 + Phenotype$DBP14_Dbp_3) / 3
Phenotype$HR <- (Phenotype$DBP15_Hr_1 + Phenotype$DBP16_Hr_2 + Phenotype$DBP17_Hr_3) / 3

#Phenotypes 2
Phenotype2 <- as.data.frame(readxl::read_xlsx("D:/PhD Thesis/Multi-Platform Proteomics Comparison/6_Other_Analysis/1_Demographics/HELIOS_DXA1_Hip_Lumbar.xlsx"))
Phenotype2 <- Phenotype2[which(Phenotype2$FREG0_PID %in% ANML_Scaling_Factors$FREG0_PID),]
Phenotype2 <- Phenotype2[which(Phenotype2$FREG14_Visit_number == "1"),]
rownames(Phenotype2) <- Phenotype2$FREG0_PID
Phenotype2 <- Phenotype2[rownames(ANML_Scaling_Factors),]
rownames(Phenotype2) <- rownames(ANML_Scaling_Factors)

Phenotype$BMD_Hip <- Phenotype2$DDH31_Total_Tscore
Phenotype$BMD_Lumbar <- Phenotype2$DDL52_Total_T_Score
rm(Phenotype2)

#Correlation
library(corrplot)
A <- cbind(ANML_Scaling_Factors[,c(9,10,11)], Phenotype[,c("FREG8_Age","FREG7_Gender","BMI","WHR","SBP","DBP","HR","DLAB15_Tc_Mmol_L","DLAB17_Hdl_Mmol_L","DLAB19_Ldl_Mmol_L","DLAB21_Trig_Mmol_L","DLAB25_Gluf_Mmol_L","DLAB27_Hba1C_Percent","DLAB81_Insulin","DLAB28_Ua_Mmol_L","DLAB30_Na_Mmol_L","DLAB31_K_Mmol_L","DLAB32_Cl_Mmol_L","DLAB33_Urea_Mmol_L","DLAB35_Creat_Umol_L","DLAB75_ALT_U_L","DLAB76_GGT_U_L","DLAB77_Bilirubin_umol_L","DLAB80_Albumin_mg_dL","DLAB83_CRP","DLAB82_VitD","BMD_Hip","BMD_Lumbar")])
colnames(A) <- c("ANML_1_5","ANML_1_200","ANML_1_20000","Age","Sex","BMI","WHR","SBP","DBP","Heart Rate","Total Cholesterol","HDL","LDL","Triglyceride","Glucose","HbA1c","Insulin","Uric Acid","Na","K","Cl","Urea","Creatinine","ALT","GGT","Bilirubin","Albumin","CRP","Vitamin D","BMD_Hip","BMD_Lumbar")
A$Sex <- ifelse(A$Sex == "M", 1, 0)

#All
corrplot(cor(A, method = "spearman", use = "pairwise.complete.obs")[c(6,1:3),-c(1:3)], method = "shade", addCoef.col = "black", tl.cex = 1.25, cl.pos = "n")

#Chinese
C <- A[which(Phenotype$FREG5_Ethnic_Group == "C"),]
corrplot(cor(C, method = "spearman", use = "pairwise.complete.obs")[c(6,1:3),-c(1:3)], method = "shade", addCoef.col = "black", tl.cex = 1.25, cl.pos = "n")

#Indian
I <- A[which(Phenotype$FREG5_Ethnic_Group == "I"),]
corrplot(cor(I, method = "spearman", use = "pairwise.complete.obs")[c(6,1:3),-c(1:3)], method = "shade", addCoef.col = "black", tl.cex = 1.25, cl.pos = "n")

#Malay
M <- A[which(Phenotype$FREG5_Ethnic_Group == "M"),]
corrplot(cor(M, method = "spearman", use = "pairwise.complete.obs")[c(6,1:3),-c(1:3)], method = "shade", addCoef.col = "black", tl.cex = 1.25, cl.pos = "n")

#Combined
corrplot(rbind(cor(A, method = "spearman", use = "pairwise.complete.obs")[c(6,1:3),-c(1:3)], rep(0,28), cor(C, method = "spearman", use = "pairwise.complete.obs")[c(6,1:3),-c(1:3)], rep(0,28), cor(I, method = "spearman", use = "pairwise.complete.obs")[c(6,1:3),-c(1:3)], rep(0,28), cor(M, method = "spearman", use = "pairwise.complete.obs")[c(6,1:3),-c(1:3)]), method = "shade", addCoef.col = "black", tl.cex = 1.25, cl.pos = "n")

#The Whole Correlation Matrix
corrplot(cor(A, method = "spearman", use = "pairwise.complete.obs"), method = "shade", addCoef.col = "black", tl.cex = 1.25, cl.pos = "n")

rm(A)
rm(C, I, M)

#Linear Regression Ethnicity

library(forestmodel)

#1:5
Scale <- data.frame(Scale = c(ANML_Scaling_Factors[which(Phenotype$FREG5_Ethnic_Group == "C"), 9], ANML_Scaling_Factors[which(Phenotype$FREG5_Ethnic_Group == "I"), 9], ANML_Scaling_Factors[which(Phenotype$FREG5_Ethnic_Group == "M"), 9]))
Scale$Ethnic_Group <- c(rep("C",167), rep("I",166), rep("M",167))
summary(lm(data = Scale, Scale ~ Ethnic_Group))
forest_model(lm(data = Scale, Scale ~ Ethnic_Group))
rm(Scale)

#1:200
Scale <- data.frame(Scale = c(ANML_Scaling_Factors[which(Phenotype$FREG5_Ethnic_Group == "C"), 10], ANML_Scaling_Factors[which(Phenotype$FREG5_Ethnic_Group == "I"), 10], ANML_Scaling_Factors[which(Phenotype$FREG5_Ethnic_Group == "M"), 10]))
Scale$Ethnic_Group <- c(rep("C",167), rep("I",166), rep("M",167))
summary(lm(data = Scale, Scale ~ Ethnic_Group))
forest_model(lm(data = Scale, Scale ~ Ethnic_Group))
rm(Scale)

#1:20000
Scale <- data.frame(Scale = c(ANML_Scaling_Factors[which(Phenotype$FREG5_Ethnic_Group == "C"), 11], ANML_Scaling_Factors[which(Phenotype$FREG5_Ethnic_Group == "I"), 11], ANML_Scaling_Factors[which(Phenotype$FREG5_Ethnic_Group == "M"), 11]))
Scale$Ethnic_Group <- c(rep("C",167), rep("I",166), rep("M",167))
summary(lm(data = Scale, Scale ~ Ethnic_Group))
forest_model(lm(data = Scale, Scale ~ Ethnic_Group))
rm(Scale)


#Machine Learning

#1:5
A <- cbind(ANML_Scaling_Factors[,9], Phenotype[,c("FREG8_Age","FREG7_Gender","FREG5_Ethnic_Group","BMI","WHR","SBP","DBP","HR","DLAB15_Tc_Mmol_L","DLAB17_Hdl_Mmol_L","DLAB19_Ldl_Mmol_L","DLAB21_Trig_Mmol_L","DLAB25_Gluf_Mmol_L","DLAB27_Hba1C_Percent","DLAB81_Insulin","DLAB28_Ua_Mmol_L","DLAB30_Na_Mmol_L","DLAB31_K_Mmol_L","DLAB32_Cl_Mmol_L","DLAB33_Urea_Mmol_L","DLAB35_Creat_Umol_L","DLAB75_ALT_U_L","DLAB76_GGT_U_L","DLAB77_Bilirubin_umol_L","DLAB80_Albumin_mg_dL","DLAB83_CRP","DLAB82_VitD","BMD_Hip","BMD_Lumbar")])
colnames(A) <- c("ANML_1_5","Age","Sex","Ethnicity","BMI","WHR","SBP","DBP","Heart_Rate","Total_Cholesterol","HDL","LDL","Triglyceride","Glucose","HbA1c","Insulin","Uric_Acid","Na","K","Cl","Urea","Creatinine","ALT","GGT","Bilirubin","Albumin","CRP","Vitamin_D","BMD_Hip","BMD_Lumbar")
A$Sex <- as.factor(A$Sex)
A$Ethnicity <- as.factor(A$Ethnicity)
A <- na.omit(A)
A <- as.data.frame(A)

library(Boruta)
set.seed(1)
Boruta1 <- Boruta(ANML_1_5 ~., data = A, doTrace = 2, maxRuns = 50000)
Boruta1
Trace <- as.data.frame(Boruta1$ImpHistory)
Trace <- Trace[,-c(30:32)]

DF <- data.frame(Importance = unlist(as.vector(Trace)), Phenotype = rep(colnames(Trace), each = 2478))
DF <- DF[-which(DF$Importance == -Inf),]

DF$Direction <- rep(NA, nrow(DF))
DF$Direction[which(DF$Phenotype %in% c("Ethnicity","BMI","WHR","SBP","DBP","Heart_Rate","Total_Cholesterol","LDL","Triglyceride","Glucose","HbA1c","Insulin","Uric_Acid","K","ALT","GGT","Albumin","CRP","BMD_Hip","BMD_Lumbar"))] <- "Negative"
DF$Direction[which(DF$Phenotype %in% c("Age","Sex","HDL","Na","Cl","Urea","Creatinine","Bilirubin","Vitamin_D"))] <- "Positive"

library(dplyr)
DF1 <- DF %>% group_by(Phenotype) %>% summarise(Median = median(Importance))
DF1 <- arrange(DF1, Median)            

#Importance Boxplot
library(ggplot2)
ggplot(DF, aes(x = Phenotype, y = Importance, fill = Direction)) +
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(limits = DF1$Phenotype) +
  xlab("Phenotype") +
  ylab("Importance") +
  ggtitle("ANML_1_5 Prediction\nBoruta Importance") +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 10, color = "black")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 12, face = "bold"))
rm(Boruta1)
rm(A)
rm(DF)
rm(Trace)

#1:200
A <- cbind(ANML_Scaling_Factors[,10], Phenotype[,c("FREG8_Age","FREG7_Gender","FREG5_Ethnic_Group","BMI","WHR","SBP","DBP","HR","DLAB15_Tc_Mmol_L","DLAB17_Hdl_Mmol_L","DLAB19_Ldl_Mmol_L","DLAB21_Trig_Mmol_L","DLAB25_Gluf_Mmol_L","DLAB27_Hba1C_Percent","DLAB81_Insulin","DLAB28_Ua_Mmol_L","DLAB30_Na_Mmol_L","DLAB31_K_Mmol_L","DLAB32_Cl_Mmol_L","DLAB33_Urea_Mmol_L","DLAB35_Creat_Umol_L","DLAB75_ALT_U_L","DLAB76_GGT_U_L","DLAB77_Bilirubin_umol_L","DLAB80_Albumin_mg_dL","DLAB83_CRP","DLAB82_VitD","BMD_Hip","BMD_Lumbar")])
colnames(A) <- c("ANML_1_200","Age","Sex","Ethnicity","BMI","WHR","SBP","DBP","Heart_Rate","Total_Cholesterol","HDL","LDL","Triglyceride","Glucose","HbA1c","Insulin","Uric_Acid","Na","K","Cl","Urea","Creatinine","ALT","GGT","Bilirubin","Albumin","CRP","Vitamin_D","BMD_Hip","BMD_Lumbar")
A$Sex <- as.factor(A$Sex)
A$Ethnicity <- as.factor(A$Ethnicity)
A <- na.omit(A)
A <- as.data.frame(A)

library(Boruta)
set.seed(1)
Boruta1 <- Boruta(ANML_1_200 ~., data = A, doTrace = 2, maxRuns = 50000)
Boruta1
Trace <- as.data.frame(Boruta1$ImpHistory)
Trace <- Trace[,-c(30:32)]

DF <- data.frame(Importance = unlist(as.vector(Trace)), Phenotype = rep(colnames(Trace), each = 29803))
DF <- DF[-which(DF$Importance == -Inf),]

DF$Direction <- rep(NA, nrow(DF))
DF$Direction[which(DF$Phenotype %in% c("Age","Ethnicity","BMI","WHR","SBP","DBP","Heart_Rate","Total_Cholesterol","LDL","Triglyceride","Glucose","HbA1c","Insulin","Uric_Acid","K","Urea","ALT","GGT","Albumin","CRP","BMD_Hip"))] <- "Negative"
DF$Direction[which(DF$Phenotype %in% c("Sex","HDL","Na","Cl","Creatinine","Bilirubin","Vitamin_D","BMD_Lumbar"))] <- "Positive"

library(dplyr)
DF2 <- DF %>% group_by(Phenotype) %>% summarise(Median = median(Importance))
DF2 <- arrange(DF2, Median)            

#Importance Boxplot
library(ggplot2)
ggplot(DF, aes(x = Phenotype, y = Importance, fill = Direction)) +
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(limits = DF2$Phenotype) +
  xlab("Phenotype") +
  ylab("Importance") +
  ggtitle("ANML_1_200 Prediction\nBoruta Importance") +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 10, color = "black")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 12, face = "bold"))
rm(Boruta1)
rm(A)
rm(DF)
rm(Trace)

#1:20000
A <- cbind(ANML_Scaling_Factors[,11], Phenotype[,c("FREG8_Age","FREG7_Gender","FREG5_Ethnic_Group","BMI","WHR","SBP","DBP","HR","DLAB15_Tc_Mmol_L","DLAB17_Hdl_Mmol_L","DLAB19_Ldl_Mmol_L","DLAB21_Trig_Mmol_L","DLAB25_Gluf_Mmol_L","DLAB27_Hba1C_Percent","DLAB81_Insulin","DLAB28_Ua_Mmol_L","DLAB30_Na_Mmol_L","DLAB31_K_Mmol_L","DLAB32_Cl_Mmol_L","DLAB33_Urea_Mmol_L","DLAB35_Creat_Umol_L","DLAB75_ALT_U_L","DLAB76_GGT_U_L","DLAB77_Bilirubin_umol_L","DLAB80_Albumin_mg_dL","DLAB83_CRP","DLAB82_VitD","BMD_Hip","BMD_Lumbar")])
colnames(A) <- c("ANML_1_20000","Age","Sex","Ethnicity","BMI","WHR","SBP","DBP","Heart_Rate","Total_Cholesterol","HDL","LDL","Triglyceride","Glucose","HbA1c","Insulin","Uric_Acid","Na","K","Cl","Urea","Creatinine","ALT","GGT","Bilirubin","Albumin","CRP","Vitamin_D","BMD_Hip","BMD_Lumbar")
A$Sex <- as.factor(A$Sex)
A$Ethnicity <- as.factor(A$Ethnicity)
A <- na.omit(A)
A <- as.data.frame(A)

library(Boruta)
set.seed(1)
Boruta1 <- Boruta(ANML_1_20000 ~., data = A, doTrace = 2, maxRuns = 50000)
Boruta1
Trace <- as.data.frame(Boruta1$ImpHistory)
Trace <- Trace[,-c(30:32)]

DF <- data.frame(Importance = unlist(as.vector(Trace)), Phenotype = rep(colnames(Trace), each = 1722))
DF <- DF[-which(DF$Importance == -Inf),]

DF$Direction <- rep(NA, nrow(DF))
DF$Direction[which(DF$Phenotype %in% c("Age","Ethnicity","BMI","WHR","SBP","DBP","Heart_Rate","Total_Cholesterol","LDL","Triglyceride","Glucose","HbA1c","Insulin","Uric_Acid","K","Urea","ALT","GGT","Albumin","CRP","BMD_Hip"))] <- "Negative"
DF$Direction[which(DF$Phenotype %in% c("Sex","HDL","Na","Cl","Creatinine","Bilirubin","Vitamin_D","BMD_Lumbar"))] <- "Positive"

library(dplyr)
DF3 <- DF %>% group_by(Phenotype) %>% summarise(Median = median(Importance))
DF3 <- arrange(DF3, Median)            

#Importance Boxplot
library(ggplot2)
ggplot(DF, aes(x = Phenotype, y = Importance, fill = Direction)) +
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(limits = DF3$Phenotype) +
  xlab("Phenotype") +
  ylab("Importance") +
  ggtitle("ANML_1_20000 Prediction\nBoruta Importance") +
  scale_fill_manual(values = c("Positive" = "#BC3C29FF", "Negative" = "#0072B5FF")) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 10, color = "black")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 12, face = "bold"))
rm(Boruta1)
rm(A)
rm(DF)
rm(Trace)

#Summary of Ranking
DF1 <- as.data.frame(DF1)
rownames(DF1) <- DF1$Phenotype
DF1$Rank <- rev(1:nrow(DF1))

DF2 <- as.data.frame(DF2)
rownames(DF2) <- DF2$Phenotype
DF2$Rank <- rev(1:nrow(DF2))

DF3 <- as.data.frame(DF3)
rownames(DF3) <- DF3$Phenotype
DF3$Rank <- rev(1:nrow(DF3))

DF2 <- DF2[rownames(DF1),]
DF3 <- DF3[rownames(DF1),]
DF <- as.data.frame(cbind(DF1, DF2, DF3))
rm(DF1, DF2, DF3)

DF <- DF[,-c(2,4,5,7,8)]
rownames(DF) <- NULL
colnames(DF) <- c("Phenotype","Rank_1_5","Rank_1_200","Rank_1_20000")
DF$Mean_Rank <- (DF$Rank_1_5 + DF$Rank_1_200 + DF$Rank_1_20000)/3
DF$Weighted_Mean_Rank <- (DF$Rank_1_5 * 8574 + DF$Rank_1_200 * 1890 + DF$Rank_1_20000 * 211)/10675
DF <- arrange(DF, Weighted_Mean_Rank)

#Heatmap
library(pheatmap)
library(RColorBrewer)
mycolors <- brewer.pal(11, "RdYlBu")
mycolors <- mycolors[c(7:11)]
bk <- unique(c(seq(0, 30, length = 100)))
row.names(DF) <- DF$Phenotype
DF <- DF[,-1]
pheatmap(DF,
         breaks = bk,
         scale = "none",
         display_numbers = T,
         number_format = "%.0f",
         fontsize_number = 10,
         number_color = "black",
         fontsize = 12,
         border_color = "white",
         color = colorRampPalette(mycolors)(100),
         show_rownames = T, show_colnames = T,
         cellwidth = 25, cellheight = 25,
         cluster_cols = F,
         cluster_rows = F,
         main = "Importance Ranking of Phenotypes")
rm(DF)
rm(bk)
rm(mycolors)


rm(ANML_Scaling_Factors, Phenotype)
gc()
