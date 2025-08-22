#pre-ANML Simulation Replacement


#Original adat Data

#Basic Setting
Sys.setenv(LANGUAGE = "en")

#Package
library(SomaDataIO)

#Data Loading
my_adat <- read_adat("SS-2453461_v5.0_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.20240521.adat")


#Simulated Data
Simulation <- read.csv("Simulated_Protein_BMI.csv")
rownames(Simulation) <- Simulation$X
Simulation <- Simulation[,-1]

#Processed pre-ANML Data

#SomaLogic Data
Dat <- read.csv("SomaLogic_pre-anml_Preprocessing/Somalogic_Merged_All.csv")
rownames(Dat) <- Dat$UniqueID

#Remove Control Samples
Dat <- Dat[which(Dat$SampleType == "Sample"),]

#500 Samples
Index <- intersect(which(Dat$tech_rep_id == "N"), which(Dat$bio_rep_id == "N"))
Index <- union(Index, which(Dat$tech_rep == 1))
Dat <- Dat[Index,]
rm(Index)
rownames(Dat) <- Dat$FREG0_PID

#Alignment
Dat <- Dat[rownames(Simulation),]
rownames(Simulation) <- paste(Dat$PlateId, Dat$well, sep = "_")


#Value Replacement
for(i in 1:nrow(Simulation)){
  for(j in 1:ncol(Simulation)){
    Position <- strsplit(rownames(Simulation)[i], "_")
    Index <- which(my_adat$PlateId == Position[[1]][1])
    Index1 <- which(my_adat$PlatePosition == Position[[1]][2])
    Index1 <- intersect(Index, Index1)
    Index2 <- which(colnames(my_adat) == colnames(Simulation)[j])
    my_adat[Index1, Index2] <- Simulation[i, j]
  }
}
rm(i, j, Position, Index, Index1, Index2)

write_adat(my_adat, "Simulation_Data_1.adat")

rm(my_adat)


#Double Check

#Data Loading
my_adat <- read_adat("Simulation_Data_1.adat")

#Extracting
rownames(my_adat) <- paste(my_adat$PlateId, my_adat$PlatePosition, sep = "_")
my_adat <- my_adat[which(my_adat$SubjectID %in% Dat$SOMA_ID), colnames(Simulation)]
my_adat <- my_adat[rownames(Simulation),]

#Comparison
sum(abs(as.matrix(Simulation) - as.matrix(my_adat)))

rm(my_adat)
rm(Simulation)
rm(Dat)

