#Load up MCF10ADMSO level 3 dataset
#Read in the data.table that holds every spot in the cell line
#and can be used to index to the images
library("ggplot2")
library(scales)
library("data.table")
library("MEMA")
library("grid")
library("knitr")
library("gplots")
library("RColorBrewer")
library(DT)
library(d3heatmap)
library(plotly)
#library(corrplot)
library(RUVnormalize)
library(ruv)
library(Rtsne)
library(XLConnect)
library(readxl)

source("MEPLINCSFunctions.R")
source("MEP-LINCS QANorm Functions.R")

#Items to be generalized
datasetName <- "MCF10ADMSO"
cellLine <- "MCF10A"
l3C <- fread(unique(paste0("../AnnotatedData/",datasetName, "_SSF","_Level3.txt")), showProgress = FALSE)
#Work with EdU signal only
dt <- l3C[,.(Barcode,MEP,Ligand,ECMp,Well,ImageID,Spot_PA_SpotCellCount,Spot_PA_SpotCellCountLog2,Spot_PA_SpotCellCountLog2RUVLoess,Nuclei_PA_Cycle_DNA2NProportionLogit,Nuclei_PA_Cycle_DNA2NProportionLogitRUVLoess,Nuclei_PA_Gated_EdUPositiveProportionLogit,Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoess)]

#Work with col1 spots only
dtc1 <- dt[grepl("COL1",ECMp),]
#remove lowest cell count spots due to high variance
dtc1f <- dtc1[dtc1$Spot_PA_SpotCellCount>2^5,]
#Use FBS replicate wells
dtc1fr <- dtc1f[grepl("FBS",dtc1f$Ligand),]
#Model proliferation vs cell count
dtc1frM <- lm(Nuclei_PA_Gated_EdUPositiveProportionLogit~poly(Spot_PA_SpotCellCountLog2,2)+Barcode, data=dtc1fr)
summary(dtc1frM)
plot(predict(dtc1frM,newdata = data.frame(Spot_PA_SpotCellCountLog2=dtc1fr$Spot_PA_SpotCellCountLog2[order(dtc1fr$Spot_PA_SpotCellCountLog2)])))
plot(x=dtc1fr$Spot_PA_SpotCellCountLog2,y=dtc1fr$Nuclei_PA_Gated_EdUPositiveProportionLogit)

#Add backtransformed values for some logit transformed signals
dt <- dt[,Nuclei_PA_Cycle_DNA2NProportionLogitRUVLoessBacktransformed:= btLogit(Nuclei_PA_Cycle_DNA2NProportionLogitRUVLoess)]
dt <- dt[,Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoessBacktransformed:= btLogit(Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoess)]
dt <- addOmeroIDs(dt)

#Build models at the MEP,ligand and ECMp level of logit(EdU+ Proportion)~spot cell count
SCCModel <- lm(Nuclei_PA_Gated_EdUPositiveProportionLogit~poly(Spot_PA_SpotCellCountLog2,4), data=dt)
summary(SCCModel)
#MEPModel <- glm(Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoess~Spot_PA_SpotCellCountLog2+MEP, data=dt)

ligandECMpModel <- lm(Nuclei_PA_Gated_EdUPositiveProportionLogit~poly(Spot_PA_SpotCellCountLog2,4)+Ligand+ECMp, data=dt)
summary(ligandECMpModel)
#step(ligandECMpModel)
ECMpModel <- lm(Nuclei_PA_Gated_EdUPositiveProportionLogit~poly(Spot_PA_SpotCellCountLog2,4)+ECMp, data=dt)
summary(ECMpModel)
#Plot data and models
#Evaluate models with R squared or deviance
#Cluster models and look for outliers