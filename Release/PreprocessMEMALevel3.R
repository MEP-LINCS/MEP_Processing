#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 2/2017


library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(RUVnormalize)
library(ruv)
library(stringr)
library(tidyr)
library(readr)
library(dplyr)

processLevel3CommandLine <- function(x, path, k=256, verbose="FALSE"){
  if(length(x)<2) stop("There must be studyName and path arguments in the command line call.")
  studyName <- x[1]
  path <- x[2]
  if((length(x)>2)) k <- x[3]
  if((length(x)>3)) verbose <- x[4]
  list(studyName, path, k, verbose)
}

#callParams <- processLevel3CommandLine(c("MCF10A_DMSO_2", "/lincs/share/lincs_user",256,TRUE))
#callParams <- processLevel3CommandLine(c("MCF10A_Neratinib_2", "/lincs/share/lincs_user",256,TRUE))
#callParams <- processLevel3CommandLine(c("HMEC240L_SS4", "/lincs/share/lincs_user",256,TRUE))
#callParams <- processLevel3CommandLine(c("HMEC122L_SS4", "/lincs/share/lincs_user",256,TRUE))
callParams <- processLevel3CommandLine(commandArgs(trailingOnly = TRUE))
studyName <-callParams[[1]]
path <-callParams[[2]]
k <- as.integer(callParams[[3]])
verbose <- as.logical(callParams[[4]])
startTime <- Sys.time()

#Read the annotated data for all plates in the study
slDT <- getSpotLevelData(studyName, path)

signalsMinimalMetadata <- grep("_SE",grep("_CP_|_PA_|Barcode|^Well$|^Spot$|^PrintSpot$|^Ligand$|^ECMp$|^Drug$|^ArrayRow$|^ArrayColumn$|^CellLine$",colnames(slDT), value=TRUE), value=TRUE, invert=TRUE)

#RUVLoess normalize all signals
if(!k==0){
  if(verbose)  message(paste("Normalizing", studyName,"\n"))
  slDT <- slDT[!is.na(slDT$Cytoplasm_CP_AreaShape_Compactness)&!is.na(slDT$Cytoplasm_CP_AreaShape_Eccentricity),]
  nDT <- normRUVLoessResiduals(slDT[,signalsMinimalMetadata, with = FALSE], k)
  nDT$NormMethod <- "RUVLoessResiduals"
  slDT$k <- k
  slDT <- merge(slDT, nDT, by = c("BW","PrintSpot"))
} else {
  slDT$NormMethod <- "none"
  slDT$k <- k
}

#Add QA flags to the data
slDT <- QASpotLevelData(slDT, lowSpotCellCountThreshold=5,
                        lowRegionCellCountThreshold = 0.4,
                        lowWellQAThreshold = .7)

if(verbose) message(paste("Writing level 3 file to disk\n"))
fwrite(data.table(slDT), paste0(path, "/study/",studyName, "/Annotated/", studyName,"_Level3.tsv"), sep = "\t", quote=FALSE)

#Write the File Annotations for Synapse to tab-delimited file
annotations <- fread(paste0(path,"/",unique(slDT$Barcode)[1],"/Analysis/",unique(slDT$Barcode)[1],"_Level2Annotations.tsv"),header = FALSE)

write.table(c(
  CellLine = annotations$V2[annotations$V1=="CellLine"],
  Preprocess = annotations$V2[annotations$V1=="Preprocess"],
  DataType = annotations$V2[annotations$V1=="DataType"],
  Consortia = annotations$V2[annotations$V1=="Consortia"],
  Drug = annotations$V2[annotations$V1=="Drug"],
  Segmentation = annotations$V2[annotations$V1=="Segmentation"],
  StainingSet = annotations$V2[annotations$V1=="StainingSet"],
  Level = "3"),
  paste0(path,"/study/",studyName, "/Annotated/", studyName,"_","Level3Annotations.tsv"), sep = "\t",col.names = FALSE, quote=FALSE)

message(paste("Elapsed time to normalize ",studyName, Sys.time()-startTime, "\n"))


