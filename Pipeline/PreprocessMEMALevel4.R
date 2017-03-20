#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 2/2017

processLevel4CommandLine <- function(x, path, verbose="FALSE"){
  if(length(x)<2) stop("There must be studyName and path arguments in the command line call.")
  studyName <- x[1]
  path <- x[2]
  if((length(x)>2)) verbose <- x[3]
  list(studyName, path, verbose)
}
library(MEMA) #merge, annotate and normalize functions

#callParams <- processLevel4CommandLine(c("MCF10A_DMSO_2", "/lincs/share/lincs_user/study", TRUE))
#callParams <- processLevel4CommandLine(c("HMEC240L_SS4", "/lincs/share/lincs_user/study",TRUE))
#callParams <- processLevel4CommandLine(c("HMEC122L_SS4", "/lincs/share/lincs_user/study",TRUE))
#callParams <- processLevel4CommandLine(c("MCF10A_Cell_Titrate_96", "/lincs/share/lincs_user/study",TRUE))
callParams <- processLevel4CommandLine(c("MCF10A_Cell_Titrate_8", "/lincs/share/lincs_user/study",TRUE))
#callParams <- processLevel4CommandLine(commandArgs(trailingOnly = TRUE))
studyName <-callParams[[1]]
path <-callParams[[2]]
verbose <- callParams[[3]]
startTime <- Sys.time()

l3DT <- fread(paste0(path,"/",studyName,"/Annotated/",studyName,"_Level3.tsv"))
#Summarize to the MEP_Drug level
mepDT <- preprocessLevel4(l3DT,seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin"))
#Add in the barcodes for each MEP_Drug
mepDT <- addBarcodes(dt3 = l3DT, dt4 = mepDT)
# Add a QA flag for spots with few replicates
mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < 3
if(verbose) message("Writing level 4 file to disk\n")
fwrite(mepDT, paste0(path, "/",studyName, "/Annotated/", studyName,"_Level4.tsv"), sep = "\t", quote=FALSE)
# #Write the File Annotations for Synapse to tab-delimited file
# annotations <- fread(paste0(path,"/",studyName,"/Annotated/",studyName,"_Level3Annotations.tsv"),header = FALSE)
# write.table(c(
#   CellLine = annotations$V2[annotations$V1=="CellLine"],
#   Preprocess = annotations$V2[annotations$V1=="Preprocess"],
#   DataType = annotations$V2[annotations$V1=="DataType"],
#   Consortia = annotations$V2[annotations$V1=="Consortia"],
#   Drug = annotations$V2[annotations$V1=="Drug"],
#   Segmentation = annotations$V2[annotations$V1=="Segmentation"],
#   StainingSet = annotations$V2[annotations$V1=="StainingSet"],
#   Level = "4"),
#   paste0(path,"/",studyName, "/Annotated/", studyName,"_","Level4Annotations.tsv"), sep = "\t",col.names = FALSE, quote=FALSE)
if(verbose) message(paste("Elapsed time for ",studyName, "is", Sys.time()-startTime, "\n"))

