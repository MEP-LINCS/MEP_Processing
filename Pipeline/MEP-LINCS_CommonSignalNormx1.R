#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 3/2017

library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(RUVnormalize)
library(ruv)
library(stringr)
library(tidyr)
library(readr)
library(dplyr)

#Command line args are positional and gathered as a character lost
#
processCommonSignalCommandLine <- function(x, path, k=256, verbose="FALSE"){
  if(length(x)<2) stop("There must be studyName List and path arguments in the command line call.")
  studyNameList <- x[1]
  path <- x[2]
  if((length(x)>2)) k <- x[3]
  if((length(x)>3)) verbose <- x[4]
  list(studyNameList, path, k, verbose)
}

#callParams <- processCommonSignalCommandLine(c("HMEC240L_SS1,HMEC240L_SS4", "/lincs/share/lincs_user",256,TRUE))
callParams <- processCommonSignalCommandLine(commandArgs(trailingOnly = TRUE))
studyNameList <-lapply(unlist(str_split(callParams[[1]],",")), function(x){
  return(x)
})
path <-callParams[[2]]
k <- as.integer(callParams[[3]])
verbose <- as.logical(callParams[[4]])
startTime <- Sys.time()

#RUV and loess normalize the common DAPI signals
l3 <- preprocessCommonSignalsLevel3(studyNameList, path, k, verbose)
#Get the level 3 raw and normalized data
l3L <- lapply(studyNameList, function(studyName){
  dt <- fread(unique(paste0(path,"/study/",studyName,"/Annotated/",studyName,"_Level3.tsv")), showProgress = FALSE)
  return(dt)
})
#Make a list of datatables containing the staining set specific columns with BW and Spot
sssL <- lapply(l3L, function(dt, l3){
  specificColumns <- c("BW","Spot",setdiff(colnames(dt),colnames(l3)))
  dt[,specificColumns, with=FALSE]
},l3=l3)
sssDT <- rbindlist(sssL, fill=TRUE)
l3C <-sssDT[l3,on=c("BW","Spot")]

#Write the normalized level3  data to disk
message("Writing level 3 combined data to disk\n")
studyNameSSC <- gsub("_.*","_SSC",studyNameList[[1]])
write.table(l3C, paste0(path,"/study/",studyNameSSC,"/Annotated/",studyNameSSC,"_Level3.tsv"), sep = "\t",row.names = FALSE, quote=FALSE)

#Summarize to the MEP_Drug_Drug1Conc level
mepDT <- preprocessLevel4(l3C[,grep("Endpoint|StainingSet|395nm|488nm|555nm|640nm|750nm|1An$|2An$|3An$|LigandSet",colnames(l3C),value=TRUE,invert=TRUE), with=FALSE],seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin"))

#Add in the barcodes for each MEP_Drug
mepDT <- addBarcodes(dt3 = l3C, dt4 = mepDT)

# Add a QA flag for spots with few replicates
mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < 3

message("Writing level 4 file to disk\n")
write.table(mepDT, paste0(path,"/study/",studyNameSSC,"/Annotated/",studyNameSSC,"_Level4.tsv"), sep = "\t",row.names = FALSE, quote=FALSE)

message("Running time:",Sys.time()-startTime)
