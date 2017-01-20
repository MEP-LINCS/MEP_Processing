#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
#1/18/2017

source("MEP_LINCS/Release/MEPLINCSFunctions.R")

#' Summarize spot level data to the MEP level
#' 
#' Median summarize the spot level normalized values the most biologically 
#' interpretable raw data at each spot, calculate the standard errors and add
#' SE columns for all median summarized data
#' @param l3 The datatable of spot level data to be summarized
#' @return A datatable of MEP level, median summarized data with standard error values and 
#'  metadata
#' @export
preprocessLevel4 <- function(l3, seNames=NULL){
  #Add a count of replicates
  l3 <- l3[,Spot_PA_ReplicateCount := .N,by="Ligand,ECMp"]
  l4Names<-grep("Loess$|RUV|Norm|^Ligand$|^ECMp|Barcode|Spot_PA_SpotCellCount$|Spot_PA_ReplicateCount$", x=names(l3),value=TRUE)
  #remove the _SE values
  l4Names <- grep("_SE|NormMethod|AnnotID",l4Names, value = TRUE, invert = TRUE)
  l4Keep<-l3[,l4Names,with=FALSE]
  l4DT<-l4Keep[,lapply(.SD,numericMedian),keyby="Ligand,ECMp,Barcode"]
  #Use seNames to select the parameters that get SE values
  if(!is.null(seNames)){
    seNamesPattern<-paste(seNames,collapse="|")
    seNames <- grep(seNamesPattern,l4Names,value=TRUE)
    l4DTse <- l4Keep[,lapply(.SD,se),keyby="Ligand,ECMp,Barcode", .SDcols=seNames]
  } else{
    l4DTse <- l4Keep[,lapply(.SD,se),keyby="Ligand,ECMp,Barcode"]
  }
  
  #Add _SE to the standard error column names
  setnames(l4DTse, grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(l4DTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(l4DTse), value = TRUE, invert = TRUE),"_SE"))
  
  l3Names <- grep("Barcode|Well|CellLine|Ligand|ECM|Endpoint488|Endpoint555|Endpoint647|EndpointDAPI|ECMp|MEP|Lx|PinDiameter", colnames(l3), value=TRUE)
  #Merge back in the replicate metadata
  mDT <- l3[,l3Names,keyby="Ligand,ECMp,Barcode", with=FALSE]
  setkey(mDT,Ligand,ECMp,Barcode)
  l4DT <- mDT[l4DT, mult="first"]
  l4DT <- l4DTse[l4DT]
  l4DT <- summarizeFBS(l4DT)
  return(l4DT)
}#End of createl4

#' Summarize cell level data to the spot level
#' 
#' Median summarize the spot level values, calculate the standard errors and add
#' SE columns for all median summarized data
#' @param dt The datatable of spot level data to be summarized
#' @param seNames a character vector of dt column names that will have standard errors calculated
#' @return A datatable of spot-level, median summarized data with standard error values and 
#'  metadata
#' @export
summarizeReplicates <- function(dt, seNames=NULL){
  #Summarize cell data to medians of the spot parameters
  parameterNames<-grep(pattern="(Children|_CP_|_PA_)",x=names(dt),value=TRUE)
  
  mepDT<-dt[,lapply(.SD, numericMedian), by="Barcode,Well,Spot", .SDcols=parameterNames]
  #Use seNames to select the parameters that get SE values
  if(!is.null(seNames)){
    seNamesPattern<-paste(seNames,collapse="|")
    seNames <- grep(seNamesPattern,parameterNames,value=TRUE)
    slDTse <- cDT[,lapply(.SD,se), by="Barcode,Well,Spot", .SDcols=seNames]
  } else{
    slDTse <- cDT[,lapply(.SD,se), by="Barcode,Well,Spot"]
  }
  
  #Add _SE to the standard error column names
  setnames(slDTse, grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE),"_SE"))
  
  #Merge back in the spot and well metadata
  metadataNames <- grep("(Row|Column|PrintOrder|Block|^ID$|Array|CellLine|Ligand|Drug|Endpoint|ECMp|MEP|Well_Ligand|ECM|ImageID|Barcode|^Well$|^PrintSpot$|^Spot$|Pin|Lx)", x=colnames(cDT), value=TRUE)
  setkey(cDT,Barcode, Well,Spot)
  mDT <- cDT[,metadataNames,keyby="Barcode,Well,Spot", with=FALSE]
  slDT <- mDT[slDT, mult="first"]
  #Merge in the standard err values
  setkey(slDTse, Barcode, Well, Spot)
  slDT <- slDTse[slDT]
  #Add a count of replicates
  slDT <- slDT[,Spot_PA_ReplicateCount := .N,by="Ligand,ECMp"]
  
  #Add the loess model of the SpotCellCount on a per well basis
  slDT <- slDT[,Spot_PA_LoessSCC := loessModel(.SD, value="Spot_PA_SpotCellCount", span=.5), by="Barcode,Well"]
  
  #Add well level QA Scores to spot level data
  slDT <- slDT[,QAScore := calcQAScore(.SD, threshold=lthresh, maxNrSpot = max(cDT$ArrayRow)*max(cDT$ArrayColumn),value="Spot_PA_LoessSCC"),by="Barcode,Well"]
  return(slDT)
}


preprocessMEMALevel4 <- function(datasetName, path, verbose=FALSE){
  startTime <- Sys.time()
  writeFiles<-TRUE
  seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")
  
  library(MEMA)#merge, annotate and normalize functions
  library(data.table)#fast file reads, data merges and subsetting
  library(parallel)#use multiple cores for faster processing
  library(RUVnormalize)
  library(ruv)
  library(stringr)
  library(googlesheets)
  library(tidyr)
  library(readr)
  
  #Set a threshold for the lowSpotReplicates flag
  lowReplicateCount <- 3
  
  l3DT <- fread(paste0(path,"/",datasetName,"/Annotated/",datasetName,"_Level3.tsv"))
  
  metadataNames <- grep("(Row|Column|PrintOrder|Block|^ID$|Array|CellLine|Ligand|Drug|Endpoint|ECMp|MEP|Well_Ligand|ECM|ImageID|Barcode|^Well$|^PrintSpot$|^Spot$|Pin|Lx)", x=colnames(l3DT), value=TRUE)
  
  mepDT <- createl4(l3DT,seNames = seNames)
  
  #Level 4
  mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < lowReplicateCount
  
  #WriteData
  if(writeFiles){
    if(verbose) cat("Writing level 4 file to disk\n")
    fwrite(data.table(format(slDT, digits = 4, trim=TRUE)), paste0(path, "/",datasetName, "/Annotated/", datasetName,"_Level4.tsv"), sep = "\t", quote=FALSE)
    
    #Write the File Annotations for Synapse to tab-delimited file
    annotations <- fread(paste0(path,"/",barcodes[1],"/Analysis/",barcodes[1],"_SpotLevelAnnotations.tsv"),header = FALSE)
    
    write.table(c(
      CellLine = annotations$V2[annotations$V1=="CellLine"],
      Preprocess = annotations$V2[annotations$V1=="Preprocess"],
      DataType = annotations$V2[annotations$V1=="DataType"],
      Consortia = annotations$V2[annotations$V1=="Consortia"],
      Drug = annotations$V2[annotations$V1=="Drug"],
      Segmentation = annotations$V2[annotations$V1=="Segmentation"],
      StainingSet = annotations$V2[annotations$V1=="StainingSet"],
      Level = "3"),
      paste0(path,"/",datasetName, "/Annotated/", datasetName,"_","Level4Annotations.tsv"), sep = "\t",col.names = FALSE, quote=FALSE)
  }
  cat("Elapsed time:", Sys.time()-startTime, "\n")
}

path <- commandArgs(trailingOnly = TRUE)[1]
datasetName <- commandArgs(trailingOnly = TRUE)[2]
res <- preprocessMEMALevel4(datasetName, path, verbose=TRUE)

