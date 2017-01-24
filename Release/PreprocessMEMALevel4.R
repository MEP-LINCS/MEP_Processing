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
preprocessLevel4 <- function(dt, seNames=NULL){
  #Add a count of replicates
  dt <- dt[,Spot_PA_ReplicateCount := .N,by="Ligand,ECMp,Drug,CellLine"]
  rawSignalNames <- grep("_SE",grep("Log2|Logit|_PA_|Intensity|AreaShape",colnames(dt), value=TRUE), value=TRUE, invert=TRUE)
  l4Signals<- dt[,lapply(.SD, numericMedian), by="Ligand,ECMp,Drug,CellLine", .SDcols=rawSignalNames]
  
  #Use seNames to select the parameters that get SE values
  if(!is.null(seNames)){
    seNamesPattern<-paste(seNames,collapse="|")
    seSignalNames <- grep(seNamesPattern,rawSignalNames,value=TRUE)
    l4Ses <- dt[,lapply(.SD,MEMA:::se),keyby="Ligand,ECMp,Drug,CellLine", .SDcols=seSignalNames]
  } else{
    l4Ses <- dt[,lapply(.SD,MEMA:::se),keyby="Ligand,ECMp,Drug,CellLine"]
  }
  
  #Add _SE to the standard error column names
  setnames(l4Ses, grep("Ligand|ECMp|Drug|CellLine",colnames(l4Ses), value = TRUE, invert = TRUE), paste0(grep("Ligand|ECMp|Drug|CellLine",colnames(l4Ses), value = TRUE, invert = TRUE),"_SE"))
  
  #Merge back in the replicate metadata
  metadataNames <- grep("_SE|Barcode|^BW$|ArrayRow|ArrayColumn|^Well$|^Spot$|^PrintSpot$|^Well_Ligand$|ImageID|QA_|ECMSet",colnames(dt), value=TRUE,invert=TRUE) %>%
    setdiff(rawSignalNames)
  mdDT <- unique(dt[,metadataNames, with=FALSE])
  setkey(l4Signals,Ligand,ECMp,Drug,CellLine)
  setkey(l4Ses,Ligand,ECMp,Drug,CellLine)
  setkey(mdDT,Ligand,ECMp,Drug,CellLine)
  l4DT <- merge(mdDT,merge(l4Signals,l4Ses))
  return(l4DT)
}#End of preprocesslevel4

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
  library(tidyr)
  
  #Set a threshold for the lowSpotReplicates flag
  lowReplicateCount <- 3
  
  l3DT <- fread(paste0(path,"/",datasetName,"/Annotated/",datasetName,"_Level3.tsv"))
  
  mepDT <- preprocessLevel4(l3DT,seNames = seNames)
  
  mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < lowReplicateCount
  #Add in the barcodes the MEPs came from
  barcodesList <- lapply(unique(mepDT$MEP_Drug), function(m){
    Barcodes=unique(l3DT$Barcode[l3DT$MEP_Drug==m])
  })
  bmdDT <- data.table(Barcode = barcodesList, MEP_Drug=unique(mepDT$MEP_Drug))
  mepDT <- merge(bmdDT,mepDT,by="MEP_Drug")
  
  #WriteData
  if(writeFiles){
    if(verbose) cat("Writing level 4 file to disk\n")
    fwrite(data.table(format(mepDT, digits = 4, trim=TRUE)), paste0(path, "/",datasetName, "/Annotated/", datasetName,"_Level4.tsv"), sep = "\t", quote=FALSE)
    
    #Write the File Annotations for Synapse to tab-delimited file
    annotations <- fread(paste0(path,"/",datasetName,"/Annotated/",datasetName,"_Level3Annotations.tsv"),header = FALSE)
    
    write.table(c(
      CellLine = annotations$V2[annotations$V1=="CellLine"],
      Preprocess = annotations$V2[annotations$V1=="Preprocess"],
      DataType = annotations$V2[annotations$V1=="DataType"],
      Consortia = annotations$V2[annotations$V1=="Consortia"],
      Drug = annotations$V2[annotations$V1=="Drug"],
      Segmentation = annotations$V2[annotations$V1=="Segmentation"],
      StainingSet = annotations$V2[annotations$V1=="StainingSet"],
      Level = "4"),
      paste0(path,"/",datasetName, "/Annotated/", datasetName,"_","Level4Annotations.tsv"), sep = "\t",col.names = FALSE, quote=FALSE)
  }
  cat("Elapsed time:", Sys.time()-startTime, "\n")
}

path <- commandArgs(trailingOnly = TRUE)[1]
datasetName <- commandArgs(trailingOnly = TRUE)[2]
res <- preprocessMEMALevel4(datasetName, path, verbose=TRUE)

