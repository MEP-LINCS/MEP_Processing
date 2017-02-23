#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
#1/18/2017

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

