#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 9/8/2016

##Introduction
library("parallel")#use multiple cores for faster processing
source("MEP_LINCS/Release/MEPLINCSFunctions.R")

#' Create an M matrix for the RUV normalization
#' 
#' Each row is for a unit to be normalized and
#' each column is a unique replicate type
#' Each row will have a 1 to indicate the replicate type
#' all other values will be 0
#' @param dt The datatable to be normalized that has columns Barcode, Well and 
#' those in the replicateCols parameter.
#' @param replicateCols A character vector of column names in dt that define a replicate well.
#' @return A matrix suitable for defining the dataset structure in RUV normalization
#' @export
createRUVM <- function(dt,replicateCols=c("CellLine","Ligand","Drug"))
{
  #Add a column that defines what makes a well a replicate
  dt$ReplicateID <- do.call(paste, c(dt[,replicateCols, with=FALSE], sep="_"))
  #Add a similar column that binds in the barcode and well locations
  dt$UnitID <- do.call(paste, c(dt[,c("Barcode","Well",replicateCols), with=FALSE], sep="_"))
  #Set up the M Matrix to denote replicate ligand wells
  nrUnits <- length(unique(dt$UnitID[dt$SignalType=="Signal"]))
  nrReplicateIDs <- length(unique(dt$ReplicateID[dt$SignalType=="Signal"]))
  M <-matrix(0, nrow = nrUnits, ncol = nrReplicateIDs)
  rownames(M) <- unique(dt$UnitID[dt$SignalType=="Signal"])
  colnames(M) <- unique(dt$ReplicateID[dt$SignalType=="Signal"])
  rownames(M) <- gsub("[|]","pipe",rownames(M))
  colnames(M) <- gsub("[|]","pipe",colnames(M))
  #Indicate the replicate ligands
  for(replicate in colnames(M)){
    #Put a 1 in the rownames that contain the column name
    M[grepl(replicate,rownames(M)),colnames(M)==replicate] <- 1
  }
  rownames(M) <- gsub("pipe","|",rownames(M))
  colnames(M) <- gsub("pipe","|",colnames(M))
  return(M)
}

#' Apply RUV normalization on a signal and its residuals
#' 
#' Assumes there are signal values in the first half of each row
#' and residuals in the second half
#' @export
RUVIIIArrayWithResiduals <- function(k, Y, M, cIdx, signalName, verboseDisplay=FALSE){
  YRUVIII <- RUVIII(Y, M, cIdx, k)
  nY <- YRUVIII[["newY"]]
  #Remove residuals
  nY <- nY[,1:(ncol(nY)/2)]
  #melt matrix to have Spot and Ligand columns
  nYm <- melt(nY, varnames=c("BW","PrintSpot"), value.name=signalName)
  nYm$BW <- as.character(nYm$BW)
  if(verboseDisplay){
    return(list(nYm=data.table(nYm), fullAlpha=YRUVIII[["fullalpha"]], W=RUVIII[["W"]]))
  }
  return(data.table(nYm))
}

#' Generate a matrix of signals and residuals
#' 
#' This function reorganizes a data.table into a matrix suitable for
#' RUV normalization.
#' 
#'@param dt a data.table with columns named BW and PrintSpot followed by one signal column 
#'@return A numeric matrix with BW rows and two sets of columns. The second
#'set of columns are the residuals from the medians of each column and have
#'"_Rel" appended to their names.
#' @export
signalResidualMatrix <- function(dt){
  signalName <- colnames(dt)[ncol(dt)]
  if(grepl("Logit", signalName)){
    fill <- log2(.01/(1-.01))
  } else if(grepl("Log", signalName)){
    fill <- log2(.001)
  } else {
    fill <- 0
  }
  
  dts <- data.table(dcast(dt[dt$SignalType=="Signal",], BW~PrintSpot, value.var=signalName, fill=fill, na.rm=TRUE))
  dtr <- data.table(dcast(dt[dt$SignalType=="Residual",], BW~PrintSpot, value.var=signalName, fill=fill, na.rm=TRUE))
  rowNames <- dts$BW
  dts <- dts[,BW := NULL]
  dtr <- dtr[,BW:=NULL]
  setnames(dtr,colnames(dtr),paste0(colnames(dtr),"_Rel"))
  dtsr <- cbind(dts,dtr)
  srm <- matrix(unlist(dtsr),nrow=nrow(dtsr))
  rownames(srm) <- rowNames
  colnames(srm) <- colnames(dtsr)
  return(srm)
}

#' Loess normalize values within an array
#'@export
loessNormArray <- function(dt){
  #Identify the Signal name
  signalName <- unique(dt$SignalName)
  setnames(dt,signalName,"Value")
  #Get the median of the replicates within the array
  dt <- dt[,mel := median(Value), by=c("BW","ECMp","Drug")]
  #Get the residuals from the spot median
  dt <- dt[,Residual := Value-mel]
  #Subtract the loess model of each array's residuals from the signal
  dt <- dt[, ValueLoess:= loessNorm(Value,Residual,ArrayRow,ArrayColumn), by="BW"]
  setnames(dt,"ValueLoess", paste0(signalName,"Loess"))
  BW <- "BW"
  PrintSpot <- "PrintSpot"
  dt <- dt[,c(paste0(signalName,"Loess"),BW,PrintSpot), with=FALSE]
  setkey(dt,BW,PrintSpot)
  return(dt)
}

#' Loess normalize an array using the spatial residuals
loessNorm <- function(Value,Residual,ArrayRow,ArrayColumn){
  dt <-data.table(Value=Value,Residual=Residual,ArrayRow=ArrayRow,ArrayColumn=ArrayColumn)
  lm <- loess(Residual~ArrayRow+ArrayColumn, dt, span=.7)
  dt$ResidualLoess<-predict(lm)
  dt <- dt[,ValueLoess := Value-ResidualLoess]
  return(ValueLoess = dt$ValueLoess)
}

#'Apply RUV and Loess Normalization to the signals in a dataset
#' @export
normRUVLoessResiduals <- function(dt, k){
  setkey(dt,CellLine,Barcode,Well,Ligand,Drug,ECMp)
  metadataNames <- "^CellLine$|Barcode|^Well$|^Spot$|^PrintSpot$|ArrayRow|ArrayColumn|^ECMp$|^Ligand$|^Drug$"
  signalNames <- grep(metadataNames,colnames(dt),invert=TRUE, value=TRUE)
  
  #Add residuals from subtracting the biological medians from each value
  residuals <- dt[,lapply(.SD,calcResidual), by="CellLine,Barcode,Well,Ligand,Drug,ECMp", .SDcols=signalNames]
  #Add within array location metadata
  residuals$Spot <- as.integer(dt$Spot)
  residuals$PrintSpot <- as.integer(dt$PrintSpot)
  residuals$ArrayRow <- dt$ArrayRow
  residuals$ArrayColumn <- dt$ArrayColumn
  #Create a signal type
  dt$SignalType <- "Signal"
  residuals$SignalType <- "Residual"
  srDT <- rbind(dt,residuals)
  
  #Add to carry metadata into matrices
  srDT$BWLD <- paste(srDT$Barcode, srDT$Well, srDT$Ligand,  srDT$Drug, sep="_") 
  
  #Create the M matrix which denotes replicates
  M <- createRUVM(srDT, replicateCols=c("CellLine","Ligand","Drug"))
  
  #Add BWL but need to check if this can be only BW
  srDT$BW <- paste(srDT$Barcode, srDT$Well, sep="_") 
  #Make a list of matrices that hold signal and residual values
  srmList <- lapply(signalNames, function(signalName, dt){
    srm <- signalResidualMatrix(dt[,.SD, .SDcols=c("BW", "PrintSpot", "SignalType", signalName)])
    return(srm)
  },dt=srDT)
  names(srmList) <- signalNames
  
  #Make a list of matrices of RUV normalized signal values
  srmRUVList <- lapply(names(srmList), function(srmName, srmList, M, k){
    Y <- srmList[[srmName]]
    #Hardcode in identification of residuals as the controls
    resStart <- ncol(Y)/2+1
    cIdx=resStart:ncol(Y)
    nY <- RUVIIIArrayWithResiduals(k, Y, M, cIdx, srmName) #Normalize the spot level data
    nY$SignalName <- paste0(srmName,"RUV")
    setnames(nY,srmName,paste0(srmName,"RUV"))
    #nY[[srmName]] <- as.vector(Y[,1:(resStart-1)]) #Add back in the raw signal (may not be needed)
    return(nY)
  }, srmList=srmList, M=M, k=k)
  
  #Reannotate with ECMp, Drug, ArrayRow and ArrayColumn as needed for loess normalization
  ECMpDT <- unique(srDT[,list(Well,PrintSpot,Spot,ECMp,Drug, ArrayRow,ArrayColumn)])
  srmERUVList <- lapply(srmRUVList, function(dt,ECMpDT){
    dt$Well <- gsub(".*_","",dt$BW)
    setkey(dt,Well,PrintSpot)
    setkey(ECMpDT,Well,PrintSpot)
    dtECMpDT <- merge(dt,ECMpDT)
    return(dtECMpDT)
  },ECMpDT=ECMpDT)
  
  #Add Loess normalized values for each RUV normalized signal
  RUVLoessList <- lapply(srmERUVList, function(dt){
    dtRUVLoess <- loessNormArray(dt)
  })
  
  #Combine the normalized signals and metadata
  signalDT <- Reduce(merge,RUVLoessList)

  return(signalDT)
}



preprocessMEMALevel3 <- function(datasetName, path, k= 256, verbose=FALSE){
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
  
  #Set a threshold for the lowSpotCellCount flag
  lowSpotCellCountThreshold <- 5
  
  #Set a threshold for the lowRegionCellCount flag
  lowRegionCellCountThreshold <- .4
  
  #Set a threshold for lowWellQA flag
  lowWellQAThreshold <- .7
  
  #Set a threshold for the lowSpotReplicates flag
  lowReplicateCount <- 3
  
  #Download Plate Tracker google sheet
  ptObj <- gs_url("https://docs.google.com/spreadsheets/d/1QefmAsK2B_no3iL-epx198pP_HVJegOlku_nliiy3eg/pubhtml")
  ptData <- gs_download(ptObj,to="PlateTracker.csv", overwrite = TRUE, verbose=FALSE) %>%
    read_csv()
  #ptData <- read_csv(gs_download(ptObj,to="PlateTracker.csv", overwrite = TRUE, verbose=FALSE))
  barcodes <- str_split(ptData[["Plate IDs"]][ptData[["Expt Name"]]==datasetName], ",") %>%
    unlist()
  
  spotDTL <- mclapply(barcodes, function(barcode, path){
    sd <- fread(paste0(path,"/",barcode,"/Analysis/",barcode,"_SpotLevel.tsv"))
  }, path=path, mc.cores=detectCores())
  slDT <- rbindlist(spotDTL)
  slDT$BW <- paste(slDT$Barcode,slDT$Well,sep="_")
  rm(spotDTL)
  gc()
  
  slDT <- slDT[!grepl("fiducial|Fiducial|gelatin|blank|air|PBS",slDT$ECMp),]
  #Debug
  slDT <- slDT[,Cytoplasm_CP_AreaShape_MaximumRadiusLog2 :=NULL]
  
  rawSignalsMinimalMetadataRXP <- grep("_SE",grep("Log2|Logit|Barcode|^Well$|^Spot$|^PrintSpot$|^Ligand$|^ECMp$|^Drug$|^ArrayRow$|^ArrayColumn$|^CellLine$",colnames(slDT), value=TRUE), value=TRUE, invert=TRUE)
  
  #Normalize each feature, pass with location and content metadata
  if(verbose) cat("Normalizing\n")
  if(!k==0){
    #Normalize each feature, pass with location and content metadata
    if(verbose)  cat("Normalizing", datasetName,"\n")
    #Debug: Add 96 well, drug and non-FBS replicate awareness to normRUVLoessResiduals
    nDT <- normRUVLoessResiduals(slDT[,rawSignalsMinimalMetadataRXP, with = FALSE], k)
    nDT$NormMethod <- "RUVLoessResiduals"
    
    #nDT has normalized RUVLoess values
    #Merge the normalized data with the original data
    setkey(nDT,BW,PrintSpot)
    setkey(slDT,BW,PrintSpot)
    #merge in the raw data to the transformed and RUVLoess 
    slDT <- merge(slDT, nDT)
  } else {
    slDT$NormMethod <- "none"
    slDT$k <- k
  }
  
  #setkey(slDT, Barcode, Well, Spot, PrintSpot, ArrayRow,ArrayColumn,ECMp,Ligand)
  
  #Label FBS with their plate index to keep separate
  slDT$Ligand[grepl("FBS",slDT$Ligand)] <- paste0(slDT$Ligand[grepl("FBS",slDT$Ligand)],"_P",match(slDT$Barcode[grepl("FBS",slDT$Ligand)], unique(slDT$Barcode)))
  
  #Write QA flags into appropriate data levels
  #Low cell count spots
  slDT$QA_LowSpotCellCount <- slDT$Spot_PA_SpotCellCount < lowSpotCellCountThreshold
  
  #Low quality DAPI
  slDT$QA_LowDAPIQuality <- FALSE
  
  #Flag spots below automatically loess QA threshold
  slDT$QA_LowRegionCellCount <- slDT$Spot_PA_LoessSCC < lowRegionCellCountThreshold
  
  #Flag wells below automatically calculated QA threshold
  slDT$QA_LowWellQA <- FALSE
  slDT$QA_LowWellQA[slDT$QAScore < lowWellQAThreshold] <- TRUE
  
  #WriteData
  if(writeFiles){
    if(verbose) cat("Writing level 3 file to disk\n")
    fwrite(data.table(format(slDT, digits = 4, trim=TRUE)), paste0(path, "/",datasetName, "/Annotated/", datasetName,"_Level3.tsv"), sep = "\t", quote=FALSE)
    
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
      paste0(path,"/",datasetName, "/Annotated/", datasetName,"_","Level3Annotations.tsv"), sep = "\t",col.names = FALSE, quote=FALSE)
  }
  cat("Elapsed time:", Sys.time()-startTime, "\n")
}

path <- commandArgs(trailingOnly = TRUE)[1]
datasetName <- commandArgs(trailingOnly = TRUE)[2]
res <- preprocessMEMALevel3(datasetName, path, verbose=TRUE)

