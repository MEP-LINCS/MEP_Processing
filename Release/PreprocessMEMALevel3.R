#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 2/1/2017

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
  barcodes <- str_split(ptData[["Plate IDs"]][ptData[["Study Name"]]==datasetName], ",") %>%
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

