#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 9/8/2016

##Introduction
library("parallel")#use multiple cores for faster processing
source("MEP_LINCS/Release/MEPLINCSFunctions.R")

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
  rm(spotDTL)
  gc()
  
  slDT <- slDT[!grepl("fiducial|Fiducial|gelatin|blank|air|PBS",slDT$ECMp),]
  #Debug
  slDT <- slDT[,Cytoplasm_CP_AreaShape_MaximumRadiusLog2 :=NULL]
  
  metadataNames <- "ObjectNumber|^Row$|^Column$|Block|^ID$|PrintOrder|Depositions|CellLine|Endpoint|WellIndex|Center|LigandAnnotID|ECMpPK|LigandPK|MEP|ECM[[:digit:]]|Ligand[[:digit:]]|Set|Well_Ligand|ImageID|Sparse|Wedge|OuterCell|Spot_PA_Perimeter|Nuclei_PA_Cycle_State|_SE|ReplicateCount|SCC|QAScore|Lx|PinDiameter"
  
  rawSignalNames <- paste(grep("_SE",gsub("Log2|Logit","",grep("Log",colnames(slDT),value=TRUE)), value=TRUE,invert=TRUE),collapse="$|^")
  #Save the un-normalized parameters to merge in later
  mdDT <- slDT[,grep(paste0(rawSignalNames,"|",metadataNames,"|Barcode|Well|^Spot$|ArrayRow|ArrayColumn|^ECMp$|^Ligand$", collapse="|"),colnames(slDT),value=TRUE), with = FALSE]
  #Identify parameters to be normalized
  signalsWithMetadata <- grep(metadataNames,grep("Log|Barcode|Well|^Spot$|^PrintSpot$|ArrayRow|ArrayColumn|^ECMp$|^Ligand$",colnames(slDT), value=TRUE),value=TRUE,invert=TRUE)
  #Normalize each feature, pass with location and content metadata
  if(verbose) cat("Normalizing\n")
  if(!k==0){
    #Normalize each feature, pass with location and content metadata
    if(verbose)  cat("Normalizing", datasetName,"\n")

    nDT <- normRUVLoessResiduals(slDT[,signalsWithMetadata, with = FALSE], k)
    nDT$NormMethod <- "RUVLoessResiduals"
    
    #nDT has transformed and RUVLoess values
    #Merge the normalized data with its metadata
    setkey(nDT,Barcode,Well,Spot,ArrayRow,ArrayColumn,ECMp,Ligand,MEP)
    setkey(mdDT,Barcode,Well,Spot,ArrayRow,ArrayColumn,ECMp,Ligand,MEP)
    #merge in the raw data to the transformed and RUVLoess 
    slDT <- merge(nDT,mdDT)
  } else {
    slDT$NormMethod <- "none"
    slDT$k <- k
  }
  setkey(slDT, Barcode, Well, Spot, PrintSpot, ArrayRow,ArrayColumn,ECMp,Ligand)
  # slDT <- merge(slDT[,signalsWithMetadata, with = FALSE], nmdDT, by=c("Barcode", "Well", "Spot", "PrintSpot", "ArrayRow","ArrayColumn","ECMp","Ligand"))
  
  #Label FBS with their plate index to keep separate
  slDT$Ligand[grepl("FBS",slDT$Ligand)] <- paste0(slDT$Ligand[grepl("FBS",slDT$Ligand)],"_P",match(slDT$Barcode[grepl("FBS",slDT$Ligand)], unique(slDT$Barcode)))
  
  #The spot level data is median summarized to the replicate level and is stored as Level 4 data and metadata.
  
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

