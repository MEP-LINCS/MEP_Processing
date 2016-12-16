#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 9/8/2016

##Introduction
library("parallel")#use multiple cores for faster processing
source("MEP_LINCS/Release/MEPLINCSFunctions.R")

preprocessMEPLINCSL3L4 <- function(ssDataset, verbose=FALSE){
  startTime <- Sys.time()
  datasetName<-ssDataset[["datasetName"]]
  ss<-ssDataset[["ss"]]
  drug<-ssDataset[["drug"]]
  cellLine<-ssDataset[["cellLine"]]
  k<-as.integer(ssDataset[["k"]])
  analysisVersion<-ssDataset[["analysisVersion"]]
  rawDataVersion<-ssDataset[["rawDataVersion"]]
  limitBarcodes<-ssDataset[["limitBarcodes"]]
  mergeOmeroIDs<-as.logical(ssDataset[["mergeOmeroIDs"]])
  calcAdjacency<-as.logical(ssDataset[["calcAdjacency"]])
  writeFiles<-as.logical(ssDataset[["writeFiles"]])
  useAnnotMetadata<-as.logical(ssDataset[["useAnnotMetadata"]])
  
  seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")
  
  library(limma)#read GAL file and strsplit2
  library(MEMA)#merge, annotate and normalize functions
  library(data.table)#fast file reads, data merges and subsetting
  library(parallel)#use multiple cores for faster processing
  library(RUVnormalize)
  library(ruv)
  library("jsonlite")#Reading in json files
  library(stringr)
  
  #Rules-based classifier thresholds for perimeter cells
  neighborsThresh <- 0.4 #Gates sparse cells on a spot
  wedgeAngs <- 20 #Size in degrees of spot wedges used in perimeter gating
  outerThresh <- 0.5 #Defines out cells used in perimeter gating
  neighborhoodNucleiRadii <- 7 #Defines the neighborhood annulus
  
  #Filter out debris based on nuclear area
  nuclearAreaThresh <- 50
  nuclearAreaHiThresh <- 4000
  
  #Only process a curated set of the data
  curatedOnly <- TRUE
  curatedCols <- "ImageNumber|ObjectNumber|AreaShape|_MedianIntensity_|_IntegratedIntensity_|_Center_|_PA_|Texture"
  
  #Do not normalized to Spot level
  normToSpot <- TRUE

  #Set a threshold for the lowSpotCellCount flag
  lowSpotCellCountThreshold <- 5
  
  #Set a threshold for the lowRegionCellCount flag
  lowRegionCellCountThreshold <- .4
  
  #Set a threshold for the loess well level QA Scores
  lthresh <- 0.6
  
  #Set a threshold for lowWellQA flag
  lowWellQAThreshold <- .7
  
  #Set a threshold for the lowSpotReplicates flag
  lowReplicateCount <- 3
  
  datasetBarcodes <- readWorksheetFromFile("DatasetManifest.xlsx", sheet=1)
  
  fileNames <- rbindlist(apply(datasetBarcodes[datasetBarcodes$DatasetName==datasetName,], 1, getMEMADataFileNames))
  
  slDT <- fread(paste0( "MEP_LINCS/AnnotatedData/", datasetName,"_",ss,"_","SpotLevel.txt"))
  
  slDT <- slDT[!grepl("fiducial|Fiducial|gelatin|blank|air|PBS",slDT$ECMp),]
  
  metadataNames <- "ObjectNumber|^Row$|^Column$|Block|^ID$|PrintOrder|Depositions|CellLine|Endpoint|WellIndex|Center|LigandAnnotID|ECMpPK|LigandPK|MEP|ECM[[:digit:]]|Ligand[[:digit:]]|Set|Well_Ligand|ImageID|Sparse|Wedge|OuterCell|Spot_PA_Perimeter|Nuclei_PA_Cycle_State|_SE|ReplicateCount|SCC|QAScore|Lx|PinDiameter"
  
  rawSignalNames <- paste(grep("_SE",gsub("Log2|Logit","",grep("Log",colnames(slDT),value=TRUE)), value=TRUE,invert=TRUE),collapse="$|^")
  #Save the un-normalized parameters to merge in later
  mdDT <- slDT[,grep(paste0(rawSignalNames,"|",metadataNames,"|Barcode|Well|^Spot$|ArrayRow|ArrayColumn|^ECMp$|^Ligand$", collapse="|"),colnames(slDT),value=TRUE), with = FALSE]
  #Identify parameters to be normalized
  signalsWithMetadata <- grep(metadataNames,grep("Log|Barcode|Well|^Spot$|^PrintSpot$|ArrayRow|ArrayColumn|^ECMp$|^Ligand$",colnames(slDT), value=TRUE),value=TRUE,invert=TRUE)
  #Normalize each feature, pass with location and content metadata
  if(verbose) {
    cat("Normalizing\n")
    #save(slDT, file="slDT.RData")
  }
  
  if(!k==0){
    #Normalize each feature, pass with location and content metadata
    if(verbose) {
      cat("Normalizing\n")
      #save(slDT, file="slDT.RData")
    }
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
  
  #Level4Data
  if(verbose) {
    cat("Creating level 4 data\n")
    #save(slDT,file="slDT.RData")
  }
  mepDT <- createl4(slDT,seNames = seNames)
  
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
  
  #Level 4
  mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < lowReplicateCount
  #WriteData
  if(writeFiles){
    if(verbose) cat("Writing level 3 file to disk\n")
    
    fwrite(data.table(format(slDT, digits = 4, trim=TRUE)), paste0( "MEP_LINCS/AnnotatedData/", datasetName,"_",ss,"_","Level3.txt"), sep = "\t", quote=FALSE)
    
    if(verbose) cat("Writing level 4 file to disk\n")
    fwrite(data.table(format(mepDT, digits = 4, trim=TRUE)), paste0( "MEP_LINCS/AnnotatedData/", datasetName,"_",ss,"_","Level4.txt"), sep = "\t", quote=FALSE)
    
  }
  cat("Elapsed time:", Sys.time()-startTime)
}

PC3df <- data.frame(datasetName=c("PC3_SS1","PC3_SS2","PC3_SS3","PC3_SS2noH3"),
                    cellLine=rep(c("PC3"), 4),
                    ss=c("SS1", "SS2","SS3","SS2noH3"),
                    drug=c("none"),
                    analysisVersion="av1.6",
                    rawDataVersion=c("v2","v2.1","v2.1", "v1"),
                    limitBarcodes=8,
                    k=7,
                    calcAdjacency=TRUE,
                    writeFiles = TRUE,
                    mergeOmeroIDs = TRUE,
                    useAnnotMetadata=TRUE,
                    stringsAsFactors=FALSE)

MCF7df <- data.frame(datasetName=c("MCF7_SS1","MCF7_SS2","MCF7_SS3"),
                     cellLine=rep(c("MCF7"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug=c("none"),
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     limitBarcodes=8,
                     k=7,
                     calcAdjacency=TRUE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     useAnnotMetadata=TRUE,
                     stringsAsFactors=FALSE)

YAPCdf <- data.frame(datasetName=c("YAPC_SS1","YAPC_SS2","YAPC_SS3"),
                     cellLine=rep(c("YAPC"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug=c("none"),
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     limitBarcodes=8,
                     k=7,
                     calcAdjacency=TRUE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     useAnnotMetadata=TRUE,
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(datasetName=c("MCF10A_SS1","MCF10A_SS2","MCF10A_SS3"),
                       cellLine="MCF10A",
                       ss=c("SS1","SS2","SS3"),
                       drug=c("none"),
                       analysisVersion="av1.7",
                       rawDataVersion="v2",
                       limitBarcodes=c(8,8,8),
                       k=c(135,135,135),
                       calcAdjacency=TRUE,
                       writeFiles = TRUE,
                       mergeOmeroIDs = TRUE,
                       useAnnotMetadata=TRUE,
                       stringsAsFactors=FALSE)

watsonMEMAs <- data.frame(datasetName=c("HCC1954_DMSO","HCC1954_Lapatinib","AU565_DMSO","AU565_Lapatinib"),
                          cellLine=c("HCC1954","HCC1954","AU565","AU565"),
                          ss=c("SS6"),
                          drug=c("DMSO","Lapatinib"),
                          analysisVersion="av1.6",
                          rawDataVersion="v2",
                          limitBarcodes=c(8,8,8,8),
                          k=c(7,7,7,7),
                          calcAdjacency=TRUE,
                          writeFiles = TRUE,
                          mergeOmeroIDs = TRUE,
                          useAnnotMetadata=FALSE,
                          stringsAsFactors=FALSE)

qualPlates <- data.frame(datasetName="MCF10A_Qual",
                         cellLine=c("MCF10A"),
                         ss=c("SS0"),
                         drug=c("none"),
                         analysisVersion="av1.6",
                         rawDataVersion="v2",
                         limitBarcodes=c(4),
                         k=c(0),
                         calcAdjacency=FALSE,
                         writeFiles = TRUE,
                         mergeOmeroIDs = TRUE,
                         useAnnotMetadata=FALSE,
                         stringsAsFactors=FALSE)

ctrlPlates <- data.frame(datasetName="HMEC_Ctrl",
                         cellLine=c("HMEC122L"),
                         ss=c("SS0"),
                         drug=c("none"),
                         analysisVersion="av1.6",
                         rawDataVersion="v2",
                         limitBarcodes=c(1),
                         k=c(0),
                         calcAdjacency=FALSE,
                         writeFiles = TRUE,
                         mergeOmeroIDs = TRUE,
                         useAnnotMetadata=FALSE,
                         stringsAsFactors=FALSE)

HMEC240L <- data.frame(datasetName=c("HMEC240L_SS1","HMEC240L_SS4"),
                       cellLine=c("HMEC240L"),
                       ss=c("SS1","SS4"),
                       drug=c("none"),
                       analysisVersion="av1.7",
                       rawDataVersion="v2",
                       limitBarcodes=c(8,8),
                       k=c(64,64),
                       calcAdjacency=TRUE,
                       writeFiles = TRUE,
                       mergeOmeroIDs = TRUE,
                       useAnnotMetadata=TRUE,
                       stringsAsFactors=FALSE)

HMEC122L <- data.frame(datasetName=c("HMEC122L_SS1","HMEC122L_SS4"),
                       cellLine=c("HMEC122L"),
                       ss=c("SS1","SS4"),
                       drug=c("none"),
                       analysisVersion="av1.7",
                       rawDataVersion="v2",
                       limitBarcodes=c(8,8),
                       k=c(64,64),
                       calcAdjacency=TRUE,
                       writeFiles = TRUE,
                       mergeOmeroIDs = TRUE,
                       useAnnotMetadata=TRUE,
                       stringsAsFactors=FALSE)
ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf,watsonMEMAs,qualPlates, ctrlPlates, HMEC240L, HMEC122L)

tcDataSet <- data.frame(datasetName=c("MCF10A_TC"),
                        cellLine=c("MCF10A"),
                        ss=c("SS4"),
                        drug=c("none"),
                        analysisVersion="av1.7",
                        rawDataVersion="v2",
                        limitBarcodes=c(3),
                        k=c(64),
                        calcAdjacency=TRUE,
                        writeFiles = TRUE,
                        mergeOmeroIDs = TRUE,
                        useAnnotMetadata=FALSE,
                        stringsAsFactors=FALSE)

Bornstein <- data.frame(datasetName=c("BornsteinOSC","BornsteinCal27"),
                        cellLine=c("OSC","Cal27"),
                        ss=c("SSA"),
                        drug=c("radiation"),
                        analysisVersion="av1.7",
                        rawDataVersion="v2",
                        limitBarcodes=c(2),
                        k=c(64),
                        calcAdjacency=TRUE,
                        writeFiles = TRUE,
                        mergeOmeroIDs = TRUE,
                        useAnnotMetadata=FALSE,
                        stringsAsFactors=FALSE)

Vertex <- data.frame(datasetName=c("Vertex1", "Vertex2"),
                     cellLine=c("LCSC-311"),
                     ss=c("SSE"),
                     drug=c("none"),
                     analysisVersion="av1.7",
                     rawDataVersion="v2",
                     limitBarcodes=c(8),
                     k=c(64),
                     calcAdjacency=FALSE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     useAnnotMetadata=FALSE,
                     stringsAsFactors=FALSE)

Baylor <- data.frame(datasetName=c("Baylor1", "Baylor2","Baylor3", "Baylor4","Baylor5", "Baylor6"),
                     cellLine=c("LM2", "LM2","LM2","SUM159","SUM159","SUM159"),
                     ss=c("SSD"),
                     drug=c("unknown"),
                     analysisVersion="av1.7",
                     rawDataVersion="v2",
                     limitBarcodes=c(2,2,1,2,2,1),
                     k=c(64),
                     calcAdjacency=FALSE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = FALSE,
                     useAnnotMetadata=FALSE,
                     stringsAsFactors=FALSE)

MLDDataSet <- data.frame(datasetName=c("MCF10ANeratinib","MCF10ADMSO","MCF10AVorinostat","MCF10ATrametinib"),
                         cellLine=c("MCF10A"),
                         ss=c("SSF"),
                         drug=c("Neratinib","DMSO","Vorinostat","Ttametinib"),
                         analysisVersion="av1.7",
                         rawDataVersion="v2",
                         limitBarcodes=c(8),
                         k=c(64),
                         calcAdjacency=TRUE,
                         writeFiles = TRUE,
                         mergeOmeroIDs = TRUE,
                         useAnnotMetadata=TRUE,
                         stringsAsFactors=FALSE)

validations <- data.frame(datasetName=c("MCF10AHighRep1","MCF10AHighRep3"),
                          cellLine=c("MCF10A"),
                          ss=c("SS4"),
                          drug=c("none"),
                          analysisVersion="av1.7",
                          rawDataVersion="v2",
                          limitBarcodes=c(4),
                          k=c(0),
                          calcAdjacency=TRUE,
                          writeFiles = TRUE,
                          mergeOmeroIDs = TRUE,
                          useAnnotMetadata=FALSE,
                          stringsAsFactors=FALSE)
library(XLConnect)
library(data.table)

tmp <- apply(Bornstein, 1, preprocessMEPLINCSL3L4, verbose=FALSE)

