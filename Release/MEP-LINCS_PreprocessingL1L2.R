#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 11/18/2015

##Introduction

#   The MEP-LINCs dataset contains imaging data from a Nikon automated microscope that is analyzed with a CellProfiler pipeline.
# 
# Part of this preprocessing of the dataset will be deprecated when the merging of the data and metadata happens within the CellProfiler part of the pipeline. For now, the metadata about the ECM proteins is read from the GAL file and the metadata about the wells (cell line, stains and ligands) is read from Excel spreadsheets.

library("parallel")#use multiple cores for faster processing
source("MEP_LINCS/Release/MEPLINCSFunctions.R")
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.3"))


processan2omero <- function (fileNames) {
  rbindlist(mclapply(fileNames, function(fn){
    #Process each file separately
    dt <- fread(fn)
    
    #Rename to preprocessing pipeline variable names
    setnames(dt,"OSpot","Spot")
    setnames(dt,"PlateID","Barcode")
    setnames(dt,"395nm","EndpointDAPI")
    setnames(dt,"488nm","Endpoint488")
    setnames(dt,"555nm","Endpoint555")
    setnames(dt,"640nm","Endpoint647")
    setnames(dt,"750nm","Endpoint750")
    #Shorten Annot names
    dt$CellLine <- gsub("_.*","",dt$CellLine)
    dt$ECMp <- dt$ECM1
    dt <- shortenHA(dt)
    dt$ECMp <-gsub("_.*","",dt$ECMp)
    dt$Ligand <-gsub("_.*","",dt$Ligand1)
    dt$EndpointDAPI <-gsub("_.*","",dt$EndpointDAPI)
    dt$Endpoint488 <-gsub("_.*","",dt$Endpoint488)
    dt$Endpoint555 <-gsub("_.*","",dt$Endpoint555)
    dt$Endpoint647 <-gsub("_.*","",dt$Endpoint647)
    #Add a WellSpace spot index that recognizes the arrays are rotated 180 degrees
    dt$PrintSpot <- dt$Spot
    nrArrayRows <- max(dt$ArrayRow)
    nrArrayColumns <- max(dt$ArrayColumn)
    dt$PrintSpot[grepl("B", dt$Well)] <- (nrArrayRows*nrArrayColumns+1)-dt$PrintSpot[grepl("B", dt$Well)]
    return(dt)
  }, mc.cores=max(4, detectCores())))
}

spotNorm <- function(x){
  return(x/median(x,na.rm=TRUE))
}

gateOnQuantile <- function(x,probs){
  gatedClass <- integer(length(x))
  gatedClass[x>quantile(x,probs=probs,na.rm=TRUE)]<-1
  return(gatedClass)
}

preprocessMEPLINCSL1Spot <- function(ssDataset, verbose=FALSE){
  startTime<- Sys.time()
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
  
  #QA flags are used to enable analyses that require minimum cell and
  #replicate counts
  
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
  
  #Select the correct path based on the system that's executing the code
  if(grepl("omero|eppec",Sys.info()[["nodename"]])){
    fileNames$Path <- fileNames$Path
  } else if(Sys.info()[["nodename"]]=="CLSBA715"){
    fileNames$Path <- fileNames$LocalPath
  } else {
    stop("Executing on unknown system and need path to input files")
  }
  
  ##Summary
  # This script prepares cell-level data and metadata for the MEP LINCs Analysis Pipeline. 
  # 
  # In the code, the variable ss determines which staining set (SS1, SS2, SS3 and others) to merge and the variable cellLine determines the cell line (PC3,MCF7, etc).
  # 
  # The merging assumes that the actual, physical B row wells (B01-B04) have been printed upside-down. That is, rotated 180 degrees resulting in the spot 1, 1 being in the lower right corner instead of the upper left corner. The metadata is matched to the actual printed orientation.
  
  #Use metadata from an2omero files or 
  #directly from GAL and xlsx files
  if(useAnnotMetadata){
    #readJSONMetadata
    fns <- fileNames$Path[fileNames$Type=="metadata"&grepl("an2omero",fileNames$Path)]
    metadata <- rbindlist(mclapply(fns[1:limitBarcodes], processan2omero, mc.cores = detectCores()))
  } else {
    # Read and clean spotmetadata
    
    #Read in the spot metadata from the gal file
    if(!length(fileNames$Path[fileNames$Type=="gal"])==1) stop("There must be 1 gal file in the dataset folders")
    smd <- readSpotMetadata(fileNames$Path[fileNames$Type=="gal"])
    
    #Relabel the column Name to ECMp
    setnames(smd, "Name", "ECMp")
    
    #Add the print order and deposition number to the metadata
    if(!length(fileNames$Path[fileNames$Type=="xml"])==1) stop("There must be 1 xml log file in the dataset folders")
    ldf <- readLogData(fileNames$Path[fileNames$Type=="xml"])
    
    spotMetadata <- merge(smd,ldf, all=TRUE)
    setkey(spotMetadata,Spot)
    #Make a rotated version of the spot metadata to match the print orientation
    spotMetadata180 <- rotateMetadata(spotMetadata)
  }
  
  # The raw data from all wells in all plates in the dataset are read in and merged with their spot and well metadata. The number of nuclei at each spot are counted and a loess model of the spot cell count is added. Then all intensity values are normalized.
  # 
  # Next, the data is filtered to remove objects with a nuclear area less than nuclearAreaThresh pixels or more than nuclearAreaHiThresh pixels.
  
  #merge_normalize_QA, echo=FALSE}
  #The next steps are to bring in the well metadata, the print order and the CP data
  if(length(fileNames$Path[fileNames$Type=="Raw"])==0) stop("No raw data files found")
  cellDataFiles <- fileNames$Path[fileNames$Type=="Raw"]
  
  barcodes <- gsub("rescan|reDAPI","",unique(fileNames$Barcode)[1:limitBarcodes])
  if(verbose) cat(paste("Reading and annotating cell level data for",cellLine,ss)," \n")
  expDTList <- lapply(barcodes, function(barcode){
    plateDataFiles <- fileNames$Path[grepl(barcode,fileNames$Barcode)&
                                       grepl("Raw",fileNames$Type)]
    splits <- strsplit2(split = "_",plateDataFiles)
    wells <- unique(splits[,dim(splits)[2]-1])
    wellDataList <- lapply(wells,function(well){
      if(verbose) cat(paste("Reading and annotating well data for",well,"of plate",barcode)," \n")
      
      wellDataFiles <- grep(well,plateDataFiles,value = TRUE)
      nucleiDataFile <- grep("Nuclei",wellDataFiles,value=TRUE,
                             ignore.case = TRUE)
      if (ss %in% c("SS1","SS3","SS4", "SS6", "SSD","SSF")){
        cellsDataFile <- grep("Cell",wellDataFiles,value=TRUE,
                              ignore.case = TRUE)
        cytoplasmDataFile <- grep("Cytoplasm",wellDataFiles,value=TRUE,
                                  ignore.case = TRUE)
      }
      #Read in csv data
      nuclei <- convertColumnNames(fread(nucleiDataFile))
      if (curatedOnly) nuclei <- nuclei[,grep(curatedCols,colnames(nuclei)), with=FALSE]
      #Remove selected columns
      nuclei <- nuclei[,grep("Euler",colnames(nuclei),invert=TRUE), with=FALSE]
      nuclei <- nuclei[,grep("AreaShape|Intensity|Number",colnames(nuclei)), with=FALSE]
      setkey(nuclei,CP_ImageNumber,CP_ObjectNumber)
      if (ss %in% c("SS1","SS3","SS4", "SS6", "SSD","SSF")){
        cells <- convertColumnNames(fread(cellsDataFile))
        if (curatedOnly) cells <- cells[,grep(curatedCols,colnames(cells)), with=FALSE]
        #Remove selected columns
        cells <- cells[,grep("Euler",colnames(cells),invert=TRUE), with=FALSE]
        cells <- cells[,grep("AreaShape|Intensity|Number",colnames(cells)), with=FALSE]
        setkey(cells,CP_ImageNumber,CP_ObjectNumber)
        cytoplasm <- convertColumnNames(fread(cytoplasmDataFile))
        if (curatedOnly) cytoplasm <- cytoplasm[,grep(curatedCols,colnames(cytoplasm)), with=FALSE]
        #Remove selected columns
        cytoplasm <- cytoplasm[,grep("Euler",colnames(cytoplasm),invert=TRUE), with=FALSE]        
        cytoplasm <- cytoplasm[,grep("AreaShape|Intensity|Number",colnames(cytoplasm)), with=FALSE]
        setkey(cytoplasm,CP_ImageNumber,CP_ObjectNumber)
      }
      
      #Add the data location as a prefix in the column names
      setnames(nuclei,paste0("Nuclei_",colnames(nuclei)))
      if (ss %in% c("SS1","SS3","SS4","SS6", "SSD","SSF")){
        setnames(cells,paste0("Cells_",colnames(cells)))
        setnames(cytoplasm,paste0("Cytoplasm_",colnames(cytoplasm)))
      }
      
      #Merge the cells, cytoplasm and nuclei data
      if (ss %in% c("SS1","SS3","SS4","SS6","SSD","SSF")){
        setkey(cells,Cells_CP_ImageNumber,Cells_CP_ObjectNumber)
        setkey(cytoplasm,Cytoplasm_CP_ImageNumber,Cytoplasm_CP_ObjectNumber)
        setkey(nuclei,Nuclei_CP_ImageNumber,Nuclei_CP_ObjectNumber)
        
        DT <- cells[cytoplasm[nuclei]]
        setnames(DT,"Cells_CP_ImageNumber","Spot")
        setnames(DT,"Cells_CP_ObjectNumber","ObjectNumber")
      } else {
        DT <- nuclei
        setnames(DT,"Nuclei_CP_ImageNumber","Spot")
        setnames(DT,"Nuclei_CP_ObjectNumber","ObjectNumber")
      }
      
      #Add the well name and barcode as parameters
      DT <- DT[,Well := well]
      DT <- DT[,Barcode := barcode]
      
      if (useAnnotMetadata) {
        DT <- merge(DT,metadata,by = c("Barcode","Well","Spot"))
      } else {
        #Merge the data with its metadata based on the row it's in
        m <- regexpr("[[:alpha:]]",well)
        row <- regmatches(well,m)
        setkey(DT,Spot)
        DT <- switch(row, A = merge(DT,spotMetadata,all.x=TRUE),
                     B = merge(DT,spotMetadata180,all.x=TRUE))
      }
      return(DT)
    })#Revert to apply when debugging
    #})
    #Create the cell data.table with spot metadata for the plate 
    pcDT <- rbindlist(wellDataList, fill = TRUE)
    rm(wellDataList)
    gc()
    
    if(!useAnnotMetadata){
      #Read the well metadata from a multi-sheet Excel file
      #wellMetadata <- data.table(readMetadata(paste0(filePath,"/Metadata/",gsub("reDAPI","",barcode),".xlsx")), key="Well")
      #Debug reDAPI plates with new file structure
      wellMetadata <- data.table(readMetadata(fileNames$Path[fileNames$Barcode==barcode&fileNames$Type=="metadata"]), key="Well")
      #merge well metadata with the data and spot metadata
      pcDT <- merge(pcDT,wellMetadata,by = "Well")
      #Add a WellSpace spot index that recognizes the arrays are rotated 180 degrees
      nrArrayRows <- max(pcDT$ArrayRow)
      nrArrayColumns <- max(pcDT$ArrayColumn)
      pcDT$PrintSpot <- pcDT$Spot
      pcDT$PrintSpot[grepl("B", pcDT$Well)] <- (nrArrayRows*nrArrayColumns+1)-pcDT$PrintSpot[grepl("B", pcDT$Well)]
    }
    
    if(mergeOmeroIDs){
      imageURLFile <- fileNames$Path[grepl(barcode,fileNames$Barcode)&grepl("imageID",fileNames$Type)]
      #Read in and merge the Omero URLs
      omeroIndex <- fread(fileNames$Path[grepl(barcode,fileNames$Barcode)&grepl("imageID",fileNames$Type)])[,list(WellName,Row,Column,ImageID)]
      m <- regexpr("Well[[:digit:]]",omeroIndex$WellName)
      wellNames <- regmatches(omeroIndex$WellName,m)
      omeroIndex$Well <- sapply(gsub("Well","",wellNames,""),FUN=switch,
                                "1"="A01",
                                "2"="A02",
                                "3"="A03",
                                "4"="A04",
                                "5"="B01",
                                "6"="B02",
                                "7"="B03",
                                "8"="B04")
      setnames(omeroIndex,"Row","ArrayRow")
      setnames(omeroIndex,"Column","ArrayColumn")
      omeroIndex <- omeroIndex[,WellName:=NULL]
      pcDT <- merge(pcDT,omeroIndex,by=c("Well","ArrayRow","ArrayColumn"))
    }
    
    if(any(grepl("Nuclei_CP_AreaShape_Area",colnames(pcDT)))){
      pcDT <- pcDT[pcDT$Nuclei_CP_AreaShape_Area > nuclearAreaThresh,]
      pcDT <- pcDT[pcDT$Nuclei_CP_AreaShape_Area < nuclearAreaHiThresh,]
    }
    
    #Change Edu back to EdU
    if(any(grepl("Edu",colnames(pcDT)))){
      edUNames <- grep("Edu",colnames(pcDT),value=TRUE)
      setnames(pcDT,edUNames,gsub("Edu","EdU",edUNames))
    }
    
    #log transform all intensity and areaShape values
    intensityNames <- grep("Intensity",colnames(pcDT), value=TRUE)
    scaledInts <- pcDT[,intensityNames, with=FALSE]*2^16
    pcDT <- cbind(pcDT[,!intensityNames, with=FALSE],scaledInts)
    transformNames <- grep("_Center_|_Eccentricity|_Orientation",grep("Intensity|AreaShape",colnames(pcDT), value=TRUE, ignore.case = TRUE), value=TRUE, invert=TRUE)
    dtLog <- pcDT[,lapply(.SD,boundedLog2),.SDcols=transformNames]
    setnames(dtLog,colnames(dtLog),paste0(colnames(dtLog),"Log2"))
    pcDT <- cbind(pcDT,dtLog)
    
    #Count the cells at each spot
    pcDT <- pcDT[,Spot_PA_SpotCellCount := .N,by="Barcode,Well,Spot"]
    pcDT <- pcDT[,Spot_PA_SpotCellCountLog2 := boundedLog2(Spot_PA_SpotCellCount)]
    
    if(any(grepl("Nuclei_CP_AreaShape_Center",colnames(pcDT)))){
      #Add the local polar coordinates and Neighbor Count
      pcDT <- pcDT[,Nuclei_PA_Centered_X :=  Nuclei_CP_AreaShape_Center_X-median(Nuclei_CP_AreaShape_Center_X)]
      pcDT <- pcDT[,Nuclei_PA_Centered_Y :=  Nuclei_CP_AreaShape_Center_Y-median(Nuclei_CP_AreaShape_Center_Y)]
      pcDT <- pcDT[, Nuclei_PA_AreaShape_Center_R := sqrt(Nuclei_PA_Centered_X^2 + Nuclei_PA_Centered_Y^2)]
      pcDT <- pcDT[, Nuclei_PA_AreaShape_Center_Theta := calcTheta(Nuclei_PA_Centered_X, Nuclei_PA_Centered_Y)]
    }
    #Set 2N and 4N DNA status
    pcDT <- pcDT[,Nuclei_PA_Cycle_State := gateOnlocalMinima(Nuclei_CP_Intensity_IntegratedIntensity_Dapi), by="Barcode,Well"]
    
    pcDT <- pcDT[,Nuclei_PA_Cycle_DNA2NProportion := calc2NProportion(Nuclei_PA_Cycle_State),by="Barcode,Well,Spot"]
    pcDT <- pcDT[,Nuclei_PA_Cycle_DNA4NProportion := 1-Nuclei_PA_Cycle_DNA2NProportion]
    
    #Logit transform DNA Proportions
    #logit(p) = log[p/(1-p)]
    if(any(grepl("Nuclei_PA_Cycle_DNA2NProportion",colnames(pcDT)))){
      pcDT <- pcDT[,Nuclei_PA_Cycle_DNA2NProportionLogit := boundedLogit(Nuclei_PA_Cycle_DNA2NProportion)]
    }
    
    if(any(grepl("Nuclei_PA_Cycle_DNA4NProportion",colnames(pcDT)))){
      pcDT <- pcDT[,Nuclei_PA_Cycle_DNA4NProportionLogit := boundedLogit(Nuclei_PA_Cycle_DNA4NProportion)]
    }
    
    #logit transform eccentricity
    if(any(grepl("Nuclei_CP_AreaShape_Eccentricity",colnames(pcDT)))){
      pcDT <- pcDT[,Nuclei_CP_AreaShape_EccentricityLogit := boundedLogit(Nuclei_CP_AreaShape_Eccentricity)]
    }
    if(normToSpot){
      #Add spot level normalizations for selected intensities
      intensityNamesAll <- grep("_CP_Intensity_Median",colnames(pcDT), value=TRUE)
      intensityNames <- grep("Norm",intensityNamesAll,invert=TRUE,value=TRUE)
      for(intensityName in intensityNames){
        #Median normalize the median intensity at each spot
        setnames(pcDT,intensityName,"value")
        pcDT <- pcDT[,paste0(intensityName,"_SpotNorm") := spotNorm(value),by="Barcode,Well,Spot"]
        setnames(pcDT,"value",intensityName)
      }
    }
    if(!useAnnotMetadata){
      #Create short display names, then replace where not unique
      #Use entire AnnotID for ligands with same uniprot IDs
      # pcDT$Ligand[grepl("NRG1",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "NRG1",annotIDs = pcDT$Ligand[grepl("NRG1",pcDT$Ligand)])
      # pcDT$Ligand[grepl("TGFB1",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "TGFB1",annotIDs = pcDT$Ligand[grepl("TGFB1",pcDT$Ligand)])
      # pcDT$Ligand[grepl("CXCL12",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "CXCL12",annotIDs = pcDT$Ligand[grepl("CXCL12",pcDT$Ligand)])
      # pcDT$Ligand[grepl("IGF1",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "IGF1",annotIDs = pcDT$Ligand[grepl("IGF1",pcDT$Ligand)])
    }
    #Add MEP and convenience labels for wells and ligands
    pcDT <- pcDT[,MEP:=paste(ECMp,Ligand,sep = "_")]
    pcDT <- pcDT[,Well_Ligand:=paste(Well,Ligand,sep = "_")]
    
    # Eliminate Variations in the Endpoint metadata
    endpointNames <- grep("End",colnames(pcDT), value=TRUE)
    endpointWL <- regmatches(endpointNames,regexpr("[[:digit:]]{3}|DAPI",endpointNames))
    setnames(pcDT,endpointNames,paste0("Endpoint",endpointWL))
    
    
    #Create staining set specific derived parameters
    if (grepl("SS2|SS4|SS6|SSA|SSD|SSE|SSF",ss)){
      #Use the entire plate to set the autogate threshold if there's no control well
      if(grepl("FBS",unique(pcDT$Ligand))){
        pcDT <- pcDT[,Nuclei_PA_Gated_EdUPositive := kmeansCluster(.SD,value =  "Nuclei_CP_Intensity_MedianIntensity_EdU",ctrlLigand = "FBS"), by="Barcode"]
      } else {
        pcDT <- pcDT[,Nuclei_PA_Gated_EdUPositive := kmeansCluster(.SD,value =  "Nuclei_CP_Intensity_MedianIntensity_EdU",ctrlLigand = "."), by="Barcode"]
      }
      
      #Modify the gate threshold to be 3 sigma from the EdU- median
      EdUDT <- pcDT[Nuclei_PA_Gated_EdUPositive==0,.(EdUMedian=median(Nuclei_CP_Intensity_MedianIntensity_EdULog2, na.rm=TRUE),EdUSD=sd(Nuclei_CP_Intensity_MedianIntensity_EdULog2, na.rm=TRUE)),by="Barcode,Well"]
      EdUDT <- EdUDT[,Median3SD := EdUMedian+3*EdUSD]
      pcDT <- merge(pcDT,EdUDT,by=c("Barcode","Well"))
      #Overide autogate for AU565 cell line in Lapatinib MEMAs
      pcDT$Nuclei_PA_Gated_EdUPositive[pcDT$Nuclei_CP_Intensity_MedianIntensity_EdULog2>pcDT$Median3SD]<-1
      if(grepl("AU565_Lapatinib|AU565_DMSO",datasetName)){
        pcDT$Nuclei_PA_Gated_EdUPositive <- 0
        pcDT$Nuclei_PA_Gated_EdUPositive[pcDT$Nuclei_CP_Intensity_MedianIntensity_EdULog2>10] <- 1 
      }
      #Calculate the EdU Positive Percent at each spot
      pcDT <- pcDT[,Nuclei_PA_Gated_EdUPositiveProportion := sum(Nuclei_PA_Gated_EdUPositive)/length(ObjectNumber),by="Barcode,Well,Spot"]
      #Logit transform EdUPositiveProportion
      #logit(p) = log[p/(1-p)]
      pcDT <-pcDT[,Nuclei_PA_Gated_EdUPositiveProportionLogit:= boundedLogit(Nuclei_PA_Gated_EdUPositiveProportion)]
    } 
    if (grepl("SS3|SS4",ss)){
      #Calculate a lineage ratio of luminal/basal or KRT19/KRT5
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT5]
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatioLog2 := boundedLog2(Cytoplasm_PA_Intensity_LineageRatio)]
      #Gate the KRT signals using kmeans clustering and calculate the spot proportions
      
      if(grepl("HMEC",cellLine)) {
        pcDT <- pcDT[,Cytoplasm_PA_Gated_KRT5Positive := gateOnQuantile(Cytoplasm_CP_Intensity_MedianIntensity_KRT5Log2,probs=.02),by="Barcode,Well"]
      } else {
        pcDT <- pcDT[,Cytoplasm_PA_Gated_KRT5Positive := kmeansCluster(.SD,value =  "Cytoplasm_CP_Intensity_MedianIntensity_KRT5",ctrlLigand = "."), by="Barcode"]
      }
      pcDT <- pcDT[,Cytoplasm_PA_Gated_KRT19Positive := kmeansCluster(.SD,value =  "Cytoplasm_CP_Intensity_MedianIntensity_KRT19",ctrlLigand = "."), by="Barcode"]
      #Determine the class of each cell based on KRT5 and KRT19 class
      #0 double negative
      #1 KRT5+, KRT19-
      #2 KRT5-, KRT19+
      #3 KRT5+, KRT19+
      pcDT <- pcDT[,Cytoplasm_PA_Gated_KRTClass := Cytoplasm_PA_Gated_KRT5Positive+2*Cytoplasm_PA_Gated_KRT19Positive]
      #Calculate the positive proportion at each spot
      pcDT <- pcDT[,Cytoplasm_PA_Gated_KRT5PositiveProportion := sum(Cytoplasm_PA_Gated_KRT5Positive)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cytoplasm_PA_Gated_KRT5PositiveProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5PositiveProportion)]
      pcDT <- pcDT[,Cytoplasm_PA_Gated_KRT19PositiveProportion := sum(Cytoplasm_PA_Gated_KRT19Positive)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cytoplasm_PA_Gated_KRT19PositiveProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT19PositiveProportion)]
      pcDT <- pcDT[,Cytoplasm_PA_Gated_KRT5KRT19NegativeProportion := sum(Cytoplasm_PA_Gated_KRTClass==0)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cytoplasm_PA_Gated_KRT5KRT19NegativeProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5KRT19NegativeProportion)]
      pcDT <- pcDT[,Cytoplasm_PA_Gated_KRT5PositiveKRT19NegativeProportion := sum(Cytoplasm_PA_Gated_KRTClass==1)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cytoplasm_PA_Gated_KRT5PositiveKRT19NegativeProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5PositiveKRT19NegativeProportion)]
      pcDT <- pcDT[,Cytoplasm_PA_Gated_KRT5NegativeKRT19PositiveProportion := sum(Cytoplasm_PA_Gated_KRTClass==2)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cytoplasm_PA_Gated_KRT5NegativeKRT19PositiveProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5NegativeKRT19PositiveProportion)]
      pcDT <- pcDT[,Cytoplasm_PA_Gated_KRT5KRT19PositiveProportion := sum(Cytoplasm_PA_Gated_KRTClass==3)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cytoplasm_PA_Gated_KRT5KRT19PositiveProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5KRT19PositiveProportion)]
    }
    if (grepl("SS4",ss)){
      #Determine the class of each cell based on KRT19 and EdU class
      #0 double negative
      #1 KRT19+, EdU-
      #2 KRT19-, EdU+
      #3 KRT19+, EdU+
      pcDT <- pcDT[,Cells_PA_Gated_EdUKRT19Class := Cytoplasm_PA_Gated_KRT19Positive+2*Nuclei_PA_Gated_EdUPositive]
      #Calculate gating proportions for EdU and KRT19
      pcDT <- pcDT[,Cells_PA_Gated_EdUKRT19NegativeProportion := sum(Cells_PA_Gated_EdUKRT19Class==0)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cells_PA_Gated_EdUKRT19NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_EdUKRT19NegativeProportion)]
      pcDT <- pcDT[,Cells_PA_Gated_EdUPositiveKRT19NegativeProportion := sum(Cells_PA_Gated_EdUKRT19Class==2)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cells_PA_Gated_EdUPositiveKRT19NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_EdUPositiveKRT19NegativeProportion)]
      pcDT <- pcDT[,Cells_PA_Gated_EdUNegativeKRT19PositiveProportion := sum(Cells_PA_Gated_EdUKRT19Class==1)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cells_PA_Gated_EdUNegativeKRT19PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_EdUNegativeKRT19PositiveProportion)]
      pcDT <- pcDT[,Cells_PA_Gated_EdUKRT19PositiveProportion := sum(Cells_PA_Gated_EdUKRT19Class==3)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cells_PA_Gated_EdUKRT19PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_EdUKRT19PositiveProportion)]
    }
    
    if (grepl("SS6",ss)){
      #Calculate a lineage ratio of luminal/basal or KRT19/KRT14
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT14]
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatioLog2 := boundedLog2(Cytoplasm_PA_Intensity_LineageRatio)]
    }
    if (grepl("SSF",ss)){
      #Determine the class of each cell based on KRT5 and EdU class
      #0 double negative
      #1 KRT5+, EdU-
      #2 KRT5-, EdU+
      #3 KRT5+, EdU+
      pcDT <- pcDT[,Cytoplasm_PA_Gated_KRT5Positive := gateOnQuantile(Cytoplasm_CP_Intensity_MedianIntensity_KRT5Log2,probs=.02),by="Barcode,Well"]
      pcDT <- pcDT[,Cells_PA_Gated_EdUKRT5Class := Cytoplasm_PA_Gated_KRT5Positive+2*Nuclei_PA_Gated_EdUPositive]
      #Calculate gating proportions for EdU and KRT19
      pcDT <- pcDT[,Cells_PA_Gated_EdUKRT5NegativeProportion := sum(Cells_PA_Gated_EdUKRT5Class==0)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cells_PA_Gated_EdUKRT5NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_EdUKRT5NegativeProportion)]
      pcDT <- pcDT[,Cells_PA_Gated_EdUPositiveKRT5NegativeProportion := sum(Cells_PA_Gated_EdUKRT5Class==2)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cells_PA_Gated_EdUPositiveKRT5NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_EdUPositiveKRT5NegativeProportion)]
      pcDT <- pcDT[,Cells_PA_Gated_EdUNegativeKRT5PositiveProportion := sum(Cells_PA_Gated_EdUKRT5Class==1)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cells_PA_Gated_EdUNegativeKRT5PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_EdUNegativeKRT5PositiveProportion)]
      pcDT <- pcDT[,Cells_PA_Gated_EdUKRT5PositiveProportion := sum(Cells_PA_Gated_EdUKRT5Class==3)/length(ObjectNumber),by="Barcode,Well,Spot"]
      pcDT <-pcDT[,Cells_PA_Gated_EdUKRT19PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_EdUKRT5PositiveProportion)]
    }
    #Move adjacency processing to within parallel code
    if(calcAdjacency){
      if(verbose){
        cat("Calculating adjacency data\n")
        #save(cDT, file="cDT.RData")
        #load(file = "~/Documents/MEP-LINCS/cDT.RData")
      } 
      
      densityRadius <- sqrt(median(pcDT$Nuclei_CP_AreaShape_Area, na.rm = TRUE)/pi)
      
      #Count the number of neighboring cells
      pcDT <- pcDT[,Nuclei_PA_AreaShape_Neighbors := cellNeighbors(.SD, radius = densityRadius*neighborhoodNucleiRadii), by = "Barcode,Well,Spot"]
      
      #Rules for classifying perimeter cells
      pcDT <- pcDT[,Spot_PA_Sparse := Nuclei_PA_AreaShape_Neighbors < neighborsThresh]
      
      #Add a local wedge ID to each cell based on conversations with Michel Nederlof
      pcDT <- pcDT[,Spot_PA_Wedge:=ceiling(Nuclei_PA_AreaShape_Center_Theta/wedgeAngs)]
      
      #Define the perimeter cell if it exists in each wedge
      #Classify cells as outer if they have a radial position greater than a thresh
      pcDT <- pcDT[,Spot_PA_OuterCell := labelOuterCells(Nuclei_PA_AreaShape_Center_R, thresh=outerThresh),by="Barcode,Well,Spot"]
      
      #Require a perimeter cell not be in a sparse region
      denseOuterDT <- pcDT[!pcDT$Spot_PA_Sparse  & pcDT$Spot_PA_OuterCell]
      denseOuterDT <- denseOuterDT[,Spot_PA_Perimeter := findPerimeterCell(.SD) ,by="Barcode,Well,Spot,Spot_PA_Wedge"]
      setkey(pcDT,Barcode,Well,Spot,ObjectNumber)
      setkey(denseOuterDT,Barcode,Well,Spot,ObjectNumber)
      pcDT <- denseOuterDT[,list(Barcode,Well,Spot,ObjectNumber,Spot_PA_Perimeter)][pcDT]
      pcDT$Spot_PA_Perimeter[is.na(pcDT$Spot_PA_Perimeter)] <- FALSE
      
    }
    ##Remove nuclear objects that dont'have cell and cytoplasm data
    if(any(grepl("SS1|SS3|SS6|SSD|SSF",ss))) pcDT <- pcDT[!is.na(pcDT$Cells_CP_AreaShape_MajorAxisLength),]
    return(pcDT)
    #}, mc.cores=max(4, detectCores()))#Revert to apply when debugging
  })
  cDT <- rbindlist(expDTList, fill = TRUE)
  rm(expDTList)
  gc()
  
  #Add the pin diameter metadata in microns
  if(grepl("MCF7|PC3|YAPC",cellLine)){
    cDT$PinDiameter <- 180
  } else {
    cDT$PinDiameter <- 350
  }
  
  # The cell level raw data and metadata is saved as Level 1 data. 
  if(writeFiles){
    #Write out cDT without normalized values as level 1 dataset
    level1Names <- grep("Norm|RUV|Loess$",colnames(cDT),value=TRUE,invert=TRUE)
    if(verbose) cat("Writing level 1 file to disk\n")
    writeTime<-Sys.time()
    
    fwrite(data.table(format(cDT[,level1Names, with=FALSE], digits = 4, trim=TRUE)), paste0( "MEP_LINCS/AnnotatedData/", datasetName,"_",ss,"_","Level1.txt"), sep = "\t", quote=FALSE)
    cat("Write time:", Sys.time()-writeTime,"\n")
    
    #### SpotLevel ####
    #The cell-level data is median summarized to the spot level and saved to disk
    if(verbose) cat("Creating spot level data\n")
    slDT <- createl3(cDT, lthresh, seNames = seNames)
    rm(cDT)
    gc()
    fwrite(data.table(format(slDT, digits = 4, trim=TRUE)), paste0( "MEP_LINCS/AnnotatedData/", datasetName,"_",ss,"_","SpotLevel.txt"), sep = "\t", quote=FALSE)
    
    #Write the pipeline parameters to  tab-delimited file
    write.table(c(
      ss=ss,
      cellLine = cellLine,
      analysisVersion = analysisVersion,
      rawDataVersion = rawDataVersion,
      neighborhoodNucleiRadii = neighborhoodNucleiRadii,
      neighborsThresh = neighborsThresh,
      wedgeAngs = wedgeAngs,
      outerThresh = outerThresh,
      nuclearAreaThresh = nuclearAreaThresh,
      nuclearAreaHiThresh = nuclearAreaHiThresh,
      curatedOnly = curatedOnly,
      curatedCols = curatedCols,
      writeFiles = writeFiles,
      limitBarcodes = limitBarcodes,
      k = k,
      normToSpot = normToSpot,
      lowSpotCellCountThreshold = lowSpotCellCountThreshold,
      lowRegionCellCountThreshold = lowRegionCellCountThreshold,
      lowWellQAThreshold = lowWellQAThreshold,
      lowReplicateCount =lowReplicateCount,
      lthresh = lthresh
    ),
    paste0("MEP_LINCS/AnnotatedData/", datasetName,"_",ss,"_","PipelineParameters.txt"), sep = "\t",col.names = FALSE, quote=FALSE)
  }
  cat("Elapsed time:", Sys.time()-startTime, "\n")
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
                       k=c(135),
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
                          k=c(7),
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
                       k=c(64),
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
                       k=c(64),
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
                     calcAdjacency=TRUE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     useAnnotMetadata=FALSE,
                     stringsAsFactors=FALSE)

Baylor <- data.frame(datasetName=c("Baylor1", "Baylor2"),
                     cellLine=c("Baylor1", "Baylor2"),
                     ss=c("SSD"),
                     drug=c("unknown"),
                     analysisVersion="av1.7",
                     rawDataVersion="v2",
                     limitBarcodes=c(5),
                     k=c(64),
                     calcAdjacency=FALSE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = FALSE,
                     useAnnotMetadata=FALSE,
                     stringsAsFactors=FALSE)

MLDDataSet <- data.frame(datasetName=c("MCF10ANeratinib","MCF10ADMSO"),
                         cellLine=c("MCF10A"),
                         ss=c("SSF"),
                         drug=c("Neratinib","DMSO"),
                         analysisVersion="av1.7",
                         rawDataVersion="v2",
                         limitBarcodes=c(8),
                         k=c(64),
                         calcAdjacency=TRUE,
                         writeFiles = TRUE,
                         mergeOmeroIDs = TRUE,
                         useAnnotMetadata=TRUE,
                         stringsAsFactors=FALSE)



library(XLConnect)
library(data.table)

tmp <- apply(Bornstein, 1, preprocessMEPLINCSL1Spot, verbose=TRUE)

