#!/bin/bash

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
# 1/16/17


#Debug dcoument this library stucture
.libPaths(c("/home/users/dane/R/x86_64-pc-linux-gnu-library/3.3","~/R/x86_64-redhat-linux-gnu-library/3.3"))
library("parallel")#use multiple cores for faster processing
source("MEP_LINCS/Release/MEPLINCSFunctions.R")
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.3"))


processan2omero <- function (fileNames) {
  rbindlist(lapply(fileNames, function(fn){
    #Process each file separately
    dt <- fread(fn,header = TRUE)
    
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
    #  }, mc.cores=max(4, detectCores())))
  }))
  
}

spotNorm <- function(x){
  return(x/median(x,na.rm=TRUE))
}

gateOnQuantile <- function(x,probs){
  gatedClass <- integer(length(x))
  gatedClass[x>quantile(x,probs=probs,na.rm=TRUE)]<-1
  return(gatedClass)
}

preprocessMEMACell <- function(barcodePath, verbose=FALSE){
  barcode <- gsub(".*/","",barcodePath)
  path <- gsub(barcode,"",barcodePath)
  if (verbose) cat("Processing plate:",barcode,"at",path,"\n")
  functionStartTime<- Sys.time()
  startTime<- Sys.time()
  
  analysisVersion<-"v1.7"
  rawDataVersion<-"v2"
  limitBarcodes<- NULL
  mergeOmeroIDs<-TRUE
  calcAdjacency<-TRUE
  writeFiles<-TRUE
  useAnnotMetadata<-TRUE
  seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")
  
  #library(limma)#read GAL file and strsplit2
  library(MEMA)#merge, annotate and normalize functions
  library(data.table)#fast file reads, data merges and subsetting
  library(parallel)#use multiple cores for faster processing
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
  
  #Use metadata from an2omero files
  if(useAnnotMetadata){
    metadata <- processan2omero(paste0(barcodePath,"/Analysis/",barcode,"_an2omero.csv"))
  } else {
    stop("Non-An! metadata not supported in this script")
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
  
  # Next, the data is filtered to remove objects with a nuclear area less than nuclearAreaThresh pixels or more than nuclearAreaHiThresh pixels.
  
  #The next steps are to bring in the well metadata, the print order and the CP data
  cellDataFilePaths <- dir(paste0(barcodePath,"/Analysis/",rawDataVersion), full.names = TRUE)
  if(length(cellDataFilePaths)==0) stop("No raw data files found")
  dataBWInfo <- data.table(Path=cellDataFilePaths,
                           Well=gsub("_","",str_extract(dir(paste0(barcodePath,"/Analysis/",rawDataVersion)),"_.*_")),
                           Location=str_extract(cellDataFilePaths,"Nuclei|Cytoplasm|Cells|Image"))
  
  startTime <- Sys.time()
  debugLimiter <- 8
  expDTList <- mclapply(unique(dataBWInfo$Well)[1:debugLimiter], function(well){
    if(verbose) cat(paste("Reading and annotating data for",barcode, well,"\n"))
    nuclei <- convertColumnNames(fread(dataBWInfo$Path[grepl("Nuclei",dataBWInfo$Location)&grepl(well,dataBWInfo$Well)]))
    if (curatedOnly) nuclei <- nuclei[,grep(curatedCols,colnames(nuclei)), with=FALSE]
    setnames(nuclei,paste0("Nuclei_",colnames(nuclei)))
    setnames(nuclei,"Nuclei_CP_ImageNumber","Spot")
    setnames(nuclei,"Nuclei_CP_ObjectNumber","ObjectNumber")
    setkey(nuclei,Spot,ObjectNumber) 
    
    if(any(grepl("Cells",dataBWInfo$Location)&grepl(well,dataBWInfo$Well))){
      cells <- convertColumnNames(fread(dataBWInfo$Path[grepl("Cells",dataBWInfo$Location)&grepl(well,dataBWInfo$Well)]))
      if (curatedOnly) cells <- cells[,grep(curatedCols,colnames(cells)), with=FALSE]
      setnames(cells,paste0("Cells_",colnames(cells)))
      setnames(cells,"Cells_CP_ImageNumber","Spot")
      setnames(cells,"Cells_CP_ObjectNumber","ObjectNumber")
      setkey(cells,Spot,ObjectNumber)
    } 
    
    if(any(grepl("Cytoplasm",dataBWInfo$Location)&grepl(well,dataBWInfo$Well))){
      cytoplasm <- convertColumnNames(fread(dataBWInfo$Path[grepl("Cytoplasm",dataBWInfo$Location)&grepl(well,dataBWInfo$Well)]))
      if (curatedOnly) cytoplasm <- cytoplasm[,grep(curatedCols,colnames(cytoplasm)), with=FALSE]
      setnames(cytoplasm,paste0("Cytoplasm_",colnames(cytoplasm)))
      setnames(cytoplasm,"Cytoplasm_CP_ImageNumber","Spot")
      setnames(cytoplasm,"Cytoplasm_CP_ObjectNumber","ObjectNumber")
      setkey(cytoplasm,Spot,ObjectNumber)
    } 
    
    #Merge the data from the different locations if it exists
    if(exists("cells")&exists("cytoplasm")) {
      dt <- cells[cytoplasm[nuclei]]
    } else {
      dt <- nuclei
    }
    #Add the well name and barcode as parameters
    dt <- dt[,Well := well]
    dt <- dt[,Barcode := barcode]
    #Remove problematic features
    dt <- dt[,grep("Euler",colnames(dt),invert=TRUE), with=FALSE]
    if (useAnnotMetadata) {
      dtm <- merge(dt,metadata,by = c("Barcode","Well","Spot"))
      dtm$PrintSpot <- dtm$Spot
      
    } else {
      #Merge the data with its metadata based on the row it's in
      m <- regexpr("[[:alpha:]]",well)
      row <- regmatches(well,m)
      setkey(DT,Spot)
      dtm <- switch(row, A = merge(DT,spotMetadata,all.x=TRUE),
                    B = merge(DT,spotMetadata180,all.x=TRUE))
      #Add a WellSpace spot index that recognizes the arrays are rotated 180 degrees
      nrArrayRows <- max(dtm$ArrayRow)
      nrArrayColumns <- max(dtm$ArrayColumn)
      dtm$PrintSpot <- dtm$Spot
      dtm$PrintSpot[grepl("B", dtm$Well)] <- (nrArrayRows*nrArrayColumns+1)-dtm$PrintSpot[grepl("B", dtm$Well)]
    }
    
    if(any(grepl("Nuclei_CP_AreaShape_Area",colnames(dtm)))){
      dtm <- dtm[dtm$Nuclei_CP_AreaShape_Area > nuclearAreaThresh,]
      dtm <- dtm[dtm$Nuclei_CP_AreaShape_Area < nuclearAreaHiThresh,]
    }
    
    #Change Edu back to EdU
    if(any(grepl("Edu",colnames(dtm)))){
      edUNames <- grep("Edu",colnames(dtm),value=TRUE)
      setnames(dtm,edUNames,gsub("Edu","EdU",edUNames))
    }
    
    
    #log transform all intensity and areaShape values
    intensityNames <- grep("Intensity",colnames(dtm), value=TRUE)
    scaledInts <- dtm[,intensityNames, with=FALSE]*2^16
    dtm <- cbind(dtm[,!intensityNames, with=FALSE],scaledInts)
    transformNames <- grep("_Center_|_Eccentricity|_Orientation",grep("Intensity|AreaShape",colnames(dtm), value=TRUE, ignore.case = TRUE), value=TRUE, invert=TRUE)
    dtLog <- dtm[,lapply(.SD,boundedLog2),.SDcols=transformNames]
    setnames(dtLog,colnames(dtLog),paste0(colnames(dtLog),"Log2"))
    dtm <- cbind(dtm,dtLog)
    
    #logit transform eccentricity
    if(any(grepl("Nuclei_CP_AreaShape_Eccentricity",colnames(dtm)))){
      dtm <- dtm[,Nuclei_CP_AreaShape_EccentricityLogit := boundedLogit(Nuclei_CP_AreaShape_Eccentricity)]
    }
    if(any(grepl("Nuclei_CP_AreaShape_Center",colnames(dtm)))){
      #Add the local polar coordinates and Neighbor Count
      dtm <- dtm[,Nuclei_PA_Centered_X :=  Nuclei_CP_AreaShape_Center_X-median(Nuclei_CP_AreaShape_Center_X)]
      dtm <- dtm[,Nuclei_PA_Centered_Y :=  Nuclei_CP_AreaShape_Center_Y-median(Nuclei_CP_AreaShape_Center_Y)]
      dtm <- dtm[, Nuclei_PA_AreaShape_Center_R := sqrt(Nuclei_PA_Centered_X^2 + Nuclei_PA_Centered_Y^2)]
      dtm <- dtm[, Nuclei_PA_AreaShape_Center_Theta := calcTheta(Nuclei_PA_Centered_X, Nuclei_PA_Centered_Y)]
    }
    #Add MEP and convenience labels for wells and ligands
    dtm <- dtm[,MEP:=paste(ECMp,Ligand,sep = "_")]
    dtm <- dtm[,Well_Ligand:=paste(Well,Ligand,sep = "_")]
    
    # Eliminate Variations in the Endpoint metadata
    endpointNames <- grep("End",colnames(dtm), value=TRUE)
    endpointWL <- regmatches(endpointNames,regexpr("[[:digit:]]{3}|DAPI",endpointNames))
    setnames(dtm,endpointNames,paste0("Endpoint",endpointWL))
    
    
    if(normToSpot){
      #Add spot level normalizations for selected intensities
      intensityNamesAll <- grep("_CP_Intensity_Median",colnames(dtm), value=TRUE)
      intensityNames <- grep("Norm",intensityNamesAll,invert=TRUE,value=TRUE)
      for(intensityName in intensityNames){
        #Median normalize the median intensity at each spot
        setnames(dtm,intensityName,"value")
        dtm <- dtm[,paste0(intensityName,"_SpotNorm") := spotNorm(value),by="Barcode,Well,Spot"]
        setnames(dtm,"value",intensityName)
      }
    }
    if(calcAdjacency){
      if(verbose) cat("Calculating adjacency data\n")
      
      densityRadius <- sqrt(median(dtm$Nuclei_CP_AreaShape_Area, na.rm = TRUE)/pi)
      
      #Count the number of neighboring cells
      dtm <- dtm[,Nuclei_PA_AreaShape_Neighbors := cellNeighbors(.SD, radius = densityRadius*neighborhoodNucleiRadii), by = "Barcode,Well,Spot"]
      
      #Rules for classifying perimeter cells
      dtm <- dtm[,Spot_PA_Sparse := Nuclei_PA_AreaShape_Neighbors < neighborsThresh]
      
      #Add a local wedge ID to each cell based on conversations with Michel Nederlof
      dtm <- dtm[,Spot_PA_Wedge:=ceiling(Nuclei_PA_AreaShape_Center_Theta/wedgeAngs)]
      
      #Define the perimeter cell if it exists in each wedge
      #Classify cells as outer if they have a radial position greater than a thresh
      dtm <- dtm[,Spot_PA_OuterCell := labelOuterCells(Nuclei_PA_AreaShape_Center_R, thresh=outerThresh),by="Barcode,Well,Spot"]
      
      #Require a perimeter cell not be in a sparse region
      denseOuterDT <- dtm[!dtm$Spot_PA_Sparse  & dtm$Spot_PA_OuterCell]
      denseOuterDT <- denseOuterDT[,Spot_PA_Perimeter := findPerimeterCell(.SD) ,by="Barcode,Well,Spot,Spot_PA_Wedge"]
      setkey(dtm,Barcode,Well,Spot,ObjectNumber)
      setkey(denseOuterDT,Barcode,Well,Spot,ObjectNumber)
      dtm <- denseOuterDT[,list(Barcode,Well,Spot,ObjectNumber,Spot_PA_Perimeter)][dtm]
      dtm$Spot_PA_Perimeter[is.na(dtm$Spot_PA_Perimeter)] <- FALSE
      
    }
    #Add the pin diameter metadata in microns
    if(any(grepl("MCF7|PC3|YAPC",unique(dtm$CellLine)))){
      dtm$PinDiameter <- 180
    } else {
      dtm$PinDiameter <- 350
    }
    #names(dtm) <- dataBWI
    return(dtm)
  }, mc.cores=detectCores())
  #Add names to the data.tables in the list
  names(expDTList) <- paste(barcode,unique(dataBWInfo$Well)[1:debugLimiter],sep="_")
  cat("converting list to datatable\n")
  startTime <- Sys.time()
  cDT <- rbindlist(expDTList)
  rm(expDTList)
  gc()
  cat(Sys.time()-startTime)
  #   
  #   if(!useAnnotMetadata){
  #     #Read the well metadata from a multi-sheet Excel file
  #     #wellMetadata <- data.table(readMetadata(paste0(filePath,"/Metadata/",gsub("reDAPI","",barcode),".xlsx")), key="Well")
  #     #Debug reDAPI plates with new file structure
  #     wellMetadata <- data.table(readMetadata(fileNames$Path[fileNames$Barcode==barcode&fileNames$Type=="metadata"]), key="Well")
  #     #merge well metadata with the data and spot metadata
  #     pcDT <- merge(pcDT,wellMetadata,by = "Well")
  #     #Add a WellSpace spot index that recognizes the arrays are rotated 180 degrees
  #     nrArrayRows <- max(pcDT$ArrayRow)
  #     nrArrayColumns <- max(pcDT$ArrayColumn)
  #     pcDT$PrintSpot <- pcDT$Spot
  #     pcDT$PrintSpot[grepl("B", pcDT$Well)] <- (nrArrayRows*nrArrayColumns+1)-pcDT$PrintSpot[grepl("B", pcDT$Well)]
  #   }
  #   
  if(mergeOmeroIDs){
    #Read in and merge the Omero URLs
    omeroIndex <- fread(paste0(barcodePath,"/Analysis/",barcode,"_imageIDs.tsv"))[,list(WellName,Row,Column,ImageID)]
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
    cDT <- merge(cDT,omeroIndex,by=c("Well","ArrayRow","ArrayColumn"))
  }
  
  #Set 2N and 4N DNA status
  cDT <- cDT[,Nuclei_PA_Cycle_State := gateOnlocalMinima(Nuclei_CP_Intensity_IntegratedIntensity_Dapi)]
  
  
  if(!useAnnotMetadata){
    #Create short display names, then replace where not unique
    #Use entire AnnotID for ligands with same uniprot IDs
    # cDT$Ligand[grepl("NRG1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "NRG1",annotIDs = cDT$Ligand[grepl("NRG1",cDT$Ligand)])
    # cDT$Ligand[grepl("TGFB1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "TGFB1",annotIDs = cDT$Ligand[grepl("TGFB1",cDT$Ligand)])
    # cDT$Ligand[grepl("CXCL12",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "CXCL12",annotIDs = cDT$Ligand[grepl("CXCL12",cDT$Ligand)])
    # cDT$Ligand[grepl("IGF1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "IGF1",annotIDs = cDT$Ligand[grepl("IGF1",cDT$Ligand)])
  }
  
  #Create staining set specific derived parameters
  if("Nuclei_CP_Intensity_MedianIntensity_EdU" %in% colnames(cDT)){
    #Use the entire plate to set the autogate threshold if there's no control well
    if(any(grepl("FBS",unique(cDT$Ligand)))){
      cDT <- cDT[,Nuclei_PA_Gated_EdUPositive := kmeansCluster(.SD,value =  "Nuclei_CP_Intensity_MedianIntensity_EdU",ctrlLigand = "FBS"), by="Barcode"]
    } else {
      cDT <- cDT[,Nuclei_PA_Gated_EdUPositive := kmeansCluster(.SD,value =  "Nuclei_CP_Intensity_MedianIntensity_EdU",ctrlLigand = "."), by="Barcode"]
    }
    
    #Modify the gate threshold to be 3 sigma from the EdU- median
    EdUDT <- cDT[Nuclei_PA_Gated_EdUPositive==0,.(EdUMedian=median(Nuclei_CP_Intensity_MedianIntensity_EdULog2, na.rm=TRUE),EdUSD=sd(Nuclei_CP_Intensity_MedianIntensity_EdULog2, na.rm=TRUE)),by="Barcode,Well"]
    EdUDT <- EdUDT[,Median3SD := EdUMedian+3*EdUSD]
    cDT <- merge(cDT,EdUDT,by=c("Barcode","Well"))
    cDT$Nuclei_PA_Gated_EdUPositive[cDT$Nuclei_CP_Intensity_MedianIntensity_EdULog2>cDT$Median3SD]<-1
    #Overide autogate for AU565 cell line in Lapatinib MEMAs
    # if(grepl("AU565_Lapatinib|AU565_DMSO",datasetName)){
    #   cDT$Nuclei_PA_Gated_EdUPositive <- 0
    #   cDT$Nuclei_PA_Gated_EdUPositive[cDT$Nuclei_CP_Intensity_MedianIntensity_EdULog2>10] <- 1
    # }
    
  }
  
  #Calculate a lineage ratio of luminal/basal or KRT19/KRT5
  if ("Cytoplasm_CP_Intensity_MedianIntensity_KRT19" %in% colnames(cDT)&
      "Cytoplasm_CP_Intensity_MedianIntensity_KRT5" %in% colnames(cDT)){
    cDT <- cDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT5]
    cDT <- cDT[,Cytoplasm_PA_Intensity_LineageRatioLog2 := boundedLog2(Cytoplasm_PA_Intensity_LineageRatio)]
    #Gate the KRT signals using kmeans clustering
    
    if(grepl("HMEC",unique(cDT$CellLine))) {
      cDT <- cDT[,Cytoplasm_PA_Gated_KRT5Positive := gateOnQuantile(Cytoplasm_CP_Intensity_MedianIntensity_KRT5Log2,probs=.02),by="Barcode,Well"]
    } else {
      cDT <- cDT[,Cytoplasm_PA_Gated_KRT5Positive := kmeansCluster(.SD,value =  "Cytoplasm_CP_Intensity_MedianIntensity_KRT5",ctrlLigand = "."), by="Barcode"]
    }
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT19Positive := kmeansCluster(.SD,value =  "Cytoplasm_CP_Intensity_MedianIntensity_KRT19",ctrlLigand = "."), by="Barcode"]
    #Determine the class of each cell based on KRT5 and KRT19 class
    #0 double negative
    #1 KRT5+, KRT19-
    #2 KRT5-, KRT19+
    #3 KRT5+, KRT19+
    cDT <- cDT[,Cytoplasm_PA_Gated_KRTClass := Cytoplasm_PA_Gated_KRT5Positive+2*Cytoplasm_PA_Gated_KRT19Positive]
    
  }
  if("Nuclei_PA_Gated_EdUPositive" %in% colnames(cDT)&
     "Cytoplasm_PA_Gated_KRT19Positive" %in% colnames(cDT)){
    #Determine the class of each cell based on KRT19 and EdU class
    #0 double negative
    #1 KRT19+, EdU-
    #2 KRT19-, EdU+
    #3 KRT19+, EdU+
    cDT <- cDT[,Cells_PA_Gated_EdUKRT19Class := Cytoplasm_PA_Gated_KRT19Positive+2*Nuclei_PA_Gated_EdUPositive]
    
  }
  
  if ("Cytoplasm_CP_Intensity_MedianIntensity_KRT14" %in% colnames(cDT)&
      "Cytoplasm_CP_Intensity_MedianIntensity_KRT19" %in% colnames(cDT)){
    #Calculate a lineage ratio of luminal/basal or KRT19/KRT14
    cDT <- cDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT14]
    cDT <- cDT[,Cytoplasm_PA_Intensity_LineageRatioLog2 := boundedLog2(Cytoplasm_PA_Intensity_LineageRatio)]
  }
  if ("Cytoplasm_CP_Intensity_MedianIntensity_KRT5Log2" %in% colnames(cDT)&
      "Nuclei_PA_Gated_EdUPositive" %in% colnames(cDT)){
    #Determine the class of each cell based on KRT5 and EdU class
    #0 double negative
    #1 KRT5+, EdU-
    #2 KRT5-, EdU+
    #3 KRT5+, EdU+
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT5Positive := gateOnQuantile(Cytoplasm_CP_Intensity_MedianIntensity_KRT5Log2,probs=.02),by="Barcode,Well"]
    cDT <- cDT[,Cells_PA_Gated_EdUKRT5Class := Cytoplasm_PA_Gated_KRT5Positive+2*Nuclei_PA_Gated_EdUPositive]
  }
  
  ##Remove nuclear objects that dont'have cell and cytoplasm data
  if("Cells_CP_AreaShape_MajorAxisLength" %in% colnames(cDT)) cDT <- cDT[!is.na(cDT$Cells_CP_AreaShape_MajorAxisLength),]
  
  # The cell level raw data and metadata is saved as Level 1 data. 
  if(writeFiles){
    #Write out cDT without normalized values as level 1 dataset
    if(verbose) cat("Writing full level 1 file to disk\n")
    writeTime<-Sys.time()
    set.seed(42)
    fwrite(data.table(format(cDT, digits = 4, trim=TRUE)), paste0(barcodePath, "/Analysis/", barcode,"_","Level1.tsv"), sep = "\t", quote=FALSE)
    cat("Write time:", Sys.time()-writeTime,"\n")
    
    #Write the pipeline parameters to  tab-delimited file
    write.table(c(
      cellLine = unique(cDT$CellLine),
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
      normToSpot = normToSpot,
      lowSpotCellCountThreshold = lowSpotCellCountThreshold,
      lowRegionCellCountThreshold = lowRegionCellCountThreshold,
      lowWellQAThreshold = lowWellQAThreshold,
      lowReplicateCount =lowReplicateCount,
      lthresh = lthresh
    ),
    paste0(barcodePath, "/Analysis/", barcode,"_","Level1PipelineParameters.txt"), sep = "\t",col.names = FALSE, quote=FALSE)
  }
  cat("Elapsed time:", Sys.time()-functionStartTime, "\n")
}


preprocessMEMASpot <- function(barcodePath, verbose=FALSE){
  barcode <- gsub(".*/","",barcodePath)
  path <- gsub(barcode,"",barcodePath)
  if (verbose) cat("Summarizing cell to spot data for plate",barcode,"at",barcodePath,"\n")
  functionStartTime<- Sys.time()
  startTime<- Sys.time()
  writeFiles<-TRUE
  seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")
  
  #library(limma)#read GAL file and strsplit2
  library(MEMA)#merge, annotate and normalize functions
  library(data.table)#fast file reads, data merges and subsetting
  library(parallel)#use multiple cores for faster processing
  library(stringr)
  
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
  
  #Read in the plate's cell level data
  cDT <- fread(paste0(barcodePath,"/Analysis/",barcode,"_Level1.tsv"))
  
  #Calculate the proportions in the gates
  calcProportion <- function(x){
    sum(x)/length(x)
  }
  
  #Count the cells at each spot at the cell level as needed by createl3
  cDT <- cDT[,Spot_PA_SpotCellCount := .N,by="Barcode,Well,Spot"]
  cDT <- cDT[,Spot_PA_SpotCellCountLog2 := boundedLog2(Spot_PA_SpotCellCount)]
  
  #Calculate proportions for gated signals
  gatedSignals <- grep("Positive|High",colnames(cDT), value=TRUE)
  proportions <- cDT[,lapply(.SD, calcProportion),by="Barcode,Well,Spot", .SDcols=gatedSignals]
  setnames(proportions,
           grep("Gated",colnames(proportions),value=TRUE),
           paste0(grep("Gated",colnames(proportions),value=TRUE),"Proportion"))
  #Calculate logits of proportions
  proportionSignals <- grep("Proportion",colnames(proportions), value=TRUE)
  proportionLogits <- proportions[,lapply(.SD, boundedLogit),by="Barcode,Well,Spot", .SDcols=proportionSignals]
  setnames(proportionLogits,
           grep("Proportion",colnames(proportionLogits),value=TRUE),
           paste0(grep("Proportion",colnames(proportionLogits),value=TRUE),"Logit"))
  proportionsDT <-merge(proportions,proportionLogits)
  
  #Create proportions and logits of mutlivariate gates
  if ("Cytoplasm_PA_Gated_KRTClass" %in% colnames(cDT)){
    #Determine the class of each cell based on KRT5 and KRT19 class
    #0 double negative
    #1 KRT5+, KRT19-
    #2 KRT5-, KRT19+
    #3 KRT5+, KRT19+
    #Calculate gating proportions for EdU and KRT19
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT5RT19NegativeProportion := sum(Cytoplasm_PA_Gated_KRTClass==0)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cytoplasm_PA_Gated_KRT5RT19NegativeProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5RT19NegativeProportion)]
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT5NegativeKRT19PositiveProportion := sum(Cytoplasm_PA_Gated_KRTClass==2)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_KRT5NegativeKRT19PositiveProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5NegativeKRT19PositiveProportion)]
    cDT <- cDT[,Cells_PA_Gated_KRT5PositiveKRT19NegativeProportion := sum(Cytoplasm_PA_Gated_KRTClass==1)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_KRT5PositiveKRT19NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_KRT5PositiveKRT19NegativeProportion)]
    cDT <- cDT[,Cells_PA_Gated_KRT5PositivedKRT19PositiveProportion := sum(Cytoplasm_PA_Gated_KRTClass==3)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_Cells_PA_Gated_KRT5PositivedKRT19PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_KRT5PositivedKRT19PositiveProportion)]
  }
  
  if ("Cells_PA_Gated_EdUKRT5Class" %in% colnames(cDT)){
    #Determine the class of each cell based on KRT5 and EdU class
    #0 double negative
    #1 KRT5+, EdU-
    #2 KRT5-, EdU+
    #3 KRT5+, EdU+
    #Calculate gating proportions for EdU and KRT5
    cDT <- cDT[,Cells_PA_Gated_EdUKRT5NegativeProportion := sum(Cells_PA_Gated_EdUKRT5Class==0)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUKRT5NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_EdUKRT5NegativeProportion)]
    cDT <- cDT[,Cells_PA_Gated_EdUPositiveKRT5NegativeProportion := sum(Cells_PA_Gated_EdUKRT5Class==2)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUPositiveKRT5NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_EdUPositiveKRT5NegativeProportion)]
    cDT <- cDT[,Cells_PA_Gated_EdUNegativeKRT5PositiveProportion := sum(Cells_PA_Gated_EdUKRT5Class==1)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUNegativeKRT5PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_EdUNegativeKRT5PositiveProportion)]
    cDT <- cDT[,Cells_PA_Gated_EdUKRT5PositiveProportion := sum(Cells_PA_Gated_EdUKRT5Class==3)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUKRT5PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_EdUKRT5PositiveProportion)]
  }
  #median summaerize the rest of the signals
  signalDT <- createl3(cDT,lthresh, seNames)
  spotDT <- merge(signalDT,proportionsDT)
  
  #Summarize 
  # The cell level raw data and metadata is saved as Level 1 data. 
  if(writeFiles){
    #Write out cDT without normalized values as level 1 dataset
    if(verbose) cat("Writing spot level data to disk\n")
    writeTime<-Sys.time()
    fwrite(data.table(format(spotDT, digits = 4, trim=TRUE)), paste0(barcodePath, "/Analysis/", barcode,"_","SpotLevel.tsv"), sep = "\t", quote=FALSE)
    cat("Write time:", Sys.time()-writeTime,"\n")
  }
  cat("Elapsed time:", Sys.time()-functionStartTime, "\n")
}

barcodePath <-commandArgs(trailingOnly = TRUE)
res <- preprocessMEMACell(barcodePath, verbose=TRUE)
res <- preprocessMEMASpot(barcodePath, verbose=TRUE)

