#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 11/18/2015

##Introduction

#   The MEP-LINCs dataset contains imaging data from a Nikon automated microscope that is analyzed with a CellProfiler pipeline.
# 
# Part of this preprocessing of the dataset will be deprecated when the merging of the data and metadata happens within the CellProfiler part of the pipeline. For now, the metadata about the ECM proteins is read from the GAL file and the metadata about the wells (cell line, stains and ligands) is read from Excel spreadsheets.

library("parallel")#use multiple cores for faster processing
source("MEP_LINCS/MEP_LINCS/MEPLINCSFunctions.R")

getFactors <- function (factors) {
  sl <- lapply(factors, function(factor){
    the_set <- factor$the_set
    if(is.null((the_set))) the_set <-"NoSet"
    content <- factor$content
    contentType <- switch(str_split(the_set,"")[[1]][1],N="ECMpBase",L="Ligand",E="ECMp","UnknownContentType")
    data.table(Content = content, Set = the_set, ContentType = contentType)
  })
  sDT <- rbindlist(sl)
  sDT$ContentShortName <- names(factors)
  return(sDT)
}

getParameters <- function(parameters) {
  sl <- lapply(parameters, function(parameter){
    ss <- parameter$the_set
    pContent <- parameter$content
    data.table(ParameterContent = pContent, StainingSet = ss)
  })
  sDT <- rbindlist(sl)
  splits <- str_split_fixed(sDT$ParameterContent,"-",n=5)
  sDT$StainType <- splits[,1]
  sDT$Endpoint <- ""
  sDT$Endpoint[sDT$StainType %in% c("cstain", "antibody1")] <-strsplit2(splits[sDT$StainType %in% c("cstain", "antibody1"),2], "_")[,1]
  sDT$Channel <- ""
  sDT$Channel[sDT$StainType %in% c("cstain")] <-splits[sDT$StainType %in% c("cstain"),3]
  sDT$Channel[sDT$StainType %in% c("antibody2")] <-splits[sDT$StainType %in% c("antibody2"),4]
  sDT$Animal <- ""
  if(any(sDT$StainType %in% c("antibody1","antibody2"))){
    sDT$Animal[sDT$StainType %in% c("antibody1","antibody2")] <-strsplit2(splits[sDT$StainType %in% c("antibody1","antibody2"),3], "_")[,1]
    #Match antibody1 animal to antibody2 animal
    sDTA1 <- sDT[sDT$StainType=="antibody1",list(Endpoint,Animal)]
    sDTA2 <- sDT[sDT$StainType=="antibody2",list(Channel,Animal)]
    ch <- merge(sDTA1,sDTA2, by="Animal")
    #then copy antibody2 channel into antibody 1 channel
    sDT <- merge(sDT,ch,by="Endpoint", all=TRUE)
    sDT$Channel.y[is.na(sDT$Channel.y)] <- ""
    sDT$Channel <- paste0(sDT$Channel.x,sDT$Channel.y)
    sDT <- sDT[,list(Endpoint,StainingSet,StainType,Channel,Animal.x)]
    sDT <- sDT[!sDT$Endpoint==""]
    setnames(sDT,"Animal.x","Animal")
  }
  return(sDT)
}

getSample <- function (samples) {
  sl <- lapply(samples, function(sample){
    data.table(CellLine = gsub(".*-","",gsub("_.*","",sample$content)), Passage = sample$fraction, CellSeedCount = sample$value)
  })
  sDT <- rbindlist(sl)
}

processJSON <- function (fileNames) {
  rbindlist(lapply(fileNames, function(fn){
    #Process each file separately
    plateMetadata <- fromJSON(fn)
    
    #Store and remove the welltype data
    welltype <- plateMetadata$welltype
    
    #Get nr rows and columns and use to calculate spot index
    columnParms <- as.integer(str_split_fixed(str_split_fixed(welltype,"[|]",2)[1,2],"_",3))
    nrArrayColumns <- columnParms[2]*columnParms[3]
    rowParms <- as.integer(str_split_fixed(str_split_fixed(welltype,"[|]",2)[1,1],"_",3))
    nrArrayRows <- rowParms[2]*rowParms[3]
    nrArraySpots <- nrArrayRows*nrArrayColumns
    
    #Get the barcode
    barcode <- plateMetadata$assayrun
    
    #Delete all of the items at the end of the json file
    plateMetadata$annot_id <- NULL
    plateMetadata$assayrun <- NULL
    plateMetadata$assaytype <- NULL
    plateMetadata$label <- NULL
    plateMetadata$welltype <- NULL
    
    #Get the well name and spot metadata
    spotmMdList <- lapply(plateMetadata, function(spotMetadata){
      #Get the well alphanumeric
      wellIndices <- as.integer(str_split_fixed(spotMetadata["i|i"],"[|]",2))
      wellName <- wellAN(2,4)[(wellIndices[1]-1)*4+wellIndices[2]]
      #Read the spot specific metadata
      spotFactors <- getFactors(spotMetadata$content$factor)
      #Get the spot information
      ArrayPositions <- str_split_fixed(spotMetadata["ii|ii"],"[|]",2)
      #browser()
      m <- regexpr(".*_",spotFactors$ContentShortName[spotFactors$ContentType=="ECMp"])
      ECMp <- gsub("_","",regmatches(spotFactors$ContentShortName[spotFactors$ContentType=="ECMp"], m))
      mDT <- data.table(Barcode=barcode,
                        Well = wellName,
                        ArrayRow = as.integer(str_split_fixed(ArrayPositions[1,1],"_",2)[1,2]),
                        ArrayColumn = as.integer(str_split_fixed(ArrayPositions[1,2],"_",2)[1,2]),
                        ECMp = ECMp)
      return(mDT)
    })
    spotMdDT <- rbindlist(spotmMdList)
    
    #Read well metadata from the first spot in each well
    wellMdList <- lapply(plateMetadata[as.integer(names(plateMetadata)) %in% seq(1,length(plateMetadata),by=nrArraySpots)], function(spotMetadata){
      #Get the well alphanumeric
      wellIndices <- as.integer(str_split_fixed(spotMetadata["i|i"],"[|]",2))
      wellName <- wellAN(2,4)[(wellIndices[1]-1)*4+wellIndices[2]]
      #Get protein info from factors
      spotFactors <- getFactors(spotMetadata$content$factor)
      #Get stain info from parameters
      spotParameters <- getParameters(spotMetadata$content$parameter)
      #Get cell line info from sample
      cellLineParameters <- getSample(spotMetadata$content$sample)
      m <- regexpr(".*_",spotFactors$ContentShortName[spotFactors$ContentType=="Ligand"])
      ligand <- gsub("_","",regmatches(spotFactors$ContentShortName[spotFactors$ContentType=="Ligand"], m))
      EndpointDAPI <- spotParameters$Endpoint[grepl("DAPI",spotParameters$Channel)]
      Endpoint488 <- spotParameters$Endpoint[grepl("488",spotParameters$Channel)]
      if(length(Endpoint488)==0) Endpoint488<- "none"
      
      Endpoint555 <- spotParameters$Endpoint[grepl("555|Orange",spotParameters$Channel)]
      if(length(Endpoint555)==0) Endpoint555<- "none"
      
      Endpoint647 <- spotParameters$Endpoint[grepl("647|Red",spotParameters$Channel)]
      if(length(Endpoint647)==0) Endpoint647<- "none"
      
      mDT <- data.table(Barcode=barcode,
                        CellLine = cellLineParameters$CellLine,
                        Well = wellName,
                        Ligand = ligand,
                        EndpointDAPI = EndpointDAPI,
                        Endpoint488 = Endpoint488,
                        #StainType488 = spotParameters$StainType[grepl("488",spotParameters$Channel)],
                        #Animal488 = spotParameters$Animal[grepl("488",spotParameters$Channel)],
                        Endpoint555 = Endpoint555,
                        #StainType555 = spotParameters$StainType[grepl("555|Orange",spotParameters$Channel)],
                        #Animal555 = spotParameters$Animal[grepl("555|Orange",spotParameters$Channel)],
                        Endpoint647 = Endpoint647)
      #StainType647 = spotParameters$StainType[grepl("647|Red",spotParameters$Channel)],
      #Animal647 = spotParameters$Animal[grepl("647|Red",spotParameters$Channel)]
    })
    wellMdDT <- rbindlist(wellMdList)
    mdDT <- merge(spotMdDT,wellMdDT,by=c("Barcode","Well"))
    mdDT$Spot <- as.integer(nrArrayColumns*(mdDT$ArrayRow-1)+mdDT$ArrayColumn)
    #Add a WellSpace spot index that recognizes the arrays are rotated 180 degrees
    mdDT$PrintSpot <- mdDT$Spot
    mdDT$PrintSpot[grepl("B", mdDT$Well)] <- (nrArrayRows*nrArrayColumns+1)-mdDT$PrintSpot[grepl("B", mdDT$Well)]
    return(mdDT)
  }))
}


preprocessMEPLINCSL1Spot <- function(ssDataset, verbose=FALSE){
  startTime<- Sys.time()
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
  useJSONMetadata<-as.logical(ssDataset[["useJSONMetadata"]])
  
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
  
  fileNames <- rbindlist(apply(datasetBarcodes[datasetBarcodes$CellLine==cellLine&datasetBarcodes$StainingSet==ss&datasetBarcodes$Drug==drug&datasetBarcodes$Version==rawDataVersion,], 1, getMEMADataFileNames))
  
  ##Summary
  # This script prepares cell-level data and metadata for the MEP LINCs Analysis Pipeline. 
  # 
  # In the code, the variable ss determines which staining set (SS1, SS2 or SS3) to merge and the variable cellLine determines the cell line (PC3,MCF7, etc).
  # 
  # The merging assumes that the actual, physical B row wells (B01-B04) have been printed upside-down. That is, rotated 180 degrees resulting in the spot 1, 1 being in the lower right corner instead of the upper left corner. The metadata is matched to the actual printed orientation.
  
  #Use metadata from Annot JSON files or 
  #directly from GAL and xlsx files
  if(useJSONMetadata){
    #readJSONMetadata
    fns <- fileNames$Path[fileNames$Type=="metadata"]
    metadata <- rbindlist(mclapply(fns, processJSON, mc.cores = detectCores()))
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
  
  barcodes <- unique(fileNames$Barcode)[1:limitBarcodes]
  if(verbose) cat(paste("Reading and annotating cell level data for",cellLine,ss)," \n")
  expDTList <- mclapply(barcodes, function(barcode){
    plateDataFiles <- fileNames$Path[grepl(barcode,fileNames$Barcode)&
                                       grepl("Raw",fileNames$Type)]
    wells <- unique(strsplit2(split = "_",plateDataFiles)[,2])
    wellDataList <- lapply(wells,function(well){
      
      wellDataFiles <- grep(well,plateDataFiles,value = TRUE)
      nucleiDataFile <- grep("Nuclei",wellDataFiles,value=TRUE,
                             ignore.case = TRUE)
      if (ss %in% c("SS1","SS3", "SS6")){
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
      if (ss %in% c("SS1","SS3", "SS6")){
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
      if (ss %in% c("SS1","SS3","SS6")){
        setnames(cells,paste0("Cells_",colnames(cells)))
        setnames(cytoplasm,paste0("Cytoplasm_",colnames(cytoplasm)))
      }
      
      #Merge the cells, cytoplasm and nuclei data
      if (ss %in% c("SS1","SS3","SS6")){
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
      
      if (useJSONMetadata) {
        DT <- merge(DT,metadata,by = c("Barcode","Well","Spot"))
      } else {
        #Merge the data with its metadata based on the row it's in
        m <- regexpr("[[:alpha:]]",well)
        row <- regmatches(well,m)
        setkey(DT,Spot)
        DT <- switch(row, A = merge(DT,spotMetadata,all=TRUE),
                     B = merge(DT,spotMetadata180,all=TRUE))
      }
      return(DT)
    })
    #Create the cell data.table with spot metadata for the plate 
    pcDT <- rbindlist(wellDataList, fill = TRUE)
    #rm(wellDataList)
    
    if(!useJSONMetadata){
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
      #imageURLFiles <- grep("reDAPI", allImageFiles, invert=TRUE, value=TRUE)
      imageURLFile <- fileNames$Path[fileNames$Barcode==barcode&fileNames$Type=="imageID"]
      #Read in and merge the Omero URLs
      omeroIndex <- fread(fileNames$Path[fileNames$Barcode==barcode&fileNames$Type=="imageID"])[,list(WellName,Row,Column,ImageID)]
      omeroIndex$Well <- sapply(gsub("Well","",strsplit2(omeroIndex$WellName,"_")[,2],""),FUN=switch,
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
    
    #Add spot level normalizations for selected intensities
    if(normToSpot){
      intensityNamesAll <- grep("_CP_Intensity_Median",colnames(pcDT), value=TRUE)
      intensityNames <- grep("Norm",intensityNamesAll,invert=TRUE,value=TRUE)
      for(intensityName in intensityNames){
        #Median normalize the median intensity at each spot
        pcDT <- pcDT[,paste0(intensityName,"_SpotNorm") := medianNorm(.SD,intensityName),by="Barcode,Well,Spot"]
      }
    }
    if(!useJSONMetadata){
      #Create short display names, then replace where not unique
      #Use entire AnnotID for ligands with same uniprot IDs
      pcDT$Ligand[grepl("NRG1",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "NRG1",annotIDs = pcDT$Ligand[grepl("NRG1",pcDT$Ligand)])
      pcDT$Ligand[grepl("TGFB1",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "TGFB1",annotIDs = pcDT$Ligand[grepl("TGFB1",pcDT$Ligand)])
      pcDT$Ligand[grepl("CXCL12",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "CXCL12",annotIDs = pcDT$Ligand[grepl("CXCL12",pcDT$Ligand)])
      pcDT$Ligand[grepl("IGF1",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "IGF1",annotIDs = pcDT$Ligand[grepl("IGF1",pcDT$Ligand)])
    }
    pcDT$MEP <- paste(pcDT$ECMp,pcDT$Ligand,sep = "_")
    
    #Add a convenience label for wells and ligands
    pcDT$Well_Ligand <- paste(pcDT$Well,pcDT$Ligand,sep = "_")
    
    # Eliminate Variations in the Endpoint metadata
    endpointNames <- grep("End",colnames(pcDT), value=TRUE)
    endpointWL <- regmatches(endpointNames,regexpr("[[:digit:]]{3}|DAPI",endpointNames))
    setnames(pcDT,endpointNames,paste0("Endpoint",endpointWL))
    
    
    #Create staining set specific derived parameters
    if(grepl("SS0|SS1",ss)){
    } else if (grepl("SS2|SS6",ss)){
      pcDT <- pcDT[,Nuclei_PA_Gated_EdUPositive := kmeansCluster(.SD,value =  "Nuclei_CP_Intensity_MedianIntensity_EdU",ctrlLigand = "FBS"), by="Barcode"]
      #pcDT <- pcDT[,Nuclei_PA_Gated_EduPositive := gateOnlocalMinima(Nuclei_CP_Intensity_MedianIntensity_Edu)-1, by="Barcode,Well"]
      #Calculate the EdU Positive Percent at each spot
      pcDT <- pcDT[,Nuclei_PA_Gated_EdUPositiveProportion := sum(Nuclei_PA_Gated_EdUPositive)/length(ObjectNumber),by="Barcode,Well,Spot"]
      #Logit transform EdUPositiveProportion
      #logit(p) = log[p/(1-p)]
      pcDT <-pcDT[,Nuclei_PA_Gated_EdUPositiveProportionLogit:= boundedLogit(Nuclei_PA_Gated_EdUPositiveProportion)]
      
      
    } else if (grepl("SS3",ss)){
      #Calculate a lineage ratio of luminal/basal or KRT19/KRT5
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT5]
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatioLog2 := boundedLog2(Cytoplasm_PA_Intensity_LineageRatio)]
      
    } else if (grepl("SS6",ss)){
      #Calculate a lineage ratio of luminal/basal or KRT19/KRT14
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT14]
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatioLog2 := boundedLog2(Cytoplasm_PA_Intensity_LineageRatio)]
      
    } else stop("Invalid ss parameter")
    
    #Move adjacency processing to within parallel code
    if(calcAdjacency){
      #cDTList <- mclapply(barcodes, function(barcode, dt){
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
    if(any(grepl("SS1|SS3|SS6",ss))) pcDT <- pcDT[!is.na(pcDT$Cells_CP_AreaShape_MajorAxisLength),]
    return(pcDT)
  }, mc.cores=max(4, detectCores()*.75))#Revert to apply when debugging
  #})
  cDT <- rbindlist(expDTList, fill = TRUE)
  
  rm(expDTList)
  gc()
  
  
 
  # The cell level raw data and metadata is saved as Level 1 data. 
  
  if(writeFiles){
    #Write out cDT without normalized values as level 1 dataset
    level1Names <- grep("Norm|RUV3|Loess$",colnames(cDT),value=TRUE,invert=TRUE)
    if(verbose) cat("Writing level 1 file to disk\n")
    writeTime<-Sys.time()
    fwrite(cDT[,level1Names, with=FALSE], paste0( "./AnnotatedData/", unique(cDT$CellLine),"_",ss,"_",rawDataVersion,"_",analysisVersion,"_","Level1.txt"), sep = "\t", quote=FALSE)
    cat("Write time:", Sys.time()-writeTime,"\n")
    
    #### SpotLevel ####
    #The cell-level data is median summarized to the spot level and saved to disk
    if(verbose) cat("Creating spot level data\n")
    slDT <- createl3(cDT, lthresh, seNames = seNames)
    rm(cDT)
    gc()
    fwrite(slDT, paste0( "./AnnotatedData/", unique(slDT$CellLine),"_",ss,"_",rawDataVersion,"_",analysisVersion,"_","SpotLevel.txt"), sep = "\t", quote=FALSE)
    
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
    paste0("./AnnotatedData/", cellLine,"_",ss,"_",analysisVersion,"_","PipelineParameters.txt"), sep = "\t",col.names = FALSE, quote=FALSE)
  }
  cat("Elapsed time:", Sys.time()-startTime)
}

PC3df <- data.frame(cellLine=rep(c("PC3"), 4),
                    ss=c("SS1", "SS2","SS3","SS2noH3"),
                    drug=c("none"),
                    analysisVersion="av1.6",
                    rawDataVersion=c("v2","v2.1","v2.1", "v1"),
                    limitBarcodes=8,
                    k=7,
                    calcAdjacency=TRUE,
                    writeFiles = TRUE,
                    mergeOmeroIDs = TRUE,
                    useJSONMetadata=TRUE,
                    stringsAsFactors=FALSE)

MCF7df <- data.frame(cellLine=rep(c("MCF7"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug=c("none"),
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     limitBarcodes=8,
                     k=7,
                     calcAdjacency=TRUE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     useJSONMetadata=TRUE,
                     stringsAsFactors=FALSE)

YAPCdf <- data.frame(cellLine=rep(c("YAPC"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug=c("none"),
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     limitBarcodes=8,
                     k=7,
                     calcAdjacency=TRUE,
                     writeFiles = TRUE,
                     mergeOmeroIDs = TRUE,
                     useJSONMetadata=TRUE,
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(cellLine="MCF10A",
                       ss=c("SS1","SS2","SS3"),
                       drug=c("none"),
                       analysisVersion="av1.6",
                       rawDataVersion="v2",
                       limitBarcodes=c(8,8,8),
                       k=c(7,7,7),
                       calcAdjacency=TRUE,
                       writeFiles = TRUE,
                       mergeOmeroIDs = TRUE,
                       useJSONMetadata=TRUE,
                       stringsAsFactors=FALSE)

watsonMEMAs <- data.frame(cellLine=c("HCC1954","HCC1954","AU565","AU565"),
                          ss=c("SS6"),
                          drug=c("DMSO","Lapatinib"),
                          analysisVersion="av1.6",
                          rawDataVersion="v2",
                          limitBarcodes=c(8,2,2,2),
                          k=c(7,1,1,1),
                          calcAdjacency=TRUE,
                          writeFiles = TRUE,
                          mergeOmeroIDs = TRUE,
                          useJSONMetadata=FALSE,
                          stringsAsFactors=FALSE)

ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf,watsonMEMAs)

library(XLConnect)
library(data.table)

tmp <- apply(ssDatasets[c(16),], 1, preprocessMEPLINCSL1Spot, verbose=TRUE)

