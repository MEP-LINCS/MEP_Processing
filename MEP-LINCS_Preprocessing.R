
#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 11/18/2015

##Introduction

#   The MEP-LINCs dataset contains imaging data from a Nikon automated microscope that is analyzed with a CellProfiler pipeline.
# 
# Part of this preprocessing of the dataset will be deprecated when the merging of the data and metadata happens within the CellProfiler part of the pipeline. For now, the metadata about the ECM proteins is read from the GAL file and the metadata about the wells (cell line, stains and ligands) is read from Excel spreadsheets.

library("parallel")#use multiple cores for faster processing

#source("MEPLINCSFunctions.R")

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
    return(mdDT)
  }))
}


preprocessMEPLINCS <- function(ssDataset, verbose=FALSE){
  ss<-ssDataset[["ss"]]
  cellLine<-ssDataset[["cellLine"]]
  k<-ssDataset[["k"]]
  analysisVersion<-ssDataset[["analysisVersion"]]
  rawDataVersion<-ssDataset[["rawDataVersion"]]
  limitBarcodes<-ssDataset[["limitBarcodes"]]
  mergeOmeroIDs<-as.logical(ssDataset[["mergeOmeroIDs"]])
  calcAdjacency<-as.logical(ssDataset[["calcAdjacency"]])
  writeFiles<-as.logical(ssDataset[["writeFiles"]])
  useJSONMetadata<-as.logical(ssDataset[["useJSONMetadata"]])
  
  seNames=c("DNA2N","SpotCellCount","Edu","MitoTracker","KRT","Lineage","Fibrillarin")
  
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
    fns <- dir(paste0("./",cellLine,"/",ss,"/Metadata"),pattern = ".json",full.names = TRUE)
    metadata <- rbindlist(mclapply(fns, processJSON, mc.cores = detectCores()))
  } else {
    # Read and clean spotmetadata
    
    #Read in the spot metadata from the gal file
    smd <- readSpotMetadata(dir(paste0("./",cellLine,"/",ss,"/Metadata/"),pattern = ".gal",full.names = TRUE))
    #Relabel the column Name to ECMpAnnotID
    setnames(smd, "Name", "ECMpAnnotID")
    
    #Add the print order and deposition number to the metadata
    ldf <- readLogData(dir(paste0("./",cellLine,"/",ss,"/Metadata"),pattern = ".xml",full.names = TRUE))
    spotMetadata <- merge(smd,ldf, all=TRUE)
    setkey(spotMetadata,Spot)
    #Make a rotated version of the spot metadata to match the print orientation
    spotMetadata180 <- rotateMetadata(spotMetadata)
    ARowMetadata <- data.table(spotMetadata,Well=rep(c("A01", "A02","A03","A04"),each=nrow(spotMetadata)))
    BRowMetadata <- data.table(spotMetadata180,Well=rep(c("B01", "B02","B03","B04"),each=nrow(spotMetadata180)))
  }
  
  # The raw data from all wells in all plates in the dataset are read in and merged with their spot and well metadata. The number of nuclei at each spot are counted and a loess model of the spot cell count is added. Then all intensity values are normalized.
  # 
  # Next, the data is filtered to remove objects with a nuclear area less than nuclearAreaThresh pixels or more than nuclearAreaHiThresh pixels.
  
  #merge_normalize_QA, echo=FALSE}
  #The next steps are to bring in the well metadata, the print order and the CP data
  
  cellDataFiles <- dir(paste0("./",cellLine,"/", ss,"/RawData/",rawDataVersion),full.names = TRUE)
  splits <- strsplit2(strsplit2(cellDataFiles,split = "_")[,1],"/")
  
  barcodes <- unique(splits[,ncol(splits)])[1:limitBarcodes]
  #if(rawDataVersion=="v1.1") barcodes <- gsub("reDAPI","",barcodes)
  if(verbose) cat(paste("Reading and annotating cell level data for",cellLine,ss)," \n")
  expDTList <- lapply(barcodes, function(barcode){
    plateDataFiles <- grep(barcode,cellDataFiles,value = TRUE)
    wells <- unique(strsplit2(split = "_",plateDataFiles)[,2])
    wellDataList <- lapply(wells,function(well){
      
      wellDataFiles <- grep(well,plateDataFiles,value = TRUE)
      nucleiDataFile <- grep("Nuclei",wellDataFiles,value=TRUE,
                             ignore.case = TRUE)
      if (ss %in% c("SS1","SS3")){
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
      
      setkey(nuclei,CP_ImageNumber,CP_ObjectNumber)
      if (ss %in% c("SS1","SS3")){
        cells <- convertColumnNames(fread(cellsDataFile))
        if (curatedOnly) cells <- cells[,grep(curatedCols,colnames(cells)), with=FALSE]
        cells <- cells[,grep("Euler",colnames(cells),invert=TRUE), with=FALSE]
        setkey(cells,CP_ImageNumber,CP_ObjectNumber)
        cytoplasm <- convertColumnNames(fread(cytoplasmDataFile))
        if (curatedOnly) cytoplasm <- cytoplasm[,grep(curatedCols,colnames(cytoplasm)), with=FALSE]
        cytoplasm <- cytoplasm[,grep("Euler",colnames(cytoplasm),invert=TRUE), with=FALSE]
        setkey(cytoplasm,CP_ImageNumber,CP_ObjectNumber)
      }
      
      #Add the data location as a prefix in the column names
      setnames(nuclei,paste0("Nuclei_",colnames(nuclei)))
      if (ss %in% c("SS1","SS3")){
        setnames(cells,paste0("Cells_",colnames(cells)))
        setnames(cytoplasm,paste0("Cytoplasm_",colnames(cytoplasm)))
      }
      #Merge the cells, cytoplasm and nuclei data
      if (ss %in% c("SS1","SS3")){
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
    rm(wellDataList)
    
    if(!useJSONMetadata){
      #Read the well metadata from a multi-sheet Excel file
      wellMetadata <- data.table(readMetadata(paste0("./",cellLine,"/",ss,"/Metadata/",gsub("reDAPI","",barcode),".xlsx")), key="Well")
      
      #merge well metadata with the data and spot metadata
      pcDT <- merge(pcDT,wellMetadata,by = "Well")
    }
    
    
    if(mergeOmeroIDs){
      if(any(grepl("v1.1|v2.1",rawDataVersion))){
        imageURLFiles <- grep("reDAPI_imageIDs",dir(paste0("./",cellLine,"/", ss,"/Metadata/"),full.names = TRUE), value=TRUE)
      } else {
        allImageFiles <- grep("imageIDs",dir(paste0("./",cellLine,"/", ss,"/Metadata/"),full.names = TRUE), value=TRUE)
        imageURLFiles <- grep("reDAPI", allImageFiles, invert=TRUE, value=TRUE)
      }
      
      #Read in and merge the Omero URLs
      omeroIndex <- fread(grep(barcode, imageURLFiles, value=TRUE))[,list(WellName,Row,Column,ImageID)]
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
    #Count the cells at each spot
    pcDT <- pcDT[,Spot_PA_SpotCellCount := .N,by="Barcode,Well,Spot"]
    pcDT <- pcDT[,Spot_PA_SpotCellCountLog2 := log2(Spot_PA_SpotCellCount)]
    #log transform all intensity values
    #pcDT <- logIntensities(pcDT)
    intensityNames <- grep("Intensity",colnames(pcDT), value=TRUE, ignore.case = TRUE)
    
    
    dtLog <- pcDT[,lapply(.SD,boundedLog2),.SDcols=intensityNames]
    setnames(dtLog,colnames(dtLog),paste0(colnames(dtLog),"Log2"))
    pcDT <- cbind(pcDT,dtLog)
    
    if(any(grepl("Nuclei_CP_AreaShape_Center",colnames(pcDT)))){
      #Add the local polar coordinates and Neighbor Count
      pcDT <- pcDT[,Nuclei_PA_Centered_X :=  Nuclei_CP_AreaShape_Center_X-median(Nuclei_CP_AreaShape_Center_X)]
      pcDT <- pcDT[,Nuclei_PA_Centered_Y :=  Nuclei_CP_AreaShape_Center_Y-median(Nuclei_CP_AreaShape_Center_Y)]
      pcDT <- pcDT[, Nuclei_PA_AreaShape_Center_R := sqrt(Nuclei_PA_Centered_X^2 + Nuclei_PA_Centered_Y^2)]
      pcDT <- pcDT[, Nuclei_PA_AreaShape_Center_Theta := calcTheta(Nuclei_PA_Centered_X, Nuclei_PA_Centered_Y)]
    }
    #Set 2N and 4N DNA status
    #pcDT <- pcDT[,Nuclei_PA_Cycle_State := kmeansDNACluster(Nuclei_CP_Intensity_IntegratedIntensity_Dapi), by="Barcode,Well"]
    pcDT <- pcDT[,Nuclei_PA_Cycle_State := gateOnlocalMinima(Nuclei_CP_Intensity_IntegratedIntensity_Dapi), by="Barcode,Well"]
    
    pcDT <- pcDT[,Nuclei_PA_Cycle_DNA2NProportion := calc2NProportion(Nuclei_PA_Cycle_State),by="Barcode,Well,Spot"]
    pcDT$Nuclei_PA_Cycle_DNA4NProportion <- 1-pcDT$Nuclei_PA_Cycle_DNA2NProportion
    
    #Logit transform DNA Proportions
    #logit(p) = log[p/(1-p)]
    if(any(grepl("Nuclei_PA_Cycle_DNA2NProportion",colnames(pcDT)))){
      DNA2NImpute <- pcDT$Nuclei_PA_Cycle_DNA2NProportion
      DNA2NImpute[DNA2NImpute==0] <- .01
      DNA2NImpute[DNA2NImpute==1] <- .99
      pcDT$Nuclei_PA_Cycle_DNA2NProportionLogit <- log2(DNA2NImpute/(1-DNA2NImpute))
    }
    
    if(any(grepl("Nuclei_PA_Cycle_DNA4NProportion",colnames(pcDT)))){
      DNA4NImpute <- pcDT$Nuclei_PA_Cycle_DNA4NProportion
      DNA4NImpute[DNA4NImpute==0] <- .01
      DNA4NImpute[DNA4NImpute==1] <- .99
      pcDT$Nuclei_PA_Cycle_DNA4NProportionLogit <- log2(DNA4NImpute/(1-DNA4NImpute))
    }
    
    #logit transform eccentricity
    if(any(grepl("Nuclei_CP_AreaShape_Eccentricity",colnames(pcDT)))){
      EccImpute <- pcDT$Nuclei_CP_AreaShape_Eccentricity
      EccImpute[EccImpute==0] <- .01
      EccImpute[EccImpute==1] <- .99
      pcDT$Nuclei_PA_AreaShape_EccentricityLogit <- log2(EccImpute/(1-EccImpute))
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
      pcDT$ECMp <- gsub("_.*","",pcDT$ECMpAnnotID)
      pcDT$Ligand <- gsub("_.*","",pcDT$LigandAnnotID)
      
      #Use entire AnnotID for ligands with same uniprot IDs
      pcDT$Ligand[grepl("NRG1",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "NRG1",annotIDs = pcDT$LigandAnnotID[grepl("NRG1",pcDT$Ligand)])
      pcDT$Ligand[grepl("TGFB1",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "TGFB1",annotIDs = pcDT$LigandAnnotID[grepl("TGFB1",pcDT$Ligand)])
      pcDT$Ligand[grepl("CXCL12",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "CXCL12",annotIDs = pcDT$LigandAnnotID[grepl("CXCL12",pcDT$Ligand)])
      pcDT$Ligand[grepl("IGF1",pcDT$Ligand)] <- simplifyLigandAnnotID(ligand = "IGF1",annotIDs = pcDT$LigandAnnotID[grepl("IGF1",pcDT$Ligand)])
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
    } else if (grepl("SS2",ss)){
      pcDT <- pcDT[,Nuclei_PA_Gated_EduPositive := kmeansCluster(.SD,value =  "Nuclei_CP_Intensity_MedianIntensity_Edu",ctrlLigand = "FBS"), by="Barcode"]
      #pcDT <- pcDT[,Nuclei_PA_Gated_EduPositive := gateOnlocalMinima(Nuclei_CP_Intensity_MedianIntensity_Edu)-1, by="Barcode,Well"]
      #Calculate the EdU Positive Percent at each spot
      pcDT <- pcDT[,Nuclei_PA_Gated_EduPositiveProportion := sum(Nuclei_PA_Gated_EduPositive)/length(ObjectNumber),by="Barcode,Well,Spot"]
      #Logit transform EduPositiveProportion
      #logit(p) = log[p/(1-p)]
      EdUppImpute <- pcDT$Nuclei_PA_Gated_EduPositiveProportion
      EdUppImpute[EdUppImpute==0] <- .01
      EdUppImpute[EdUppImpute==1] <- .99
      pcDT$Nuclei_PA_Gated_EduPositiveLogit <- log2(EdUppImpute/(1-EdUppImpute))
      
    } else if (grepl("SS3",ss)){
      #Calculate a lineage ratio of luminal/basal or KRT19/KRT5
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT5]
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatioLog2 := log2(Cytoplasm_PA_Intensity_LineageRatio)]
      
    } else stop("Invalid ss parameter")
    return(pcDT)
    
  #}, mc.cores=detectCores())#Use all cores in production
  })#Revert to apply when debugging
  
  cDT <- rbindlist(expDTList, fill = TRUE)
  rm(expDTList)
  
  #Change to mclapply
  if(calcAdjacency){
    #cDTList <- mclapply(barcodes, function(barcode, dt){
    if(verbose){
      cat("Calculating adjacency data\n")
      #save(cDT, file="cDT.RData")
    } 
    
    densityRadius <- sqrt(median(cDT$Nuclei_CP_AreaShape_Area)/pi)
    
    #Count the number of neighboring cells
    cDT <- cDT[,Nuclei_PA_AreaShape_Neighbors := cellNeighbors(.SD, radius = densityRadius*neighborhoodNucleiRadii), by = "Barcode,Well,Spot"]
    
    #Rules for classifying perimeter cells
    cDT <- cDT[,Spot_PA_Sparse := Nuclei_PA_AreaShape_Neighbors < neighborsThresh]
    
    #Add a local wedge ID to each cell based on conversations with Michel Nederlof
    cDT <- cDT[,Spot_PA_Wedge:=ceiling(Nuclei_PA_AreaShape_Center_Theta/wedgeAngs)]
    
    #Define the perimeter cell if it exists in each wedge
    #Classify cells as outer if they have a radial position greater than a thresh
    cDT <- cDT[,Spot_PA_OuterCell := labelOuterCells(Nuclei_PA_AreaShape_Center_R, thresh=outerThresh),by="Barcode,Well,Spot"]
    
    #Require a perimeter cell not be in a sparse region
    denseOuterDT <- cDT[!cDT$Spot_PA_Sparse  & cDT$Spot_PA_OuterCell]
    denseOuterDT <- denseOuterDT[,Spot_PA_Perimeter := findPerimeterCell(.SD) ,by="Barcode,Well,Spot,Spot_PA_Wedge"]
    setkey(cDT,Barcode,Well,Spot,ObjectNumber)
    setkey(denseOuterDT,Barcode,Well,Spot,ObjectNumber)
    cDT <- denseOuterDT[,list(Barcode,Well,Spot,ObjectNumber,Spot_PA_Perimeter)][cDT]
    cDT$Spot_PA_Perimeter[is.na(cDT$Spot_PA_Perimeter)] <- FALSE
    
  }
  
  # After merging the metadata with the cell-level data, several types of derived parameters are added. These include:
  #   
  #   The origin of coordinate system is placed at the median X and Y of each spot and the local cartesian and polar coordinates are added to the dataset.
  # 
  # The number of nuclei within three nuclear radii of `r densityRadius ` around each nuclei is counted and stored as a neighbor count parameter. The neighbor count value is thresholded to classify each cell as Sparse or not.The distance from the local origin is used to classify each cell as an OuterCell or not. The Sparse, OutCell and Wedge classifications are used to classify each cell as a Perimeter cell or not. 
  # 
  # For staining set 2, each cell is classified as EdU+ or EdU-. The threshold for EdU+ is based on kmeans threshold of the mean EdU intensity from the control well of each plate.
  # 
  # The intensity values are normalized at each spot so that spot-level variations can be analyzed.
  # 
  # The cell level raw data and metadata is saved as Level 1 data. naiveReplicateRUV does not create a level 2 dataset.
  
  
  #The cell-level data is median summarized to the spot level and then normalized. The spot level data and metadata are saved as Level 3 data.
  ##Remove nuclear objects that dont'have cell and cytoplasm data
  if(verbose) cat("Creating level 3 data\n")
  
  if(any(grepl("SS1|SS3",ss))) cDT <- cDT[!is.na(cDT$Cells_CP_AreaShape_MajorAxisLength),]
  #### Level3 ####
  slDT <- createl3(cDT, lthresh,seNames = seNames)
  
  if(writeFiles){
    #Write out cDT without normalized values as level 1 dataset
    level1Names <- grep("Norm|RUV3|Loess$",colnames(cDT),value=TRUE,invert=TRUE)
    if(verbose) cat("Writing level 1 file to disk\n")
    fwrite(cDT[,level1Names, with=FALSE], file.path = paste0("./",cellLine,"/", ss, "/AnnotatedData/", unique(cDT$CellLine),"_",ss,"_",rawDataVersion,"_",analysisVersion,"_","Level1.txt"),sep="\t", verbose=TRUE)
    #write.table(format(cDT[,level1Names, with=FALSE], digits=4, trim=TRUE), paste0("./",cellLine,"/", ss, "/AnnotatedData/", unique(cDT$CellLine),"_",ss,"_",rawDataVersion,"_",analysisVersion,"_","Level1.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
    
    normParmameterNames <- grep("Norm|RUV3|Loess$",colnames(cDT), value=TRUE)
    rawParameterNames <- gsub("_?[[:alnum:]]*?Norm$|_?[[:alnum:]]*?RUV3|_?[[:alnum:]]*?Loess$", "", normParmameterNames)
    metadataNormNames <- colnames(cDT)[!colnames(cDT) %in% rawParameterNames]
    #Paste back in the QA and selected raw data
    
    level2Names <- c(metadataNormNames,
                     grep("Nuclei_CP_Intensity_MedianIntensity_Dapi$|Cytoplasm_CP_Intensity_MedianIntensity_Actin$|Cytoplasm_CP_Intensity_MedianIntensity_CellMask$|Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker$|Nuclei_CP_Intensity_MedianIntensity_H3$|Nuclei_CP_Intensity_MedianIntensity_Fibrillarin$|Nuclei_CP_Intensity_MedianIntensity_Edu$|Cytoplasm_CP_Intensity_MedianIntensity_KRT5$|Cytoplasm_CP_Intensity_MedianIntensity_KRT19$|Spot_PA_SpotCellCount$", colnames(cDT), value = TRUE))
    
    #Write out cDT with normalized values as level 2 dataset
    #if(verbose) cat("Writing level 2 file to disk\n")
   # write.table(format(cDT[,level2Names, with = FALSE], digits=4, trim=TRUE), paste0("./",cellLine,"/", ss, "/AnnotatedData/", unique(cDT$CellLine),"_",ss,"_",rawDataVersion,"_",analysisVersion,"_","Level2.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
    rm(cDT)
  }  
  
  #save(slDT,file="slDT.RData")
  slDT <- slDT[!grepl("fiducial|Fiducial|gelatin|blank|air|PBS",slDT$ECMp),]
  
  metadataNames <- "ObjectNumber|^Row$|^Column$|Block|^ID$|PrintOrder|Depositions|CellLine|Endpoint|WellIndex|Center|ECMpAnnotID|LigandAnnotID|ECMpPK|LigandPK|MEP|Well_Ligand|ImageID|Sparse|Wedge|OuterCell|Spot_PA_Perimeter|Nuclei_PA_Cycle_State|_SE|ReplicateCount|SCC|QAScore"
  
  #Save the un-normalized parameters to merge in later
  mdDT <- slDT[,grep(paste(metadataNames,"Barcode|Well|^Spot$|ArrayRow|ArrayColumn|^ECMp$|^Ligand$",sep="|"),colnames(slDT),value=TRUE), with = FALSE]
  #Identify parameters to be normalized
  signalsWithMetadata <- grep(metadataNames,colnames(slDT),value=TRUE,invert=TRUE)
  #Normalize each feature, pass with location and content metadata
  if(verbose) {
    cat("Normalizing\n")
    #save(slDT, file="slDT.RData")
  }
  nDT <- normRUV3LoessResiduals(slDT[,signalsWithMetadata, with = FALSE], k)
  nDT$NormMethod <- "RUV3LoessResiduals"
  #Merge the normalized data with its metadata
  setkey(nDT,Barcode,Well,Spot,ArrayRow,ArrayColumn,ECMp,Ligand,MEP)
  setkey(mdDT,Barcode,Well,Spot,ArrayRow,ArrayColumn,ECMp,Ligand,MEP)
  nmdDT <- merge(nDT,mdDT)
  #merge spot level normalized and raw data
  setkey(slDT, Barcode, Well, Spot,ArrayRow,ArrayColumn,ECMp,Ligand)
  slDT <- merge(slDT[,signalsWithMetadata, with = FALSE], nmdDT)
  
  #Label FBS with their plate index to keep separate
  slDT$Ligand[grepl("FBS",slDT$Ligand)] <- paste0(slDT$Ligand[grepl("FBS",slDT$Ligand)],"_P",match(slDT$Barcode[grepl("FBS",slDT$Ligand)], barcodes))
  #Add QAScore and Spot_PA_LoessSCC to cell level data
  #setkey(cDT,Barcode, Well, Spot)
  
  #cDT <- cDT[slDT[,list(Barcode, Well, Spot, QAScore, Spot_PA_LoessSCC)]]
  #The spot level data is median summarized to the replicate level and is stored as Level 4 data and metadata.
  
  #Level4Data
  if(verbose) {
    cat("Creating level 4 data\n")
    #save(slDT,file="slDT.RData")
  }
  mepDT <- createl4(slDT,seNames = seNames)
  
  #Write QA flags into appropriate data levels
  #Low cell count spots
  #cDT$QA_LowSpotCellCount <- cDT$Spot_PA_SpotCellCount < lowSpotCellCountThreshold
  slDT$QA_LowSpotCellCount <- slDT$Spot_PA_SpotCellCount < lowSpotCellCountThreshold
  
  #Low quality DAPI
  #cDT$QA_LowDAPIQuality <- FALSE
  slDT$QA_LowDAPIQuality <- FALSE
  
  #Flag spots below automatically loess QA threshold
  #cDT$QA_LowRegionCellCount <- cDT$Spot_PA_LoessSCC < lowRegionCellCountThreshold
  slDT$QA_LowRegionCellCount <- slDT$Spot_PA_LoessSCC < lowRegionCellCountThreshold
  
  #Flag wells below automatically calculated QA threshold
  slDT$QA_LowWellQA <- FALSE
  slDT$QA_LowWellQA[slDT$QAScore < lowWellQAThreshold] <- TRUE
  #cDT$QA_LowWellQA <- FALSE
  #cDT$QA_LowWellQA[cDT$QAScore < lowWellQAThreshold] <- TRUE
  
  #Level 4
  mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < lowReplicateCount
  #WriteData
  if(writeFiles){
    if(verbose) cat("Writing level 3 file to disk\n")
    write.table(format(slDT, digits = 4, trim=TRUE), paste0("./",cellLine,"/", ss, "/AnnotatedData/", unique(slDT$CellLine),"_",ss,"_",rawDataVersion, "_",analysisVersion,"_","Level3.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
    
    if(verbose) cat("Writing level 4 file to disk\n")
    write.table(format(mepDT, digits = 4, trim=TRUE), paste0("./",cellLine,"/",ss, "/AnnotatedData/", unique(slDT$CellLine),"_",ss,"_",rawDataVersion,"_",analysisVersion,"_","Level4.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
    
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
    paste0("./",cellLine,"/",ss, "/AnnotatedData/", cellLine,"_",ss,"_",analysisVersion,"_","PipelineParameters.txt"), sep = "\t",col.names = FALSE, quote=FALSE)
  }
}

PC3df <- data.frame(cellLine=rep(c("PC3"), 4),
                    ss=c("SS1", "SS2","SS3","SS2noH3"),
                    analysisVersion="av1.5",
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
                     analysisVersion="av1.5",
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
                     analysisVersion="av1.5",
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
                       analysisVersion="av1.5",
                       rawDataVersion="v2",
                       limitBarcodes=c(8,8,5),
                       k=c(7,7,4),
                       calcAdjacency=FALSE,
                       writeFiles = TRUE,
                       mergeOmeroIDs = TRUE,
                       useJSONMetadata=TRUE,
                       stringsAsFactors=FALSE)

ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf)

tmp <- apply(ssDatasets[c(12),], 1, preprocessMEPLINCS, verbose=TRUE)
