#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 11/18/2015


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
  
  #Use metadata from an2omero files or 
  #directly from GAL and xlsx files
  if(useAnnotMetadata){
    fns <- fileNames$Path[fileNames$Type=="metadata"&grepl("an2omero",fileNames$Path)]
    metadata <- rbindlist(mclapply(fns[1:limitBarcodes], processan2omero,mc.cores=detectCores()))
  } else {
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
  if(length(fileNames$Path[fileNames$Type=="Raw"])==0) stop("No raw data files found")
  cellDataFiles <- fileNames$Path[fileNames$Type=="Raw"]
  fileNames$DataBW <- paste(fileNames$Barcode,fileNames$Well,sep="_")
  
  dataBWs <- unique(fileNames$DataBW) %>%
    grep("NA",.,value=TRUE,invert=TRUE)
  
  dataBWInfo <- lapply(dataBWs,function(dataBW){
    #browser()
    dplyr::filter(fileNames,DataBW==dataBW)
  })
  
  startTime <- Sys.time()
  debugLimiter <- 64
  expDTList <- mclapply(dataBWInfo[1:debugLimiter], function(dataBWI){
    if(verbose) cat(paste("Reading and annotating data for",unique(dataBWI$DataBW)),"\n")
    nuclei <- convertColumnNames(fread(dataBWI$Path[grepl("Nuclei",dataBWI$Location)]))
    if (curatedOnly) nuclei <- nuclei[,grep(curatedCols,colnames(nuclei)), with=FALSE]
    setnames(nuclei,paste0("Nuclei_",colnames(nuclei)))
    setnames(nuclei,"Nuclei_CP_ImageNumber","Spot")
    setnames(nuclei,"Nuclei_CP_ObjectNumber","ObjectNumber")
    setkey(nuclei,Spot,ObjectNumber) 

    
    if(any(grepl("Cells",dataBWI$Location))){
      cells <- convertColumnNames(fread(dataBWI$Path[grepl("Cells",dataBWI$Location)]))
      if (curatedOnly) cells <- cells[,grep(curatedCols,colnames(cells)), with=FALSE]
      setnames(cells,paste0("Cells_",colnames(cells)))
      setnames(cells,"Cells_CP_ImageNumber","Spot")
      setnames(cells,"Cells_CP_ObjectNumber","ObjectNumber")
      setkey(cells,Spot,ObjectNumber)
    } 
    
    if(any(grepl("Cytoplasm",dataBWI$Location))){
      cytoplasm <- convertColumnNames(fread(dataBWI$Path[grepl("Cytoplasm",dataBWI$Location)]))
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
    dt <- dt[,Well := unique(dataBWI$Well)]
    dt <- dt[,Barcode := unique(dataBWI$Barcode)]
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
    #Count the cells at each spot
    dtm <- dtm[,Spot_PA_SpotCellCount := .N,by="Barcode,Well,Spot"]
    dtm <- dtm[,Spot_PA_SpotCellCountLog2 := boundedLog2(Spot_PA_SpotCellCount)]
    
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
  names(expDTList) <- dataBWs[1:debugLimiter]
  cat("converting list to datatable")
  cDT <- rbindlist(expDTList)
  Sys.time()-startTime
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
    barcodes <- unique(cDT$Barcodes)
    for(barcode in unique(cDT$Barcode)){
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
      cDT <- merge(cDT,omeroIndex,by=c("Well","ArrayRow","ArrayColumn"))
    }
  }

  
  #Set 2N and 4N DNA status
  cDT <- cDT[,Nuclei_PA_Cycle_State := gateOnlocalMinima(Nuclei_CP_Intensity_IntegratedIntensity_Dapi), by="Barcode,Well"]
  
  cDT <- cDT[,Nuclei_PA_Cycle_DNA2NProportion := calc2NProportion(Nuclei_PA_Cycle_State),by="Barcode,Well,Spot"]
  cDT <- cDT[,Nuclei_PA_Cycle_DNA4NProportion := 1-Nuclei_PA_Cycle_DNA2NProportion]
  
  #Logit transform DNA Proportions
  #logit(p) = log[p/(1-p)]
  if(any(grepl("Nuclei_PA_Cycle_DNA2NProportion",colnames(cDT)))){
    cDT <- cDT[,Nuclei_PA_Cycle_DNA2NProportionLogit := boundedLogit(Nuclei_PA_Cycle_DNA2NProportion)]
  }
  
  if(any(grepl("Nuclei_PA_Cycle_DNA4NProportion",colnames(cDT)))){
    cDT <- cDT[,Nuclei_PA_Cycle_DNA4NProportionLogit := boundedLogit(Nuclei_PA_Cycle_DNA4NProportion)]
  }
  
  
  if(!useAnnotMetadata){
    #Create short display names, then replace where not unique
    #Use entire AnnotID for ligands with same uniprot IDs
    # cDT$Ligand[grepl("NRG1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "NRG1",annotIDs = cDT$Ligand[grepl("NRG1",cDT$Ligand)])
    # cDT$Ligand[grepl("TGFB1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "TGFB1",annotIDs = cDT$Ligand[grepl("TGFB1",cDT$Ligand)])
    # cDT$Ligand[grepl("CXCL12",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "CXCL12",annotIDs = cDT$Ligand[grepl("CXCL12",cDT$Ligand)])
    # cDT$Ligand[grepl("IGF1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "IGF1",annotIDs = cDT$Ligand[grepl("IGF1",cDT$Ligand)])
  }
  
  
  #Create staining set specific derived parameters
  if (grepl("SS2|SS4|SS6|SSA|SSD|SSE|SSF",ss)){
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
    #Overide autogate for AU565 cell line in Lapatinib MEMAs
    cDT$Nuclei_PA_Gated_EdUPositive[cDT$Nuclei_CP_Intensity_MedianIntensity_EdULog2>cDT$Median3SD]<-1
    if(grepl("AU565_Lapatinib|AU565_DMSO",datasetName)){
      cDT$Nuclei_PA_Gated_EdUPositive <- 0
      cDT$Nuclei_PA_Gated_EdUPositive[cDT$Nuclei_CP_Intensity_MedianIntensity_EdULog2>10] <- 1
    }
    #Calculate the EdU Positive Percent at each spot
    cDT <- cDT[,Nuclei_PA_Gated_EdUPositiveProportion := sum(Nuclei_PA_Gated_EdUPositive)/length(ObjectNumber),by="Barcode,Well,Spot"]
    #Logit transform EdUPositiveProportion
    #logit(p) = log[p/(1-p)]
    cDT <-cDT[,Nuclei_PA_Gated_EdUPositiveProportionLogit:= boundedLogit(Nuclei_PA_Gated_EdUPositiveProportion)]
  }
  if (grepl("SS3|SS4",ss)){
    #Calculate a lineage ratio of luminal/basal or KRT19/KRT5
    cDT <- cDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT5]
    cDT <- cDT[,Cytoplasm_PA_Intensity_LineageRatioLog2 := boundedLog2(Cytoplasm_PA_Intensity_LineageRatio)]
    #Gate the KRT signals using kmeans clustering and calculate the spot proportions
    
    if(grepl("HMEC",cellLine)) {
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
    #Calculate the positive proportion at each spot
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT5PositiveProportion := sum(Cytoplasm_PA_Gated_KRT5Positive)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cytoplasm_PA_Gated_KRT5PositiveProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5PositiveProportion)]
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT19PositiveProportion := sum(Cytoplasm_PA_Gated_KRT19Positive)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cytoplasm_PA_Gated_KRT19PositiveProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT19PositiveProportion)]
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT5KRT19NegativeProportion := sum(Cytoplasm_PA_Gated_KRTClass==0)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cytoplasm_PA_Gated_KRT5KRT19NegativeProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5KRT19NegativeProportion)]
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT5PositiveKRT19NegativeProportion := sum(Cytoplasm_PA_Gated_KRTClass==1)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cytoplasm_PA_Gated_KRT5PositiveKRT19NegativeProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5PositiveKRT19NegativeProportion)]
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT5NegativeKRT19PositiveProportion := sum(Cytoplasm_PA_Gated_KRTClass==2)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cytoplasm_PA_Gated_KRT5NegativeKRT19PositiveProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5NegativeKRT19PositiveProportion)]
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT5KRT19PositiveProportion := sum(Cytoplasm_PA_Gated_KRTClass==3)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cytoplasm_PA_Gated_KRT5KRT19PositiveProportionLogit:= boundedLogit(Cytoplasm_PA_Gated_KRT5KRT19PositiveProportion)]
  }
  if (grepl("SS4",ss)){
    #Determine the class of each cell based on KRT19 and EdU class
    #0 double negative
    #1 KRT19+, EdU-
    #2 KRT19-, EdU+
    #3 KRT19+, EdU+
    cDT <- cDT[,Cells_PA_Gated_EdUKRT19Class := Cytoplasm_PA_Gated_KRT19Positive+2*Nuclei_PA_Gated_EdUPositive]
    #Calculate gating proportions for EdU and KRT19
    cDT <- cDT[,Cells_PA_Gated_EdUKRT19NegativeProportion := sum(Cells_PA_Gated_EdUKRT19Class==0)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUKRT19NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_EdUKRT19NegativeProportion)]
    cDT <- cDT[,Cells_PA_Gated_EdUPositiveKRT19NegativeProportion := sum(Cells_PA_Gated_EdUKRT19Class==2)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUPositiveKRT19NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_EdUPositiveKRT19NegativeProportion)]
    cDT <- cDT[,Cells_PA_Gated_EdUNegativeKRT19PositiveProportion := sum(Cells_PA_Gated_EdUKRT19Class==1)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUNegativeKRT19PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_EdUNegativeKRT19PositiveProportion)]
    cDT <- cDT[,Cells_PA_Gated_EdUKRT19PositiveProportion := sum(Cells_PA_Gated_EdUKRT19Class==3)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUKRT19PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_EdUKRT19PositiveProportion)]
  }
  
  if (grepl("SS6",ss)){
    #Calculate a lineage ratio of luminal/basal or KRT19/KRT14
    cDT <- cDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT14]
    cDT <- cDT[,Cytoplasm_PA_Intensity_LineageRatioLog2 := boundedLog2(Cytoplasm_PA_Intensity_LineageRatio)]
  }
  if (grepl("SSF",ss)){
    #Determine the class of each cell based on KRT5 and EdU class
    #0 double negative
    #1 KRT5+, EdU-
    #2 KRT5-, EdU+
    #3 KRT5+, EdU+
    cDT <- cDT[,Cytoplasm_PA_Gated_KRT5Positive := gateOnQuantile(Cytoplasm_CP_Intensity_MedianIntensity_KRT5Log2,probs=.02),by="Barcode,Well"]
    cDT <- cDT[,Cells_PA_Gated_EdUKRT5Class := Cytoplasm_PA_Gated_KRT5Positive+2*Nuclei_PA_Gated_EdUPositive]
    #Calculate gating proportions for EdU and KRT19
    cDT <- cDT[,Cells_PA_Gated_EdUKRT5NegativeProportion := sum(Cells_PA_Gated_EdUKRT5Class==0)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUKRT5NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_EdUKRT5NegativeProportion)]
    cDT <- cDT[,Cells_PA_Gated_EdUPositiveKRT5NegativeProportion := sum(Cells_PA_Gated_EdUKRT5Class==2)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUPositiveKRT5NegativeProportionLogit:= boundedLogit(Cells_PA_Gated_EdUPositiveKRT5NegativeProportion)]
    cDT <- cDT[,Cells_PA_Gated_EdUNegativeKRT5PositiveProportion := sum(Cells_PA_Gated_EdUKRT5Class==1)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUNegativeKRT5PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_EdUNegativeKRT5PositiveProportion)]
    cDT <- cDT[,Cells_PA_Gated_EdUKRT5PositiveProportion := sum(Cells_PA_Gated_EdUKRT5Class==3)/length(ObjectNumber),by="Barcode,Well,Spot"]
    cDT <-cDT[,Cells_PA_Gated_EdUKRT19PositiveProportionLogit:= boundedLogit(Cells_PA_Gated_EdUKRT5PositiveProportion)]
  }
  
  ##Remove nuclear objects that dont'have cell and cytoplasm data
  if(any(grepl("SS1|SS3|SS6|SSD|SSF",ss))) cDT <- cDT[!is.na(cDT$Cells_CP_AreaShape_MajorAxisLength),]
  
  # The cell level raw data and metadata is saved as Level 1 data. 
  if(writeFiles){
    #Write out cDT without normalized values as level 1 dataset
    level1Names <- grep("Norm|RUV|Loess$",colnames(cDT),value=TRUE,invert=TRUE)
    if(verbose) cat("Writing level 1 file to disk\n")
    writeTime<-Sys.time()
    set.seed(42)

    fwrite(data.table(format(cDT[sample(x = 1:nrow(cDT),size = .1*nrow(cDT),replace = FALSE),level1Names, with=FALSE], digits = 4, trim=TRUE)), paste0( "MEP_LINCS/AnnotatedData/", datasetName,"_",ss,"_","Level1.txt"), sep = "\t", quote=FALSE)

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

MLDDataSet <- data.frame(datasetName=c("MCF10ANeratinib2","MCF10ADMSO2","MCF10AVorinostat","MCF10ATrametinib"),
                         cellLine=c("MCF10A"),
                         ss=c("SSF"),
                         drug=c("Neratinib","DMSO","Vorinostat","Trametinib"),
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
                          k=c(64),
                          calcAdjacency=TRUE,
                          writeFiles = TRUE,
                          mergeOmeroIDs = TRUE,
                          useAnnotMetadata=FALSE,
                          stringsAsFactors=FALSE)



library(XLConnect)
library(data.table)

tmp <- apply(MLDDataSet[2,], 1, preprocessMEPLINCSL1Spot, verbose=TRUE)


