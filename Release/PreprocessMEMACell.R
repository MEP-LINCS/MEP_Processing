#!/bin/bash

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
# 1/16/17

simplifyLigandAnnotID <- function (ligand, annotIDs) 
{
  if (length(annotIDs)) {
    splits <- limma::strsplit2(annotIDs, split = "_")
    us <- apply(splits, 1, function(x) {
      if (x[length(x)] == "") {
        if (length(x) < 2) 
          stop("There are not enough substrings in a ligandAnnotID")
        u <- x[length(x) - 1]
      }
      else {
        u <- x[length(x)]
      }
      return(u)
    })
    ligands <- paste(ligand, us, sep = "_")
  }
  else ligands <- annotIDs
  return(ligands)
}

.libPaths("/home/users/dane/R/x86_64-redhat-linux-gnu-library/3.3")

preprocessMEMACell <- function(barcodePath, verbose=FALSE){
  barcode <- gsub(".*/","",barcodePath)
  path <- gsub(barcode,"",barcodePath)
  if (verbose) cat("Processing plate:",barcode,"at",path,"\n")
  functionStartTime<- Sys.time()
  startTime<- Sys.time()
  
  analysisVersion<-"v1.8"
  rawDataVersion<-"v1"
  limitBarcodes<- NULL
  mergeOmeroIDs<-TRUE
  calcAdjacency<-TRUE
  writeFiles<-TRUE
  useAnnotMetadata<-FALSE
  
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
    fn <- dir(paste0(barcodePath,"/Analysis"),pattern = "xml",full.names = TRUE)
    if(!length(fn)==1) stop(paste("There must be 1 xml file in the",barcode, "Analysis folder"))
    ldf <- readLogData(fn)
    #Read the well metadata from a multi-sheet Excel file
    fn <- dir(paste0(barcodePath,"/Analysis"),pattern = "xlsx",full.names = TRUE)
    if(!length(fn)==1) stop(paste("There must be 1 xlsx metadata file in the",barcode, "Analysis folder"))
    wellMetadata <- data.table(readMetadata(fn), key="Well")
    MEMA8Well <- setequal(unique(wellMetadata$Well),c("A01","A02","A03","A04","B01","B02","B03","B04"))
    MEMA96Well <- length(unique(wellMetadata$Well))==96
    #Read in the spot metadata from the gal file
    fn <- dir(paste0(barcodePath,"/Analysis"),pattern = "gal",full.names = TRUE)
    if(!length(fn)==1) stop(paste("There must be 1 gal file in the",barcode, "Analysis folder"))
    if(MEMA8Well){ #Read and prepare 8 well metadata
      smd <- readSpotMetadata(fn)
      setnames(smd, "Name", "ECMp")
      spotMetadata <- merge(smd,ldf, all=TRUE)
      setkey(spotMetadata,Spot)
      #Make a rotated version of the spot metadata for 8 well plates
      spotMetadata180 <- rotateMetadata(spotMetadata)
    } else if (MEMA96Well) {
      #Use a 12 well gal file 8 times in a 96 well plate
      smd <- readSpotMetadata(fn)
      setnames(smd, "Name", "ECMp")
      #The ArrayRow and ArrayColumn indices are oriented with A01 well in the upper left
      #These are coordinates within each well
      #The gal file coordinates are rotated 90 degrees ccw from the array coordinates
      nrArrayRow <- max(smd$Column)
      nrArrayColumn <- max(smd$Row)
      smd <- smd[,ArrayRow := (nrArrayColumn+1-Column)]
      smd <- smd[,ArrayColumn := Row]
      #Assign Spots as sequential indices within each array
      smd <- smd[,Spot := ArrayColumn+(ArrayRow-1)*nrArrayRow]
      #PrintSpots are the same as spot in these plates
      smd <- smd[,PrintSpot := Spot]
      #Copy the 12 pin spot metadata 8 times to fill the 96 well plate
      smdPlate <- rbindlist(lapply(1:(length(unique(wellMetadata$Well))/length(unique(smd$Block))),function(x){
        smd <- cbind(smd,PlatePrintIndex = x)
      }))
      ###Assign PrintHead rows and columns in order to define wells
      #PlatePrint coordinates are for the print head sectors in the 96 wellplate
      nrPlatePrintCol <- 2
      nrPlatePrintRow <- max(smdPlate$PlatePrintIndex)/nrPlatePrintCol
      smdPlate$PlatePrintRow <- ceiling(smdPlate$PlatePrintIndex/nrPlatePrintCol)
      smdPlate$PlatePrintCol <- ((smdPlate$PlatePrintIndex-1) %% nrPlatePrintCol)+1
      #PrintHead coords are the rows and columns of pins/blocks within the print head
      nrPrintHeadCol <- 2
      nrPrintHeadRow <- max(smdPlate$Block)/nrPrintHeadCol
      smdPlate$PrintHeadRow <- ceiling(smdPlate$Block/nrPrintHeadCol)
      smdPlate$PrintHeadCol <- ((smdPlate$Block-1) %% nrPrintHeadCol)+1
      #Use PlatePrint and PrintHead to determine WellIndex
      smdPlate$PlateRow <- (smdPlate$PlatePrintRow-1)*nrPrintHeadCol+(3-smdPlate$PrintHeadCol)
      smdPlate$PlateCol <- smdPlate$PlateCol <- (smdPlate$PlatePrintCol-1)*nrPrintHeadRow+smdPlate$PrintHeadRow
      nrPlateRow <- max(smdPlate$PlateRow)
      nrPlateCol <- max(smdPlate$PlateCol)
      smdPlate <- smdPlate[,WellIndex := (PlateRow-1)*nrPlateCol+PlateCol]
      smdPlate <-smdPlate[,Well:=wellAN(nrPlateRow, nrPlateCol)[smdPlate$WellIndex]]
      #Cleanup working columns used to create Well values
      smdPlate <- smdPlate[,PlatePrintRow :=NULL]
      smdPlate <- smdPlate[,PlatePrintCol :=NULL]
      smdPlate <- smdPlate[,PrintHeadRow :=NULL]
      smdPlate <- smdPlate[,PrintHeadCol :=NULL]
      smdPlate <- smdPlate[,WellIndex :=NULL]
    } else {
      stop("Only 8 well and 96 well plates are supported")
    }
  }
  #Gather filenames of raw data
  cellDataFilePaths <- dir(paste0(barcodePath,"/Analysis/",rawDataVersion), full.names = TRUE)
  if(length(cellDataFilePaths)==0) stop("No raw data files found")
  dataBWInfo <- data.table(Path=cellDataFilePaths,
                           Well=gsub("_","",str_extract(dir(paste0(barcodePath,"/Analysis/",rawDataVersion)),"_.*_")),
                           Location=str_extract(cellDataFilePaths,"Nuclei|Cytoplasm|Cells|Image"))
  startTime <- Sys.time()
  #Determine which pipeline created the data
  CPPipeline <- "Nuclei" %in% dataBWInfo$Location
  ICPipeline <- any(grepl("96well",dataBWInfo$Path))
  #Gather data from either CP or INCell
  if(CPPipeline) {
    dtL <- mclapply(unique(dataBWInfo$Well), function(well){
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
      return(dt)
    },mc.cores=detectCores())
  } else if(ICPipeline) { #read and convert INCell data to CP format
    if(verbose) cat(paste("Reading and annotating INCell data for",barcode,"\n"))
    #Read and combine the 2 header rows after the summary information
    hdrRows <- read.csv(cellDataFilePaths,skip = 18, nrows=2, header=FALSE,stringsAsFactors = FALSE)
    hdr <- sub("^_","",paste(hdrRows[1,],hdrRows[2,],sep="_"))
    #Read the cell level and spot summary data
    df <- read.csv(cellDataFilePaths,skip = 20, header=FALSE, stringsAsFactors = FALSE)
    #remove the spot summary data and convert to a data.table
    if(any(which(df$V1==""))) {
      df <- df[-(min(which(df$V1=="")):nrow(df)),]
    }
    dt <- data.table(df)
    #Name the columns
    setnames(dt,names(dt), hdr)
    if("NA_NA" %in% colnames(dt)) dt <- dt[,NA_NA := NULL]
    #Create a spot column from the field value
    dt$Spot <- gsub(".*fld ","",dt$Well)
    dt$Spot <- as.numeric(gsub(")","",dt$Spot))
    #Convert well names to alphanumeric with 2 digit columns
    wellRow <- str_match(dt$Well,"[:alpha:]")
    wells <- str_match(dt$Well,"[:digit:][:digit:]?") %>%
      as.numeric() %>%
      sprintf("%02d",.) %>%
      paste0(wellRow,.)
    
    #Convert all columns besides the Well to numeric values
    dt <- dt[,lapply(.SD, as.numeric),.SDcols = grep("Well",colnames(dt),value=TRUE,invert=TRUE)]
    dt$Well <- wells
    # dt$PlateRow <- wellRow
    # dt <- dt[,PlateCol := as.numeric(gsub("[[:alpha:]]","",dt$Well))]
    #Convert INCell names to CP versions
    dt <- convertColumnNames(dt)
    #Assume first nuclear channel is DAPI
    #Create a list of lists with IC names and corresponding CP names if available
    ICtoCPNames <- list(
      list("Nuclei_NucIntensity","Nuclei_CP_Intensity_MedianIntensity_Dapi"),
      list("Nuclei_NucArea","Nuclei_CP_AreaShape_Area"),
      list("Nuclei_NuccgX","Nuclei_CP_AreaShape_Center_X"),
      list("Nuclei_NuccgY","Nuclei_CP_AreaShape_Center_Y"),
      list("Nuclei_NucElongation","Nuclei_CP_AreaShape_Eccentricity"),
      list("Nuclei_NucCellIntensity","Cell_CP_Intensity_MedianIntensity_Dapi"),
      list("Nuclei_IxANuc","Nuclei_CP_Intensity_IntegratedIntensity_Dapi"),
      list("Cells_CellIntensity",paste0("Cytoplasm_CP_Intensity_MedianIntensity_",unique(wellMetadata$EndPoint488))),
      list("Reference1_CellIntensity",paste0("Cytoplasm_CP_Intensity_MedianIntensity_",unique(wellMetadata$EndPoint555))),
      list("Reference2_NucIntensity",paste0("Nuclei_CP_Intensity_MedianIntensity_",unique(wellMetadata$EndPoint647))),
      list("Cell","ObjectNumber")
    )
    #Change names within dt to match downstream CP pipeline
    foo <- sapply(ICtoCPNames,function(x) {
      if(x[[1]] %in% colnames(dt)) setnames(dt,x[[1]],x[[2]])
    })
    rm(foo)
    
    dt$Barcode <- barcode
    
    dtL <-list(dt)
  } else {
    stop("Only CP and IC pipelines are supported")
  }
  
  expDTList <- mclapply(dtL, function(dt){
    #Remove problematic features
    dt <- dt[,grep("Euler",colnames(dt),invert=TRUE), with=FALSE]
    if(useAnnotMetadata) {
      dtm <- merge(dt,metadata,by = c("Barcode","Well","Spot"))
      dtm$PrintSpot <- dtm$Spot
      
    } else {
      #For 8 well plates, merge the data with its metadata based on the row it's in
      
      if(MEMA8Well) {
        m <- regexpr("[[:alpha:]]",well)
        row <- regmatches(well,m)
        setkey(dt,Spot)
        dtm <- switch(row, A = merge(dt,spotMetadata,all.x=TRUE),
                      B = merge(dt,spotMetadata180,all.x=TRUE))
        #Add a WellSpace spot index that recognizes the arrays are rotated 180 degrees
        nrArrayRows <- max(dtm$ArrayRow)
        nrArrayColumns <- max(dtm$ArrayColumn)
        dtm$PrintSpot <- dtm$Spot
        dtm$PrintSpot[grepl("B", dtm$Well)] <- (nrArrayRows*nrArrayColumns+1)-dtm$PrintSpot[grepl("B", dtm$Well)]
        #merge well metadata with the data and spot metadata
        dtm <- merge(dtm,wellMetadata,by = "Well")
        #Add a WellSpace spot index that recognizes the arrays are rotated 180 degrees
        nrArrayRows <- max(dtm$ArrayRow)
        nrArrayColumns <- max(dtm$ArrayColumn)
        dtm$PrintSpot <- dtm$Spot
        dtm$PrintSpot[grepl("B", dtm$Well)] <- (nrArrayRows*nrArrayColumns+1)-dtm$PrintSpot[grepl("B", dtm$Well)]
        
      } else if(MEMA96Well) {
        #merge the data with the spot metadata for the plate
        dtm <- merge(dt,smdPlate,by = c("Well","Spot"))
      }
      dtm <- merge(wellMetadata,dtm,by="Well")
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
    #Scale CP pipelin values
    if(CPPipeline) {
      intensityNames <- grep("Intensity",colnames(dtm), value=TRUE)
      scaledInts <- dtm[,intensityNames, with=FALSE]*2^16
      dtm <- cbind(dtm[,!intensityNames, with=FALSE],scaledInts)
    }
    #log transform all intensity and areaShape values
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
    dtm <- dtm[,MEP_Drug:=paste(MEP,Drug,sep = "_")]
    
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
    return(dtm)
  }, mc.cores=detectCores())
  
  #Add names to the data.tables in the list
  #names(expDTList) <- paste(barcode,unique(dataBWInfo$Well),sep="_")
  cDT <- rbindlist(expDTList)
  rm(expDTList,dtL)
  gc()
  
  if(mergeOmeroIDs){
    #Read in and merge the Omero URLs
    omeroIndex <- fread(paste0(barcodePath,"/Analysis/",barcode,"_imageIDs.tsv"))[,list(WellName,Row,Column,ImageID)]
    if(MEMA8Well){
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
    } else if(MEMA96Well){
      #Convert well names to alphanumeric with 2 digit columns
      wellRow <- str_match(omeroIndex$WellName,"[:alpha:]-") %>%
        str_replace("-","")
      omeroIndex$Well <- str_match(omeroIndex$WellName,"-[:digit:]*") %>%
        str_replace("-","") %>%
        as.numeric() %>%
        sprintf("%02d",.) %>%
        paste0(wellRow,.)
      setnames(omeroIndex,"Row","ArrayRow")
      setnames(omeroIndex,"Column","ArrayColumn")
      omeroIndex <- omeroIndex[,WellName:=NULL]
      cDT <- merge(cDT,omeroIndex,by=c("Well","ArrayRow","ArrayColumn"))
    } else {
      stop("Only 8 and 96 well plates can be merged with omero IDs")
    }
  }
  
  if (verbose) cat("Gating cell level data for plate",barcode,"\n")  
  #Set 2N and 4N DNA status
  cDT <- cDT[,Nuclei_PA_Cycle_State := gateOnlocalMinima(Nuclei_CP_Intensity_IntegratedIntensity_Dapi)]
  
  if(!useAnnotMetadata){
    #Create short display names, then replace where not unique
    #Use entire AnnotID for ligands with same uniprot IDs
    cDT$Ligand[grepl("NRG1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "NRG1",annotIDs = cDT$Ligand[grepl("NRG1",cDT$Ligand)])
    cDT$Ligand[grepl("TGFB1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "TGFB1",annotIDs = cDT$Ligand[grepl("TGFB1",cDT$Ligand)])
    cDT$Ligand[grepl("CXCL12",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "CXCL12",annotIDs = cDT$Ligand[grepl("CXCL12",cDT$Ligand)])
    cDT$Ligand[grepl("IGF1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "IGF1",annotIDs = cDT$Ligand[grepl("IGF1",cDT$Ligand)])
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
    # formatTime<-Sys.time()
    # if(verbose) cat("Formatting",barcode,"level 1 data\n")
    # dt <- data.table(format(cDT, digits = 4, trim=TRUE))
    # if(verbose) cat("Format time for ",barcode,":", Sys.time()-writeTime,"\n")
    if(verbose) cat("Writing",barcode,"level 1 full file to disk\n")
    writeTime<-Sys.time()
    fwrite(cDT, paste0(barcodePath, "/Analysis/", barcode,"_","Level1.tsv"), sep = "\t", quote=FALSE)
    cat("Write time:", Sys.time()-writeTime,"\n")
    if(verbose) cat("Writing",barcode,"level 1 subset file to disk\n")
    writeTime<-Sys.time()
    set.seed(42)
    fwrite(cDT[sample(x=1:nrow(cDT),size = .1*nrow(cDT),replace=FALSE),], paste0(barcodePath, "/Analysis/", barcode,"_","Level1Subset.tsv"), sep = "\t", quote=FALSE)
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
    
    #Write the File Annotations for Synapse to tab-delimited file
    write.table(c(
      CellLine = unique(cDT$CellLine),
      Preprocess = analysisVersion,
      DataType = "Quanititative Imaging",
      Consortia = "MEP-LINCS",
      Drug = unique(cDT$Drug),
      Segmentation = rawDataVersion,
      StainingSet = gsub("Layout.*","",unique(cDT$StainingSet)),
      Level = "1"
    ),
    paste0(barcodePath, "/Analysis/", barcode,"_","Level1Annotations.tsv"), sep = "\t",col.names = FALSE, quote=FALSE)
  }
  cat("Elapsed time:", Sys.time()-functionStartTime, "\n")
}

barcodePath <-commandArgs(trailingOnly = TRUE)
#barcodePath <- "/lincs/share/lincs_user/LI8X00771"
#barcodePath <- "/lincs/share/lincs_user/lincs96well/LI9V01612"
res <- preprocessMEMACell(barcodePath, verbose=TRUE)

