#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
# 2/1/2017
######
# #' Merge in Omero IDs to link images to t
# if(mergeOmeroIDs){
#   #Read in and merge the Omero URLs
#   omeroIndex <- fread(paste0(barcodePath,"/Analysis/",barcode,"_imageIDs.tsv"))[,list(WellName,Row,Column,ImageID)]
#   if(MEMA8Well){
#     m <- regexpr("Well[[:digit:]]",omeroIndex$WellName)
#     wellNames <- regmatches(omeroIndex$WellName,m)
#     omeroIndex$Well <- sapply(gsub("Well","",wellNames,""),FUN=switch,
#                               "1"="A01",
#                               "2"="A02",
#                               "3"="A03",
#                               "4"="A04",
#                               "5"="B01",
#                               "6"="B02",
#                               "7"="B03",
#                               "8"="B04")
#     setnames(omeroIndex,"Row","ArrayRow")
#     setnames(omeroIndex,"Column","ArrayColumn")
#     omeroIndex <- omeroIndex[,WellName:=NULL]
#     cDT <- merge(cDT,omeroIndex,by=c("Well","ArrayRow","ArrayColumn"))
#   } else if(MEMA96Well){
#     #Convert well names to alphanumeric with 2 digit columns
#     wellRow <- str_match(omeroIndex$WellName,"[:alpha:]-") %>%
#       str_replace("-","")
#     omeroIndex$Well <- str_match(omeroIndex$WellName,"-[:digit:]*") %>%
#       str_replace("-","") %>%
#       as.numeric() %>%
#       sprintf("%02d",.) %>%
#       paste0(wellRow,.)
#     setnames(omeroIndex,"Row","ArrayRow")
#     setnames(omeroIndex,"Column","ArrayColumn")
#     omeroIndex <- omeroIndex[,WellName:=NULL]
#     cDT <- merge(cDT,omeroIndex,by=c("Well","ArrayRow","ArrayColumn"))
#   } else {
#     stop("Only 8 and 96 well plates can be merged with omero IDs")
#   }
# }

######
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
  
  #Set a threshold for the loess well level QA Scores
  lthresh <- 0.6
  
  #Read in the plate's cell level data and annotations
  cDT <- fread(paste0(barcodePath,"/Analysis/",barcode,"_Level1.tsv"))
  annotations <- fread(paste0(barcodePath,"/Analysis/",barcode,"_Level1Annotations.tsv"),header = FALSE)
  
  #Count the cells at each spot at the cell level as needed by createl3
  cDT <- cDT[,Spot_PA_SpotCellCount := .N,by="Barcode,Well,Spot"]
  cDT <- cDT[,Spot_PA_SpotCellCountLog2 := boundedLog2(Spot_PA_SpotCellCount)]
  
  #Calculate DNA proportions based on cell cycle state
  cDT <- cDT[,Nuclei_PA_Cycle_DNA2NProportion := calc2NProportion(Nuclei_PA_Cycle_State),by="Barcode,Well,Spot"]
  cDT <- cDT[,Nuclei_PA_Cycle_DNA2NProportionLogit := boundedLogit(Nuclei_PA_Cycle_DNA2NProportion)]
  cDT <- cDT[,Nuclei_PA_Cycle_DNA4NProportion := 1-Nuclei_PA_Cycle_DNA2NProportion]
  cDT <- cDT[,Nuclei_PA_Cycle_DNA4NProportionLogit := boundedLogit(Nuclei_PA_Cycle_DNA4NProportion)]
  
  #If present, Calculate proportions for gated signals
  gatedSignals <- grep("Positive|High",colnames(cDT), value=TRUE)
  if(length(gatedSignals)>0){
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
  }
  
  
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
  
  #median summarize the rest of the signals
  signalDT <- summarizeToSpot(cDT,lthresh, seNames)
  if(exists("proportionsDT")) {
    spotDT <- merge(signalDT,proportionsDT)
  } else {
    spotDT <- signalDT
  }
  
  #Summarize 
  # The cell level raw data and metadata is saved as Level 1 data. 
  if(writeFiles){
    #Write out cDT without normalized values as level 1 dataset
    if(verbose) cat("Writing spot level data to disk\n")
    writeTime<-Sys.time()
    fwrite(data.table(format(spotDT, digits = 4, trim=TRUE)), paste0(barcodePath, "/Analysis/", barcode,"_","SpotLevel.tsv"), sep = "\t", quote=FALSE)
    cat("Write time:", Sys.time()-writeTime,"\n")
    
    #Write the File Annotations for Synapse to tab-delimited file
    write.table(c(
      CellLine = annotations$V2[annotations$V1=="CellLine"],
      Preprocess = annotations$V2[annotations$V1=="Preprocess"],
      DataType = annotations$V2[annotations$V1=="DataType"],
      Consortia = annotations$V2[annotations$V1=="Consortia"],
      Drug = annotations$V2[annotations$V1=="Drug"],
      Segmentation = annotations$V2[annotations$V1=="Segmentation"],
      StainingSet = annotations$V2[annotations$V1=="StainingSet"],
      Level = "Spot"),
      paste0(barcodePath, "/Analysis/", barcode,"_","SpotLevelAnnotations.tsv"), sep = "\t",col.names = FALSE, quote=FALSE)
  }
  cat("Elapsed time:", Sys.time()-functionStartTime, "\n")
}

#Debug purposes:  barcodePath <- "/lincs/share/lincs_user/LI8X00850"
#Debug purposes:  barcodePath <- "/lincs/share/lincs_user/lincs96well/LI9V01612"

barcodePath <-commandArgs(trailingOnly = TRUE)
res <- preprocessMEMASpot(barcodePath, verbose=TRUE)
