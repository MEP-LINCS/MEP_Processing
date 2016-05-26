library("jsonlite")#Reading in json files
library(data.table)
library(stringr)
library(parallel)
library(limma)


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
    content <- sample$content
    data.table(CellLine = sample$content, Passage = sample$fraction, CellSeedCount = sample$value)
  })
  sDT <- rbindlist(sl)
}

processJSON <- function (fileNames) {
  rbindlist(lapply(fileNames, function(fn){
    #Process each file separately
    plateMetadata <- fromJSON(fn)
    #Store and remove the welltype data
    welltype <- plateMetadata$welltype
    #Delete all of the items at the end of the json file
    plateMetadata$annot_id <- NULL
    plateMetadata$assayrun <- NULL
    plateMetadata$assaytype <- NULL
    plateMetadata$label <- NULL
    plateMetadata$welltype <- NULL
    #Get the barcode from the name, TODO change when barcode is in the data
    barcode <- str_replace_all(str_extract(fn,"[-][[:alnum:]]*_"),"[-_]*","")
    #Get the factor list for each spot
    #browser()
    mdList <- lapply(plateMetadata, function(spotMetadata){
      #Get protein info from factors
      spotFactors <- getFactors(spotMetadata$content$factor)
      #Get stain info from parameters
      spotParameters <- getParameters(spotMetadata$content$parameter)
      #Get cell line info from sample
      cellLineParameters <- getSample(spotMetadata$content$sample)
      #Get the well and spot information
      ArrayPositions <- str_split_fixed(spotMetadata["ii|ii"],"[|]",2)
      #cast the metadata into a single row for the spot
      wellRow <-as.integer(str_split_fixed(spotMetadata["i|i"],"[|]",2))[1]
      wellColumn <- as.integer(str_split_fixed(spotMetadata["i|i"],"[|]",2))[2]
      #browser()
      mDT <- data.table(CellLine = cellLineParameters$CellLine,
                        Well = wellAN(2,4)[(wellRow-1)*4+wellColumn],
                        ArrayRow = str_split_fixed(ArrayPositions[1,1],"_",2)[1,2],
                        ArrayColumn = str_split_fixed(ArrayPositions[1,2],"_",2)[1,2],
                        Passage = cellLineParameters$Passage,
                        CellSeedCount = cellLineParameters$CellSeedCount,
                        ECMp = spotFactors$ContentShortName[spotFactors$ContentType=="ECMp"],
                        ECMpPK = spotFactors$Content[spotFactors$ContentType=="ECMp"],
                        ECMpSet = spotFactors$Set[spotFactors$ContentType=="ECMp"],
                        Ligand = spotFactors$ContentShortName[spotFactors$ContentType=="Ligand"],
                        LigandPK = spotFactors$Content[spotFactors$ContentType=="Ligand"],
                        LigandSet = spotFactors$Set[spotFactors$ContentType=="Ligand"],
                        EndpointDAPI = spotParameters$Endpoint[grepl("DAPI",spotParameters$Channel)],
                        Endpoint488 = spotParameters$Endpoint[grepl("488",spotParameters$Channel)],
                        StainType488 = spotParameters$StainType[grepl("488",spotParameters$Channel)],
                        Animal488 = spotParameters$Animal[grepl("488",spotParameters$Channel)],
                        Endpoint555 = spotParameters$Endpoint[grepl("555|Orange",spotParameters$Channel)],
                        StainType555 = spotParameters$StainType[grepl("555|Orange",spotParameters$Channel)],
                        Animal555 = spotParameters$Animal[grepl("555|Orange",spotParameters$Channel)],
                        Endpoint647 = spotParameters$Endpoint[grepl("647|Red",spotParameters$Channel)],
                        StainType647 = spotParameters$StainType[grepl("647|Red",spotParameters$Channel)],
                        Animal647 = spotParameters$Animal[grepl("647|Red",spotParameters$Channel)])

    })
    mdDT <- rbindlist(mdList)
    mdDT$ECMpBase <- unique(mdDT$ECMp[grepl("^COL1_Own$",mdDT$ECMp)])
    return(mdDT)
  }))
}

fileNames <- dir("./json", pattern=".json", full.names = TRUE)

metadata <- rbindlist(mclapply(fileNames, processJSON, mc.cores=detectCores()))
