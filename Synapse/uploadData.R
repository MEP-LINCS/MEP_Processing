library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(githubr)
library(limma)
library(data.table)

# Take row of data frame and remove file name
# Convert to a list to use as Synapse annotations
toAnnotationList <- function(x) {
  if("Level3ID" %in% colnames(x)) x <- select(x,-Level3ID)
  if("Used3" %in% colnames(x)) x <- select(x,-Used3)
  if("Used4" %in% colnames(x)) x <- select(x,-Used4)
  if("Level1ID" %in% colnames(x)) x <- select(x,-Level1ID)
  as.list(x %>% select(-c(filename)))
}

uploadToSynapseLevel1 <- function(x, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName=activityName,
                  forceVersion=FALSE,
                  executed=scriptLink)
  obj
}

# Upload to Synapse and set annotations
uploadToSynapseLevel3 <- function(x, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName=activityName,
                  used=c(x$Level1ID),
                  contentType="text/tab-separated-values",
                  forceVersion=FALSE,
                  executed=scriptLink)
  obj
}

uploadToSynapseLevel4 <- function(x, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName=activityName,
                  used=c(x$Level3ID),
                  contentType="text/tab-separated-values",
                  forceVersion=FALSE,
                  executed=scriptLink)
  obj
}

uploadToSynapseSSCData <- function(x, used, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName=activityName,
                  used=used,
                  contentType="text/tab-separated-values",
                  forceVersion=FALSE,
                  executed=scriptLink)
  obj
}

uploadFileLevel1 <- function(x){
  dataDir <- paste("../AnnotatedData")
  # Take file names and turn into basic annotation set
  # Replace this with a better way to get basic annotations from 
  # a standardized source
  dataFile <- data.frame(CellLine=x$CellLine,
                         StainingSet=x$StainingSet,
                         Drug=x$Drug,
                         Preprocess=x$Preprocess,
                         Segmentation=x$Segmentation,
                         DataType=dataType,
                         Consortia="MEP-LINCS",
                         Level=x$Level,
                         filename=x$FilePath,
                         stringsAsFactors = FALSE)
  uploadToSynapseLevel1(dataFile, parentId=synapseAnnotatedDataDir)
}

uploadFileLevel3 <- function(x){
  dataDir <- paste("../AnnotatedData")
  # Take file names and turn into basic annotation set
  # Replace this with a better way to get basic annotations from 
  # a standardized source
  dataFile <- data.frame(CellLine=x$CellLine,
                         StainingSet=x$StainingSet,
                         Drug=x$Drug,
                         Preprocess=x$Preprocess,
                         Segmentation=x$Segmentation,
                         DataType=dataType,
                         Consortia="MEP-LINCS",
                         Level=x$Level,
                         filename=x$FilePath,
                         Level1ID=x$Level1ID,
                         stringsAsFactors = FALSE)
  uploadToSynapseLevel3(dataFile, parentId=synapseAnnotatedDataDir)
}

uploadFileLevel4 <- function(x){
  dataDir <- paste("../AnnotatedData")
  # Take file names and turn into basic annotation set
  # Replace this with a better way to get basic annotations from 
  # a standardized source
  dataFile <- data.frame(CellLine=x$CellLine,
                         StainingSet=x$StainingSet,
                         Drug=x$Drug,
                         Preprocess=x$Preprocess,
                         Segmentation=x$Segmentation,
                         DataType=dataType,
                         Consortia="MEP-LINCS",
                         Level=x$Level,
                         filename=x$FilePath,
                         Level3ID=x$Level3ID,
                         stringsAsFactors = FALSE)
  if("1" %in% colnames(x)) dataFile$Used1 <- x[["1"]]
  if("3" %in% colnames(x)) dataFile$Used3 <- x[["3"]]
  if("4" %in% colnames(x)) dataFile$Used4 <- x[["4"]]
  uploadToSynapseLevel4(dataFile, parentId=synapseAnnotatedDataDir)
}

uploadFileSSCData <- function(x){
  dataDir <- paste("../AnnotatedData")
  # Take file names and turn into basic annotation set
  # Replace this with a better way to get basic annotations from 
  # a standardized source
  dataFile <- unique(data.frame(CellLine=x$CellLine,
                         StainingSet="SSC",
                         Drug=x$Drug,
                         Preprocess=x$Preprocess,
                         Segmentation=x$Segmentation,
                         DataType=dataType,
                         Consortia="MEP-LINCS",
                         Level=x$Level,
                         filename=x$FilePath,
                         stringsAsFactors = FALSE))
  #Separate used from annotations on the datafile
  x <- filter(x,!grepl("SSC",StainingSet))
  if(unique(x$Level==3)) used <- c(x[["3"]])
  if(unique(x$Level==4)) used <- c(x[["4"]])
  uploadToSynapseSSCData(dataFile, used=used, parentId=synapseAnnotatedDataDir)
}

synapseLogin()
synapseRawDataDir <- "syn6167751"
synapseAnnotatedDataDir <- "syn5713302"
synapseReportDir <- "syn4939350"
synapseMetadataDir <- "syn5007108"

#Get github repo link
repo <- getRepo("MEP-LINCS/MEP_LINCS", ref="branch", refName="master")

#Select datasets to upload
datasetAnns <- rbind(
  data.frame(CellLine=rep(c("PC3"), 4),
             StainingSet=c("SS1", "SS2","SS3","SS2noH3"),
             Drug="none",
             Preprocess="av1.4",
             Segmentation=c("v2","v2.1","v2.1", "v1"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF7"), 3),
             StainingSet=c("SS1", "SS2","SS3"),
             Drug="none",
             Preprocess="av1.4",
             Segmentation=c("v2","v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("YAPC"), 3),
             StainingSet=c("SS1","SS2","SS3"),
             Drug="none",
             Preprocess="av1.4",
             Segmentation=c("v2","v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF10A"), 4),
             StainingSet=c("SS1","SS2","SS3","SSC"),
             Drug="none",
             Preprocess="av1.7",
             Segmentation="v2",
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("HMEC240L"), 3),
             StainingSet=c("SS1","SS4","SSC"),
             Drug="none",
             Preprocess="av1.7",
             Segmentation="v2",
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("HMEC122L"), 3),
             StainingSet=c("SS1","SS4","SSC"),
             Drug="none",
             Preprocess="av1.7",
             Segmentation="v2",
             stringsAsFactors=FALSE))

#Begin data or report specific code
scriptName <- "PreprocessingL1L2"
dataType <- "Quantitative Imaging"
activityName <- "Annotate data"

#Get GitHub link for next uploads
scriptLink <- getPermlink(repo, paste0("Release/MEP-LINCS_",scriptName,".R"))

#Get the filepaths for the zipped level 1 data on the graylab server
filePaths <- grep("[LA]_SS[1234C]_Level1.zip",dir("../AnnotatedData",full.names = TRUE), value=TRUE)

#Get annotations from filename
splits <- strsplit2(filePaths,"_")
serverFiles <- data.frame(FilePath=filePaths,StainingSet=splits[,2], CellLine=gsub("../AnnotatedData/","",splits[,1]),Level=as.numeric(gsub("Level|.txt|.zip","",splits[,3])), stringsAsFactors = FALSE)

ssDatasets <- merge(datasetAnns,serverFiles, by=c("CellLine","StainingSet"))

#Upload the level 1 files without a provenance to lower level files
res <- dlply(ssDatasets, c("CellLine","StainingSet"), uploadFileLevel1)

#Prepare to upload level 3 data
scriptName <- "PreprocessingL3L4"
dataType <- "Quantitative Imaging"
activityName <- "Annotate data"
scriptLink <- getPermlink(repo, paste0("Release/MEP-LINCS_",scriptName,".R"))

#Get the filepaths for the level 3 data on the graylab server
filePaths <- grep("[LA]_SS[1234]_Level3",dir("../AnnotatedData",full.names = TRUE), value=TRUE)
#Get annotations from filename
splits <- strsplit2(filePaths,"_")
serverFiles <- data.frame(FilePath=filePaths,StainingSet=splits[,2], CellLine=gsub("../AnnotatedData/","",splits[,1]),Level=as.numeric(gsub("Level|.txt","",splits[,3])), stringsAsFactors = FALSE)

ssDatasets <- merge(datasetAnns,serverFiles, by=c("CellLine","StainingSet"))

#Get the level 1 synIDs for the provenance
synIDs <- synQuery(paste("select id, name, Level, CellLine, StainingSet from file where parentId=='",  synapseAnnotatedDataDir,"'"))
colnames(synIDs)<-gsub("file.","",colnames(synIDs))
synIDs <- filter(synIDs,!grepl("none|All",name))
synIDs <- select(synIDs, -name)
synIDs <- filter(synIDs, Level=="1")
synIDs <- select(synIDs, -Level)

#Merge in the synIDs of the level 1 data
ssDatasets <- merge(ssDatasets,synIDs)
ssDatasets <- rename(ssDatasets,Level1ID=id)

#Upload the level 3 files without a provenance to lower level files
res <- dlply(ssDatasets, c("CellLine","StainingSet"), uploadFileLevel3)

#Get the filepaths for the level 4 data on the graylab server
filePaths <- grep("[LA]_SS[1234]_Level4",dir("../AnnotatedData",full.names = TRUE), value=TRUE)
splits <- strsplit2(filePaths,"_")
serverFiles <- data.frame(FilePath=filePaths,StainingSet=splits[,2], CellLine=gsub("../AnnotatedData/","",splits[,1]),Level=as.numeric(gsub("Level|.txt","",splits[,3])), stringsAsFactors = FALSE)
ssDatasets <- merge(datasetAnns,serverFiles, by=c("CellLine","StainingSet"))

synIDs <- synQuery(paste("select id, name, Level, CellLine, StainingSet from file where parentId=='",  synapseAnnotatedDataDir,"'"))
colnames(synIDs)<-gsub("file.","",colnames(synIDs))
synIDs <- filter(synIDs,!grepl("none|All",name))
synIDs <- select(synIDs, -name)
synIDs <- filter(synIDs, Level=="3")
synIDs <- select(synIDs, -Level)

#Merge in the synIDs of the level 3 data
ssDatasets <- merge(ssDatasets,synIDs)
ssDatasets <- rename(ssDatasets,Level3ID=id)

#Begin the level 4 uploads
res <- dlply(ssDatasets, c("CellLine","StainingSet"), uploadFileLevel4)

#Get the info for the level 3 SSC data
scriptName <- "CommonSignalNormx1"
dataType <- "Quantitative Imaging"
activityName <- "Normalize and summarize data"
scriptLink <- getPermlink(repo, paste0("Release/MEP-LINCS_",scriptName,".R"))
filePaths <- grep("[LA]_SSC_Level3",dir("../AnnotatedData",full.names = TRUE), value=TRUE)
splits <- strsplit2(filePaths,"_")
serverFiles <- data.frame(FilePath=filePaths,StainingSet=splits[,2], CellLine=gsub("../AnnotatedData/","",splits[,1]),Level=as.numeric(gsub("Level|.txt","",splits[,3])), stringsAsFactors = FALSE)
ssDatasets <- merge(datasetAnns,serverFiles, by=c("CellLine","StainingSet"))

synIDs <- synQuery(paste("select id, name, Level, CellLine, StainingSet from file where parentId=='",  synapseAnnotatedDataDir,"'"))
colnames(synIDs)<-gsub("file.","",colnames(synIDs))
synIDs <- filter(synIDs,!grepl("none|All",name))
synIDs <- select(synIDs, -name)
synIDs <- filter(synIDs, grepl("1|3|4|C",Level))
synIDs <- dcast(synIDs,CellLine+StainingSet~Level,value.var="id")

#Merge in the synIDs of the level 3 data
ssDatasets <- select(ssDatasets, -StainingSet)
ssDatasets <- merge(ssDatasets,synIDs,by=c("CellLine"))

#Begin the level 3 SSC data uploads
res <- dlply(ssDatasets, c("CellLine"), uploadFileSSCData)

#Upload level 4 SSC data
filePaths <- grep("[LA]_SSC_Level4",dir("../AnnotatedData",full.names = TRUE), value=TRUE)
splits <- strsplit2(filePaths,"_")
serverFiles <- data.frame(FilePath=filePaths,StainingSet=splits[,2], CellLine=gsub("../AnnotatedData/","",splits[,1]),Level=as.numeric(gsub("Level|.txt","",splits[,3])), stringsAsFactors = FALSE)
ssDatasets <- merge(datasetAnns,serverFiles, by=c("CellLine","StainingSet"))

synIDs <- synQuery(paste("select id, name, Level, CellLine, StainingSet from file where parentId=='",  synapseAnnotatedDataDir,"'"))
colnames(synIDs)<-gsub("file.","",colnames(synIDs))
synIDs <- filter(synIDs,!grepl("none|All",name))
synIDs <- select(synIDs, -name)
synIDs <- filter(synIDs, grepl("1|3|4|C",Level))
synIDs <- dcast(synIDs,CellLine+StainingSet~Level,value.var="id")

#Merge in the synIDs of the level 3 data
ssDatasets <- select(ssDatasets, -StainingSet)
ssDatasets <- merge(ssDatasets,synIDs,by=c("CellLine"))

#Begin the level 4 SSC data uploads
res <- dlply(ssDatasets, c("CellLine"), uploadFileSSCData)
