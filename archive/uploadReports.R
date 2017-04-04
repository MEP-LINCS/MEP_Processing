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
  if("Used1" %in% colnames(x)) x <- select(x,-Used1)
  if("Used3" %in% colnames(x)) x <- select(x,-Used3)
  if("Used4" %in% colnames(x)) x <- select(x,-Used4)
  if("Level1SynID" %in% colnames(x)) x <- select(x,-Level1SynID)
  as.list(x %>% select(-c(filename)))
}

synapseLogin()
synapseAnnotatedDataDir <- "syn5713302"
synapseReportDir <- "syn4939350"
synapseMetadataDir <- "syn7213947"

#Get github repo link
repo <- getRepo("MEP-LINCS/MEP_LINCS", ref="branch", refName="master")

#Select datasets to upload
datasetAnns <- rbind(
  data.frame(CellLine=rep(c("MCF10A"), 4),
             StainingSet=c("SS1","SS2","SS3","SSC"),
             Drug="none",
             Preprocess="v1.8",
             Segmentation="v2",
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("HMEC240L"), 3),
             StainingSet=c("SS1","SS4","SSC"),
             Drug="none",
             Preprocess="v1.8",
             Segmentation="v2",
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("HMEC122L"), 3),
             StainingSet=c("SS1","SS4","SSC"),
             Drug="none",
             Preprocess="v1.8",
             Segmentation="v2",
             stringsAsFactors=FALSE))

#Upload the reports
uploadToSynapseReports <- function(x, used, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName=activityName,
                  forceVersion=FALSE,
                  used=used,
                  executed=scriptLink)
  obj
}

uploadReports <- function(x){
  dataDir <- paste("../QAReports/")
  dataFile <- unique(data.frame(CellLine=x$CellLine,
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
  if(unique(x$Level=="QA_Cell")) used <- c(x[["1"]])
  if(unique(x$Level=="QA_SpotMEP")) used <- c(x[["3"]],x[["4"]])
  if(grepl("Analysis",unique(x$Level))) used <- c(x[["3"]],x[["4"]])
  uploadToSynapseReports(dataFile, used=used, parentId=synapseReportDir)
}

dataType <- "Quantitative Imaging"
activityName <- "QA data"

#Get the filepaths for the QA reports on the graylab server
filePaths <- grep("QA_Cell_.*_SS..html",dir("../QAReports/",full.names = TRUE), value=TRUE)
splits <- strsplit2(filePaths,"_")

serverFiles <- data.frame(FilePath=filePaths,StainingSet=gsub(".html","",splits[,5]), CellLine=splits[,4],Level=paste(splits[,2],splits[,3],sep="_"), stringsAsFactors = FALSE)

ssDatasets <- merge(datasetAnns,serverFiles, by=c("CellLine","StainingSet"))

synIDs <- synQuery(paste("select id, name, Level, CellLine, StainingSet from file where parentId=='",  synapseAnnotatedDataDir,"'"))
colnames(synIDs)<-gsub("file.","",colnames(synIDs))
synIDs <- filter(synIDs,!grepl("none|All",name))
synIDs <- select(synIDs, -name)
synIDs <- filter(synIDs, grepl("1|3|4|C",Level))
synIDs <- dcast(synIDs,CellLine+StainingSet~Level,value.var="id")
ssDatasets <- merge(synIDs,ssDatasets,by=c("CellLine","StainingSet"))

reportType <- "QACellLevel"
scriptLink <- getPermlink(repo, paste0("Release/MEP-LINCS_",reportType,".Rmd"))

res <- dlply(ssDatasets, c("CellLine","StainingSet"), uploadReports)


#Get the filepaths for the QA Spot reports on the graylab server
filePaths <- grep("QA_SpotMEP_.*_SS..html",dir("../QAReports/",full.names = TRUE), value=TRUE)
splits <- strsplit2(filePaths,"_")

serverFiles <- data.frame(FilePath=filePaths,StainingSet=gsub(".html","",splits[,5]), CellLine=splits[,4],Level=paste(splits[,2],splits[,3],sep="_"), stringsAsFactors = FALSE)

ssDatasets <- merge(datasetAnns,serverFiles, by=c("CellLine","StainingSet"))
synIDs <- synQuery(paste("select id, name, Level, CellLine, StainingSet from file where parentId=='",  synapseAnnotatedDataDir,"'"))
colnames(synIDs)<-gsub("file.","",colnames(synIDs))
synIDs <- filter(synIDs,!grepl("none|All",name))
synIDs <- select(synIDs, -name)
synIDs <- filter(synIDs, grepl("1|3|4|C",Level))
synIDs <- dcast(synIDs,CellLine+StainingSet~Level,value.var="id")
ssDatasets <- merge(synIDs,ssDatasets,by=c("CellLine","StainingSet"))

reportType <- "QASpotMEPLevel"
scriptLink <- getPermlink(repo, paste0("Release/MEP-LINCS_",reportType,".Rmd"))

res <- dlply(ssDatasets, c("CellLine","StainingSet"), uploadReports)

#Upload Analysis Files
#Get the filepaths for the analysis reports on the graylab server
activityName <- "Analyze Staining Set Data"
filePaths <- grep("Analysis_.*_SS.[_]?.html",dir("../AnalysisReports/",full.names = TRUE), value=TRUE)
splits <- gsub(".html","",strsplit2(filePaths,"_"))

serverFiles <- data.frame(FilePath=filePaths,StainingSet=splits[,4], CellLine=splits[,3],Level=splits[,2], stringsAsFactors = FALSE)
ssDatasets <- merge(datasetAnns,serverFiles, by=c("CellLine","StainingSet"))
synIDs <- synQuery(paste("select id, name, Level, CellLine, StainingSet from file where parentId=='",  synapseAnnotatedDataDir,"'"))
colnames(synIDs)<-gsub("file.","",colnames(synIDs))
synIDs <- filter(synIDs,!grepl("none|All",name))
synIDs <- select(synIDs, -name)
synIDs <- filter(synIDs, grepl("3|4",Level))
synIDs <- filter(synIDs, grepl("1|2|3|4",StainingSet))
synIDs <- dcast(synIDs,CellLine+StainingSet~Level,value.var="id")

#Merge in the synIDs of the level 3 data
ssDatasets <- merge(ssDatasets,synIDs,by=c("CellLine","StainingSet"))
reportType <- "Analysis"
scriptLink <- getPermlink(repo, paste0("Release/MEP-LINCS_",reportType,".Rmd"))
res <- dlply(ssDatasets, c("CellLine","StainingSet"), uploadReports)

#Get the filepaths for the SB analysis reports on the graylab server
#and used data on Synapse
activityName <- "Analyze Cell Line Data"
filePaths <- grep("AnalysisSB_.*_SS.[_]?.html",dir("../AnalysisReports/",full.names = TRUE), value=TRUE)
splits <- gsub(".html","",strsplit2(filePaths,"_"))

serverFiles <- data.frame(FilePath=filePaths,StainingSet=splits[,4], CellLine=splits[,3],Level=splits[,2], stringsAsFactors = FALSE)
ssDatasets <- merge(datasetAnns,serverFiles, by=c("CellLine","StainingSet"))
ssDatasets <- select(ssDatasets, -StainingSet)
synIDs <- synQuery(paste("select id, name, Level, CellLine, StainingSet from file where parentId=='",  synapseAnnotatedDataDir,"'"))
colnames(synIDs)<-gsub("file.","",colnames(synIDs))
synIDs <- filter(synIDs,!grepl("none|All",name))
synIDs <- select(synIDs, -name)
synIDs <- filter(synIDs, grepl("3|4|C",Level))
synIDs <- dcast(synIDs,CellLine+StainingSet~Level,value.var="id")
ssDatasets <- merge(synIDs,ssDatasets,by=c("CellLine"))
reportType <- "AnalysisSB"
scriptLink <- getPermlink(repo, paste0("Release/MEP-LINCS_",reportType,"_SSC.Rmd"))
res <- dlply(ssDatasets, c("CellLine"), uploadReports)

