#setup and Synpase login
library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(githubr)
library(limma)
library(data.table)
synapseLogin()

synapseAnnotatedDataDir <- "syn5713302"
synapseReportDir <- "syn4939350"
synapseMetadataDir <- "syn7213947"

#functions
getFilePaths <- function(x){
  path <-"/lincs/share/lincs_user/study"
  filePath <- dir(paste0(path,"/",x[["Study"]],"/Reports"), full.names = TRUE) %>%
    grep(x[["Report"]],.,value=TRUE)
  return(filePath)
}

#Define study and report level
#Select datasets to upload
studies <- data.frame(Study = rep(c("MCF10A_SS1",
                                        "MCF10A_SS2",
                                        "MCF10A_SS3",
                                        "HMEC122L_SS1",
                                        "HMEC122L_SS4",
                                        "HMEC240L_SS1",
                                        "HMEC240L_SS4"
),each=3),
Report = c("Cell", "SpotMep", "Analysis"),
stringsAsFactors=FALSE)

#Get file paths from lincs server
studies$FileName <- unlist(apply(studies, 1, getFilePaths))

#get permlink from GitHub
repo <- getRepo("MEP-LINCS/MEP_LINCS", ref="branch", refName="develop")

uploadReports <- function(x){
  path <- "/lincs/share/lincs_user/study"
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



#Get files used from Synpase
synIDs <- synQuery(paste("select id, name, Level, CellLine, StainingSet, Study from file where parentId=='",  synapseAnnotatedDataDir,"'"))
colnames(synIDs)<-gsub("file.","",colnames(synIDs))
synIDs <- filter(synIDs,!grepl("none|All",name))
synIDs <- select(synIDs, -name)
synIDs <- filter(synIDs, grepl("3|4|C",Level))
synIDs <- dcast(synIDs,Study~Level,value.var="id")
ssDatasets <- merge(synIDs,studies,by=c("Study"))
###Debug need to have level 3 and 4 files before proceeding
#Add annotations, used files and executed script link

#Upload file with annotations, used files and executed script