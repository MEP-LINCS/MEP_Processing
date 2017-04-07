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
synapsePreReleaseView <- "syn7494072"

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
Report = c("Cell", "SpotMEP", "Analysis"),
Activity = c("QA Cell Level", "QA Spot and MEP Level", "Analyze Replicates"),
stringsAsFactors=FALSE)
studies <- rbind(studies, data.frame(Study=c("MCF10A_SSC","HMEC122L_SSC","HMEC240L_SSC"),
                                     Report="SSC",
                                     Activity="Analyze Cell Line"))

#Get file paths from lincs server
studies$FileName <- unlist(apply(studies, 1, getFilePaths))

#get permlink from GitHub
repo <- getRepo("MEP-LINCS/MEP_LINCS", ref="branch", refName="develop")

# Convert to a list to use as Synapse annotations
toAnnotationList <- function(x) {
  as.list(x %>% select(-c(filename, Activity)))
}
# Take row of data frame with filename and annots
# Upload to Synapse and set annotations
uploadToSynapse <- function(x, used, scriptLink, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName=x$Activity,
                  used=used,
                  forceVersion=FALSE,
                  executed=scriptLink)
  obj
}

uploadReports <- function(x){
  dataFile <- unique(data.frame(Study=x$Study,
                                StainingSet=x$StainingSet,
                                CellLine=x$CellLine,
                                Drug=x$Drug,
                                Preprocess=x$Preprocess,
                                Segmentation=x$Segmentation,
                                DataType=x$DataType,
                                Consortia=x$Consortia,
                                Level=x$Level,
                                Activity=x$Activity,
                                filename=x$FileName,
                                stringsAsFactors = FALSE))
  #Separate used from annotations on the datafile
  if(unique(x$Level=="QA_Cell")){
    used <- unlist(x[["1"]])
    scriptLink <- getPermlink(repo, "Reports/MEP-LINCS_QACellLevel.Rmd")
  } 
  if(unique(x$Level=="QA_SpotMEP")){
    used <- c(x[["3"]],x[["4"]])
    scriptLink <- getPermlink(repo, "Reports/MEP-LINCS_QASpotMEPLevel.Rmd")
  } 
  if(grepl("Analysis",unique(x$Level))){
    used <- c(x[["3"]],x[["4"]])
    scriptLink <- getPermlink(repo, "Reports/MEP-LINCS_Analysis.Rmd")
  }
  if(grepl("SSC",unique(x$Level))){
    used <- c(x[["3"]],x[["4"]])
    scriptLink <- getPermlink(repo, "Reports/MEP-LINCS_AnalysisSB_SSC.Rmd")
  }
  uploadToSynapse(dataFile, used=used, scriptLink=scriptLink, parentId=synapseReportDir)
}

#Begin assembling file locations and annotations
#Get all files in th eprerelease are on Synpase
synIDs <- synTableQuery(sprintf('SELECT id,Level,StainingSet,Study,CellLine,Drug,Preprocess,Segmentation,DataType,Consortia from %s', synapsePreReleaseView))
synIDs <- data.frame(synIDs@values)
#Remove entities without STudy annotation
synIDs <- synIDs[!is.na(synIDs$Study),]
synIDs34 <- filter(synIDs, grepl("3|4",Level))
used <- dcast(synIDs34,Study~Level,value.var="id")
#For QA_Cell levels, add a level 1 column that's a list of the 8 barcodes in each staining set
#find the barcodes for each study in used
getLevel1Barcodes <- function(x){
  barcodes <- synIDs$id[synIDs$Study==x$Study&
                          synIDs$Level=="1"&
                          synIDs$DataType=="Quantitative Imaging"]
}
level1used <- rbind(dlply(used, "Study", getLevel1Barcodes))
#Add the list of barcodes as a column labeld 1 into the used datatable
used[["1"]] <- t(level1used)
#Add annotations, used files and executed script link
studies$Study <- tolower(studies$Study)
reports <- merge(used,studies,by=c("Study"))
reports$Level <-reports$Report
#Creat a level for each report
reports$Level[reports$Level %in% c("Cell","SpotMEP")] <- paste0("QA_",reports$Level[reports$Level %in% c("Cell","SpotMEP")])
reportsAnnotated <- unique(merge(reports,select(synIDs34,-c(id,Level)),by="Study"))
#Upload file with annotations, used files and executed script
res <- dlply(reportsAnnotated,"FileName", uploadReports)

