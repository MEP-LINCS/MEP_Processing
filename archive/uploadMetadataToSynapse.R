library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(githubr)
library(limma)
library(XLConnect)
library(data.table)

synapseLogin()
synapseRawDataDir <- "syn8565039"
synapseAnnotatedDataDir <- "syn5713302"
synapseReportDir <- "syn4939350"
synapseMetadataDir <- "syn7213947"
synapsePreReleaseView <- "syn7494072"
synapseAn2Study <- "syn8440875"

#Select datasets to upload
datasetAnns <- rbind(
  data.frame(CellLine=rep(c("MCF10A"), 4),
             StainingSet=c("SS1","SS2","SS3","SSC"),
             Drug="none",
             Preprocess="av1.8",
             Segmentation="v2",
             Study=c("mcf10a_ss1","mcf10a_ss2","mcf10a_ss3","mcf10a_ssc"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("HMEC240L"), 3),
             StainingSet=c("SS1","SS4","SSC"),
             Drug="none",
             Preprocess="av1.8",
             Segmentation="v2",
             Study=c("hmec240l_ss1","hmec240l_ss4","hmec240l_ssc"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("HMEC122L"), 3),
             StainingSet=c("SS1","SS4","SSC"),
             Drug="none",
             Preprocess="av1.8",
             Segmentation="v2",
             Study=c("hmec122l_ss1","hmec122l_ss4","hmec122l_ssc"),
             stringsAsFactors=FALSE)
)

datasetAnns$Consortia <- "MEP-LINCS"
datasetAnns$DataType <- "Metadata"
datasetAnns$Level <- "Metadata"
datasetAnns$Activity <- "Annotate Data"


#Get barcodes in study from Synapse file, deparse and convert to a data.table
study2Barcodes <- fread(synGet(synapseAn2Study)@filePath, header = TRUE, stringsAsFactors = FALSE)

barcodeLists <- apply(study2Barcodes,1,function(x){
  barcodes <- str_split(x[["Barcode"]],",")
  names(barcodes) <- x[["StudyName"]]
  return(barcodes)
})
barcodesDT <- lapply(barcodeLists, function(x){
  data.table(Study=names(x), Barcode=unlist(x))
})
studyDT <- rbindlist(barcodesDT)

#Add in the dataset annotations
ssDatasets <- merge(datasetAnns,studyDT, by="Study")

#Add the file paths on the lincs server
ssDatasets$FileName <-  paste0("/lincs/share/lincs_user/",ssDatasets$Barcode,"/Analysis/",ssDatasets$Barcode,"_an2omero.csv")

#get permlink from GitHub
repo <- getRepo("MEP-LINCS/MEP_Processing", ref="branch", refName="develop")

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

uploadFiles <- function(x){
  dataFile <- unique(data.frame(Study=x$Study,
                                StainingSet=x$StainingSet,
                                CellLine=x$CellLine,
                                Drug=x$Drug,
                                DataType=x$DataType,
                                Consortia=x$Consortia,
                                Level=x$Level,
                                Activity=x$Activity,
                                filename=x$FileName,
                                stringsAsFactors = FALSE))
  #Separate used from annotations on the datafile

  uploadToSynapse(dataFile, used=NULL, scriptLink=NULL, parentId=synapseMetadataDir)
}

#Upload file with annotations, used files and executed script
res <- dlply(ssDatasets,"Barcode", uploadFiles)

