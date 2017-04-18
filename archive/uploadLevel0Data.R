library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(githubr)
library(limma)
library(XLConnect)
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

uploadToSynapseLevel0 <- function(x, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  
  obj <- synStore(obj, 
                  activityName=activityName,
                  forceVersion=FALSE)
  obj
}

uploadFileLevel0 <- function(x){
  dataFile <- data.frame(CellLine=x$CellLine,
                         StainingSet=x$StainingSet,
                         Drug=x$Drug,
                         Preprocess=x$Preprocess,
                         Segmentation=x$Segmentation,
                         Barcode=x$Barcode,
                         Well=x$Well,
                         Location=x$Location,
                         DataType="Quantitative Imaging",
                         Consortia="MEP-LINCS",
                         Level="0",
                         Study=x$Study,
                         filename=x$Path,
                         stringsAsFactors = FALSE)
  uploadToSynapseLevel0(dataFile, parentId=synapseRawDataDir)
}

synapseLogin()
synapseRawDataDir <- "syn9611756" #Pre-releaseYr3
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
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF10A"), 1),
             StainingSet=c("SSL"),
             Drug="none",
             Preprocess="av1.8",
             Segmentation="v2",
             Study=c("mcf10a_mema_v2"),
             stringsAsFactors=FALSE)
)

#Define annotations and provenance values
activityName <- "Analyze Images"

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

#Get raw data files from the local server
fnames <- apply(ssDatasets[ssDatasets$Barcode %in% c("LI8X00801","LI8X00802","LI8X00803","LI8X00804","LI8X00805","LI8X00806","LI8X00807","LI8X00808"),], 1, function(x){
  data.table(Study=x[["Study"]],
             CellLine=x[["CellLine"]],
             StainingSet =x[["StainingSet"]],
             Drug = x[["Drug"]],
             Preprocess = x[["Preprocess"]],
             Segmentation = x[["Segmentation"]],
             Barcode =x[["Barcode"]],
             Path=dir(paste0("/lincs/share/lincs_user/",x[["Barcode"]],"/Analysis/",x[["Segmentation"]]),full.names = TRUE))
})

fnamesDT <- rbindlist(fnames)
#Pull well and locations out of filenames
fnamesDT$Well <- str_extract(fnamesDT$Path,"_[[:alpha:]][[:digit:]]{2}_") %>%
  str_replace_all("_","")
fnamesDT$Location <- str_extract(fnamesDT$Path,"Image|Cells|Nuclei|Cytoplasm")
#Upload the level o files without a provenance to lower level files
res <- dlply(fnamesDT, c("Path"), uploadFileLevel0)

