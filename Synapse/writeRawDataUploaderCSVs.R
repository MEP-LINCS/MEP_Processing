library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(rGithubClient)

# Take row of data frame and remove file name
# Convert to a list to use as Synapse annotations
toAnnotationList <- function(x) {
  as.list(x %>% select(-filename))
}

# Take row of data frame with filename and annots
# Upload to Synapse and set annotations
writeCSV <- function(x) {
  write.csv(x, file=paste0(paste(unique(x$CellLine),unique(x$StainingSet),unique(x$Segmentation),sep="_"),".csv"), row.names=FALSE)

}

#Select datasets to upload
ssDatasets <- rbind(
  data.frame(CellLine=rep(c("PC3"), 4),
             StainingSet=c("SS1", "SS2","SS3","SS2noH3"),
             Segmentation=c("v2","v2.1","v2.1", "v1"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF7"), 3),
             StainingSet=c("SS1", "SS2","SS3"),
             Segmentation=c("v2","v2","v2"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("YAPC"), 3),
             StainingSet=c("SS1","SS2","SS3"),
             Segmentation=c("v2"),
             stringsAsFactors=FALSE),
  data.frame(CellLine=rep(c("MCF10A"), 2),
             StainingSet=c("SS1","SS3"),
             Segmentation=c("v2","v2"),
             stringsAsFactors=FALSE)
)


getPaths <- function(x){
  dataDir <- paste("/Users/dane/Documents/MEP-LINCS",x$CellLine,x$StainingSet,"RawData",x$Segmentation, sep = "/")
  # Take file names and turn into basic annotation set
  # Replace this with a better way to get basic annotations from 
  # a standardized source
  dataFiles <- data.frame(filename=list.files(path=dataDir), stringsAsFactors = FALSE) %>%
    mutate(Level=0,
           Consortia="MEP-LINCS",
           DataType="Quantitative Imaging",
           CellLine=x$CellLine,
           StainingSet=x$StainingSet,
           Segmentation=x$Segmentation) %>% 
    mutate(basename=str_replace(filename, "\\.csv", "")) %>% 
    separate(basename, c("Barcode", "Well", "Location"))

  writeCSV(dataFiles)
  
}

dataFiles <- do.call(rbind, dlply(ssDatasets[7,], c("CellLine","StainingSet"), getPaths))


