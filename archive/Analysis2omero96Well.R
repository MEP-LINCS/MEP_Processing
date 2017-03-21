library(data.table)
library(magrittr)
library(readxl)

#Read in curated signals with descriptions
f <- read_excel("FeatureDescriptionsLevel4.xls")

#Read in the dataset manifest
dsm <- read_excel("~/GrayLabData/dane/MEP-LINCS/DatasetManifest.xlsx")

#Get a list of dataset file names that will get new analysis2omero files
parentPath <- dsm$Path[grepl("lincs96well",dsm$DatasetName)]
dirNames <- dir(parentPath, full.names = TRUE)
#Get the full paths to the level3 files in each barcode directory
l3FileNames <- grep("Level3",lapply(dirNames, function(x) dir(paste0(x,"/Analysis"), pattern="Level3", full.names=TRUE)),value = TRUE)

tmpl <- lapply(l3FileNames, function(fn){
  l3 <- fread(fn,key = "Barcode")
  barcodes <- unique(l3$Barcode)
  for (barcode in barcodes){
    dt <- l3[barcode,f$Binding_IC[f$Binding_IC %in% colnames(l3)], with=FALSE]
    dt <- setnames(dt,colnames(dt),f$Name[f$Binding_IC %in% colnames(dt)])
    #Only display metadata and raw values in omero
    dt <- dt[,grep("normal",colnames(dt), value=TRUE, invert=TRUE), with=FALSE]
    fwrite(dt,file=paste0(parentPath,"/",barcode,"/Analysis/",barcode,"_analysis2omero.tsv"),sep = "\t")
  }
})
