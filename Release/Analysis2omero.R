library(data.table)
library(magrittr)
library(readxl)

#Read in curated signals with descriptions
f <- read_excel("FeatureDescriptionsLevel4.xls")

#Read in the dataset manifest
dsm <- read_excel("~/GrayLabData/dane/MEP-LINCS/DatasetManifest.xlsx")

#Get a list of dataset file names that will get new analysis2omero files
parentPath <- "/graylab/share/dane/MEP-LINCS/MEP_LINCS/AnnotatedData/"
l3FileNames <- dir(parentPath, pattern="Level3", full.names = TRUE) %>% 
  grep("Trametinib",., value=TRUE)
l3FileNames <- l3FileNames[!grepl("SSC|SSR|_CS_",l3FileNames)]

tmpl <- lapply(l3FileNames, function(fn){
  datasetName <- gsub("_.*","",gsub(".*/","",fn))
  l3 <- fread(fn,key = "Barcode")
  dataPath <- dsm$Path[grepl(datasetName,dsm$DatasetName)]
  barcodes <- unique(l3$Barcode)
  for (barcode in barcodes){
    dt <- l3[barcode,f$Binding_CP[f$Binding_CP %in% colnames(l3)], with=FALSE]
    dt <- setnames(dt,colnames(dt),f$Name[f$Binding_CP %in% colnames(dt)])
    #Only display metadata and raw values in omero
    dt <- dt[,grep("normal",colnames(dt), value=TRUE, invert=TRUE), with=FALSE]
    fwrite(dt,file=paste0(dataPath,barcode,"/Analysis/",barcode,"_analysis2omero.tsv"),sep = "\t")
  }
})
