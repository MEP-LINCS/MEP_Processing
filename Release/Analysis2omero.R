library(data.table)
library(magrittr)
library(readxl)

#Read in curated signals with descriptions
f <- read_excel("FeatureDescriptionsOmero.xls")

#Get a list of barcodes that will get new annotate2omero files
parentPath <- "/graylab/share/dane/MEP-LINCS/MEP_LINCS/AnnotatedData/"
l3FileNames <- dir(parentPath, pattern="Level3", full.names = TRUE) %>% 
  grep("MCF10A|HMEC|MCF7|PC3|YAPC",., value=TRUE) %>%
  grep("DMSO|Neratinib",., value=TRUE)
l3FileNames <- l3FileNames[!grepl("SSC|SSR|_CS_",l3FileNames)]

tmpl <- lapply(l3FileNames, function(fn){
  l3 <- fread(fn,key = "Barcode")
  barcodes <- unique(l3$Barcode)
  for (barcode in barcodes){
    dt <- l3[barcode,f$Binding[f$Binding %in% colnames(l3)], with=FALSE]
    dtd <- copy(dt)
    dtd <- setnames(dtd,colnames(dtd),f$Name[f$Binding %in% colnames(dtd)])
    fwrite(dtd,file=paste0("/lincs/share/lincs_user/",barcode,"/Analysis/",barcode,"_analysis2omero.tsv"),sep = "\t")
  }
})
