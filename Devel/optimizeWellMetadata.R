#optimze reading the well metadata
#Mark Dane 2/2016

readMetadata <- function (xlsFile) 
{
  sheetList <- sapply(gdata::sheetNames(path.expand(xlsFile)), 
                      gdata::read.xls, xls = path.expand(xlsFile), simplify = FALSE, 
                      stringsAsFactors = TRUE, check.names = FALSE, row.names = "Row/Column")
  nrRows <- dim(sheetList[[1]])[1]
  nrCols <- as.numeric(max(colnames(sheetList[[1]])))
  nrWells = nrRows * nrCols
  sheetDF <- data.frame(lapply(sheetList, function(df, nrCols) {
    dfM <- matrix(t(df[, 1:nrCols]), byrow = TRUE)
  }, nrCols = nrCols), WellIndex = 1:nrWells, Well = wellAN(nrRows, 
                                                            nrCols), check.names = TRUE, stringsAsFactors = FALSE)
  return(sheetDF)
}


getMetadataFileNames <- function(path){
  mdFiles <- data.frame(FilePaths = dir(path, full.names = TRUE), stringsAsFactors = FALSE)
  mdFiles$Type <- "WellMetadata"
  mdFiles$Type[grepl(".gal",mdFiles$FilePaths,ignore.case = TRUE)] <- "GAL"
  mdFiles$Type[grepl(".xml",mdFiles$FilePaths,ignore.case = TRUE)] <- "PrintLog"
  mdFiles$Type[grepl("imageIDs",mdFiles$FilePaths,ignore.case = TRUE)] <- "ImageIDs"
  return(mdFiles)
}

mdf <- getMetadataFileNames("PC3/SS1/Metadata")
wmd <- readMetadata(fileName)

wmDTL <- mclapply(mdf$FilePaths[mdf$Type=="WellMetadata"], function(fileName){
  dt <- data.table(readMetadata(fileName), key="Well")
  dt <- dt[,Barcode := sub(".*/","",sub(".xlsx","",fileName))]
  return(convertColumnNames(dt))
}, mc.cores = detectCores())
wmDT <- rbindlist(wmDTL)
