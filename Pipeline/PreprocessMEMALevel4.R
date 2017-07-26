#!/usr/bin/env Rscript

#author: "Mark Dane"
# 2017

library(MEMA)
library(githubr)
suppressPackageStartupMessages(library(optparse))

# Get the command line arguments and options
# returns a list with options and args elements
getCommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output"),
    make_option(c("-i", "--inputPath"), type="character", default=NULL, metavar="PATH",
                help="Path to local input data directory or Synapse ID for a File View."),
    make_option(c("--synapseStore"), type="character", default=NULL, metavar="SYNAPSEID",
                help="Store output file in Synapse directory (provide Synapse ID of Folder to store).")
  )
  parser <- OptionParser(usage = "%prog [options] STUDY OUTPUTFILE", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}


#Specify the command line options
cl <- getCommandLineArgs()
studyName <- cl$args[1]
ofname <- cl$args[2]

verbose <- cl$options$verbose
if(file.exists(cl$options$inputPath)){
  useSynapse <- FALSE
} else {
  useSynapse <- TRUE
}

startTime <- Sys.time()
if(useSynapse){
  suppressPackageStartupMessages(library(synapseClient))
  suppressMessages(synapseLogin())
  
  level <- "3"
  
  levelQuery <- sprintf('SELECT id,Segmentation,Preprocess,Study,DataType,Consortia,StainingSet,CellLine,Drug from %s WHERE Study="%s" AND Level="%s"',
                        cl$options$inputPath, studyName, level)
  levelRes <- synTableQuery(levelQuery)
  
  if (nrow(levelRes@values) > 1) {
    stop(sprintf("Found more than one Level %s file for Study %s", level, studyName))
  }
  
  dataPath <- getFileLocation(synGet(levelRes@values$id))
} else {
  dataPath <- cl$options$inputPath
}

l3DT <- fread(dataPath)

#Summarize to the MEP_Drug level
mepDT <- preprocessLevel4(l3DT,seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin"))
#Add in the barcodes for each MEP_Drug
mepDT <- addBarcodes(dt3 = l3DT, dt4 = mepDT)
# Add a QA flag for spots with few replicates
mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < 3

#reduce the numeric values to 4 significant digits
shorten <- function(x){
  if(class(x)=="numeric") x <- signif(x,4)
  return(x)
}

for (j in colnames(mepDT)) data.table::set(mepDT, j = j, value = shorten(mepDT[[j]]))

fwrite(mepDT, file=ofname, sep = "\t", quote=FALSE)
if(verbose) message("Writing level 4 file to disk\n")
if(!is.null(cl$options$synapseStore)){
  if(verbose) message(sprintf("Writing to Synapse Folder %s", cl$options$synapseStore))
  #get permlink from GitHub
  scriptLink <- "https://github.com/MEP-LINCS/MEP_Processing/"
  repo <- try(getRepo("MEP-LINCS/MEP_Processing", ref="branch", refName="master"),silent = TRUE)
  if(!class(repo)=="try-error" ) scriptLink <- getPermlink(repo, "Pipeline/PreprocessMEMALevel4.R")
  synFile <- File(ofname, parentId=cl$options$synapseStore)
  synSetAnnotations(synFile) <- list(CellLine = levelRes@values$CellLine,
                                     Study = levelRes@values$Study,
                                     Preprocess = levelRes@values$Preprocess,
                                     DataType = levelRes@values$DataType,
                                     Consortia = levelRes@values$Consortia,
                                     Drug = levelRes@values$Drug,
                                     Segmentation = levelRes@values$Segmentation,
                                     StainingSet = levelRes@values$StainingSet,
                                     Level = "4")
  
  synFile <- synStore(synFile,
                      used=c(levelRes@values$id),
                      executed=scriptLink,
                      forceVersion=FALSE)
}
if(verbose) message(paste("Elapsed time for ",studyName, "is", Sys.time()-startTime, "\n"))

