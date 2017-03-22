#!/usr/bin/env Rscript

#author: "Mark Dane"
# 2/2017


library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(RUVnormalize)
library(ruv)
library(stringr)
library(tidyr)
library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

#' Get the command line arguments and options
#' @return A list with options and args elements
#'@export
getL3CommandLineArgs <- function(){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
                help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false",
                dest="verbose", help="Print little output"),
    make_option(c("-s", "--Synapse"), action="store_true", default=TRUE,
                help="use Synpase for file accesses  [default]"),
    make_option(c("-l", "--localServer"), action="store_false",
                dest="Synapse", help="use a local server for file access"),
    make_option(c("-w", "--writeFiles"), action="store_true", default=TRUE,
                help="write output files to disk  [default]"),
    make_option(c("-n", "--noWriteFiles"), action="store_false",
                dest="writeFiles", help="do not write output to disk"),
    make_option(c("-k", "--k"), type="integer", default=256,
                help="Number of factors to use in RUV normalization [default %default]",
                metavar="number")
  )
  parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
  arguments <- parse_args(parser, positional_arguments = 2)
}

#Specify the command line options
###Debug
cl <-list(options=list(k=256,
                       verbose=TRUE,
                       Synapse=FALSE,
                       writeFiles=FALSE),
          args=c("HMEC122L_SS1", "/lincs/share/lincs_user"))
####

startTime <- Sys.time()
#Get and decode command line arguments and options
cl <- getL3CommandLineArgs()
studyName <- cl$args[1]
path <- cl$args[2]

opt <- cl$options
k <- opt$k
verbose <- opt$verbose
useSynapse <- opt$Synapse
writeFiles <- opt$writeFiles


#Read the annotated data for all plates in the study
if(useSynapse){
  stop("Synapse not supported for level 3 data yet ")
} else {
  slDT <- getSpotLevelData(studyName, path)
}

signalsMinimalMetadata <- grep("_SE",grep("_CP_|_PA_|Barcode|^Well$|^Spot$|^PrintSpot$|^Ligand$|^ECMp$|^Drug$|^Drug1Conc$|^ArrayRow$|^ArrayColumn$|^CellLine$",colnames(slDT), value=TRUE), value=TRUE, invert=TRUE)

#RUVLoess normalize all signals
if(!k==0){
  if(verbose)  message(paste("Normalizing", studyName,"\n"))
  #slDT <- slDT[!is.na(slDT$Cytoplasm_CP_AreaShape_Compactness)&!is.na(slDT$Cytoplasm_CP_AreaShape_Eccentricity),]
  nDT <- normRUVLoessResiduals(slDT[,signalsMinimalMetadata, with = FALSE], k)
  nDT$NormMethod <- "RUVLoessResiduals"
  slDT$k <- k
  slDT <- merge(slDT, nDT, by = c("BW","PrintSpot"))
} else {
  slDT$NormMethod <- "none"
  slDT$k <- k
}

#Add QA flags to the data
slDT <- QASpotLevelData(slDT, lowSpotCellCountThreshold=5,
                        lowRegionCellCountThreshold = 0.4,
                        lowWellQAThreshold = .7)

if(writeFiles){
  if(verbose) message(paste("Writing level 3 file to disk\n"))
  if(useSynapse){
    stop("Synapse not supported for level 3 data yet ")
  } else {
    fwrite(data.table(slDT), paste0(path, "/study/",studyName, "/Annotated/", studyName,"_Level3.tsv"), sep = "\t", quote=FALSE)
  }

}
message(paste("Elapsed time to normalize ",studyName, Sys.time()-startTime, "\n"))


