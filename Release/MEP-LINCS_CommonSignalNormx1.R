library("rmarkdown")
source("MEP_LINCS/Release/MEPLINCSFunctions.R")

PC3df <- data.frame(cellLine=rep(c("PC3"), 3),
                    k=135L,
                    ss=c("SS1", "SS2","SS3"),
                    drug="none",
                    analysisVersion="av1.6",
                    rawDataVersion=c("v2","v2.1","v2.1"),
                    stringsAsFactors=FALSE)

MCF7df <- data.frame(cellLine=rep(c("MCF7"), 3),
                     k=135L,
                     ss=c("SS1", "SS2","SS3"),
                     drug="none",
                     analysisVersion=c("av1.6"),
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

YAPCdf <- data.frame(cellLine=rep(c("YAPC"), 3),
                     k=135L,
                     ss=c("SS1", "SS2","SS3"),
                     drug="none",
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(cellLine=rep(c("MCF10A"), 3),
                       k=135L,
                       ss=c("SS1","SS2","SS3"),
                       drug="none",
                       analysisVersion="av1.7",
                       rawDataVersion=c("v2","v2","v2"),
                       stringsAsFactors=FALSE)

HMEC240L <- data.frame(cellLine=rep(c("HMEC240L"), 2),
                       k=135L,
                       ss=c("SS1","SS4"),
                       drug="none",
                       analysisVersion="av1.7",
                       rawDataVersion=c("v2","v2"),
                       stringsAsFactors=FALSE)

HMEC122L <- data.frame(cellLine=rep(c("HMEC122L"), 2),
                       k=135L,
                       ss=c("SS1","SS4"),
                       drug="none",
                       analysisVersion="av1.7",
                       rawDataVersion=c("v2","v2"),
                       stringsAsFactors=FALSE)

ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf,HMEC240L,HMEC122L)


MLDDataSet <- data.frame(datasetName=c("MCF10A_Neratinib"),
                         cellLine=c("MCF10A"),
                         ss=c("SSF"),
                         drug=c("Neratinib"),
                         analysisVersion="av1.7",
                         rawDataVersion="v2",
                         stringsAsFactors=FALSE)
startTime <- Sys.time()
x <- MLDDataSet
datasetName <- x[["datasetName"]]
cellLine <- unique(x[["cellLine"]])
analysisVersion <- unique(x[["analysisVersion"]])
k <- unique(x[["k"]])

#RUV and loess normalize the common DAPI signals
l3R <- preprocessCommonSignals1x(x=x, k=k, verbose=TRUE)

#Bind the staining set specific values, leaving NAs
l3C <- rbindlist(lapply(unique(l3R$StainingSet), function(ss){
  dt <- fread(unique(paste0("MEP_LINCS/AnnotatedData/",x[["cellLine"]],"_",ss,"_Level3.txt")), showProgress = FALSE)
  dt <- addOmeroIDs(dt)
  return(dt)
}),use.names=TRUE,fill=TRUE)

#Write the normalized replicate signal data to disk
cat("Writing level 3 common signal data to disk\n")
write.table(format(l3R, digits = 4, trim=TRUE), paste0("MEP_LINCS/AnnotatedData/", cellLine,"_SSR_","Level3.txt"), sep = "\t",row.names = FALSE, quote=FALSE)

cat("Writing level 3 combined data to disk\n")
write.table(format(l3C, digits = 4, trim=TRUE), paste0("MEP_LINCS/AnnotatedData/", cellLine,"_SSC_","Level3.txt"), sep = "\t",row.names = FALSE, quote=FALSE)

#Median summarize the l3 replicates, 
l3RCS <- level4CommonSignals(l3R)

#mean summarize replicates
l4R <- meanSummarizel4(l3RCS)

l3RMetadata <- unique(l3R[,grep("MEP|Lx|PinDiameter",colnames(l3R),value=TRUE), with=FALSE])
l3RMetadata$MEP <- gsub("FBS_P.","FBS",l3RMetadata$MEP)
l3RMetadata <- unique(l3RMetadata)
#Merge back in the LINCS IDs and pin diameter
l4R <- merge(l4R,l3RMetadata,by="MEP")

#Read the level 4 staining set data and remove the common signals
level4DataL <- lapply(unique(l3R$StainingSet), function(ss){
  dt <-fread(unique(paste0("MEP_LINCS/AnnotatedData/",x[["cellLine"]],"_",ss,"_Level4.txt")), showProgress = FALSE)
  dt <- dt[,grep("MEP|MitoTracker|EdU|Fibrillarin|KRT",colnames(dt), value=TRUE), with=FALSE]
  return(dt)
})

#Merge the 2 or 3 staining sets by MEP
if(length(level4DataL)==2){
  level4Data <- merge(level4DataL[[1]],
                      level4DataL[[2]],
                      by = "MEP")
} else if(length(level4DataL)==3) {
  level4Data <- merge(merge(level4DataL[[1]],
                            level4DataL[[2]],
                            by = "MEP"),
                      level4DataL[[3]],
                      by = "MEP")
} else(stop("There must be 2 or 3 staining sets in the level 4 data"))

#Combine the common signal replicates with the staining set signals
l4C <- merge(l4R, level4Data, by = "MEP")

cat("Writing level 4 file to disk\n")
write.table(format(l4R, digits = 4, trim=TRUE), paste0("MEP_LINCS/AnnotatedData/", cellLine,"_SSR_","Level4.txt"), sep = "\t",row.names = FALSE, quote=FALSE)

write.table(format(l4C, digits = 4, trim=TRUE), paste0("MEP_LINCS/AnnotatedData/", cellLine,"_SSC_","Level4.txt"), sep = "\t",row.names = FALSE, quote=FALSE)

cat("Running time:",Sys.time()-startTime)
