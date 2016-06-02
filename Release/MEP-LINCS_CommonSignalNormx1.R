library("rmarkdown")
source("MEP_LINCS/Release/MEPLINCSFunctions.R")

PC3df <- data.frame(cellLine=rep(c("PC3"), 3),
                    ss=c("SS1", "SS2","SS3"),
                    drug="none",
                    analysisVersion="av1.6",
                    rawDataVersion=c("v2","v2.1","v2.1"),
                    stringsAsFactors=FALSE)

MCF7df <- data.frame(cellLine=rep(c("MCF7"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug="none",
                     analysisVersion=c("av1.6"),
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

YAPCdf <- data.frame(cellLine=rep(c("YAPC"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     drug="none",
                     analysisVersion="av1.6",
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(cellLine=rep(c("MCF10A"), 3),
                       ss=c("SS1","SS2","SS3"),
                       drug="none",
                       analysisVersion="av1.6",
                       rawDataVersion=c("v2","v2","v2"),
                       stringsAsFactors=FALSE)

ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf)

x <- ssDatasets[c(10:12),]
cellLine <- unique(x[["cellLine"]])
analysisVersion <- unique(x[["analysisVersion"]])

l3n <- preprocessCommonSignals1x(x=x, k=135L, verbose=TRUE)

l4n <- level4CommonSignals(l3n)

#WriteData
cat("Writing level 3 file to disk\n")
write.table(format(l3n, digits = 4, trim=TRUE), paste0("MEP_LINCS/AnnotatedData/", cellLine,"_CS_",unique(x[["drug"]]),"_",analysisVersion,"_","Level3.txt"), sep = "\t",row.names = FALSE, quote=FALSE)

cat("Writing level 4 file to disk\n")
write.table(format(l4n, digits = 4, trim=TRUE),paste0("MEP_LINCS/AnnotatedData/", cellLine,"_CS_",unique(x[["drug"]]),"_",analysisVersion,"_","Level4.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
