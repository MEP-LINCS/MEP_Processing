#MEP-LINCS Analysis functions
#Mark Dane 9/2015

#Preprocessing Functions

#Functions to create or expose in MEMA
calcTheta <- function(x,y) {
  z <- x + 1i * y
  res <- 90 - Arg(z) / pi * 180
  res %% 360
}

spotCellDensities<- function (spot, radius = (max(spot$Cells_Location_Center_X) - min(spot$Cells_Location_Center_X))/5) 
{
  distMatrix <- as.matrix(dist(spot[, list(Cells_Location_Center_X, Cells_Location_Center_Y)]))
  count <- apply(distMatrix, 2, function(x) {
    sum(x <= radius) - 1
  })
  cellDensity <- count/(pi * radius^2)
  return(cellDensity)
}

cellNeighbors<- function (spot, radius = (max(spot$Nuclei_CP_AreaShape_Center_X) - min(spot$Nuclei_CP_AreaShape_Center_X))/5) 
{
  distMatrix <- as.matrix(dist(spot[, list(Nuclei_CP_AreaShape_Center_X, Nuclei_CP_AreaShape_Center_Y)]))
  count <- apply(distMatrix, 2, function(x) {
    sum(x <= radius) - 1
  })
  return(count)
}

medianNorm <- function(DT, value){
  normedValues <- DT[, value, with = FALSE]/median(unlist(DT[, value, with = FALSE]), na.rm=TRUE)
}

normDataset <- function(dt){
  parmNormedList <- lapply(grep("_CP_|_QI_|_PA_|SpotCellCount|Lineage",colnames(dt),value = TRUE), function(parm){
    dt <- dt[,paste0(parm,"_MedNorm") := normWellsWithinPlate(.SD, value=parm, baseECM = ".*",baseGF = "FBS"), by="Barcode"]
    #     parmNormed <- pcDT[,normWellsWithinPlate(.SD, value=parm, baseECM = ".*",baseGF = "HighSerum"), by="Barcode"]
    #     parmNormed <- parmNormed[,Barcode := NULL]
    return(dt)
  })
  return(parmNormedList[[length(parmNormedList)]])
}

calcGroupRatios <- function(x,group,signal){
  #browser()
  medianInGroup <- median(x[[signal]][x[[group]]], na.rm=TRUE)
  medianOutGroup <- median(x[[signal]][!x[[group]]], na.rm=TRUE)
  return(medianInGroup/medianOutGroup)
} 

scaleToMedians <- function(x, normBase){
  #browser()
  if(!length(x) == length(normBase)) stop("vector to be normalized and base must be the same length")
  xn <- as.numeric(x)/normBase
  return(xn)
}

create8WellPseudoImage <- function(DT, pr, prDisplay){
  highThresh = .998
  #move outliers to maximum displayed value
  DT[[pr]][DT[[pr]]>=quantile(DT[[pr]],probs = highThresh,na.rm=TRUE)] <- as.integer(quantile(DT[[pr]],probs = highThresh,na.rm=TRUE))
  p <- ggplot(DT, aes_string(x="ArrayColumn", y="ArrayRow",colour=pr))+
    geom_point(size=1)+
    scale_y_reverse()+   scale_x_continuous(breaks= c(min(DT$ArrayColumn),round(mean(c(min(DT$ArrayColumn),max(DT$ArrayColumn)))),max(DT$ArrayColumn)))+
    scale_colour_gradient(low = "white", high = "red")+
    guides(colour = guide_legend(prDisplay, keywidth = .5, keyheight = .5))+
    ggtitle(paste("\n\n",prDisplay,"for",unique(DT$CellLine), "cells in plate",unique(DT$Barcode)))+
    xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(.8)),axis.title.x = element_text(size=rel(.5)),plot.title = element_text(size = rel(.5)),legend.text=element_text(size = rel(.4)),legend.title=element_text(size = rel(.3)))+
    facet_wrap(~Well, ncol=4)
}

create8WellHistograms <- function(DT, pr, prDisplay, binwidth = diff(quantile(DT[[pr]],probs = c(0,.98),na.rm=TRUE))/50, upperProb = .99, ncol = 4) {
  
  p <- ggplot(DT, aes_string(x=pr))+
    geom_histogram(binwidth = binwidth)+
    scale_x_continuous(limits = quantile(DT[[pr]],probs = c(0,upperProb),na.rm=TRUE))+
    ggtitle(paste("\n\n",prDisplay,"in",unique(DT$CellLine), "cells in plate",unique(DT$Barcode)))+
    xlab(prDisplay)+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(.8)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), plot.title = element_text(size = rel(.5)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.3)))+
    facet_wrap(~Well, ncol=ncol)
}

ubHeatmap <- function(DT, title = NULL, cols = plateCol, activeThresh = .95) {
  #browser()
  DT$Barcode <- as.factor(DT$Barcode)
  #Get medians of high serum numeric features
  fvDTHS <- DT[grepl("FBS", DT$MEP)]
  hsMedians <- data.frame(t(as.matrix(apply(fvDTHS[,grep("MEP|Barcode", colnames(fvDTHS),invert=TRUE),with=FALSE],2,median, na.rm = TRUE))),MEP="FBS", Barcode = NA, stringsAsFactors = FALSE)
  #Replace all FBS rows with one row of medians as the last row
  DT<- rbind(DT[!grepl("FBS", DT$MEP)],hsMedians)
#   
#   #Normalize the feature vectors to the medians
#   #to equally weight all features
#   fvDTNorm <- rbindlist(apply(DT[,grep("MEP|Barcode",colnames(DT),invert=TRUE), with=FALSE], 1, scaleToMedians, normBase = hsMedians[,grep("MEP|Barcode",colnames(hsMedians),invert=TRUE)]))
  
  #Calculate the dist matrix with euclidean method
  dmm <- as.matrix(dist(DT[,grep("MEP|Barcode",colnames(DT),invert=TRUE), with=FALSE]), labels=TRUE)
  #Extract the distance to the high serum medians
  distHS <- dmm[which(DT$MEP == "FBS"),]
  #Name the distance values
  names(distHS) <- DT$MEP
  
  #Select the most active by distance from the median control fv
  dmmThresh <- quantile(distHS,probs = activeThresh)
  activeMEPs <- distHS[distHS>dmmThresh]
  #browser()
  #Create an active MEP subset matrix of the normalized data
  activeFV <- DT[DT$MEP %in% names(activeMEPs)]
  #Remove MEP and Barcode columns and convert to a matrix
  activeFVM <- as.matrix(activeFV[,grep("MEP|Barcode",colnames(activeFV),invert=TRUE), with=FALSE])
  
  #This assignment of names retains the order after matrix coercion
  rownames(activeFVM) <- names(activeMEPs)
  
  #Cluster the active MEPs, scaling the inputs
  #plot.new()
  heatmap.2(activeFVM,scale="column", col = bluered, trace = "none", cexRow=.5, cexCol=.9, key=FALSE, main = paste(title), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
            lwid=c(1.5,0.2,2.5,2.5),mar=c(20,5), RowSideColors=cols[activeFV$Barcode], colRow = cols[activeFV$Barcode])
  
  return(activeFV)
}

heatmapNoBC <- function(DT, title = NULL, cols = plateCol, activeThresh = .95) {

  activeFV <- DT
  #Remove MEP and Barcode columns and convert to a matrix
  activeFVM <- as.matrix(activeFV[,grep("MEP|Barcode",colnames(activeFV),invert=TRUE), with=FALSE])
  
  #This assignment of names retains the order after matrix coercion
  rownames(activeFVM) <- names(DT$MEP)
  
  #Cluster the active MEPs, scaling the inputs
  #plot.new()
  heatmap.2(activeFVM,scale="column", col = bluered, trace = "none", cexRow=.5, cexCol=.9, key=FALSE, main = paste(title), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
            lwid=c(1.5,0.2,2.5,2.5),mar=c(20,5), RowSideColors=cols[activeFV$Barcode], colRow = cols[activeFV$Barcode])
  
  return(activeFV)
}

findPerimeterCell <- function(x){
  #browser()
  if(!nrow(x)==0){
    perimeterLogicals <- vector(length=nrow(x))
    perimeterLogicals[which.max(x$Nuclei_PA_AreaShape_Center_R)] <- TRUE
  }
  return(perimeterLogicals)
}

labelOuterCells <- function(x, thresh=.75){
  #browser()
  outerLogicals <- NULL
  if(!length(x)==0){
    outerLogicals <- x>quantile(x,probs = thresh, na.rm=TRUE)
  }
  return(outerLogicals)
}


kmeansDNACluster <- function (x, centers = 2) 
{
  #browser()
  x <- data.frame(x)
  xkmeans <- kmeans(x, centers = centers)
  #Swap cluster IDs to make sure cluster 2 has higher values
  if(centers==2){
    if(xkmeans$centers[1] > xkmeans$centers[2]){
      tmp <- xkmeans$cluster == 1
      xkmeans$cluster[xkmeans$cluster == 2] <- 1L
      xkmeans$cluster[tmp] <- 2L
    }
  }
  return(xkmeans$cluster)
}

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

simplifyLigandAnnotID <- function(ligand,annotIDs){
  if(length(annotIDs)){
  splits <- strsplit2(annotIDs, split = "_")
  ligands <- paste(ligand,splits[,ncol(splits)], sep = "_")
  } else ligands <- annotIDs
  return(ligands)
}

normRZSWellsWithinPlate <- function(DT, value, baseECM, baseGF) {
  if(!c("ECMpAnnotID") %in% colnames(DT)) stop(paste("DT must contain a ECMpAnnotID column."))
  if(!c(value) %in% colnames(DT)) stop(paste("DT must contain a", value, "column."))
  if("LigandAnnotID" %in% colnames(DT)){
    valueMedian <- median(unlist(DT[(grepl(baseECM, DT$ECMpAnnotID)  & grepl(baseGF,DT$LigandAnnotID)),value, with=FALSE]), na.rm = TRUE)
    valueMAD <- mad(unlist(DT[(grepl(baseECM, DT$ECMpAnnotID)  & grepl(baseGF,DT$LigandAnnotID)),value, with=FALSE]), na.rm = TRUE)
  } else if (c("Growth.Factors") %in% colnames(DT)) {
    valueMedian <- median(unlist(DT[(grepl(baseECM, DT$ECMpAnnotID)  & grepl(baseGF,DT$Growth.Factors)),value, with=FALSE]), na.rm = TRUE)
    valueMAD <- mad(unlist(DT[(grepl(baseECM, DT$ECMpAnnotID)  & grepl(baseGF,DT$Growth.Factors)),value, with=FALSE]), na.rm = TRUE)
  } else stop (paste("DT must contain a Growth.Factors or LigandAnnotID column."))
  normedValues <- (DT[,value,with=FALSE]-valueMedian)/valueMAD
  return(normedValues)
}

normRZSDataset <- function(dt){
  parmNormedList <- lapply(grep("_CP_|_QI_|_PA_|SpotCellCount|Lineage",colnames(dt),value = TRUE), function(parm){
    dt <- dt[,paste0(parm,"_RZSNorm") := normRZSWellsWithinPlate(.SD, value=parm, baseECM = ".*",baseGF = "FBS"), by="Barcode"]
    return(dt)
  })
  return(parmNormedList[[length(parmNormedList)]])
}


calc2NProportion <- function(x){
  #x numeric vector of cycle states with values of 1 for 2n and 2 for 4N
  #return proportion of cells in x that are in 2N
  if(!length(x)) stop("Calculating 2N/4N proportion on an empty group")
  if(sum(x==1)) {
    proportion2N <- sum(x==1)/length(x)
  } else proportion2N <- 0
  return(proportion2N)
}
