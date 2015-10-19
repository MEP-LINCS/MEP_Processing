#MEP-LINCS Analysis functions
#Mark Dane 10/2015

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
  browser()
  DT$Barcode <- as.factor(DT$Barcode)
  #Get medians of high serum numeric features
  fvDTHS <- DT[grepl("FBS", DT$MEP)]
  hsMedians <- data.frame(t(as.matrix(apply(fvDTHS[,grep("MEP|Barcode", colnames(fvDTHS),invert=TRUE),with=FALSE],2,median, na.rm = TRUE))),MEP="FBS", Barcode = NA, stringsAsFactors = FALSE)
  #Replace all FBS rows with one row of medians as the last row
  DT<- rbind(DT[!grepl("FBS", DT$MEP)],hsMedians)

  #Calculate the dist matrix with euclidean method
  dmm <- as.matrix(dist(DT[,grep("MEP|Barcode",colnames(DT),invert=TRUE), with=FALSE]), labels=TRUE)
  #Extract the distance to the high serum medians
  distHS <- dmm[which(DT$MEP == "FBS"),]
  #Name the distance values
  names(distHS) <- DT$MEP
  
  #Select the most active by distance from the median control fv
  dmmThresh <- quantile(distHS,probs = activeThresh, na.rm = TRUE)
  activeMEPs <- distHS[distHS>dmmThresh]
  #browser()
  #Create an active MEP subset matrix of the normalized data
  activeFV <- DT[DT$MEP %in% names(activeMEPs)]
  #Remove MEP and Barcode columns and convert to a matrix
  activeFVM <- as.matrix(activeFV[,grep("MEP|Barcode",colnames(activeFV),invert=TRUE), with=FALSE])
  
  #This assignment of names retains the order after matrix coercion
  rownames(activeFVM) <- activeFV$MEP
  
  #Cluster the active MEPs, scaling the inputs
  heatmap.2(activeFVM,scale="column", col = bluered, trace = "none", cexRow=.5, cexCol=.9, key=FALSE, main = paste(title), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1,5.0), lwid=c(.5,0.2,2.5,1.5),
            mar=c(25,1),
            RowSideColors=cols[activeFV$Barcode], colRow = cols[activeFV$Barcode], na.rm = TRUE)
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
    #Find the last non-empty substring
    us <- apply(splits,1,function(x){
      if(x[length(x)]==""){
        if (length(x)< 2) stop("There are not enough substrings in a ligandANnotID")
        u <- x[length(x)-1]
      } else{
        u <- x[length(x)]
      }
      return(u)
    })
    ligands <- paste(ligand,us, sep = "_")
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
  if(sum(grepl("[^12]",x)))stop("Invalid cycle state passed to calc2NProportion")
  if(sum(x==1)) {
    proportion2N <- sum(x==1)/length(x)
  } else proportion2N <- 0
  return(proportion2N)
}

plotTotalDAPI <- function(l1, barcodes){
  for (barcode in barcodes){
    mDT <- l1[l1$Barcode == barcode]
    mDT <- mDT[mDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi > quantile(mDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi, probs=.01, na.rm=TRUE) & mDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi < quantile(mDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi,probs=.98, na.rm=TRUE)]
    mDT <- mDT[,DNAThresh := min(Nuclei_CP_Intensity_IntegratedIntensity_Dapi[Nuclei_PA_Cycle_State==2]), by="Well"]
    p <- ggplot(mDT, aes(x=Nuclei_CP_Intensity_IntegratedIntensity_Dapi))+geom_bar(binwidth = 2)+
      geom_vline(data = mDT, aes(xintercept = DNAThresh), colour = "blue")+
      facet_wrap(~Ligand, nrow=2, scales="free_x")+
      #xlim(0,quantile(mDT$TotalIntensityDAPI,probs=.98, na.rm=TRUE))+
      ggtitle(paste("\n\n","Total DAPI Signal,",barcode))+
      ylab("Count")+xlab("Total Intensity DAPI")+
      theme(strip.text = element_text(size = 5))
    suppressWarnings(print(p))
  }
}


plotSCCHeatmapsQAHistograms <- function(l3, barcodes){
  for (barcode in barcodes){
    DT <-l3[l3$Barcode==barcode,]
    #Remove the fiducial entries
    setkey(DT,ECMp)
    DT <- DT[!"fiducial"]
    DT <- DT[!"blank"]
    
    p <- create8WellPseudoImage(DT, pr = "Spot_PA_SpotCellCount",prDisplay = "Spot Cell Count")
    suppressWarnings(print(p))
    
    wellScores <- unique(DT[,list(Well, QAScore=sprintf("%.2f",QAScore))])
    
    p <- ggplot(DT, aes(x=Spot_PA_LoessSCC))+
      geom_histogram(binwidth=.04)+
      geom_vline(xintercept=lthresh, colour="blue")+
      geom_text(data=wellScores, aes(label=paste0("QA\n",QAScore)), x = 2, y = 30, size = rel(3), colour="red")+
      ggtitle(paste("\n\n\n\n","QA on Loess Model of Spot Cell Count for",unique(DT$CellLine), "cells in plate",unique(DT$Barcode)))+xlab("Normalized Spot Cell Count")+xlim(0,3)+
      theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(.5)), axis.title.x = element_text( size=rel(.5)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), axis.title.y = element_text( size=rel(.5)), plot.title = element_text(size = rel(.5)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.3)))+
      facet_wrap(~Well, ncol=4)
    suppressWarnings(print(p))
  }
}

RZScore <- function(x){
  xMedian <- median(x, na.rm=TRUE)
  xMad <-mad(x, na.rm=TRUE)
  if(xMad == 0){ zscores <- NA
  } else zscores <- (x-xMedian)/xMad
  return(zscores)
}

se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

filterl4 <- function(dt,lowQALigands){
  #Remove failed QA wells
  l4QA<- dt[!dt$Ligand %in% lowQALigands]
  
  setkey(l4QA, "ECMp")
  l4QA <- l4QA[!"blank"]
  l4QA <- l4QA[!"fiducial"]
  l4QA <- l4QA[,grep("Center_X|Center_Y|Center_Theta",colnames(l4QA),value = TRUE, invert = TRUE), with = FALSE]
  
  #Define features for clustering
  fv <- paste("^Barcode","MEP",
              "Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker_RZSNorm",
              "Nuclei_CP_AreaShape_Area_RZSNorm",
              "Nuclei_CP_AreaShape_Eccentricity_RZSNorm",
              "Nuclei_CP_AreaShape_Perimeter_RZSNorm",
              "Nuclei_CP_Intensity_MedianIntensity_Dapi_RZSNorm",
              "Spot_PA_SpotCellCount_RZSNorm",
              "Nuclei_PA_AreaShape_Neighbors_RZSNorm",
              "Nuclei_PA_Cycle_2NProportion_RZSNorm$",
              "Nuclei_CP_Intensity_MedianIntensity_Edu_RZSNorm",
              "Nuclei_PA_Gated_EduPositiveProportion_RZSNorm_RZSNorm",
              "Cytoplasm_CP_Intensity_MedianIntensity_KRT19_RZSNorm",
              "Cytoplasm_CP_Intensity_MedianIntensity_KRT5_RZSNorm",
              "Cytoplasm_PA_Intensity_LineageRatio_RZSNorm$",
              sep="$|^")
  
  fv <- grep(fv, colnames(l4QA), value = TRUE)
  #Create numeric feature vectors datatable
  fvDT <- l4QA[,fv,with = FALSE]
  return(fvDT)
}

plotSCCRobustZScores <- function(dt, thresh = 3){
  #Filter our FBS MEPs then plot spot cell count robust Z scores
  #browser()
  dt <- dt[!grepl("FBS",dt$MEP)]
  p <- ggplot(dt, aes(x=Spot_PA_SpotCellCount_RZSNorm_RobustZ))+geom_bar(binwidth = .1)+
    geom_vline(xintercept = c(-thresh,thresh), colour = "blue")+
    ggtitle(paste("\n\n","MEP Normalized Spot Cell Count Robust Z Scores Distribution"))+
    ylab("Count")+xlab("Normalized Spot Cell Count Robust Z Scores")+
    theme(strip.text = element_text(size = 5))
  suppressWarnings(print(p))
  
}


combineSSs <- function(SSs){
  #browser()
  l4List <- lapply(SSs, function(ss){
    l4 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine,"_",ss,"_Level4.txt"), showProgress = FALSE)
    setkey(l4,"ECMp")
    l4 <- l4[!"fiducial"]
    l4 <- l4[!"blank"]
    l4$SS <- ss
    return(l4)
  })
  
  l4SS1 <- l4List[[1]]
  l4SS2 <- l4List[[2]]
  l4SS3 <- l4List[[3]]
  
  setkey(l4SS1,"LigandAnnotID","ECMpAnnotID")
  setkey(l4SS2,"LigandAnnotID","ECMpAnnotID")
  setkey(l4SS3,"LigandAnnotID","ECMpAnnotID")
  
  #Bind the data
  DT <- data.table(l4SS1, l4SS2, l4SS3, check.names = TRUE)
}


integrateSSs <- function(SSs, cellLine = "PC3"){
  #browser()
  l4List <- lapply(SSs, function(ss){
    l4 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine,"_",ss,"_Level4.txt"), showProgress = FALSE)
    setkey(l4,"ECMp")
    l4 <- l4[!"fiducial"]
    l4 <- l4[!"blank"]
    setkey(l4, "MEP")
    l4$SS <- ss
    return(l4)
  })
  
  l4SS1 <- l4List[[1]]
  l4SS2 <- l4List[[2]]
  l4SS3 <- l4List[[3]]
  
  setkey(l4SS1,"LigandAnnotID","ECMpAnnotID")
  setkey(l4SS2,"LigandAnnotID","ECMpAnnotID")
  setkey(l4SS3,"LigandAnnotID","ECMpAnnotID")
  
  #Bind the data using the common MEPs
  DT <- data.table(l4SS1, l4SS2, l4SS3, check.names = TRUE)
  
  #Median summarize the FBS rows
  setkey(DT,"MEP")
  DTFBS <- DT[grepl("FBS", DT$MEP)]
  #Get the medians of each numeric parameter
  parms <- colnames(DTFBS)[unlist(lapply(DTFBS,class)) %in% c("numeric","integer")]
  FBSMedians <- data.frame(t(as.matrix(apply(DTFBS[, parms,with=FALSE],2,median))),MEP="FBS", stringsAsFactors = FALSE)
  
  #Merge the metadata back in with the data
  metadata <- colnames(DTFBS)[!unlist(lapply(DTFBS,class)) %in% c("numeric","integer")]
  FBSMetadata <- DTFBS[, metadata, with = FALSE]
  FBSMetadata$MEP <- "FBS"
  FBSMetadata$MEP.1 <- "FBS"
  FBSMetadata$MEP.2 <- "FBS"
  FBSMetadata$ECMp <- NA
  FBSMetadata$ECMp.1 <- NA
  FBSMetadata$ECMp.2 <- NA
  FBSMetadata$ECMpAnnotID <- NA
  FBSMetadata$ECMpAnnotID.1 <- NA
  FBSMetadata$ECMpAnnotID.2 <- NA
  FBSMetadata$Well <- NA
  FBSMetadata$Well.1 <- NA
  FBSMetadata$Well.2 <- NA
  FBSMetadata$Barcode <- NA
  FBSMetadata$Barcode.1 <- NA
  FBSMetadata$Barcode.2 <- NA
  
  FBSMetadata <- unique(FBSMetadata)
  
  FBSMisOrdered <- cbind(FBSMetadata[,MEP:=NULL],FBSMedians)
  
  #Replace all FBS rows with one row of medians as the last row
  DT1FBS<- rbind(DT[!grepl("FBS", DT$MEP)],FBSMisOrdered,use.names=TRUE)
  
}

createl3 <- function(cDT, lthresh = lthresh){
  #Summarize cell data to medians of the spot parameters
  parameterNames<-grep(pattern="(Children|_CP_|_PA_|Barcode|^Spot$|^Well$)",x=names(cDT),value=TRUE)
  
  #Remove any spot-normalized and cell level parameters
  parameterNames <- grep("SpotNorm|^Nuclei_PA_Gated_EduPositive$|^Nuclei_PA_Gated_EduPositive_RZSNorm$",parameterNames,value=TRUE,invert=TRUE)
  
  #Remove any raw parameters
  parameterNames <- grep("Barcode|^Spot$|^Well$|Norm|Nuclei_CP_Intensity_MedianIntensity_Dapi$|Cytoplasm_CP_Intensity_MedianIntensity_Actin$|Cytoplasm_CP_Intensity_MedianIntensity_CellMask$|Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker$|Nuclei_CP_Intensity_MedianIntensity_H3$|Nuclei_CP_Intensity_MedianIntensity_Fibrillarin$|Nuclei_CP_Intensity_MedianIntensity_Edu$|Cytoplasm_CP_Intensity_MedianIntensity_KRT5$|Cytoplasm_CP_Intensity_MedianIntensity_KRT19$|Spot_PA_SpotCellCount$", parameterNames, value = TRUE)
  
  cDTParameters<-cDT[,parameterNames,with=FALSE]
  
  slDT<-cDTParameters[,lapply(.SD,numericMedian),keyby="Barcode,Well,Spot"]
  slDTse <- cDTParameters[,lapply(.SD,se),keyby="Barcode,Well,Spot"]
  
  #Add _SE to the standard error column names
  setnames(slDTse, grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE),"_SE"))
  
  #Merge back in the spot and well metadata
  #TODO: Convert the logic to not name the metadata
  metadataNames <- grep("(Row|Column|PrintOrder|Block|^ID$|Array|CellLine|Ligand|Endpoint|ECMp|MEP|Barcode|^Well$|^Spot$)", x=colnames(cDT), value=TRUE)
  setkey(cDT,Barcode, Well,Spot)
  mDT <- cDT[,metadataNames,keyby="Barcode,Well,Spot", with=FALSE]
  slDT <- mDT[slDT, mult="first"]
  #Merge in the standard errr values
  slDT <- slDTse[slDT]
  #Add a count of replicates
  slDT <- slDT[,Spot_PA_ReplicateCount := .N,by="LigandAnnotID,ECMpAnnotID"]
  
  #Add the loess model of the SpotCellCount on a per well basis
  slDT <- slDT[,Spot_PA_LoessSCC := loessModel(.SD, value="Spot_PA_SpotCellCount", span=.5), by="Barcode,Well"]
  
  #Add well level QA Scores to spot level data
  slDT <- slDT[,QAScore := calcQAScore(.SD, threshold=lthresh, maxNrSpot = max(cDT$ArrayRow)*max(cDT$ArrayColumn),value="Spot_PA_LoessSCC"),by="Barcode,Well"]
}#End of create l3

createl4 <- function(l3){
  #Add a count of replicates
  l3 <- l3[,Spot_PA_ReplicateCount := .N,by="LigandAnnotID,ECMpAnnotID"]
  l4Names<-grep("Norm|LigandAnnotID|ECMpAnnotID|Barcode|Spot_PA_SpotCellCount$|Spot_PA_ReplicateCount$", x=names(l3),value=TRUE)
  #remove the _SE values
  l4Names <- grep("_SE",l4Names, value = TRUE, invert = TRUE)
  l4Keep<-l3[,l4Names,with=FALSE]
  l4DT<-l4Keep[,lapply(.SD,numericMedian),keyby="LigandAnnotID,ECMpAnnotID,Barcode"]
  
  l4DTse <- l4Keep[,lapply(.SD,se),keyby="LigandAnnotID,ECMpAnnotID,Barcode"]
  #Add _SE to the standard error column names
  setnames(l4DTse, grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(l4DTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(l4DTse), value = TRUE, invert = TRUE),"_SE"))
  
  #Merge back in the replicate metadata
  mDT <- l3[,list(Well,CellLine,Ligand,Endpoint488,Endpoint555,Endpoint647,EndpointDAPI,ECMp,MEP),keyby="LigandAnnotID,ECMpAnnotID,Barcode"]
  l4DT <- mDT[l4DT, mult="first"]
  l4DT <- l4DTse[l4DT]
  
  return(l4DT)
}#End of createl4


preprocessMEPLINCS <- function(ss, cellLine, limitBarcodes = 8){
  source("MEPLINCSFunctions.R")
  library("limma")#read GAL file and strsplit2
  library("MEMA")#merge, annotate and normalize functions
  library("data.table")#fast file reads, data merges and subsetting
  library("parallel")#use multiple cores for faster processing
  #Select a staining set
  #ss <- "SS3"
  #Select a CellLine
  #cellLine <- "YAPC"
  #select analysis version
  analysisVersion <- "v1"
  
  #Rules-based classifier thresholds for perimeter cells
  neighborsThresh <- 0.4 #Gates sparse cells on a spot
  wedgeAngs <- 20 #Size in degrees of spot wedges used in perimeter gating
  outerThresh <- 0.5 #Defines out cells used in perimeter gating
  neighborhoodNucleiRadii <- 7 #Defines the neighborhood annulus
  
  #Filter out debris based on nuclear area
  nuclearAreaThresh <- 50
  nuclearAreaHiThresh <- 4000
  
  #Only process a curated set of the data
  curatedOnly <- TRUE
  curatedCols <- "ImageNumber|ObjectNumber|_Area$|_Eccentricity|_Perimeter|_MedianIntensity_|_IntegratedIntensity_|_Center_|_PA_"
  
  #Flag to control files updates
  writeFiles <- TRUE
  
  #Do not normalized to Spot level
  normToSpot <- TRUE
  
  #QA flags are used to enable analyses that require minimum cell and
  #replicate counts
  
  #Set a threshold for the lowSpotCellCount flag
  lowSpotCellCountThreshold <- 5
  
  #Set a threshold for the lowRegionCellCount flag
  lowRegionCellCountThreshold <- .4
  
  #Set a threshold for the loess well level QA Scores
  lthresh <- 0.6
  
  #Set a threshold for lowWellQA flag
  lowWellQAThreshold <- .7
  
  #Set a threshold for the lowSpotReplicates flag
  lowReplicateCount <- 3
  
  ##Summary
  # This script prepares cell-level data and metadata for the MEP LINCs Analysis Pipeline. 
  # 
  # In the code, the variable ss determines which staining set (SS1, SS2 or SS3) to merge and the variable cellLine determines the cell line (PC3,MCF7, etc). All .txt data files in the "./RawData" folder will be merged with the well (xlsx) and log (XML) data from the "./Metadata" folder.
  # 
  # The merging assumes that the actual, physical B row wells (B01-B04) have been printed upside-down. That is, rotated 180 degrees resulting in the spot 1, 1 being in the lower right corner instead of the upper left corner. The metadata is matched to the actual printed orientation.
  
  # Read and clean spotmetadata
  
  #Read in the spot metadata from the gal file
  smd <- readSpotMetadata(paste0("./",cellLine,"/",ss,"/Metadata/20150515_LI8X001_v1.2.gal"))
  #Relabel the column Name to ECMpAnnotID
  setnames(smd, "Name", "ECMpAnnotID")
  
  #Add the print order and deposition number to the metadata
  ldf <- readLogData(paste0("./",cellLine,"/",ss,"/Metadata/20150512-112336.xml"))
  spotMetadata <- merge(smd,ldf, all=TRUE)
  setkey(spotMetadata,Spot)
  #Make a rotated version of the spot metadata to match the print orientation
  spotMetadata180 <- rotateMetadata(spotMetadata)
  ARowMetadata <- data.table(spotMetadata,Well=rep(c("A01", "A02","A03","A04"),each=nrow(spotMetadata)))
  BRowMetadata <- data.table(spotMetadata180,Well=rep(c("B01", "B02","B03","B04"),each=nrow(spotMetadata180)))
  
  # The well metadata describes the cell line, ligands and staining endpoints that are all added on a per well basis. There is one mutlisheet .xlsx file for each plate. Each filename is the plate's barcode.
  
  
  # The raw data from all wells in all plates in the dataset are read in and merged with their spot and well metadata. The number of nuclei at each spot are counted and a loess model of the spot cell count is added. Then all intensity values are normalized through dividing them by the median intensity value of the control well in the same plate.
  # 
  # Next, the data is filtered to remove objects with a nuclear area less than nuclearAreaThresh pixels or more than nuclearAreaHiThresh pixels.
  
  #merge_normalize_QA, echo=FALSE}
  #The next steps are to bring in the well metadata, the print order and the CP data
  
  cellDataFiles <- dir(paste0("./",cellLine,"/", ss,"/RawData/",analysisVersion),full.names = TRUE)
  splits <- strsplit2(strsplit2(cellDataFiles,split = "_")[,1],"/")
  
  if(limitBarcodes) {
    barcodes <- unique(splits[,ncol(splits)])[1:limitBarcodes] 
  } else barcodes <- unique(splits[,ncol(splits)])
  
  expDTList <- mlapply(barcodes, function(barcode){
    #browser()
    plateDataFiles <- grep(barcode,cellDataFiles,value = TRUE)
    wells <- unique(strsplit2(split = "_",plateDataFiles)[,2])
    wellDataList <- lapply(wells,function(well){
      #browser()
      wellDataFiles <- grep(well,plateDataFiles,value = TRUE)
      #imageDataFile <- grep("Image",wellDataFiles,value=TRUE,
      #                     ignore.case = TRUE)
      nucleiDataFile <- grep("Nuclei",wellDataFiles,value=TRUE,
                             ignore.case = TRUE)
      if (ss %in% c("SS1","SS3")){
        cellsDataFile <- grep("Cell",wellDataFiles,value=TRUE,
                              ignore.case = TRUE)
        cytoplasmDataFile <- grep("Cytoplasm",wellDataFiles,value=TRUE,
                                  ignore.case = TRUE)
      }
      #Read in csv data
      #image <- convertColumnNames(fread(imageDataFile))
      #setkey(image,CP_ImageNumber)
      nuclei <- convertColumnNames(fread(nucleiDataFile))
      if (curatedOnly) nuclei <- nuclei[,grep(curatedCols,colnames(nuclei)), with=FALSE]
      setkey(nuclei,CP_ImageNumber,CP_ObjectNumber)
      if (ss %in% c("SS1","SS3")){
        cells <- convertColumnNames(fread(cellsDataFile))
        if (curatedOnly) cells <- cells[,grep(curatedCols,colnames(cells)), with=FALSE]
        setkey(cells,CP_ImageNumber,CP_ObjectNumber)
        cytoplasm <- convertColumnNames(fread(cytoplasmDataFile))
        if (curatedOnly) cytoplasm <- cytoplasm[,grep(curatedCols,colnames(cytoplasm)), with=FALSE]
        setkey(cytoplasm,CP_ImageNumber,CP_ObjectNumber)
      }
      
      #Add the data location as a prefix in the column names
      setnames(nuclei,paste0("Nuclei_",colnames(nuclei)))
      if (ss %in% c("SS1","SS3")){
        setnames(cells,paste0("Cells_",colnames(cells)))
        setnames(cytoplasm,paste0("Cytoplasm_",colnames(cytoplasm)))
      }
      
      #Merge the cells, cytoplasm and nuclei data
      if (ss %in% c("SS1","SS3")){
        setkey(cells,Cells_CP_ImageNumber,Cells_CP_ObjectNumber)
        setkey(cytoplasm,Cytoplasm_CP_ImageNumber,Cytoplasm_CP_ObjectNumber)
        setkey(nuclei,Nuclei_CP_ImageNumber,Nuclei_CP_ObjectNumber)
        
        DT <- cells[cytoplasm[nuclei]]
        setnames(DT,"Cells_CP_ImageNumber","Spot")
        setnames(DT,"Cells_CP_ObjectNumber","ObjectNumber")
      } else {
        DT <- nuclei
        setnames(DT,"Nuclei_CP_ImageNumber","Spot")
        setnames(DT,"Nuclei_CP_ObjectNumber","ObjectNumber")
      }
      
      #Add the well name as a parameter
      DT <- DT[,Well := well]
      
      #Merge the data with its metadata based on the row it's in
      m <- regexpr("[[:alpha:]]",well)
      row <- regmatches(well,m)
      setkey(DT,Spot)
      DT <- switch(row, A = merge(DT,spotMetadata,all=TRUE),
                   B = merge(DT,spotMetadata180,all=TRUE))
      
      return(DT)
    })
    #browser()
    #Create the cell data.table with spot metadata for the plate 
    pcDT <- rbindlist(wellDataList, fill = TRUE)
    #Read the well metadata from a multi-sheet Excel file
    wellMetadata <- data.table(readMetadata(paste0("./",cellLine,"/",
                                                   ss,"/Metadata/",barcode,".xlsx")), key="Well")
    
    #merge well metadata with the data and spot metadata
    pcDT <- merge(pcDT,wellMetadata,by = "Well")
    pcDT <- pcDT[,Barcode := barcode]
    #Count the cells at each spot
    pcDT<-pcDT[,Spot_PA_SpotCellCount := .N,by="Barcode,Well,Spot"]
    
    pcDT <- pcDT[pcDT$Nuclei_CP_AreaShape_Area > nuclearAreaThresh,]
    pcDT <- pcDT[pcDT$Nuclei_CP_AreaShape_Area < nuclearAreaHiThresh,]
    
    #Add the local polar coordinates and Neighbor Count
    pcDT <- pcDT[,Nuclei_PA_Centered_X :=  Nuclei_CP_AreaShape_Center_X-median(Nuclei_CP_AreaShape_Center_X)]
    pcDT <- pcDT[,Nuclei_PA_Centered_Y :=  Nuclei_CP_AreaShape_Center_Y-median(Nuclei_CP_AreaShape_Center_Y)]
    pcDT <- pcDT[, Nuclei_PA_AreaShape_Center_R := sqrt(Nuclei_PA_Centered_X^2 + Nuclei_PA_Centered_Y^2)]
    pcDT <- pcDT[, Nuclei_PA_AreaShape_Center_Theta := calcTheta(Nuclei_PA_Centered_X, Nuclei_PA_Centered_Y)]
    
    
    #Set 2N and 4N DNA status
    pcDT <- pcDT[,Nuclei_PA_Cycle_State := kmeansDNACluster(Nuclei_CP_Intensity_IntegratedIntensity_Dapi), by="Barcode,Well"]
    #Manually reset clusters for poorly classified wells
    #This is based on review of the clusters after a prior run
    pcDT$Nuclei_PA_Cycle_State[pcDT$Barcode=="LI8X00403" & pcDT$Well=="A03" & pcDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi >150] <- 2
    
    pcDT <- pcDT[,Nuclei_PA_Cycle_DNA2NProportion := calc2NProportion(Nuclei_PA_Cycle_State),by="Barcode,Well,Spot"]
    pcDT$Nuclei_PA_Cycle_DNA4NProportion <- 1-pcDT$Nuclei_PA_Cycle_DNA2NProportion
    
    #Add spot level normalizations for selected intensities
    if(normToSpot){
      intensityNamesAll <- grep("_CP_Intensity_Median",colnames(pcDT), value=TRUE)
      intensityNames <- grep("Norm",intensityNamesAll,invert=TRUE,value=TRUE)
      for(intensityName in intensityNames){
        #Median normalize the median intensity at each spot
        pcDT <- pcDT[,paste0(intensityName,"_SpotNorm") := medianNorm(.SD,intensityName),by="Barcode,Well,Spot"]
      }
    }
    
    #Create staining set specific derived parameters
    if(ss %in% c("SS1")){
      
    } else if (ss == "SS2"){
      pcDT <- pcDT[,Nuclei_PA_Gated_EduPositive := kmeansCluster(.SD, value="Nuclei_CP_Intensity_MedianIntensity_Edu", ctrlLigand = "FBS"), by="Barcode"]
      #Calculate the EdU Positive Percent at each spot
      pcDT <- pcDT[,Nuclei_PA_Gated_EduPositiveProportion := sum(Nuclei_PA_Gated_EduPositive)/length(ObjectNumber),by="Barcode,Well,Spot"]
      #Logit transform EduPositiveProportion
      #logit(p) = log[p/(1-p)]
      EdUppImpute <- pcDT$Nuclei_PA_Gated_EduPositiveProportion
      EdUppImpute[EdUppImpute==0] <- .01
      EdUppImpute[EdUppImpute==1] <- .99
      pcDT$Nuclei_PA_Gated_EduPositiveLogit <- log2(EdUppImpute/(1-EdUppImpute))
      
    } else if (ss == "SS3"){
      #Calculate a lineage ratio of luminal/basal or KRT19/KRT5
      pcDT <- pcDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT5]
      
    } else stop("Invalid ss parameter")
    #browser()
    return(pcDT)
  }, mc.cores=detectCores())
  
  cDT <- rbindlist(expDTList, fill = TRUE)
  #TODO delete unwanted columns here such as Euler Number
  densityRadius <- sqrt(median(cDT$Nuclei_CP_AreaShape_Area)/pi)
  
  #Create short display names, then replace where not unique
  cDT$ECMp <- gsub("_.*","",cDT$ECMpAnnotID)
  cDT$Ligand <- gsub("_.*","",cDT$LigandAnnotID)
  
  #Use entire AnnotID for ligands with same uniprot IDs
  cDT$Ligand[grepl("NRG1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "NRG1",annotIDs = cDT$LigandAnnotID[grepl("NRG1",cDT$Ligand)])
  cDT$Ligand[grepl("TGFB1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "TGFB1",annotIDs = cDT$LigandAnnotID[grepl("TGFB1",cDT$Ligand)])
  cDT$Ligand[grepl("CXCL12",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "CXCL12",annotIDs = cDT$LigandAnnotID[grepl("CXCL12",cDT$Ligand)])
  
  cDT$MEP <- paste(cDT$ECMp,cDT$Ligand,sep = "_")
  
  # After merging the metadata with the cell-level data, several types of derived parameters are added. These include:
  #   
  #   The origin of coordinate system is placed at the median X and Y of each spot and the local cartesian and polar coordinates are added to the dataset.
  # 
  # The number of nuclei within three nuclear radii of `r densityRadius ` around each nuclei is counted and stored as a neighbor count parameter. The neighbor count value is thresholded to classify each cell as Sparse or not.The distance from the local origin is used to classify each cell as an OuterCell or not. The Sparse, OutCell and Wedge classifications are used to classify each cell as a Perimeter cell or not. 
  # 
  # For staining set 2, each cell is classified as EdU+ or EdU-. The threshold for EdU+ is based on kmeans threshold of the mean EdU intensity from the control well of each plate.
  # 
  # The intensity values are normalized at each spot so that spot-level variations can be analyzed.
  # 
  # The cell level raw data and metadata is saved as Level 1 data. Normalized values are added to the dataset and saved as Level 2 data.
  
  #Add the local polar coordinates and Neighbor Count
  cDT <- cDT[,Nuclei_PA_Centered_X :=  Nuclei_CP_AreaShape_Center_X-median(Nuclei_CP_AreaShape_Center_X)]
  cDT <- cDT[,Nuclei_PA_Centered_Y :=  Nuclei_CP_AreaShape_Center_Y-median(Nuclei_CP_AreaShape_Center_Y)]
  cDT <- cDT[, Nuclei_PA_AreaShape_Center_R := sqrt(Nuclei_PA_Centered_X^2 + Nuclei_PA_Centered_Y^2)]
  cDT <- cDT[, Nuclei_PA_AreaShape_Center_Theta := calcTheta(Nuclei_PA_Centered_X, Nuclei_PA_Centered_Y)]
  cDT <- cDT[,Nuclei_PA_AreaShape_Neighbors := cellNeighbors(.SD, radius = densityRadius*neighborhoodNucleiRadii), by = "Barcode,Well,Spot"]
  
  #Rules for classifying perimeter cells
  cDT <- cDT[,Spot_PA_Sparse := Nuclei_PA_AreaShape_Neighbors < neighborsThresh]
  
  #Add a local wedge ID to each cell based on conversations with Michel Nederlof
  cDT <- cDT[,Spot_PA_Wedge:=ceiling(Nuclei_PA_AreaShape_Center_Theta/wedgeAngs)]
  
  #Define the perimeter cell if it exists in each wedge
  #Classify cells as outer if they have a radial position greater than a thresh
  cDT <- cDT[,Spot_PA_OuterCell := labelOuterCells(Nuclei_PA_AreaShape_Center_R, thresh=outerThresh),by="Barcode,Well,Spot"]
  
  #Require the cell not be in a sparse region
  denseOuterDT <- cDT[!cDT$Spot_PA_Sparse  & cDT$Spot_PA_OuterCell]
  denseOuterDT <- denseOuterDT[,Spot_PA_Perimeter := findPerimeterCell(.SD) ,by="Barcode,Well,Spot,Spot_PA_Wedge"]
  setkey(cDT,Barcode,Well,Spot,ObjectNumber)
  setkey(denseOuterDT,Barcode,Well,Spot,ObjectNumber)
  cDT <- denseOuterDT[,list(Barcode,Well,Spot,ObjectNumber,Spot_PA_Perimeter)][cDT]
  cDT$Spot_PA_Perimeter[is.na(cDT$Spot_PA_Perimeter)] <- FALSE
  
  #Set 2N and 4N DNA status
  cDT <- cDT[,Nuclei_PA_Cycle_State := kmeansDNACluster(Nuclei_CP_Intensity_IntegratedIntensity_Dapi), by="Barcode,Well"]
  #Manually reset clusters for poorly classified wells
  #This is based on review of the clusters after a prior run
  cDT$Nuclei_PA_Cycle_State[cDT$Barcode=="LI8X00403" & cDT$Well=="A03" & cDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi >150] <- 2
  
  cDT <- cDT[,Nuclei_PA_Cycle_2NProportion := calc2NProportion(Nuclei_PA_Cycle_State),by="Barcode,Well,Spot"]
  cDT$Nuclei_PA_Cycle_4NProportion <- 1-cDT$Nuclei_PA_Cycle_2NProportion
  
  #Add spot level normalizations for selected intensities
  if(normToSpot){
    intensityNamesAll <- grep("_CP_Intensity_Median",colnames(cDT), value=TRUE)
    intensityNames <- grep("Norm",intensityNamesAll,invert=TRUE,value=TRUE)
    for(intensityName in intensityNames){
      #Median normalize the median intensity at each spot
      cDT <- cDT[,paste0(intensityName,"_SpotNorm") := medianNorm(.SD,intensityName),by="Barcode,Well,Spot"]
    }
  }
  
  #Create staining set specific derived parameters
  if(ss %in% c("SS1")){
    
  } else if (ss == "SS2"){
    cDT <- cDT[,Nuclei_PA_Gated_EduPositive := kmeansCluster(.SD, value="Nuclei_CP_Intensity_MedianIntensity_Edu", ctrlLigand = "FBS"), by="Barcode"]
    #Calculate the EdU Positive Percent at each spot
    cDT <- cDT[,Nuclei_PA_Gated_EduPositiveProportion := sum(Nuclei_PA_Gated_EduPositive)/length(ObjectNumber),by="Barcode,Well,Spot"]
    #Logit transform EduPositiveProportion
    #logit(p) = log[p/(1-p)]
    cDT$Nuclei_PA_Gated_EduPositiveLogit <- log2(cDT$Nuclei_PA_Gated_EduPositiveProportion/(1-cDT$Nuclei_PA_Gated_EduPositiveProportion))
    
  } else if (ss == "SS3"){
    #Calculate a lineage ratio of luminal/basal or KRT19/KRT5
    cDT <- cDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT5]
    
  } else stop("Invalid ss parameter")
  
  # Eliminate Variations in the Endpoint metadata
  endpointNames <- grep("End",colnames(cDT), value=TRUE)
  endpointWL <- regmatches(endpointNames,regexpr("[[:digit:]]{3}|DAPI",endpointNames))
  setnames(cDT,endpointNames,paste0("Endpoint",endpointWL))
  
  #Identify parameters that shouldn't be normalized
  normParameters <- grep("Sparse|Wedge|OuterCell|Spot_PA_Perimeter|Nuclei_PA_Cycle_State",colnames(cDT),value=TRUE,invert=TRUE)
  
  #Save the un-normalized parameters to merge in later
  mdKeep <- cDT[,grep("Barcode|^Well$|^Spot$|ObjectNumber|Sparse|Wedge|OuterCell|Spot_PA_Perimeter|Nuclei_PA_Cycle_State",colnames(cDT),value=TRUE), with = FALSE]
  
  #Normalized each feature by dividing by the median of its plate's FBS value well
  #cDT <- normDataset(cDT[,normParameters,with=FALSE])
  #Normalize each feature by subtracting the median of its plate's FBS value
  # and divding by its plates MAD
  cDT <- normRZSDataset(cDT[,normParameters, with = FALSE])
  cDT <- merge(cDT,mdKeep)
  
  #The cell-level data is median summarized to the spot level and coefficients of variations on the replicates are calculated. The spot level data and metadata are saved as Level 3 data.
  
  #### Level3 ####
  
  
  slDT <- createl3(cDT, lthresh)
  # #Summarize cell data to medians of the spot parameters####
  # parameterNames<-grep(pattern="(Children|_CP_|_PA_|Barcode|^Spot$|^Well$)",x=names(cDT),value=TRUE)
  # 
  # #Remove any spot-normalized and cell level parameters
  # parameterNames <- grep("SpotNorm|^Nuclei_PA_Gated_EduPositive$|^Nuclei_PA_Gated_EduPositive_RZSNorm$",parameterNames,value=TRUE,invert=TRUE)
  # 
  # #Remove any raw parameters
  # parameterNames <- grep("Barcode|^Spot$|^Well$|Norm|Nuclei_CP_Intensity_MedianIntensity_Dapi$|Cytoplasm_CP_Intensity_MedianIntensity_Actin$|Cytoplasm_CP_Intensity_MedianIntensity_CellMask$|Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker$|Nuclei_CP_Intensity_MedianIntensity_H3$|Nuclei_CP_Intensity_MedianIntensity_Fibrillarin$|Nuclei_CP_Intensity_MedianIntensity_Edu$|Cytoplasm_CP_Intensity_MedianIntensity_KRT5$|Cytoplasm_CP_Intensity_MedianIntensity_KRT19$|Spot_PA_SpotCellCount$", parameterNames, value = TRUE)
  # 
  # cDTParameters<-cDT[,parameterNames,with=FALSE]
  # 
  # slDT<-cDTParameters[,lapply(.SD,numericMedian),keyby="Barcode,Well,Spot"]
  # slDTse <- cDTParameters[,lapply(.SD,se),keyby="Barcode,Well,Spot"]
  # 
  # #Add _SE to the standard error column names
  # setnames(slDTse, grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE),"_SE"))
  # 
  # #Merge back in the spot and well metadata
  # #TODO: Convert the logic to not name the metadata
  # metadataNames <- grep("(Row|Column|PrintOrder|Block|^ID$|Array|CellLine|Ligand|Endpoint|ECMp|MEP|Barcode|^Well$|^Spot$)", x=colnames(cDT), value=TRUE)
  # setkey(cDT,Barcode, Well,Spot)
  # mDT <- cDT[,metadataNames,keyby="Barcode,Well,Spot", with=FALSE]
  # slDT <- mDT[slDT, mult="first"]
  # #Merge in the standard errr values
  # slDT <- slDTse[slDT]
  # #Add a count of replicates
  # slDT <- slDT[,Spot_PA_ReplicateCount := .N,by="LigandAnnotID,ECMpAnnotID"]
  # 
  # #Add the loess model of the SpotCellCount on a per well basis
  # slDT <- slDT[,Spot_PA_LoessSCC := loessModel(.SD, value="Spot_PA_SpotCellCount", span=.5), by="Barcode,Well"]
  # 
  # #Add well level QA Scores to spot level data
  # slDT <- slDT[,QAScore := calcQAScore(.SD,threshold=lthresh,maxNrSpot = max(cDT$ArrayRow)*max(cDT$ArrayColumn),value="Spot_PA_LoessSCC"),by="Barcode,Well"]
  
  #Add QA scores to cell level data####
  setkey(cDT,Barcode, Well, Spot)
  setkey(slDT, Barcode, Well, Spot)
  cDT <- cDT[slDT[,list(Barcode, Well, Spot, QAScore, Spot_PA_LoessSCC)]]
  #The spot level data is median summarized to the replicate level and is stored as Level 4 data and metadata.
  
  #Level4Data
  mepDT <- createl4(slDT)
  
  #Write QA flags into appropriate data levels
  #Low cell count spots
  cDT$QA_LowSpotCellCount <- cDT$Spot_PA_SpotCellCount < lowSpotCellCountThreshold
  slDT$QA_LowSpotCellCount <- slDT$Spot_PA_SpotCellCount < lowSpotCellCountThreshold
  
  #Manually flag low quality DAPI wells
  cDT$QA_LowDAPIQuality <- FALSE
  cDT$QA_LowDAPIQuality[cDT$Barcode=="LI8X00420"&cDT$Well=="B01"] <- TRUE
  cDT$QA_LowDAPIQuality[cDT$Barcode=="LI8X00425"&cDT$Well=="B01"] <- TRUE
  cDT$QA_LowDAPIQuality[cDT$Barcode=="LI8X00426"] <- TRUE
  cDT$QA_LowDAPIQuality[cDT$Barcode=="LI8X00427"] <- TRUE
  
  slDT$QA_LowDAPIQuality <- FALSE
  slDT$QA_LowDAPIQuality[slDT$Barcode=="LI8X00420"& slDT$Well=="B01"] <- TRUE
  slDT$QA_LowDAPIQuality[slDT$Barcode=="LI8X00425"& slDT$Well=="B01"] <- TRUE
  slDT$QA_LowDAPIQuality[slDT$Barcode=="LI8X00426"] <- TRUE
  slDT$QA_LowDAPIQuality[slDT$Barcode=="LI8X00427"] <- TRUE
  
  
  #Flag spots below automatically loess QA threshold
  cDT$QA_LowRegionCellCount <- cDT$Spot_PA_LoessSCC < lowRegionCellCountThreshold
  slDT$QA_LowRegionCellCount <- slDT$Spot_PA_LoessSCC < lowRegionCellCountThreshold
  
  #Flag wells below automatically calculated QA threshold
  slDT$QA_LowWellQA <- FALSE
  slDT$QA_LowWellQA[slDT$QAScore < lowWellQAThreshold] <- TRUE
  cDT$QA_LowWellQA <- FALSE
  cDT$QA_LowWellQA[cDT$QAScore < lowWellQAThreshold] <- TRUE
  
  #Level 4
  mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < lowReplicateCount
  
  #WriteData
  
  if(writeFiles){
    #Write out cDT without normalized values as level 1 dataset
    level1Names <- grep("Norm",colnames(cDT),value=TRUE,invert=TRUE)
    write.table(format(cDT[,level1Names, with=FALSE], digits=4, trim=TRUE), paste0("./",cellLine,"/", ss, "/AnnotatedData/", unique(cDT$CellLine),"_",ss,"_","Level1.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
    
    normParmameterNames <- grep("Norm",colnames(cDT), value=TRUE)
    rawParameterNames <- gsub("_?[[:alnum:]]*?Norm$", "", normParmameterNames)
    metadataNormNames <- colnames(cDT)[!colnames(cDT) %in% rawParameterNames]
    #Paste back in the QA and selected raw data
    
    level2Names <- c(metadataNormNames,
                     grep("Nuclei_CP_Intensity_MedianIntensity_Dapi$|Cytoplasm_CP_Intensity_MedianIntensity_Actin$|Cytoplasm_CP_Intensity_MedianIntensity_CellMask$|Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker$|Nuclei_CP_Intensity_MedianIntensity_H3$|Nuclei_CP_Intensity_MedianIntensity_Firbillarin$|Nuclei_CP_Intensity_MedianIntensity_Edu$|Cytoplasm_CP_Intensity_MedianIntensity_KRT5$|Cytoplasm_CP_Intensity_MedianIntensity_KRT19$|Spot_PA_SpotCellCount$", colnames(cDT), value = TRUE))
    
    #Write out cDT with normalized values as level 2 dataset
    write.table(format(cDT[,level2Names, with = FALSE], digits=4, trim=TRUE), paste0("./",cellLine,"/", ss, "/AnnotatedData/", unique(cDT$CellLine),"_",ss,"_","Level2.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
    
    write.table(format(slDT, digits = 4, trim=TRUE), paste0("./",cellLine,"/", ss, "/AnnotatedData/", unique(slDT$CellLine),"_",ss,"_","Level3.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
    
    write.table(format(mepDT, digits = 4, trim=TRUE), paste0("./",cellLine,"/",ss, "/AnnotatedData/", unique(slDT$CellLine),"_",ss,"_","Level4.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
    
    #Write the pipeline parameters to  tab-delimited file
    write.table(c(
      ss=ss,
      cellLine = cellLine,
      analysisVersion = analysisVersion,
      neighborhoodNucleiRadii = neighborhoodNucleiRadii,
      neighborsThresh = neighborsThresh,
      wedgeAngs = wedgeAngs,
      outerThresh = outerThresh,
      nuclearAreaThresh = nuclearAreaThresh,
      nuclearAreaHiThresh = nuclearAreaHiThresh,
      curatedOnly = curatedOnly,
      curatedCols = curatedCols,
      writeFiles = writeFiles,
      limitBarcodes = limitBarcodes,
      normToSpot = normToSpot,
      lowSpotCellCountThreshold = lowSpotCellCountThreshold,
      lowRegionCellCountThreshold = lowRegionCellCountThreshold,
      lowWellQAThreshold = lowWellQAThreshold,
      lowReplicateCount =lowReplicateCount,
      lthresh = lthresh
    ),
    paste0("./",cellLine,"/",ss, "/AnnotatedData/", cellLine,"_",ss,"_","PipelineParameters.txt"), sep = "\t",col.names = FALSE, quote=FALSE)
  }
}