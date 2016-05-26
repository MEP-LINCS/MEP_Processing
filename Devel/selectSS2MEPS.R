
```{r global_options, include=FALSE}
knitr::opts_chunk$set(fig.width=12, fig.height=8,
                      echo=FALSE, warning=FALSE, message=FALSE)
```


```{r }
#Author: Mark Dane, copyright 2015

source("./MEPLINCSFunctions.R")

```



```{r setup}
library("ggplot2")
library("data.table")
library("MEMA")
library("grid")
library("knitr")
library("gplots")
library("RColorBrewer")
library(DT)
library(d3heatmap)
library(plotly)

#Setup colors for Barcode and text in all heatmaps
selDark2 <- colorRampPalette(brewer.pal(8,"Dark2"))
plateCol = selDark2(8)
hmcols<-colorRampPalette(c("blue","white","red"))(16)

l1 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine,"_",ss,"_",analysisVersion,"_Level1.txt"), showProgress = FALSE)
l2 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine,"_",ss,"_",analysisVersion,"_Level2.txt"), showProgress = FALSE)
l3 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine, "_",ss,"_",analysisVersion,"_Level3.txt"), showProgress = FALSE)
l4 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine, "_",ss,"_",analysisVersion,"_Level4.txt"), showProgress = FALSE)

barcodes <- sort(unique(l3$Barcode))


#Set a threshold for filtering wells on their QA score
wellQAThresh <- 0.7

#TODO: Read this from Level 3 data
lthresh <- 0.6

#Number of PCS components to use
nrPCs <- 9

#Z score threshold for extreme spot cell count
SCCZscoreThresh <- 3

#Spot cell count threshold for HF dataset
HFSCCThresh <- 20

#Replicate count threshold for HF Dataset
HFRepThresh <- 3

# https://meplincs.ohsu.edu/webclient/

if(!analysisVersion=="v1"){
  l3$OmeroDetailURL <- paste0('<a href="https://meplincs.ohsu.edu/webclient/img_detail/',l3$ImageID,'/"',' target="_blank">Omero</a>')
  l3$OmeroThumbnailURL <- paste0('<a href="https://meplincs.ohsu.edu/webclient/render_thumbnail/',l3$ImageID,'/"',' target="_blank">Omero</a>')
  l3$OmeroImageURL <- paste0('<a href="https://meplincs.ohsu.edu/webclient/render_image/',l3$ImageID,'/"',' target="_blank">Omero</a>')
}

```
#MEP-LINCS Validation Set
####date: `r Sys.Date()`
#Select MEPs for validation
#1/2015 Mark Dane

library(MEMA)
library(ggplot2)
library(plotly)
library(reshape2)
library(data.table)

getl4HF <- function(cellLine, analysisVersion){
  ss="SS2"
  
  l3 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine, "_",ss,"_",analysisVersion,"_Level3.txt"), showProgress = FALSE)
  
  #Set a threshold for filtering wells on their QA score
  wellQAThresh <- 0.7
  
  #TODO: Read this from Level 3 data
  lthresh <- 0.6
  
  #Number of PCS components to use
  nrPCs <- 9
  
  #Z score threshold for extreme spot cell count
  SCCZscoreThresh <- 3
  
  #Spot cell count threshold for HF dataset
  HFSCCThresh <- 20
  
  #Replicate count threshold for HF Dataset
  HFRepThresh <- 3
  
  #Remove the fiducial and blank data
  setkey(l3,ECMp)
  l3F <- l3[!"fiducial"]
  l3F <- l3F[!"blank"]
  setkey(l3F,"Ligand")
  l3F <- l3F[!grepl("FBS", l3F$Ligand)]
  
  #Remove failed QA wells
  l3F <- l3F[!l3F$QA_LowWellQA]
  #Filter for high spot occupancy and good quality DAPI####
  l3HF <- l3F[l3F$Spot_PA_SpotCellCount > HFSCCThresh]
  l3HF <- l3HF[!l3HF$QA_LowDAPIQuality]
  l3HF <- l3HF[!l3HF$QA_LowRegionCellCount]
  l4HF <- createl4(l3HF)
  l4HF <- l4HF[l4HF$Spot_PA_ReplicateCount >= HFRepThresh]
  return(l4HF)
}


validationMEPS <- function(cellLine, analysisVersion){
  ss="SS2"
  
  l3 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine, "_",ss,"_",analysisVersion,"_Level3.txt"), showProgress = FALSE)
  
  barcodes <- sort(unique(l3$Barcode))
  
  
  #Set a threshold for filtering wells on their QA score
  wellQAThresh <- 0.7
  
  #TODO: Read this from Level 3 data
  lthresh <- 0.6
  
  #Number of PCS components to use
  nrPCs <- 9
  
  #Z score threshold for extreme spot cell count
  SCCZscoreThresh <- 3
  
  #Spot cell count threshold for HF dataset
  HFSCCThresh <- 20
  
  #Replicate count threshold for HF Dataset
  HFRepThresh <- 3
  
  l3$OmeroDetailURL <- paste0('<a href="https://meplincs.ohsu.edu/webclient/img_detail/',l3$ImageID,'/"',' target="_blank">Omero</a>')
  l3$OmeroThumbnailURL <- paste0('<a href="https://meplincs.ohsu.edu/webclient/render_thumbnail/',l3$ImageID,'/"',' target="_blank">Omero</a>')
  l3$OmeroImageURL <- paste0('<a href="https://meplincs.ohsu.edu/webclient/render_image/',l3$ImageID,'/"',' target="_blank">Omero</a>')
  
  #Remove the fiducial and blank data
  setkey(l3,ECMp)
  l3F <- l3[!"fiducial"]
  l3F <- l3F[!"blank"]
  setkey(l3F,"Ligand")
  l3F <- l3F[!grepl("FBS", l3F$Ligand)]
  
  #Remove failed QA wells
  l3F <- l3F[!l3F$QA_LowWellQA]
  #Filter for high spot occupancy and good quality DAPI####
  l3HF <- l3F[l3F$Spot_PA_SpotCellCount > HFSCCThresh]
  l3HF <- l3HF[!l3HF$QA_LowDAPIQuality]
  l3HF <- l3HF[!l3HF$QA_LowRegionCellCount]
  l4HF <- createl4(l3HF)
  l4HF <- l4HF[l4HF$Spot_PA_ReplicateCount >= HFRepThresh]
  
  set.seed(1234)
  setkey(l4HF,Nuclei_PA_Gated_EduPositiveLogit_Norm)
  midMEPs <- round(quantile(1:nrow(l4HF), probs=c(.25,.75)))
  selectedMEPS <-l4HF[c(1:10, 
                        sample(seq(midMEPs[1], midMEPs[2], by=1), size = 5),
                        (nrow(l4HF)-9):nrow(l4HF))]
  
   p <- ggplot(l3HF, aes(Nuclei_PA_Gated_EduPositiveLogit,Nuclei_PA_Gated_EduPositiveLogit_Norm,colour=Barcode))+geom_point(size=.8, alpha=.8)+
     ggtitle("EdU Positive Normalized vs Raw (logit)")
  print(p)
  
  p <- ggplot(l3HF, aes(x=factor(Barcode), y=Nuclei_PA_Gated_EduPositiveLogit, colour=Barcode))+geom_boxplot()+
    coord_cartesian(ylim = c(-10,0))+
    ggtitle("EdU Positive Raw (logit)")
  print(p)
  
  p <- ggplot(l3HF, aes(x=factor(Barcode), y=Nuclei_PA_Gated_EduPositiveLogit_Norm, colour=Barcode))+geom_boxplot()+
    coord_cartesian(ylim = c(-10,0))+
    ggtitle("EdU Positive Normalized (logit)")
  print(p)

 return(selectedMEPS=selectedMEPS)
}

PC3MEPS <- validationMEPS(cellLine = "PC3",analysisVersion = "v1.3")
PC3MEPS$MEP[PC3MEPS$MEP=="CDH8_KNG1"&PC3MEPS$CellLine=="PC3"] <- "NID1_FGF6"
PC3MEPS$ECMp[PC3MEPS$MEP=="NID1_FGF6"] <- "NID1"
PC3MEPS$Ligand[PC3MEPS$MEP=="NID1_FGF6"] <- "FGF6"
PC3MEPS$ECMp[PC3MEPS$MEP=="NID1_FGF6"] <- "NID1"
PC3MEPS$Nuclei_PA_Gated_EduPositiveLogit_Norm[PC3MEPS$MEP=="NID1_FGF6"] <- -7.17

MCF7MEPS <- validationMEPS(cellLine = "MCF7",analysisVersion = "v1.3")
MEPS <-setkey(rbind(PC3MEPS,MCF7MEPS),MEP)


#write.csv(setkey(MEPS[,list(MEP,ECMp, Ligand, CellLine,Nuclei_PA_Gated_EduPositiveLogit_Norm,Barcode)],CellLine,Nuclei_PA_Gated_EduPositiveLogit_Norm), file = paste0("validationMEPS_v1.csv"))


p <- ggplot(l3, aes(Nuclei_PA_Gated_EduPositiveLogit,Nuclei_PA_Gated_EduPositiveLogit_Norm,colour=Barcode))+geom_point()
print(p)

#What are the EdU+ values for the 50 MEPs in the PC3 and MCF7 cell lines
#Merge the l4 HF MCF7 and PC3 datasets
l4HFDT <- rbindlist(lapply(c("PC3","MCF7"), getl4HF, analysisVersion=analysisVersion))
#Plot the rank ordered EdU+ values, colored by cellline, shape by in set
l4HFDT$MC <-paste(l4HFDT$MEP,l4HFDT$CellLine,sep="_")
l4HFDT$Selected <- FALSE
l4HFDT$Selected[l4HFDT$MEP %in% MEPS$MEP] <- TRUE

#write.csv(setkey(l4HFDT[l4HFDT$Selected,list(MEP,ECMp, Ligand, CellLine,Nuclei_PA_Gated_EduPositiveLogit_Norm,Barcode)],CellLine,Nuclei_PA_Gated_EduPositiveLogit_Norm), file = paste0("validationMEPSInBothCellLines_v1.csv"))

pdf(paste0(cellLine,"_EdUHF.pdf"))

for(cellLine in c("PC3","MCF7")){
  setkey(l4HFDT,CellLine)
  DT <- l4HFDT[cellLine]
  pdf(paste0(cellLine,"_EdUHF.pdf"))
  
    p <- ggplot(DT, aes(x =reorder(MEP, Nuclei_PA_Gated_EduPositiveLogit_Norm), y = Nuclei_PA_Gated_EduPositiveLogit_Norm, colour=Selected, shape=Selected, size=Selected))+
    #geom_errorbar(aes(ymin=Nuclei_PA_Gated_EduPositiveLogit_Norm-Nuclei_PA_Gated_EduPositiveLogit_Norm_SE, ymax=Nuclei_PA_Gated_EduPositiveLogit_Norm+Nuclei_PA_Gated_EduPositiveLogit_Norm_SE), width=.01, colour="black") +
    xlab("MEP")+ylab("Normalized, EdU+ Proportion Ratio")+
    geom_point(alpha = 1)+
    theme( axis.text.x=element_blank(), axis.ticks=element_blank(),  panel.grid.major = element_blank())+
    ggtitle(paste("Validation MEPs in",cellLine,"Ordered by EdU+ Proportion"))
  
  print(p)
  dev.off()
}

#Plot the responses in PC3 vs MCF7 for the MEPS validation set
MEPSL <- lapply(c("MCF7","PC3"),function(cellLine){
  fread(paste0("./",cellLine,"/SS2/AnnotatedData/",cellLine, "_SS2_v1.3_Level4.txt"), showProgress = FALSE)
})
MEPSDT <- rbindlist(MEPSL)
setkey(MEPSDT,MEP)
setkey(MEPS,MEP)
MEPSDT <- MEPSDT[MEPS[,list(MEP)]]
MEPSDTc <- dcast(MEPSDT,MEP~CellLine,value.var = "Nuclei_PA_Gated_EduPositiveLogit_Norm")
p <- ggplot(MEPSDTc, aes(PC3,MCF7,colour=MEP))+
  geom_point()

ggplotly(p)
