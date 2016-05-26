library("data.table")
library("ggplot2")

dtList<- lapply(c("MCF7", "PC3", "YAPC"), function(cellLine){
  cellLineList <- lapply(c("SS1", "SS2","SS3"), function(ss, cellLine){
    l1 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine,"_",ss,"_Level1.txt"), showProgress = FALSE)[,list(CellLine,Barcode,Well,Spot,ECMp,Ligand,Spot_PA_SpotCellCount,Nuclei_CP_Intensity_MedianIntensity_Dapi)]
    l1$StainingSet = ss
    return(l1)
  }, cellLine=cellLine)
  cellLineDT <- rbindlist(cellLineList)
})
dt <- rbindlist(dtList)

setkey(dt,Ligand)
dtfbs <- dt["FBS",list(SpotCellCount=.N, NuclearDapi=median(Nuclei_CP_Intensity_MedianIntensity_Dapi, na.rm=TRUE), StainingSet),by="CellLine,Barcode,Well,Spot"]
setkey(dtfbs,CellLine)

pdf("v1FBSWells.pdf", height = 4)
for(cellLine in c("MCF7", "PC3", "YAPC")){
  dt <- dtfbs[cellLine]
p <- ggplot(dt, aes(x=Barcode, y=SpotCellCount, colour=StainingSet))+
  geom_boxplot()+
  ggtitle(paste("Spot Cell Count in FBS Wells for",cellLine))+
  ylim(0,100)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(.8)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.3)))
print(p)
}

for(cellLine in c("MCF7", "PC3", "YAPC")){
  dt <- dtfbs[cellLine]
  p <- ggplot(dt, aes(x=Barcode, y=NuclearDapi, colour=StainingSet))+
    geom_boxplot()+
    ggtitle(paste("Median Nuclear Intensity in FBS Wells for",cellLine))+
    ylim(0,.25)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(.8)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.3)))
  print(p)
}
dev.off()