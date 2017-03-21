library("rmarkdown")
library("ruv")
library("data.table")

plotLEHmap <- function(dt, fill, titleExpression, limits){
  p <- ggplot(dt, aes_string(x="Ligand", y="ECMp", fill = fill))+
    geom_tile()+
    scale_fill_gradient(low="white",high="red",oob = scales::squish,
                        limits=limits)+
    ggtitle(titleExpression)+
    xlab("")+ylab("")+
    guides(fill=FALSE)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(.6)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(.5)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.7)),panel.grid.major = element_line(linetype = 'blank'),panel.background = element_rect(fill = "dimgray"))
  suppressWarnings(print(p))
}


addMarginValues <- function(dt, mValue, MValue){
  #browser()
  dt.l <- dt[,list(m.l = median(get(mValue), na.rm = TRUE),
                   M.l = median(get(MValue), na.rm= TRUE)),
             by="Ligand"]
  
  dte. <- dt[,list(me. = median(get(mValue), na.rm = TRUE),
                   Me. = median(get(MValue), na.rm= TRUE)),
             by="ECMp"]
  
  setkey(dt,Ligand)
  setkey(dt.l,Ligand)
  dtm <- merge(dt,dt.l)
  
  setkey(dtm,ECMp)
  setkey(dte.,ECMp)
  dtmm <- merge(dtm, dte.) 
  return(dtmm)
}

plotCenteredBoxes <- function(dt, value, groupBy, colourBy=NULL, titleExpression, yLimits=NULL){
  
  if(is.null(colourBy))  p <- ggplot(dt,aes_string(x=groupBy, y=value))
  else p <- ggplot(dt,aes_string(x=groupBy, y=value, colour=colourBy))
  
  p <- p + geom_boxplot()+ 
    ggtitle(titleExpression)+
    xlab("")+ylab("")+
    guides(fill=FALSE)+
    coord_cartesian(ylim=yLimits)+
    theme_grey(base_size = 12, base_family = "")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(1)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(1)))
  suppressWarnings(print(p))
}

calcResidual <- function(x){
  mel <- median(x, na.rm=TRUE)
  return(x-mel)
}

####################################

x <- data.table(CellLine=rep(c("MCF10A"),each=1),
                        StainingSet = rep(c("SS20")),
                        Signal=rep(c("SCC"), each=1),
                        Method=c("RUV3Loess"),
                        inputFileName=rep(c("../MEP-LINCS/PC3/SS2noH3/AnnotatedData/PC3_SS2noH3_v1.3_Level1.txt"), each=1))

callQASVD <- function(x){
  render(paste0("Mep-LINCS_QASVD.Rmd"),
         output_file = paste0("Mep-LINCS_SourcePlateQA_",x[["CellLine"]],"_",x[["Signal"]],".html"),
         output_format = "html_document") 
}

x <- unique(x[,list(CellLine, StainingSet, Signal, inputFileName)])

apply(x, 1, callQASVD)

