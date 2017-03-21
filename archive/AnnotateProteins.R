library(data.table)
library(XLConnect)

#Use any level 4 dataset that uses Annot metadata
l4 <- fread(paste0("../AnnotatedData/MCF10A_SS1_none_v2_av1.6_Level4.txt"), showProgress = FALSE)
#Reduce to ECMp and Ligand columns
types <- l4[,.(ECMp,Ligand)]
#Combine to get a type column
typesM <- unique(data.frame(protein=c(types$Ligand, types$ECMp),
                     Class=rep(c("Ligand","ECMp"), each=dim(types)[1]), stringsAsFactors = FALSE))
#TOdo: Clean up FBS...
typesM$protein[grepl("FBS",typesM$protein)] <- "FBS"

#Read in protein brick
proteinBrick <- data.table(readWorksheetFromFile("../protein_brick_20160615_114217_human.xlsx", sheet="protein_brick_20160615_114217_h",startRow=1, header=TRUE))
#reduce to Protein column
proteins <- unique(proteinBrick[,.(protein)])
#Remove lagging specifier
proteins$protein<- sub("_.*","",proteins$protein)

#Read in proteinset brick
proteinsetBrick <- fread("../proteinset_brick_20160615_114607_human.txt")
#reduce to Protein column
proteinsets <- unique(proteinsetBrick[,.(proteinset)])
#Remove lagging specifier
proteinsets$proteinset<- sub("_.*","",proteinsets$proteinset)
setnames(proteinsets, "proteinset", "protein")

proteins <- rbind(proteins,proteinsets,
                  data.frame(protein=c("hyaluronicacidgreaterthan500kDa","hyaluronicacidlessthan500kDa","FBS")))
#merge common proteins
proteinsA <-unique(merge(proteins,typesM,by="protein"))

write.csv(proteinsA,file="../MEP-LINCSReagentClasses.csv",row.names = FALSE)

