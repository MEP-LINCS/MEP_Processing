library("data.table")

l3 <- fread(paste0("~/Documents/ME Watson/Lapatinib MEMAs/AnnotatedData/AU565_SS6_v2_av1.6_Level1.txt"), showProgress = FALSE)
l4 <- fread(paste0("../AnnotatedData/",cellLine, "_",ss,"_",drug,"_",rawDataVersion,"_",analysisVersion,"_Level4.txt"), showProgress = FALSE)
