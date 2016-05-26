SS2StainingSet <- "SS2noH3"
PC3_SS2noH3_l4 <- fread(paste0("./PC3/SS2noH3/AnnotatedData/PC3_",SS2StainingSet,"_v1_av1.4_Level4.txt"), showProgress = FALSE)

SS2StainingSet <- "SS2"
PC3_SS2_l4 <- fread(paste0("./PC3/SS2/AnnotatedData/PC3_",SS2StainingSet,"_v2.1_av1.4_Level4.txt"), showProgress = FALSE)

annotIDs_phase1 <- unique(c(PC3_SS2_l4$ECMpAnnotID, PC3_SS2_l4$LigandAnnotID))

write.csv(annotIDs_phase1,"AnnotIDs_Phase1.csv",row.names=FALSE,quote=FALSE)

annotIDs_phase2 <- unique(c(PC3_SS2noH3_l4$ECMpAnnotID, PC3_SS2noH3_l4$LigandAnnotID))

write.csv(annotIDs_phase2,"AnnotIDs_Phase2.csv",row.names=FALSE,quote=FALSE)
