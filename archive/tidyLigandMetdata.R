library(XLConnect)
library(MEMA)

#Read in ligand metadata
ligandMetadata <- readWorksheetFromFile("Ligand-cellsource_effects_v1.xlsx", sheet="Sheet1",startRow=1, header=TRUE)
ligandMetadataDT <- convertColumnNames(data.table(ligandMetadata))
ligandMetadataDT$Ligand <- sub("_[[:alnum:]]*$","",ligandMetadataDT$protein)
ligandMetadataDT$Growth <- grepl("Growth",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Migration <- grepl("Migration",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Invasion <- grepl("Invasion",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$TherapeuticResistance <- grepl("therapeutic resistance",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Metastasis <- grepl("metastasis",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Survival <- grepl("survival",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Angiogenesis <- grepl("angiogenesis",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Proliferation <- grepl("Proliferation",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Immunomodulatory <- grepl("Immunomodulatory",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$TumorProgression <- grepl("tumor progression",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$TherapeuticResponse <- grepl("therapeutic response",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$StemCellMaintenance <- grepl("stem cell maintenance",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$TumorInitiation <- grepl("tumor initiation",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$VascularHomeostasis <- grepl("vascular homeostasis",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Inflammation <- grepl("inflammation",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$ImmuneEffects <- grepl("immune effects",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Progression <- grepl("progression",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Adhesion <- grepl("adhesion",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$TumorEffects <- grepl("tumor effects",ligandMetadataDT$Effect, ignore.case = TRUE)
ligandMetadataDT$Liver <- grepl("liver",ligandMetadataDT$TissueExpression, ignore.case = TRUE)
ligandMetadataDT$BoneMarrow <- grepl("bone marrow",ligandMetadataDT$TissueExpression, ignore.case = TRUE)
ligandMetadataDT$Brain <- grepl("brain",ligandMetadataDT$TissueExpression, ignore.case = TRUE)
ligandMetadataDT$Lung <- grepl("lung",ligandMetadataDT$TissueExpression, ignore.case = TRUE)
ligandMetadataDT$Thyroid <- grepl("thyroid",ligandMetadataDT$TissueExpression, ignore.case = TRUE)
ligandMetadataDT$LymphSystem <- grepl("lymph node|lymphatics|lymphoid",ligandMetadataDT$TissueExpression, ignore.case = TRUE)
ligandMetadataDT$Spleen <- grepl("spleen",ligandMetadataDT$TissueExpression, ignore.case = TRUE)
ligandMetadataDT$Muscle <- grepl("muscle",ligandMetadataDT$TissueExpression, ignore.case = TRUE)

#Delete the columns that had combined annotations
ligandMetadataDT <- ligandMetadataDT[,Effect :=NULL]
ligandMetadataDT <- ligandMetadataDT[,TissueExpression :=NULL]

fwrite(ligandMetadataDT,file.path = "LigandCellSourceEffects_v1.tsv",sep="\t")
     
