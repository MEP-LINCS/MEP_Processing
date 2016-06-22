library("data.table")

for(cellLine in c("AU565","HCC1954")){
  for(drug in c("DMSO","Lapatinib")){
    dt <- fread(paste0("./AnnotatedData/",cellLine,"_SS6_",drug,"_v2_av1.6_Level4.txt"), showProgress = FALSE)
    dtCurated <- (dt[,grep("Barcode|^Well$|^Spot$|Array|MEP|ECMp|^Ligand$|Cells_CP_AreaShape_AreaLog2|Nuclei_CP_AreaShape_AreaLog2|Cytoplasm_CP_Intensity_IntegratedIntensity_KRT14Log2|Cytoplasm_CP_Intensity_IntegratedIntensity_KRT19Log2|Cytoplasm_CP_Intensity_MedianIntensity_KRT14Log2|Cytoplasm_CP_Intensity_MedianIntensity_KRT19Log2|Nuclei_CP_Intensity_IntegratedIntensity_DapiLog2|Nuclei_CP_Intensity_IntegratedIntensity_EdULog2|Nuclei_CP_Intensity_MedianIntensity_DapiLog2|Nuclei_CP_Intensity_MedianIntensity_EdULog2|Nuclei_PA_Cycle_DNA2NProportionLogit|Nuclei_PA_Cycle_DNA4NProportionLogit|Nuclei_PA_Gated_EdUPositiveProportionLogit",colnames(dt), value=TRUE), with=FALSE])
    dtCurated <- dtCurated[,grep("RUV3$|_SE",colnames(dtCurated), value=TRUE, invert=TRUE), with=FALSE]
    fwrite(dtCurated, file.path = paste0("./AnnotatedData/",cellLine,"_SS6_",drug,"_v2_av1.6_Level4Curated.csv"))
  }
}



