library(data.table)
library(magrittr)

featureNames <- grep("MEP|^ECMp$|^Ligand$|Nuclei_CP_AreaShape_AreaLog2RUVLoess|Nuclei_CP_AreaShape_PerimeterLog2RUVLoess|Nuclei_CP_Intensity_IntegratedIntensity_DapiLog2RUVLoess|Nuclei_CP_Intensity_MedianIntensity_DapiLog2RUVLoess|Spot_PA_SpotCellCountLog2RUVLoess|Nuclei_PA_Cycle_DNA2NProportionLogitRUVLoess|^Spot_PA_SpotCellCountLog2$|Cytoplasm_CP_Intensity_MedianIntensity_.*Log2RUVLoess|Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoess$|Cytoplasm_PA_Gated_KRT19PositiveProportionLogitRUVLoess",colnames(dt),value=TRUE) %>%
  grep("_SE|RLE",.,value=TRUE,invert=TRUE)
shortNames <- c("")
l <- list(list(FeatureName="MEP",DisplayName="MEP",Description="The microenvironment perturbation that is the combination of an ECM protein and a ligand."),
          list("Ligand",
               "Ligand",
               "The cells were grown in this soluble ligand, growth factor or cytokine."),
          list("ECMp",
               "ECM Protein",
               "The cells were grown on this extracellular matrix protein and Collagen 1."),
          list("Well",
               "Well",
               "The alphanumeric label for the well."),
          list("Spot",
               "Spot",
               "The spot number in the MEMA array starting in the upper left corner and increasing acroos each row and then down each column."),
          list("Nuclei_CP_AreaShape_AreaLog2RUVLoess",
               "Nuclear Area (log2)",
               "The RUVLoess normalized nuclear area measured in pixels."),
          list("Nuclei_CP_AreaShape_PerimeterLog2RUVLoess",
               "Nuclear Perimeter (log2)",
               "The RUVLoess normalized nuclear perimeter measured in pixels."),
          list("Nuclei_CP_Intensity_IntegratedIntensity_DapiLog2RUVLoess",
               "Total DAPI Intensity (log2)",
               "The RUVLoess normalized sum of the pixel intensities of the DAPI channel measured in the nucleus."),
          list("Nuclei_CP_Intensity_MedianIntensity_DapiLog2RUVLoess",
               "Median DAPI Intensity (log2)",
               "The RUVLoess normalized median intensity of the DAPI pixels measured in the nucleus."),
          list("Spot_PA_SpotCellCountLog2RUVLoess",
               "Spot Cell Count (log2)",
               "The RUVLoess normalized number of nuclei segmented at each spot."),
          list("Nuclei_PA_Cycle_DNA2NProportionLogitRUVLoess",
               "DNA 2N Proportion",
               "The RUVLoess normalized proportion of cells at each spot gated in to the 2N DNA population based on total DAPI intensity."),
          list("Spot_PA_SpotCellCountLog2",
               "Spot Cell Count (raw, log2)",
               "The raw number of nuclei segmented at each spot."),
          list("Cytoplasm_CP_Intensity_MedianIntensity_KRT19Log2RUVLoess",
               "Median KRT19 Intensity (log2)",
               "The RUVLoess normalized median intensity of the KRT19 luminal marker measured in the cytoplasm."),
          list("Cytoplasm_CP_Intensity_MedianIntensity_KRT5Log2RUVLoess",
               "Median KRT5 Intensity (log2)",
               "The RUVLoess normalized median intensity of the KRT5 basal marker measured in the cytoplasm."),
          list("Cytoplasm_PA_Gated_KRT19PositiveProportionLogitRUVLoess",
               "KRT19 High Proportion",
               "The RUVLoess normalized proportion of cells at each spot in the KRT19 high gate."),
          list("Cytoplasm_CP_Intensity_MedianIntensity_MitoTrackerLog2RUVLoess",
               "Median MitoTracker Intensity (log2)",
               "The RUVLoess normalized median intensity of the MitoTracker metabolism marker measured in the cytoplasm."),
          list("Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoess",
               "EdU+ Proportion",
               "The RUVLoess normalized proportion of cells at each spot that are positive for the EdU proliferation marker.")
)
f <- rbindlist(l)

#Get a list of barcodes that will get new annotate2omero files
parentPath <- "/graylab/share/dane/MEP-LINCS/MEP_LINCS/AnnotatedData/"
l3FileNames <- dir(parentPath, pattern="Level3", full.names = TRUE) %>% 
  grep("MCF10A|HMEC|MCF7|PC3|YAPC",., value=TRUE)

tmpl <- lapply(l3FileNames[1], function(fn){
  l3 <- fread(fn,key = "Barcode")
  barcodes <- unique(l3$Barcode)
  for (barcode in barcodes[1]){
    dt <- l3[barcode,f$FeatureName[f$FeatureName %in% colnames(l3)], with=FALSE]
    dtd <- setnames(dt,colnames(dt),f$DisplayName[f$FeatureName %in% colnames(dt)])
    fwrite(dt,file=paste0("/lincs/share/lincs_user/",barcode,"/Analysis/",barcode,"_analysis2omero.tsv"),sep = "\t")
    fwrite(dtd,file=paste0("/lincs/share/lincs_user/",barcode,"/Analysis/",barcode,"_analysis2omeroDisplay.tsv"),sep = "\t")
    
  }
})
