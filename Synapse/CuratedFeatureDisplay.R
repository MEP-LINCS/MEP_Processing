library(synapseClient)
synapseLogin()
d <- synGet("syn7122655")
dt <- fread(getFileLocation(d))
curatedFeatures <- grep("MEP|^ECMp$|^Ligand$|Nuclei_CP_AreaShape_AreaLog2RUVLoess|Nuclei_CP_AreaShape_PerimeterLog2RUVLoess|Nuclei_CP_Intensity_IntegratedIntensity_DapiLog2RUVLoess|Nuclei_CP_Intensity_MedianIntensity_DapiLog2RUVLoess|Spot_PA_SpotCellCountLog2RUVLoess|Nuclei_PA_Cycle_DNA2NProportionLogitRUVLoess|^Spot_PA_SpotCellCountLog2$|Cytoplasm_CP_Intensity_MedianIntensity_.*Log2RUVLoess|Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoess$|Cytoplasm_PA_Gated_KRT19PositiveProportionLogitRUVLoess",colnames(dt),value=TRUE)
featureNames <- grep("_SE|RLE",curatedFeatures,value=TRUE,invert=TRUE)
shortNames <- c("")
l <- list(list(FeatureName="MEP",DisplayName="MEP",Description="The microenvironment perturbation that is the combination of an ECM protein and a ligand."),
          list("Ligand","Ligand","The cells were grown in this soluble ligand, growth factor or cytokine."),
          list("ECMp","ECM Protein","The cells were grown on this extracellular matrix protein and Collagen 1."),
          list("Nuclei_CP_AreaShape_AreaLog2RUVLoess", "Nuclear Area", "The nuclear area measured in pixels."),
          list("Nuclei_CP_AreaShape_PerimeterLog2RUVLoess", "Nuclear Perimeter","The nuclear perimeter measured in pixels."),
          list("Nuclei_CP_Intensity_IntegratedIntensity_DapiLog2RUVLoess","Total DAPI Intensity (Norm)","The sum of the pixel intensities of the DAPI channel measured in the nucleus."),
          list("Nuclei_CP_Intensity_MedianIntensity_DapiLog2RUVLoess","Median DAPI Intensity (Norm)","The median intensity of the DAPI pixels measured in the nucleus."),
          list("Spot_PA_SpotCellCountLog2RUVLoess","Spot Cell Count (Norm)","The normalized number of nuclei segmented at each spot."),
          list("Nuclei_PA_Cycle_DNA2NProportionLogitRUVLoess","DNA 2N Proportion","The normalized proportion of cells at each spot gated in to the 2N DNA population based on total DAPI intensity."),
          list("Spot_PA_SpotCellCountLog2","Spot Cell Count (Norm)","The raw number of nuclei segmented at each spot."),
          list("Cytoplasm_CP_Intensity_MedianIntensity_KRT19Log2RUVLoess","Median KRT19 Intensity","The median intensity of the KRT19 luminal marker measured in the cytoplasm."),
          list("Cytoplasm_CP_Intensity_MedianIntensity_KRT5Log2RUVLoess","Median KRT5 Intensity","The median intensity of the KRT5 basal marker measured in the cytoplasm."),
          list("Cytoplasm_PA_Gated_KRT19PositiveProportionLogitRUVLoess","KRT19 High Proportion","The proportion of cells at each spot in the KRT19 high gate."),
          list("Cytoplasm_CP_Intensity_MedianIntensity_MitoTrackerLog2RUVLoess","Median MitoTracker Intensity","The median intensity of the MitoTracker metabolism marker measured in the cytoplasm."),
          list("Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoess","EdU+ Proportion","The proportion of cells at each spot that are positive for the EdU proliferation marker.")
)
f <- rbindlist(l)
fileName <- "MEP_LINCS/AnnotatedData/CuratedFeatures.tsv"
fwrite(f,fileName,sep="\t")
synapseAnnotatedDataDir <- "syn5713302"

obj <- File(fileName, parentId=synapseAnnotatedDataDir)
obj <- synStore(obj,contentType="text/tab-separated-values")
obj
