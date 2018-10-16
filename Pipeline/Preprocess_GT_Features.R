#!/usr/bin/env Rscript

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
library(tidyverse)
library(MEMA)

barcodes <- c("LI8V01171", "LI8V01172", "LI8V01173a", "LI8V01174","LI8V01175", "LI8V01176", "LI8V01177", "LI8V01178", "LI8V01179")
raw_data_dir <- "/graylab/share/thibaulg/Mema/HCC1143_COL1 - Features"
known_errors <- c(3074163, 3101463)

#Replace channel names with staining set specific values
#ssT
channel_0 <- "Dapi"
channel_1 <- "KRT14"
channel_2 <- "EdU"
channel_3 <- "VIM"

#get feature files
files <- tibble(Full_filename = dir(raw_data_dir, full.names = TRUE)) %>%
  mutate(Filename = str_remove(Full_filename, ".*/"),
         Channel = str_extract(Filename, "Channel [[:digit:]]*"),
         Channel = str_replace(Channel, " ","_"),
         Channel = str_replace(Channel, "Channel_0", channel_0),
         Channel = str_replace(Channel, "Channel_1", channel_1),
         Channel = str_replace(Channel, "Channel_2", channel_2),
         Channel = str_replace(Channel, "Channel_3", channel_3),
         Compartment = str_extract(Filename, "Nuclei|CellsMinusNuclei|Cells"),
         Compartment = str_replace(Compartment, "CellsMinusNuclei", "Cytoplasm"),
         ImageID = str_extract(Filename,"[[:digit:]]*"),
         ImageID = as.integer(ImageID)) %>%
  filter(!ImageID %in% known_errors)

imageIDs <- unique(files$ImageID)[sample(1:length(unique(files$ImageID)),size = 10, replace = FALSE)]
#imageIDs <- unique(files$ImageID)

addPrefix <-function(x, compartment, channel) paste0(compartment,"_GT_",channel,"_",x)

#Read in raw data with each row a cell
#0 area cyto compartments result in NA values in Cyto features
start_time <- Sys.time()
raw_data <- lapply(imageIDs, function(imageID){
  file_set <- files %>%
    filter(ImageID == imageID)
  foo <- lapply(1:nrow(file_set), function(i){
    rd <- suppressWarnings(suppressMessages(read_delim(file_set$Full_filename[i],delim = " ",col_names = TRUE))) %>%
      select(-starts_with("X")) %>%
      select(-matches("Haralick|SZM|_PS_|_LBP_|_Rank")) %>%
      rename_all(.funs=addPrefix,compartment=file_set$Compartment[i],channel=file_set$Channel[i])
    rd$Label <- rd[[paste0(file_set$Compartment[i],"_GT_",channel=file_set$Channel[i],"_Label")]]
    rd <- rd %>%
      select(-matches(paste0(file_set$Compartment[i],"_GT_",channel=file_set$Channel[i],"_Label"))) 
  }) %>% 
    plyr::join_all(by="Label") %>%
    mutate(ImageID = imageID)
}) %>%
  bind_rows()
Sys.time()-start_time

#Read in the metadata
metadata <- lapply(barcodes, function(barcode){
  inputPath <- "/lincs/share/lincs_user"
  metadataFiles <- list(annotMetadata=paste0(inputPath,"/",barcode,"/Analysis/",barcode,"_an2omero.csv"))
  metadata <- MEMA::getMetadata(metadataFiles, TRUE)
  
  imageIDs <-  paste0(inputPath,"/",barcode,"/Analysis/",barcode,"_imageIDs.tsv") %>%
    read_tsv(col_types = cols(
      PlateID = col_character(),
      WellName = col_character(),
      Row = col_integer(),
      Column = col_integer(),
      ImageID = col_integer())) %>%
    rename(ArrayRow = Row,
           ArrayColumn = Column) %>%
    mutate(WellIndex = str_remove(WellName, ".*Well"),
           WellIndex = as.integer(WellIndex),
           Well = c("A01","A02","A03","A04","B01","B02","B03","B04")[WellIndex],
           Spot = ArrayColumn+(ArrayRow-1)*20) %>%
    select(-WellName)
  
  #merge the metadata with the ImageIDs
  full_metadata <- full_join(metadata,imageIDs, by = c("Well", "Spot", "WellIndex", "ArrayRow", "ArrayColumn"))
}) %>%
  bind_rows()

level1_data <- raw_data %>%
  inner_join(metadata, by="ImageID") %>%
  mutate(Nuclei_CP_Intensity_MedianIntensity_Dapi = Nuclei_GT_Dapi_MeanIntensity,
         Nuclei_CP_Intensity_IntegratedIntensity_Dapi = Nuclei_GT_Dapi_MeanIntensity*Nuclei_GT_Dapi_Surface,
         Nuclei_CP_Intensity_MedianIntensity_EdU = Nuclei_GT_EdU_MeanIntensity,
         Nuclei_CP_Intensity_IntegratedIntensity_EdU = Nuclei_GT_EdU_MeanIntensity*Nuclei_GT_EdU_Surface,
         Cytoplasm_CP_Intensity_MedianIntensity_KRT14 = Cytoplasm_GT_KRT14_MeanIntensity,
         Cytoplasm_CP_Intensity_IntegratedIntensity_KRT14 = Cytoplasm_GT_KRT14_MeanIntensity*Cytoplasm_GT_KRT14_Surface,
         Cytoplasm_CP_Intensity_IntegratedIntensity_VIM = Cytoplasm_GT_VIM_MeanIntensity*Cytoplasm_GT_VIM_Surface,
         Cytoplasm_CP_Intensity_MedianIntensity_VIM = Cytoplasm_GT_VIM_MeanIntensity,
         Well_Ligand = paste0(Well,"_",Ligand)) %>%
  data.table::data.table() %>%
  gateCells()

# # Filter our debris and cell clusters
# dtL <- lapply(dtL, function(dt){
#   filterObjects(dt,nuclearAreaThresh = 200, nuclearAreaHiThresh = 4000)})
# 
# # Add local XY and polar coordinates
# # dtL <- mclapply(dtL, addPolarCoords, mc.cores = detectCores())
# dtL <- lapply(dtL, addPolarCoords)
# 
# # Add spot level normalizations for median intensities
# # dtL <- mclapply(dtL,spotNormIntensities, mc.cores = detectCores())
# dtL <- lapply(dtL,spotNormIntensities)
# 
# # Add adjacency values
# # dtL <- mclapply(dtL, calcAdjacency, mc.cores = detectCores())
# dtL <- lapply(dtL, calcAdjacency)
# 
# #reduce the numeric values to 4 significant digits
# shorten <- function(x){
#   if(class(x)=="numeric") x <- signif(x,4)
#   return(x)
# }
# for (j in colnames(cDT)) data.table::set(cDT, j = j, value = shorten(cDT[[j]]))

write_csv(level1_data, "/graylab/share/dane/MEP-LINCS/DevelopmentExperiments/MEMA_Analysis/HCC1143_level1_GT_4.csv")

level2_numeric_data <- level1_data %>%
  group_by(Barcode,Well,Spot) %>%
  mutate(Spot_PA_SpotCellCount=n()) %>%
  summarise_if(is.numeric, median) 
# %>%
#   mutate_if(is.numeric, signif, digits=4)

level2_character_data <- level1_data %>%
  group_by(Barcode,Well,Spot) %>%
  summarise_if(is.character, unique)

level2_logical_data <- level1_data %>%
  group_by(Barcode,Well,Spot) %>%
  summarise_if(is.logical, unique)

level2_data <- full_join(level2_character_data,level2_numeric_data) %>%
  right_join(level2_logical_data)

write_csv(level2_data, "/graylab/share/dane/MEP-LINCS/DevelopmentExperiments/MEMA_Analysis/HCC1143_level2_GT_4.csv")

