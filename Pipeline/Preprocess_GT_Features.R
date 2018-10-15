#!/usr/bin/env Rscript

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
library(tidyverse)
library(MEMA)

barcodes <- paste0("LI8V0117",1:9)
raw_data_dir <- "/graylab/share/thibaulg/Mema/HCC1143_COL1 - Features"
known_errors <- c(3074163, 3101463)

#get feature files
files <- tibble(Full_filename = dir(raw_data_dir, full.names = TRUE)) %>%
  mutate(Filename = str_remove(Full_filename, ".*/"),
         Channel = str_extract(Filename, "Channel [[:digit:]]*"),
         Channel = str_replace(Channel, " ","_"),
         Compartment = str_extract(Filename, "Nuclei|CellsMinusNuclei|Cells"),
         ImageID = str_extract(Filename,"[[:digit:]]*"),
         ImageID = as.integer(ImageID)) %>%
  filter(!ImageID %in% known_errors)

imageIDs <- unique(files$ImageID[1:10])

addPrefix <-function(x, compartment, channel) paste0(compartment,"_",channel,"_",x)

raw_data <- lapply(imageIDs, function(imageID){
  file_set <- files %>%
    filter(ImageID == imageID)
  foo <- lapply(1:nrow(file_set), function(i){
    rd <- read_delim(file_set$Full_filename[i],delim = " ",col_names = TRUE) %>%
      select(-starts_with("X")) %>%
      rename_all(.funs=addPrefix,compartment=file_set$Compartment[i],channel=file_set$Channel[i])
    rd$Label <- paste0(file_set$Compartment[i],"_",channel=file_set$Channel[i],"_Label")
    rd <- rd %>%
      select(-matches(paste0(file_set$Compartment[i],"_",channel=file_set$Channel[i],"_Label")))
  }) %>% 
    plyr::join_all(by="Label") %>%
    mutate(ImageID = imageID)
}) %>%
  bind_rows()

#Read in the metadata
metadata <- lapply(barcodes, function(barcode){
  inputPath <- "/lincs/share/lincs_user"
  metadataFiles <- list(annotMetadata=paste0(inputPath,"/",barcode,"/Analysis/",barcode,"_an2omero.csv"))
  metadata <- getMetadata(metadataFiles, TRUE)
  
  imageIDs <-  paste0(inputPath,"/",barcode,"/Analysis/",barcode,"_imageIDs.tsv") %>%
    read_tsv() %>%
    rename(ArrayRow = Row,
           ArrayColumn = Column) %>%
    mutate(WellIndex = str_remove(WellName, ".*Well"),
           WellIndex = as.integer(WellIndex),
           Well = c("A01","A02","A03","A04","B01","B02","B03","B04")[WellIndex],
           Spot = ArrayColumn+(ArrayRow-1)*20) %>%
    select(-WellName)
  
  #merge the metadata with the ImageIDs
  full_metadata <- full_join(metadata,imageIDs)
}) %>%
  bind_rows()

level1_data <- raw_data %>%
  left_join(metadata, by="ImageID")

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

#write_csv(level1_data, "/graylab/share/dane/MEP-LINCS/DevelopmentExperiments/MEMA_Analysis/HCC1143_level1_GT.csv")

level2_numeric_data <- level1_data %>%
  group_by(Barcode,Well,Spot) %>%
  mutate(Spot_PA_SpotCellCount=n()) %>%
  summarise_if(is.numeric, median) %>%
  mutate_if(is.numeric, signif, digits=4)

level2_character_data <- level1_data %>%
  group_by(Barcode,Well,Spot) %>%
  summarise_if(is.character, unique)

level2_logical_data <- level1_data %>%
  group_by(Barcode,Well,Spot) %>%
  summarise_if(is.logical, unique)

level2_data <- full_join(level2_character_data,level2_numeric_data) %>%
  right_join(level2_logical_data)

#write_csv(level2_data, "/graylab/share/dane/MEP-LINCS/DevelopmentExperiments/MEMA_Analysis/HCC1143_level2_GT.csv")
