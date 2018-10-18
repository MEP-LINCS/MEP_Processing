#!/usr/bin/env Rscript

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
library(data.table)
library(tidyverse)
library(MEMA)
library(optparse)

#Get the path to the plateID from the command line

# Get the command line arguments and options
# returns a list with options and args elements
getCommandLineArgs <- function(){
  parser <- OptionParser(usage = "%prog [options] PLATEID")
  arguments <- parse_args(parser, positional_arguments = 1)
}

#Specify the command line options
cl <- getCommandLineArgs()
plateID <- cl$args[1]
#plateID <- "LI8V01171"
message("processing plate: ",plateID)

#get feature files
files <- tibble(Full_filename = dir(plateID, pattern = ".txt",full.names = TRUE, recursive = TRUE)) %>%
  mutate(Filename = str_remove(Full_filename, ".*/"),
         PlateID = str_remove(Filename, "_.*"),
         Well = str_extract(Filename, "_[[:alpha:]][[:digit:]]+_"),
         Well = str_remove_all(Well, "_"),
         Spot = str_remove_all(Filename, ".*_|.txt"),
         Spot = as.integer(Spot),
         Well_Spot = paste0(Well,"_",Spot),
         Compartment = str_extract(Filename,"Cells|Cytoplasm|Nuclei"),
         Marker = str_remove_all(Filename,"^([^_]*_){3}|_.*"))

well_spots <- unique(files$Well_Spot)
addPrefix <-function(x, compartment, marker) paste0(compartment,"_GT_",marker,"_",x)

#Read in raw data with each row a cell
#0 area cyto compartments result in NA values in Cyto features
start_time <- Sys.time()
raw_data <- lapply(well_spots, function(well_spot){
  file_set <- files %>%
    filter(Well_Spot == well_spot)
  foo <- lapply(1:nrow(file_set), function(i){
    rd <- suppressWarnings(suppressMessages(read_delim(file_set$Full_filename[i],delim = " ",col_names = TRUE))) %>%
      select(-starts_with("X")) %>%
      select(-matches("Haralick|_PS_|_LBP_|_Rank")) %>%
      rename_all(.funs=addPrefix,compartment=file_set$Compartment[i],marker=file_set$Marker[i])
    rd$ObjectNumber <- rd[[paste0(file_set$Compartment[i],"_GT_",marker=file_set$Marker[i],"_Label")]]
    rd <- rd %>%
      select(-matches(paste0(file_set$Compartment[i],"_GT_",marker=file_set$Marker[i],"_Label"))) 
  }) %>% 
    plyr::join_all(by="ObjectNumber") %>%
    mutate(Well_Spot = well_spot,
           Well = str_remove(Well_Spot,"_.*"),
           Spot = str_remove(Well_Spot, ".*_"),
           Spot = as.integer(Spot))
  return(foo)
}) %>%
  bind_rows()
Sys.time()-start_time

#Read in the metadata
  inputPath <- "/lincs/share/lincs_user"
  metadataFiles <- list(annotMetadata=paste0(inputPath,"/",plateID,"/Analysis/",plateID,"_an2omero.csv"))
  metadata <- MEMA::getMetadata(metadataFiles, TRUE)
  
  imageIDs <-  paste0(inputPath,"/",plateID,"/Analysis/",plateID,"_imageIDs.tsv") %>%
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

level1_data <- raw_data %>%
  inner_join(metadata, by=c("Well", "Spot")) %>%
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

writeLevelData <- function(x, level){
  #Make a v4 directory if needed
  eppec_data_path <- getwd()
  output_path <- paste0(eppec_data_path,"/",plateID,"/Analysis")
  if(!dir.exists(output_path)) dir.create(output_path)
  write_csv(x, paste0(output_path,"/",plateID,"_level",level,".csv"))
}

#Add proportions for signals with multivariate gating and non-conforming gate values
addSpotProportions(level1_data)

#Calculate proportions for binary gated signals
gatedSignals <- grep("Proportion", grep("Positive|High",colnames(level1_data), value=TRUE), value=TRUE, invert=TRUE)

if(length(gatedSignals)>0){
  proportions <- level1_data[,lapply(.SD, calcProportion),by="Barcode,Well,Spot", .SDcols=gatedSignals]
  setnames(proportions,
           grep("Gated",colnames(proportions),value=TRUE),
           paste0(grep("Gated",colnames(proportions),value=TRUE),"Proportion"))
  level1_data <- left_join(level1_data, proportions, by = c("Barcode", "Well", "Spot"))
}

writeLevelData(level1_data,1)

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

level2_data <- full_join(level2_character_data,level2_numeric_data, by = c("Barcode", "Well", "Spot")) %>%
  right_join(level2_logical_data, by = c("Barcode", "Well", "Spot"))

writeLevelData(level2_data,2)
