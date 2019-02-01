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
  arguments <- parse_args(parser, positional_arguments = 3)
}

#Debug add signif reduction
writeLevelData <- function(x, level){
  write_csv(x, paste0(output_path,plateID,"/Analysis/GT1/",plateID,"_level",level,".csv"))
}

#Specify the command line options
if(!interactive()){
  cl <- getCommandLineArgs()
  input_path <- cl$args[1]
  output_path <- cl$args[2]
  plateID <- cl$args[3]
} else {
  input_path <- "/lincs/share/lincs_user/"
  output_path <- "/lincs/share/lincs_user/"
  plateID <- "LI8V01171"
}

message("processing GT1 data in plate: ",input_path,plateID)

#get feature files
files <- tibble(Full_filename = dir(  input_path <- paste0(input_path,plateID,"/Analysis/GT1/"), pattern = "level_0.csv",full.names = TRUE)) %>%
  mutate(Filename = str_remove(Full_filename, ".*/"),
         PlateID = str_remove(Filename, "_.*"),
         Well = str_extract(Filename, "_[[:digit:]]+_"),
         Well = str_remove_all(Well, "_"),
         Well = as.integer(Well),
         Well = wellAN(2,4)[Well])

#Read in raw data with each row a cell
#0 area cyto compartments result in NA values in Cyto features
raw_data <- map(files$Full_filename,read_csv) %>%
  bind_rows

#Read in the metadata
  metadataFiles <- list(annotMetadata=paste0(input_path, plateID, "/Analysis/", plateID, "_an2omero.csv"))
  metadata <- MEMA::getMetadata(metadataFiles, TRUE) %>%
    mutate_if(is.logical, as.character)
  
  imageIDs <-  paste0(input_path, plateID,"/Analysis/", plateID,"_imageIDs.tsv") %>%
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
  inner_join(full_metadata, by=c("ImageID")) %>%
  mutate(ObjectNumber = Label,
         Nuclei_CP_Intensity_MedianIntensity_Dapi = Nuclei_GT_Dapi_MeanIntensity,
         Nuclei_CP_Intensity_IntegratedIntensity_Dapi = Nuclei_GT_Dapi_MeanIntensity*Nuclei_GT_Dapi_Surface,
         Nuclei_CP_Intensity_MedianIntensity_EdU = Cells_GT_EdU_MeanIntensity,
         Nuclei_CP_Intensity_IntegratedIntensity_EdU = Cells_GT_EdU_MeanIntensity*Cells_GT_EdU_Surface,
         Cytoplasm_CP_Intensity_MedianIntensity_KRT14 = Cytoplasm_GT_KRT14_MeanIntensity,
         Cytoplasm_CP_Intensity_IntegratedIntensity_KRT14 = Cytoplasm_GT_KRT14_MeanIntensity*Cytoplasm_GT_KRT14_Surface,
         Cytoplasm_CP_Intensity_IntegratedIntensity_VIM = Cytoplasm_GT_VIM_MeanIntensity*Cytoplasm_GT_VIM_Surface,
         Cytoplasm_CP_Intensity_MedianIntensity_VIM = Cytoplasm_GT_VIM_MeanIntensity,
         Well_Ligand = paste0(Well,"_",Ligand)) %>%
  data.table::data.table() %>%
  gateCells()

#assign sectors
assign_sectors <- function(df){
  #assign a sector  based on concentric regions of equal area
  #work in pixels where 1600 pixels = 500 um
  spot_rad <- 200*(1600/500)
  mid_rad <- sqrt(spot_rad^2/1.5)
  ctr_rad <- sqrt(mid_rad^2/2)
  plate_median_x <-  median(level1_data$Nuclei_GT_Dapi_CentroidX)
  plate_median_y <-  median(level1_data$Nuclei_GT_Dapi_CentroidY)
  spot_center_tol <- 50
  
  #convert to polar coordinates
  #let the spot center wander +/-50 pixels about the typical center
  #set sector logicals
  df_with_sectors <- df %>%
    group_by(Well, Spot) %>%
    mutate(Spot_center_x = median(Nuclei_GT_Dapi_CentroidX),
           Spot_center_x = scales::squish(Spot_center_x, c(plate_median_x-spot_center_tol, plate_median_x+spot_center_tol)),
           Spot_center_y = median(Nuclei_GT_Dapi_CentroidY),
           Spot_center_y = scales::squish(Spot_center_y, c(plate_median_y-spot_center_tol,plate_median_y+spot_center_tol)),
           Nuclei_GT_Dapi_SpotCentroidX = Nuclei_GT_Dapi_CentroidX-Spot_center_x,
           Nuclei_GT_Dapi_SpotCentroidY = Nuclei_GT_Dapi_CentroidY-Spot_center_y,
           Nuclei_GT_Dapi_SpotCentroidR = sqrt(Nuclei_GT_Dapi_SpotCentroidX^2+Nuclei_GT_Dapi_SpotCentroidY^2),
           Nuclei_GT_Dapi_SpotCentroidTheta = atan2(Nuclei_GT_Dapi_SpotCentroidY, Nuclei_GT_Dapi_SpotCentroidX)*180/pi,
           Nuclei_GT_spatial_SectorCenter = Nuclei_GT_Dapi_SpotCentroidR <= ctr_rad,
           Nuclei_GT_spatial_SectorMiddle = Nuclei_GT_Dapi_SpotCentroidR > ctr_rad&Nuclei_GT_Dapi_SpotCentroidR <= mid_rad,
           Nuclei_GT_spatial_SectorOuter = Nuclei_GT_Dapi_SpotCentroidR > mid_rad & Nuclei_GT_Dapi_SpotCentroidR <= spot_rad,
           Nuclei_GT_spatial_Quadrant1 = Nuclei_GT_Dapi_SpotCentroidTheta >= 0 & Nuclei_GT_Dapi_SpotCentroidTheta < 90,
           Nuclei_GT_spatial_Quadrant2 = Nuclei_GT_Dapi_SpotCentroidTheta >= 90,
           Nuclei_GT_spatial_Quadrant3 = Nuclei_GT_Dapi_SpotCentroidTheta <= -90,
           Nuclei_GT_spatial_Quadrant4 = Nuclei_GT_Dapi_SpotCentroidTheta < 0 & Nuclei_GT_Dapi_SpotCentroidTheta > -90
           ) %>%
  ungroup() %>%
    data.table::data.table() 
}

level1_data <- assign_sectors(level1_data)

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

writeLevelData(level1_data, 1)

level2_numeric_data <- level1_data %>%
  group_by(Well,Spot) %>%
  mutate(Spot_PA_SpotCellCount=n()) %>%
  summarise_if(is.numeric, median)  %>%
  ungroup()
#   mutate_if(is.numeric, signif, digits=4)

level2_character_data <- level1_data %>%
  group_by(Well,Spot) %>%
  summarise_if(is.character, unique) %>%
  ungroup()

# level2_logical_data <- level1_data %>%
#   group_by(Well,Spot) %>%
#   summarise_if(is.logical, unique) %>%
#   ungroup()

level2_spatial_data <- level1_data %>%
    group_by(Well, Spot) %>%
    summarise(Nuclei_GT_spatial_SectorCenterProportion = sum(Nuclei_GT_spatial_SectorCenter)/n(),
              Nuclei_GT_spatial_SectorMiddleProportion = sum(Nuclei_GT_spatial_SectorMiddle)/n(),
              Nuclei_GT_spatial_SectorOuterProportion = sum(Nuclei_GT_spatial_SectorOuter)/n(),
              Nuclei_GT_spatial_Quadrant1Proportion = sum(Nuclei_GT_spatial_Quadrant1)/n(),
              Nuclei_GT_spatial_Quadrant2Proportion = sum(Nuclei_GT_spatial_Quadrant2)/n(),
              Nuclei_GT_spatial_Quadrant3Proportion = sum(Nuclei_GT_spatial_Quadrant3)/n(),
              Nuclei_GT_spatial_Quadrant4Proportion = sum(Nuclei_GT_spatial_Quadrant4)/n()
              ) %>%
  ungroup()

level2_data <- full_join(level2_character_data, level2_numeric_data, by = c("Well", "Spot")) %>%
  #right_join(level2_logical_data, by = c("Well", "Spot")) %>%
  right_join(level2_spatial_data, by = c("Well", "Spot"))
  
writeLevelData(level2_data,2)
