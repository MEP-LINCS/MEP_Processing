#Structure files and directories for processing
library(tidyverse)

process8wellImageIDs <- function(barcode) {
  imageIDs <-  paste0(input_data_path,"/",barcode,"/Analysis/",barcode,"_imageIDs.tsv") %>%
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
}

#Execute script from mass directory
#Will copy txt files from flat organization to organized by barcode_well_compartment_marker
#Reads imageID metadata to get values based on imageID key
#Requires meatdata file in Analysis subdirectory with Well and Beacon columns

raw_data_path <- getwd()
raw_data_path <- "/graylab/share/thibaulg/Mema/HCC1143_COL1 - Features"
eppec_data_path <- "/eppec/storage/groups/heiserlab/lincs/"

#Replace channel names with staining set specific values
#ssT
channel_0 <- "Dapi"
channel_1 <- "KRT14"
channel_2 <- "EdU"
channel_3 <- "VIM"

#Image IDs flaged as known errors during upstream processing
known_errors <- c(3074163, 3101463)

#identify barcodes for getting imageID metadata
barcodes <- c("LI8V01171", "LI8V01172", "LI8V01173a", "LI8V01174","LI8V01175", "LI8V01176", "LI8V01177", "LI8V01178", "LI8V01179")
input_data_path <- "/lincs/share/lincs_user"

#Get imagedID metadata
imageID_metadata <- lapply(barcodes, process8wellImageIDs) %>%
  bind_rows()

#get feature files
files <- tibble(Full_filename = dir(raw_data_path, full.names = TRUE)) %>%
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

#join file names with metadata and create directory path
file_metadata <- left_join(files, imageID_metadata, by="ImageID") %>%
  mutate(dir_path = paste(PlateID,Well,Compartment, Channel,sep="_"),
         new_filename = paste0(dir_path,"_",Spot,".txt"))

copyMkdir <- function(full_filenames, dir_paths, new_filenames, plateIDs){
  plateID <- unique(plateIDs)
  dir_name <- unique(dir_paths)
  full_eppec_dir_name <- paste0(eppec_data_path,"/",plateID,"/Analysis/v5/",dir_name)
  if(!dir.exists(paste0(eppec_data_path,"/",plateID))) dir.create(paste0(eppec_data_path,"/",plateID))
  if(!dir.exists(paste0(eppec_data_path,"/",plateID,"/Analysis"))) dir.create(paste0(eppec_data_path,"/",plateID,"/Analysis"))
  if(!dir.exists(paste0(eppec_data_path,"/",plateID,"/Analysis/v5"))) dir.create(paste0(eppec_data_path,"/",plateID,"/Analysis/v5"))
  if(!dir.exists(full_eppec_dir_name)) dir.create(full_eppec_dir_name)
  foo <- lapply(paste0(full_filenames,"|",new_filenames), function(x){
    ffn <- str_remove(x,"[|].*")
    nfn <- paste0(eppec_data_path, plateID, "/Analysis/v5/", dir_name, "/", str_remove(x, ".*[|]"))
    file.copy(ffn, nfn)
  })
}

#create directories and copy/move files
res <- file_metadata %>%
  group_by(dir_path) %>%
  mutate(res = copyMkdir(Full_filename,dir_path,new_filename, PlateID))
