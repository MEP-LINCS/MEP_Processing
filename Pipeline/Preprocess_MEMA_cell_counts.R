#author: "Mark Dane"

library(tidyverse)

barcodes <- c("LI8X00907", "LI8V01150")

addCellCounts <- function(barcode){
  data_dir <- "/lincs/share/lincs_user/"
  l2 <- read_tsv(paste0(data_dir, barcode, "/Analysis/",barcode,"_Level2.tsv"))
  count_filename <- paste0(data_dir,barcode,"/Analysis/Count.txt")
  if(!file.exists(count_filename)) stop(paste0(count_filename, " file not found"))
  cell_counts <- read_delim(count_filename, "\"", escape_double = FALSE, trim_ws = TRUE) %>%
    mutate(ImageID = Name,
           Cell_count = X3) %>%
    select(ImageID, Cell_count) %>%
    right_join(l2)
}

l2List <- lapply(barcodes, function(barcode){
  l2 <- addCellCounts(barcode)
  write_tsv(l2, path=paste0(data_dir, barcode,"/Analysis/",barcode,"_Level2.tsv"))
})




