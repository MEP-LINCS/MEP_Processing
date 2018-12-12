#author: "Mark Dane"

library(tidyverse)

barcodes <- c("LI8X00907", "LI8V01150")

addCellCounts <- function(barcode){
  data_dir <- "/lincs/share/lincs_user/"
  l2 <- read_tsv(paste0(data_dir, barcode, "/Analysis/",barcode,"_Level2.tsv"))
  
  count_filenames <- dir(paste0(data_dir,barcode,"/Analysis/"), pattern="Count.*.txt", full.names = TRUE)
  res <- lapply(count_filenames, function(fn){
    cell_counts <- read_delim(fn, "\"", escape_double = FALSE, trim_ws = TRUE) %>%
      mutate(ImageID = Name,
             Cell_count = X3) %>%
      select(ImageID, Cell_count)
  }) %>%
    bind_rows() %>%
    right_join(l2)

}

l2List <- lapply(barcodes, function(barcode){
  l2 <- addCellCounts(barcode)
  write_tsv(l2, path=paste0(data_dir, barcode,"/Analysis/",barcode,"_Level2.tsv"))
})




