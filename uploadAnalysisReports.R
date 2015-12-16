library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(rGithubClient)

synapseLogin()

repo <- getRepo("MEP-LINCS/MEP_LINCS_Pilot", ref="branch", refName="updateReports")
thisScript <- getPermlink(repo, "uploadAnalysisReports.R")

synapseDevelAnalysisDir <- "syn5555312"


# Take row of data frame and remove file name
# Convert to a list to use as Synapse annotations
toAnnotationList <- function(x) {
  as.list(x %>% select(-filename))
}

# Take row of data frame with filename and annots
# Upload to Synapse and set annotations
uploadToSynapse <- function(x, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  obj <- synStore(obj, 
                  activityName="Upload", 
                  contentType="text/tab-separated-values",
                  executed=thisScript)
  obj
}

for(cellLine in c("MCF7","PC3","YAPC")){
  for (ss in c("SS1","SS2","SS3")){
    dataDir <- paste("/Users/dane/Documents/MEP-LINCS",cellLine,ss,"AnnotatedData", sep = "/")
    # Take file names and turn into basic annotation set
    # Replace this with a better way to get basic annotations from 
    # a standardized source
    dataFiles <- data.frame(filename=list.files(path=dataDir, pattern = ".txt", full.names = TRUE), stringsAsFactors = FALSE) %>%
      mutate(fileType="tsv",
             basename=str_replace(filename, ".*/", "")) %>% 
      mutate(basename=str_replace(basename, "\\.txt", "")) %>% 
      separate(basename, c("CellLine","StainingSet","level" ))
    dataFiles$level <- str_replace(dataFiles$level, "Level", "")
    
    res <- dlply(dataFiles[, ], .(filename), uploadToSynapse, parentId=synapseAnnotatedDataDir)
    
  }
}