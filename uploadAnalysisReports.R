library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(rGithubClient)

analysisVersion <- "v1"

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
                  executed=thisScript)
  obj
}

for(cellLine in c("MCF7","PC3","YAPC")){
  for (ss in c("SS1","SS2","SS3")){
    dataDir <- paste0("/Users/dane/Documents/MEP-LINCS/AnalysisReports")
    # Take file names and turn into basic annotation set
    # Replace this with a better way to get basic annotations from 
    # a standardized source
    dataFiles <- data.frame(filename=list.files(path=dataDir, pattern = ".html", full.names = TRUE), stringsAsFactors = FALSE) %>%
      mutate(fileType="html",
             level="AnalysisReportDevel",
             basename=str_replace(filename, ".*/MEP-LINCS_Analysis_", "")) %>% 
      mutate(basename=str_replace(basename, "\\.html", "")) %>% 
      separate(basename, c("CellLine","analysisVersion","StainingSet" ))
    
    res <- dlply(dataFiles[, ], .(filename), uploadToSynapse, parentId=synapseDevelAnalysisDir)
    
  }
}