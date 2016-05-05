# Make the pre-release manifest table
library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(tidyr)
library(knit2synapse)

library(synapseClient)
synapseLogin()

annotatedFolder <- "syn5706203"
rawFolder <- "syn5706233"
reportFolder <- "syn5007815"
metadataFolder <- "syn4997970"

query <- paste('select id,name,versionNumber,Level,Barcode,CellLine,StainingSet,Location,Well',
               'from file where parentId=="%s"')

annotFiles <- synQuery(sprintf(query, annotatedFolder), blockSize = 400)$collectAll()
rawFiles <- synQuery(sprintf(query, rawFolder), blockSize = 400)$collectAll()
reportFiles <- synQuery(sprintf(query, reportFolder), blockSize = 400)$collectAll()
metadataFiles <- synQuery(sprintf(query, metadataFolder), blockSize = 400)$collectAll()

allFiles <- rbind(annotFiles, rawFiles, reportFiles, metadataFiles)
colnames(allFiles) <- gsub(".*\\.", "", colnames(allFiles))

allFiles <- allFiles %>%
  mutate(versionNumber=as.numeric(versionNumber)) %>% 
  arrange(Level, CellLine, StainingSet) %>% 
  select(id,versionNumber,name,Level,Barcode,CellLine,StainingSet,Location,Well)

tableName <- sprintf("Release %s", format(Sys.time(), "%d-%b-%Y %H%M%S"))
tblCols <- as.tableColumns(allFiles)
schema <- TableSchema(name=tableName, columns=tblCols$tableColumns, 
                      parent="syn2862345")
tbl <- synStore(Table(tableSchema = schema, values=allFiles))

manifestTableId <- tbl@schema@properties$id

rmarkdown::render("./releaseWiki.Rmd", 
                  params = list(manifestTableId=manifestTableId),
                  clean=FALSE)

knit2synapse::knitfile2synapse(file="./releaseWiki.knit.md", 
                               owner='syn4215176', overwrite=TRUE,
                               knitmd=FALSE)
