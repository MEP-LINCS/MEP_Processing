# Make the pre-release manifest table
library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(tidyr)

library(synapseClient)
synapseLogin()

annotatedFolder <- "syn5706203"
rawFolder <- "syn5706233"
reportFolder <- "syn5007815"

query <- paste('select id,name,versionNumber,Level,Barcode,CellLine,StainingSet,Location,Well',
               'from file where parentId=="%s"')

annotFiles <- synQuery(sprintf(query, annotatedFolder), blockSize = 400)$collectAll()
rawFiles <- synQuery(sprintf(query, rawFolder), blockSize = 400)$collectAll()
reportFiles <- synQuery(sprintf(query, reportFolder), blockSize = 400)$collectAll()

allFiles <- rbind(annotFiles, rawFiles, reportFiles) 
colnames(allFiles) <- gsub(".*\\.", "", colnames(allFiles))

allFiles <- allFiles %>%
  # filter(CellLine %in% c("PC3", "MCF7"), StainingSet %in% c("SS2")) %>% 
  mutate(versionNumber=as.numeric(versionNumber)) %>% 
  arrange(Level, CellLine, StainingSet) %>% 
  select(id,versionNumber,name,Level,Barcode,CellLine,StainingSet,Location,Well)

tableName <- "Pre-release manifest"
tblCols <- as.tableColumns(allFiles)
schema <- TableSchema(name=tableName, columns=tblCols$tableColumns, 
                      parent="syn2862345")
tbl <- synStore(Table(tableSchema = schema, values=allFiles))

