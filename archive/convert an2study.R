library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(magrittr)
library(synapseClient)

miscDir <- 'syn8698827'

#Read the an2study file from the server
an2study <- read_delim("~/GrayLabData/dane/MEP-LINCS/MEP_LINCS/archive/an2study.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
if("PlateIDs" %in% colnames(an2study)) names(an2study)[names(an2study)=="PlateIDs"] <- "Barcode"
an2study <- select(an2study,StudyName,Barcode)
an2study$StudyName <- tolower(an2study$StudyName)
fwrite(an2study, "~/GrayLabData/dane/MEP-LINCS/MEP_LINCS/archive/an2study.tsv", sep="\t")

#logon to Synapse
synapseLogin()
#Write file up to an2study synID
synFile <- File("~/GrayLabData/dane/MEP-LINCS/MEP_LINCS/archive/an2study.tsv", parentId=miscDir)
res <- synStore(synFile)

