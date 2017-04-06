library(readr)
library(stringr)
an2study <- read_delim("~/GrayLabData/dane/MEP-LINCS/MEP_LINCS/archive/an2study.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
an2study$StudyName <- tolower(an2study$StudyName)
fwrite(an2study,"an2study.tsv",sep="\t")
