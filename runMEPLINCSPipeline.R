
#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 9/30/2015

##Introduction

#   The MEP-LINCs dataset contains imaging data from a Nikon automated microscope that is analyzed with a CellProfiler pipeline.
# 
# Part of this preprocessing of the dataset will be deprecated when the merging of the data and metadata happens within the CellProfiler part of the pipeline. For now, the metadata about the ECM proteins is read from the GAL file and the metadata about the wells (cell line, stains and ligands) is read from Excel spreadsheets.

source("MEPLINCSFunctions.R")
library("limma")#read GAL file and strsplit2
library("MEMA")#merge, annotate and normalize functions
library("data.table")#fast file reads, data merges and subsetting
library("parallel")#use multiple cores for faster processing

mclapply(c("SS2"), preprocessMEPLINCS, cellLine = "MCF7", mc.cores=detectCores())
