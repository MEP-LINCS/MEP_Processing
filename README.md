# Microenvironment Perturbation - Library of Integrated Network Signatures (MEP-LINCS) Processing Pipeline

This repository holds scripts and utilities to preprocess, QA, normalize and explore Microenvironment Microarray (MEMA) experiments. MEMAs are used in Oregon Health and Science's (OHSU) MEP-LINCS project. An overveiw of the project along with all data and protocols are avilable at https://www.synapse.org/#!Synapse:syn2862345/wiki/72486. These scripts are open souce R code that use the MEMA package available from https://github.com/MEP-LINCS/MEMA.

The MEMA package includes vignettes that demonstrate how to run the pipeline on a subset of the data. The documentation below details how to use the scripts in this repository to genrate a full set of MEP-LINCS datasets and reports.  

The Pipeline folder in the repo contains R and shell scripts used to generate the MEP-LINCS datasets. The Reports folder has R knitr scripts that generate MEP-LINCS QA and Analysis reports. These reports contain interactive elements that help explore the data.  
  
#### Dependencies

This pipeline uses version 1.0.0 of the `MEMA` package (https://github.com/MEP-LINCS/MEMA/releases/tag/v1.0.0).
  
#### Pipeline Preprocessing Scripts

There are five preprocessing scripts that each read in data and generate the next level. These scripts run from the command line with arguments that control their operation. Help on the scripts arguments is shown when including a -h or --help as a command line argument.

<br>

#### Raw to Level 1

The PreprocessMEMACell.R script reads in one plate's data and metadata, merges them, gates some of the signals and writes out the results as Level 1 data in a tab-delimited file. The input data comes from an image segmentation pipeline. The metadata data comes from either OHSU's Annot (An!) database or from structured Excel files. The output can be written to the Synapse website or to a local server. Each row in the output file represents data from one cell in the experiment. The columns are either metadata about the cell, intensity and morphology measurements from the segmentation pipeline or derived features computed within the script. Help for the PreprocessMEMACell.R script is:


```

Usage: PreprocessMEMACell.R [options] BARCODE OUTPUTFILE


Options:
	-v, --verbose
		Print extra output

	-e, --excelMetadata
		Get metadata from Excel files instead of from the An! database

	-i PATH, --inputPath=PATH
		Path to local input data directory or Synapse ID for a File View.

	-r RAWDATAVERSION, --rawDataVersion=RAWDATAVERSION
		Raw data version from local server [default "v2"]

	--synapseStore=SYNAPSEID
		Store output file in Synapse directory (provide Synapse ID of Folder to store).

	-h, --help
		Show this help message and exit
		
```
The required BARCODE command line argument identifies the plate to be processed. Values will be similar to LI8X00641. The OUTPUTFILE argument is the full path on the local server for the output file. On many systems an output name that starts with "/tmp/" can be used if the output is to be stored on Synapse and not retained on the local server.  

<br>

#### Level 1 to Level 2

PreprocessMEMALevel2.R reads in the Level 1 data and median summarizes it from the cell level to the spot level. Additionally, this script computes the proportions of cells in the gates. Help for the PreprocessMEMALevel2.R script is:  

```
Usage: PreprocessMEMALevel2.R [options] BARCODE OUTPUTFILE


Options:
	-v, --verbose
		Print extra output

	-i PATH, --inputPath=PATH
		Path to local input data directory or Synapse ID for a File View.

	--synapseStore=SYNAPSEID
		Store output file in Synapse directory (provide Synapse ID of Folder to store).

	-h, --help
		Show this help message and exit

```
The required command line arguments are the same as in PreprocessMEMACell.R.

<br>

#### Level 2 to Level 3

PreprocessMEMALevel3.R reads in the Level 2 data and normalizes it to create Level 3 data. The RUVLoessResidual normalization method uses the k parameter to determine how many factors of unwanted variation to remove. The raw and normalized, spot-level data are written out in tab-delimited files where each row is a spot. Help for the PreprocessMEMALevel3.R script is:  

```

Usage: PreprocessMEMALevel3.R [options] STUDY OUTPUTFILE


Options:
	-v, --verbose
		Print extra output

	-k NUMBER, --k=NUMBER
		Number of factors to use in RUV normalization [default 256]

	-i PATH, --inputPath=PATH
		Path to local input data directory or Synapse ID for a File View.

	--synapseStore=SYNAPSEID
		Store output file in Synapse directory (provide Synapse ID of Folder to store).

	-h, --help
		Show this help message and exit


```
The required STUDY command line argument identifies the MEP-LINCS study to be processed. A study is a collection of plates that will be normalized together. These values must match a value in the StudyName field of the an2study.tsv file at https://www.synapse.org/#!Synapse:syn8440875. The OUTPUTFILE argument is the full path on the local server for the output file.  

<br>

#### Level 3 to Level 4

PreprocessMEMALevel4.R reads in the Level 3 data and median summarizes the replicates to create Level 4 data. The raw and normalized, replicate-level data are written out in tab-delimited files where each row contains values from one microenvironment perturbation or MEP. MEPs are combinations of ligands and ECM proteins, possibly with drugs. Help for the PreprocessMEMALevel4.R script is:  

```

Usage: PreprocessMEMALevel4.R [options] STUDY OUTPUTFILE


Options:
	-v, --verbose
		Print extra output

	-i PATH, --inputPath=PATH
		Path to local input data directory or Synapse ID for a File View.

	--synapseStore=SYNAPSEID
		Store output file in Synapse directory (provide Synapse ID of Folder to store).

	-h, --help
		Show this help message and exit

```
The required command line arguments are the same as in PreprocessMEMALevel3.R.   

<br>

#### Combined Staining Sets

All MEP-LINCS datasets include DAPI to identify the nucleii. This generates replicates across the datasets that can be combined to cretae more robust signals. PreprocessCombinedStudies.R reads in the Level 2 data from 2 or 3 studies, normalizes the common intensity and morphology signals and then appends the unique level 3 signals from the rest of the datasets. This Level 3 data is then median summarized across the MEP replicates to create Level 4 data. The raw and normalized, Level 3 and Level 4 data are written out in tab-delimited files as staining set SSC. Help for the PreprocessCombinedStudies.R script is:  

```

Usage: PreprocessCombinedStudies.R [options] STUDY STUDY [STUDY]


Options:
	-v, --verbose
		Print extra output

	-k NUMBER, --k=NUMBER
		Number of factors to use in RUV normalization [default 256]

	-i SYNID, --inputSynID=SYNID
		Synapse ID for a File View with input data.

	--synapseStore=SYNAPSEID
		Store output file in Synapse directory (provide Synapse ID of Folder to store).

	-h, --help
		Show this help message and exit
		
```
The required command line arguments are two or three study names as described in PreprocessMEMALevel3.R.  


<br>

#### Reports
The \*.Rmd files in the Reports folder of this repo are sourced by the MEP-LINCS_ReportsWrapper.R or MEP-LINCSAnalysisStoryboardWrapper.R scripts when they are executed on a local sever containing the MEP-LINCS datasets. 


		
