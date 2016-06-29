## ---- Basic setup ----------------------------------------

# run below codes anytime opening this file
setwd("/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge");
rm(list = ls())
load(".RData")

download_path = "/media/ducdo/UUI/Bioinformatics/DMI_DreamChallenge"
downloaded_folder <- "DREAM11_subchallenge1_initial_networks"
file_names = paste(download_path,"/",downloaded_folder,"/",list.files(downloaded_folder)
                   ,sep = "")

#setwd("D:/Bioinformatics/DMI_DreamChallenge")
# source('http://depot.sagebase.org/CRAN.R')
# pkgInstall(c("synapseClient"))

## ---- Log in Synapse Server ------------------------------
library(synapseClient)
synapseLogin(username = '...',password = '...')

## ----  Sub-challenge 1  ------------------------------------

# Obtain a pointer and download the data
entity <- synGet(id='syn6169211', downloadLocation = download_path)
#cat(getFileLocation(entity), "\n");unzip(getFileLocation(entity))

ppi_1 = read.table(file = file_names[1], header = TRUE, sep = "\t")
ppi_2 = read.table(file = file_names[2], header = TRUE, sep = "\t")
signal_directed = read.table(file = file_names[3], header = TRUE, sep = "\t")
coexpr = read.table(file = file_names[4], header = TRUE, sep = "\t")
cancer = read.table(file = file_names[5], header = TRUE, sep = "\t")
homology = read.table(file = file_names[6], header = TRUE, sep = "\t")

## ----  Sub-challenge 2  ------------------------------------
entity2 <- synGet(id='syn6169210', downloadLocation = download_path)
#cat(getFileLocation(entity2), "\n");unzip(getFileLocation(entity2))
downloaded_folder2 <- "DREAM11_subchallenge2_initial_networks"
file_names2 = paste(download_path,"/",downloaded_folder2,"/",list.files(downloaded_folder2)
                   ,sep = "")
ppi_aligned_1 = read.table(file = file_names2[1], header = TRUE, sep = "\t")
ppi_aligned_2 = read.table(file = file_names2[2], header = TRUE, sep = "\t")
signal_directed_aligned = read.table(file = file_names2[3], header = TRUE, sep = "\t")
coexpr_aligned = read.table(file = file_names2[4], header = TRUE, sep = "\t")
cancer_aligned = read.table(file = file_names2[5], header = TRUE, sep = "\t")
homology_aligned = read.table(file = file_names2[6], header = TRUE, sep = "\t")
