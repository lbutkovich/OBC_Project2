# metaboAnalystR_data_processing
# Lazarina Butkovich, created 8/10/25
# This workflow performs data processing in MetaboAnalystR, in order to run the processed data through additional tools in the MetaboAnalyst web platform.
# 2 Sample Types: 
# MetaboAnalyst: use the Statistical Analysis [metadata table] (https://www.metaboanalyst.ca/MetaboAnalyst/ModuleView.xhtml)


# Clean global environment
rm(list = ls())

############################################
# Install and Load Required Packages
############################################
# # Function to install (if needed) and load packages
# install_and_load <- function(pkg, bioc = FALSE, github = NULL, quiet = TRUE) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     message(paste("Installing package:", pkg))
#     if (!is.null(github)) {
#       if (!requireNamespace("devtools", quietly = TRUE))
#         install.packages("devtools")
#       devtools::install_github(github, build = TRUE,
#                               build_vignettes = FALSE,
#                               build_manual = FALSE)
#     } else if (bioc) {
#       if (!requireNamespace("BiocManager", quietly = TRUE))
#         install.packages("BiocManager")
#       BiocManager::install(pkg)
#     } else {
#       install.packages(pkg)
#     }
#   }
#   suppressPackageStartupMessages(
#     library(pkg, character.only = TRUE)
#   )
# }

# # Set options for package installation
# options(install.packages.compile.from.source = "always")

# # List of required packages
# # MetaboAnalystR dependencies (Bioconductor)
# bioc_packages <- c(
#   "impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", 
#   "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", 
#   "siggenes", "BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea"
# )

# # CRAN packages
# cran_packages <- c(
#   "devtools", "crmn", "httr", "qs", "readxl", "ggrepel", "ellipse", 
#   "vegan", "pls", "rjson", "pheatmap", "ggplot2", "iheatmapr"
# )

# # Install and load all required packages
# message("Installing and loading Bioconductor packages...")
# invisible(lapply(bioc_packages, function(pkg) install_and_load(pkg, bioc = TRUE)))

# message("Installing and loading CRAN packages...")
# invisible(lapply(cran_packages, function(pkg) install_and_load(pkg)))

# # Install and load MetaboAnalystR from GitHub. Specified a working version of MetaboAnalystR (can update as needed)
# message("Installing and loading MetaboAnalystR...")
# devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, ref = "5aee8b4f0d27c27864198a6fd99414575d693836", build_vignettes = FALSE, build_manual =F)

# # Load required libraries explicitly for this script
# library(ggplot2)  # For violin plots


##############
# Values to Change
##############
# Folders:
main_dir <- "C:\\Users\\lazab\\Documents\\github\\Open_Bootcamp_Collective_Bioinformatics\\OBC_Project2"
output_folder <- "output"
metaboanalystR_output_folder <- "MetaboAnalystR_output"

# Create metaboanalystR_output_folder in output_folder if it doesn't exist
metaboanalystR_output_folder_dir <- paste(main_dir, output_folder, metaboanalystR_output_folder, sep = "\\")
if (!dir.exists(metaboanalystR_output_folder_dir)) {
  dir.create(metaboanalystR_output_folder_dir)
}

# Metabolite Filenames:
# Formatted metabolite intensity data in output folder, with 2 batches to be processed separately then combined.
metabolites_data_batch_1_filename <- "metabolites_input_for_analyst_batch_1.csv" 
metabolites_data_batch_2_filename <- "metabolites_input_for_analyst_batch_2.csv"
# metadata with sample type information per SampleID, in input folder
metabolites_metadata_filename <- "metabolites_metadata.xlsx" 


##############
# Organize Output Directory
##############
# Empty the MetaboAnalystR_output folder 
file.remove(list.files(metaboanalystR_output_folder_dir, full.names = TRUE))

# Set directory to MetaboAnalystR output folder so that created files get placed there
setwd(metaboanalystR_output_folder_dir)


##############
# Load MetaboAnalystR and Data Table
##############
# Load MetaboAnalystR
library(MetaboAnalystR)

# Initialize data object mSet for MetaboAnalystR
# data.type: pktable = peak intensity table
mSet <- InitDataObjects("pktable", "stat", FALSE);

# Read metabolite data for batch
# "colu" = sampleIDs in columns, unpaired
mSet <- Read.TextData(mSet, paste(main_dir, output_folder, metabolites_data_batch_2_filename, sep = "\\"), "colu", "disc");

# Sanity check data
mSet <- SanityCheckData(mSet)


##############
# Replace Missing Values
##############
# Replace missing values with 1/2 of the value of the smallest non-zero value for each feature across samples
mSet <- ReplaceMin(mSet)


##############
# Perform Data Filtering
##############
# Filter features based on QC samples (remove samples with RSD >40% in QC samples). No low-variance filter or low-abundance filter applied
mSet <- FilterVariable(mSet, "T", 40, "iqr", 0, "mean", 0)


##############
# Normalize Data
##############
mSet<-PreparePrenormData(mSet)

mSet<-Normalization(mSet, "CompNorm", "LogNorm", "AutoNorm", "IS", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 150, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 150, width=NA)
