# metaboAnalystR_data_processing
# Lazarina Butkovich, created 8/10/25
# This workflow performs data processing in MetaboAnalystR, in order to run the processed data through additional tools in the MetaboAnalyst web platform.
# 2 Sample Types: 
# MetaboAnalyst: use the Statistical Analysis [one factor] (https://www.metaboanalyst.ca/MetaboAnalyst/ModuleView.xhtml)


# Clean global environment
rm(list = ls())

#############################################
# Install and Load Required Packages
#############################################
# Function to install (if needed) and load packages
install_and_load <- function(pkg, bioc = FALSE, github = NULL, quiet = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    if (!is.null(github)) {
      if (!requireNamespace("devtools", quietly = TRUE))
        install.packages("devtools")
      devtools::install_github(github, build = TRUE,
                              build_vignettes = FALSE,
                              build_manual = FALSE)
    } else if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

# Set options for package installation
options(install.packages.compile.from.source = "always")

# List of required packages
# MetaboAnalystR dependencies (Bioconductor)
bioc_packages <- c(
  "impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", 
  "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", 
  "siggenes", "BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea"
)

# CRAN packages
cran_packages <- c(
  "devtools", "crmn", "httr", "qs", "readxl", "ggrepel", "ellipse", 
  "vegan", "pls", "rjson", "pheatmap", "ggplot2", "iheatmapr"
)

# Install and load all required packages
message("Installing and loading Bioconductor packages...")
invisible(lapply(bioc_packages, function(pkg) install_and_load(pkg, bioc = TRUE)))

message("Installing and loading CRAN packages...")
invisible(lapply(cran_packages, function(pkg) install_and_load(pkg)))

# Install and load MetaboAnalystR from GitHub
message("Installing and loading MetaboAnalystR...")
install_and_load("MetaboAnalystR", github = "xia-lab/MetaboAnalystR")

# Load required libraries explicitly for this script
library(ggplot2)  # For violin plots


##############
# Values to Change
##############
# Folders:
main_dir <- "C:\\Users\\lazab\\Documents\\github\\Open_Bootcamp_Collective_Bioinformatics\\OBC_Project2"
input_folder <- "input"
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

