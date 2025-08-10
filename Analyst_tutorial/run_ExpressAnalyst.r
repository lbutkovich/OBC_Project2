############################################################

# run_ExpressAnalyst.r
# Lazarian Butkovich
# Created 7/25/25

# Based on "Web-based multi-omics integration using the Analyst software suite" (2024) Nature Protocols. https://doi.org/10.1038/s41596-023-00950-4

# ExpressAnalyst web-based GUI: https://www.expressanalyst.ca/ExpressAnalyst/home.xhtml
# ExpressAnalystR package: https://github.com/xia-lab/ExpressAnalystR
# Note, visual analytics are largely hosted on the web version; R version is primarily for statistical analysis.

# Tutorial dataset source: "Multi-omics profiling of living human pancreatic islet donors reveals heterogeneous beta cell trajectories towards type 2 diabetes." (2021) Nat. Metab. https://doi.org/10.1038/s42255-021-00420-9

# Background:
# Donor types, in order of disease progression:
# ND = no diabetes
# IGT = impaired glucose tolerance
# T3cD = type 3c diabetes (impaired glucose response caused by pancreatic injury or inflammation)
# T2D = type 2 diabetes
# Study: 135 donors, not all omics collected from all donors
# HBA1C: Hemoglobin A1C levels measured, reflects average blood glucose levels over the past 2-3 months

############################################################

###
# Values
###
OVERALL_FOLDER <- "Analyst_tutorial"
INPUT_FOLDER <- "input"
METADATA_FILENAME <- "metadata.csv"
RAW_RNASEQ_DATA_FILENAME <- "raw_rnaseq.txt"

OUTPUT_FOLDER <- "Multiomics_protocol"
OUTPUT_SUBFOLDER <- "output_ExpressAnalyst"


###
# Install and load required packages
###
# Load pacman
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)

pacman::p_load(igraph, RColorBrewer, qs, rjson, RSQLite)

# (added) make sure Bioconductor deps for ExpressAnalystR are available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Install needed BioC deps quietly/non-interactive
BiocManager::install(
  c("metagenomeSeq", "Biobase", "SummarizedExperiment"),
  ask = FALSE,
  update = FALSE
)

# Load devtools
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

# Install ExpressAnalystR without documentation
devtools::install_github(
  "xia-lab/ExpressAnalystR",
  dependencies = TRUE,
  upgrade      = "never",
  build        = TRUE,
  build_opts   = c("--no-resave-data", "--no-manual", "--no-build-vignettes")
)

library(ExpressAnalystR)

###
# Organize Input Files and Set wd to Output Subfolder
###
maindir <- getwd()
setwd(OVERALL_FOLDER)

# Create output folder if it doesn't exist
if (!dir.exists(OUTPUT_FOLDER)) {
  dir.create(OUTPUT_FOLDER)
}

# Create subfolder for ExpressAnalyst output if it doesn't exist
if (!dir.exists(file.path(OUTPUT_FOLDER, OUTPUT_SUBFOLDER))) {
  dir.create(file.path(OUTPUT_FOLDER, OUTPUT_SUBFOLDER))
}

# Copy the input files to the Output subfolder
# Build full paths for the two required input files (they live in "./input")
input_files <- file.path(INPUT_FOLDER, c(METADATA_FILENAME, RAW_RNASEQ_DATA_FILENAME))
if (!all(file.exists(input_files))) {
  stop(
    "Input files not found: ",
    paste(input_files[!file.exists(input_files)], collapse = ", ")
  )
}
file.copy(
  from      = input_files,
  to        = file.path(OUTPUT_FOLDER, OUTPUT_SUBFOLDER),
  overwrite = TRUE
)

# Set wd to the output subfolder. This way ExpressAnalystR will save output files there automatically.
setwd(file.path(OUTPUT_FOLDER, OUTPUT_SUBFOLDER))
outputdir <- getwd()

# Build absolute paths – ExpressAnalystR is much happier with them
raw_file_path  <- normalizePath(RAW_RNASEQ_DATA_FILENAME , winslash = "/")
meta_file_path <- normalizePath(METADATA_FILENAME        , winslash = "/")

# +++++ NEW: verify that the two input files are really here +++++
missing <- c(RAW_RNASEQ_DATA_FILENAME, METADATA_FILENAME)[
  !file.exists(c(RAW_RNASEQ_DATA_FILENAME, METADATA_FILENAME))
]
if (length(missing)) {
  stop(
    "After copy step the following input file(s) are still missing in ",
    outputdir, ": ", paste(missing, collapse = ", ")
  )
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###
# Initialize R Objects
###
Init.Data(FALSE)

# Set analysis type to single gene expression matrix
SetAnalType("onedata")

# ---------- data import ----------
dataSets <- tryCatch(
  ReadTabExpressData(raw_file_path, meta_file_path, "false", "default"),
  error = function(e) {
    stop("ReadTabExpressData failed – check that the two input files are in correct tab-delimited format: ",
         conditionMessage(e))
  }
)

# ---------- gene annotation (optional fall-back) ----------
dataSets <- tryCatch(
  PerformDataAnnot(raw_file_path, "hsa", "count", "embl_gene", "sum"),
  error = function(e) {
    warning("PerformDataAnnot failed – proceeding without annotation: ",
            conditionMessage(e))
    dataSets    # return the object obtained from the previous step
  }
)

SummarizeQC(raw_file_path, "proc", 0.1)

# Determine how to set the order of disease progression
# to-do

dataSets <- PerformNormalization(raw_file_path, "RLE", 15, 4, "true", "false", "sum")
SummarizeQC(raw_file_path, "proc", 0.1)
PlotDataNsig(raw_file_path, "qc_norm_nsig80_1_", "72", "png")
PlotDataDendrogram(raw_file_path, "qc_norm_dendrogram_1_", 0.1, "72", "png")
PlotDataGini(raw_file_path, "qc_norm_gini_1_", 0.95, "72", "png")
SetSelectedMetaInfo(RAW_RNASEQ_DATA_FILENAME,"Diagnosis", "NA", F)
SetSelectedMetaInfo(RAW_RNASEQ_DATA_FILENAME,"Diagnosis", "NA", F)
SetSelectedMetaInfo(RAW_RNASEQ_DATA_FILENAME,"Diagnosis", "NA", F)
adj.vec <- c("Age", "BMI")
MultiCovariateRegression(raw_file_path, "Diagnosis", "ND", "T3D", "NA", FALSE)


