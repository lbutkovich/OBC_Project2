"""
format_inputs_for_analyst.py
Lazarina Butkovich
created 8/9/25

Takes as input supplementary datasets (.xlsx) from Zhang et al. (2025) for metabolomics and proteomics data.
Generates .xlsx files for each omics type, formatted for input into the Analyst tools for one-factor analysis (MetaboAnalyst for metabolomics, ExpressAnalyst for proteomics).

Zhang et al. (2025): https://doi.org/10.1038/s41467-024-52262-0
"""

import pandas as pd
import os
from os.path import join as pjoin

"""
Functions
"""
def create_batch_dataframe(metabolites_data_raw, sample_ids):
    """Create a batch dataframe with SampleID column and specified sample columns."""
    columns = ["ID"] + sample_ids
    batch_data = metabolites_data_raw[columns].copy()
    batch_data.rename(columns={"ID": "SampleID"}, inplace=True)
    return batch_data

def add_group_labels_row(dataframe, sample_ids, metabolites_metadata):
    """Add group labels as the first row in the dataframe."""
    # Create a mapping from SampleID to Group
    sample_to_group = dict(zip(metabolites_metadata["SampleID"], metabolites_metadata["Group"]))
    
    # Create the group labels row
    group_labels = ["Group"]  # First column label
    for sample_id in sample_ids:
        group_labels.append(sample_to_group.get(sample_id, "Unknown"))

    # Insert the group labels row as the first row (index 0)
    group_row = pd.DataFrame([group_labels], columns=dataframe.columns)
    dataframe = pd.concat([group_row, dataframe], ignore_index=True)

    return dataframe

"""
Values
"""
# Folder names:
INPUT_FOLDER = "input"
OUTPUT_FOLDER = "output"

# Input Filenames:
# Zhang et al. (2025) Supplementary Data 2 = "Raw intensity values of 326 annotated metabolites for 580 samples in discovery study cohort" 
# metabolites_raw.xlsx: "41467_2024_52262_MOESM5_ESM.xlsx" without header
# metabolites_metadata.xlsx: header info for "41467_2024_52262_MOESM5_ESM.xlsx" 
METABOLITES_DATA_RAW_FILENAME = "metabolites_raw.xlsx"
METABOLITES_METADATA_FILENAME = "metabolites_metadata.xlsx"

# Zhang et al. (2025) Supplementary Data 13: "Raw normalized intensity values of 524 annotated proteins for 183 samples"
# proteins_raw.xlsx = "41467_2024_52262_MOESM16_ESM.xlsx" without header
# proteins_metadata.xlsx: header info for "41467_2024_52262_MOESM16_ESM.xlsx"
PROTEINS_DATA_RAW_FILENAME = "proteins_raw.xlsx"

# Output Filenames:
# Metabolites
METABOLITES_OUTPUT_BATCH_1_FILENAME = "metabolites_input_for_analyst_batch_1.csv"
METABOLITES_OUTPUT_BATCH_2_FILENAME = "metabolites_input_for_analyst_batch_2.csv"
# Proteins
PROTEINS_OUTPUT_FILENAME = "proteins_input_for_analyst.csv"
PROTEINS_METADATA_FILENAME = "proteins_metadata.xlsx"


"""
Create Output Folder if it doesn't exist
"""
if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)


"""
Import Raw Data
"""
metabolites_data_raw = pd.read_excel(pjoin(INPUT_FOLDER, METABOLITES_DATA_RAW_FILENAME))
metabolites_metadata = pd.read_excel(pjoin(INPUT_FOLDER, METABOLITES_METADATA_FILENAME))

proteins_data_raw = pd.read_excel(pjoin(INPUT_FOLDER, PROTEINS_DATA_RAW_FILENAME))
proteins_metadata = pd.read_excel(pjoin(INPUT_FOLDER, PROTEINS_METADATA_FILENAME))

"""
Specify QC Groups in Metadata
"""
# For values of "-" in column "Group", replace with "QC"
metabolites_metadata = metabolites_metadata.replace({"Group": {"-": "QC"}})
proteins_metadata = proteins_metadata.replace({"Group": {"-": "QC"}})

"""
Create Output DataFrames for Metabolite Batches
"""
# Determine sampledIDs in "Batch 1" or "Batch 2" based on values of column "Batch number" in metabolites_metadata.
batch_1_sampleIDs = []
batch_2_sampleIDs = []
for index, row in metabolites_metadata.iterrows():
    if row["Batch number"] == "Batch 1":
        batch_1_sampleIDs.append(row["SampleID"])
    elif row["Batch number"] == "Batch 2":
        batch_2_sampleIDs.append(row["SampleID"])
    else:
        raise ValueError(f"Unexpected Batch number: {row['Batch number']}")

# Create batch dataframes
metabolites_output_batch_1 = create_batch_dataframe(metabolites_data_raw, batch_1_sampleIDs)
metabolites_output_batch_2 = create_batch_dataframe(metabolites_data_raw, batch_2_sampleIDs)

# Add group labels as a second row in the dataframes
metabolites_output_batch_1 = add_group_labels_row(metabolites_output_batch_1, batch_1_sampleIDs, metabolites_metadata)
metabolites_output_batch_2 = add_group_labels_row(metabolites_output_batch_2, batch_2_sampleIDs, metabolites_metadata)

# Save to output folder
metabolites_output_batch_1.to_csv(pjoin(OUTPUT_FOLDER, METABOLITES_OUTPUT_BATCH_1_FILENAME), index=False)
metabolites_output_batch_2.to_csv(pjoin(OUTPUT_FOLDER, METABOLITES_OUTPUT_BATCH_2_FILENAME), index=False)


"""
Create Output Dataframe for Proteins
"""
# Starting from the 5th column, fetch the SampleIDs in the proteins_data_raw
proteins_sample_ids = proteins_data_raw.columns[4:].tolist()
# For proteins_output df, include "Protein Accession" and all the columns for SampleIDs
proteins_output = proteins_data_raw[["Protein Accession"] + proteins_sample_ids].copy()
# Rename "Protein Accession" to "SampleID" for ExpressAnalyst formatting
proteins_output.rename(columns={"Protein Accession": "SampleID"}, inplace=True)

# Add group labels as a second row in the dataframe
proteins_output = add_group_labels_row(proteins_output, proteins_sample_ids, proteins_metadata)

# Save to output folder
proteins_output.to_csv(pjoin(OUTPUT_FOLDER, PROTEINS_OUTPUT_FILENAME), index=False)
