"""
combine_normalized_batches.py
Lazarina Butkovich
created 8/9/25

From normalized peak intensity dataframes for metabolomics data from Zhang et al. (2025), this script combines Batch 1 and Batch 2 data, based on overlap of shared rows (feature ID values in the first column) 

Zhang et al. (2025): https://doi.org/10.1038/s41467-024-52262-0
"""

import pandas as pd
import os
from os.path import join as pjoin


"""
Values
"""
# Folder names:
OUTPUT_FOLDER = "output"
INPUT_FOLDER = "input"
# INPUT_FOLDER = pjoin(OUTPUT_FOLDER, "MetaboAnalyst_manual_results")

BATCH_1_NORMALIZED_FILENAME = "data_normalized_batch_1.csv"
BATCH_2_NORMALIZED_FILENAME = "data_normalized_batch_2.csv"


"""
Import data as pd dataframes
"""
batch_1_df = pd.read_csv(pjoin(INPUT_FOLDER, BATCH_1_NORMALIZED_FILENAME), header=[0,1])
batch_2_df = pd.read_csv(pjoin(INPUT_FOLDER, BATCH_2_NORMALIZED_FILENAME), header=[0,1])

# Name the first column, top row as "SampleID" for both dfs
batch_1_df.columns = pd.MultiIndex.from_tuples([("SampleID", "Label")] + list(batch_1_df.columns[1:]))
batch_2_df.columns = pd.MultiIndex.from_tuples([("SampleID", "Label")] + list(batch_2_df.columns[1:]))


"""
Create combined_df
"""
# Initialize combined_df as having the columns of batch_1_df and batch_2_df
combined_df = pd.concat([batch_1_df.set_index(batch_1_df.columns[0]), 
                        batch_2_df.set_index(batch_2_df.columns[0])], 
                       axis=1, join='outer')

# Reset index to bring SampleID back as a column
combined_df = combined_df.reset_index()

# For rows with any missing values (ie: features that were not kept in BOTH batch 1 and batch 2), remove those rows
combined_df = combined_df.dropna()

# Save combined_df to a CSV file
combined_df.to_csv(pjoin(OUTPUT_FOLDER, "data_normalized_with_QCs.csv"), index=False)

# For downstream analysis, we want to remove QC samples (Label value is "QC")
# Remove columns with value in "Label" header row of "QC"
combined_df = combined_df.loc[:, combined_df.columns.get_level_values(1) != "QC"]

# Save combined_df to a CSV file
combined_df.to_csv(pjoin(OUTPUT_FOLDER, "data_normalized_final.csv"), index=False)

"""
Print relevant statistics for "data_normalized_final.csv"
"""
# Print how many features were kept after combining
print(f"Number of features kept after combining: {combined_df.shape[0]}")
# Print number of samples
print(f"Number of samples kept after combining: {combined_df.shape[1] - 1}")

count_pt_dvt = 0
count_control = 0
# Fetch header row for Label (2nd row in header)
label_row = combined_df.columns.get_level_values(1)
for val in label_row:
    if val == "pt-DVT":
        count_pt_dvt += 1
    elif val == "Control":
        count_control += 1

print(f"Number of samples in pT-DVT group: {count_pt_dvt}")
print(f"Number of samples in Control group: {count_control}")