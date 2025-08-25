# Multi-omics Analysis for Diagnosis of Post-trauma Blood Clots

Lazarina Butkovich, Selina Jessica, Vansika K, Nadia Tasevski

Open Bootcamp Collective, Project 2 - 2025

Study: Zhang et al. "Integrated landscape ofplasma metabolism and proteome ofpatients with post-traumatic deep vein thrombosis" Nature Communications (2025). https://doi.org/10.1038/s41467-024-52262-0

## Objectives
- For this portfolio bioinformatics project, we replicate analyses of Zhang et al.
- Metabolomics, single-omics
- Proteomics, single-omics (in progress)
- Clinical Parameters (in progress)
- Multi-omics (in progress)
- Predictive Model for pt-DVT (in progress)

## Background
- Study Design
    - Zhang et al. collected plasma metabolomics and proteomics data for **680 individuals**, with and without **deep vein thrombosis (DVT)** (aka blood clot) after trauma, referred to as **post-trauma deep vein thrombosis (pt-DVT)**.
    - DVT is a major health problem that can lead to complication, including post-thrombotic syndrome, recurrent DVT, and life-threatening pulmonary embolism. Risk of DVT is elevated post-trauma (aka physical injury).
    - Failure to diagnose pt-DVT is a leading cause of death post-trauma. Through their multi-omics analysis, Zhang et al. proposed (i) protein and metabolite biomarkers for early diagnosis of pt-DVT and (ii) potential therapeutic strategies for pt-DVT.
    - In order to develop a predictive model for pt-DVT, Zhang et al. recruited patients for a discovery group (N=580, with n_pt_DVT=252 and n_without_DVT=328 ) and a separate validation group (N=100, with n_pt_DVT=50, n_without_DVT=50).
- Metabolomics Dataset
    - After MS/MS identification, 326 metabolites were determined for metabolomic analysis in all 580 samples of the discovery cohort.
    - To process the data:
        - Process analytical batches separately (batch 1 and batch 2). By scaling by batch, the sensitivity differences between batches is corrected, and the relative abundances matrix of all peaks becomes comparable.
        - Filter to remove noise
            - Remove peaks whose RSD (relative standard deviation) is 40%< in QC samples
            - Remove peaks with more than 50% null values in a single group or all groups
        - Impute data:
            - For remaining features with missing intensity values, replace with half the minimum value for that metabolite across all samples.
        - Normalize by dividing intensities by the intensity of the internal standard (IS) (see data column "ID", value "IS" in Supplementary Data 2; we normalized by the positive ionization "IS")
        - Log-transform
        - Scale by median centered (we scaled via auto scaling - mean-center then divide by the standard deviation of each variable)
    - Use Supplementary Data 3 to directly compare metabolomics data processing of Zhang et al. to the results in this portfolio project.

## Scripts Overview for Metabolomics Analysis - Lazarina Butkovich
1. format_inputs_for_analyst.py
    - formats both metabolomic and proteomic data from edited version of Zhang et al Supplementary Datasets (Supplementary Datasets 2 and 13, respectively)
    - [Input files made available](https://drive.google.com/drive/folders/16BtKjSYnBhux7OyCesLEZcAfHZCnFDlo?usp=drive_link)
2. metaboanalystR_data_processing.r
    - Performs single-omic metabolomic analysis
    - Run batches separately then combine post-normalization
    - Run the combined, normalized data through [MetaboAnalyst](https://www.metaboanalyst.ca/MetaboAnalyst/ModuleView.xhtml) web platform, ["Statistical Analysis [one factor]"](https://www.metaboanalyst.ca/MetaboAnalyst/upload/StatUploadView.xhtml) (skip normalization steps)