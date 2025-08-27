# pt-DVT ML Replication (Discovery → External Validation)

This folder contains a reproducible pipeline that replicates the SVM-based classifier from the Nature Communications study on post-traumatic deep vein thrombosis (pt-DVT) using the paper’s discovery and external validation cohorts.

What it does

- Parses the paper’s Supplementary Excel sheets into tidy per-sample matrices with labels from the Group row.
- Removes internal standards so they aren’t used as features.
- (Optional) Covariate adjustment: regress each metabolite on user-provided covariates (e.g., age/sex/trauma_time) and use residuals.
- Ranks metabolites on discovery with Welch’s t-test; selects top-K (default 9) without leakage.
- Trains an SVM (polynomial kernel, degree 3) on discovery; evaluates once on validation (ROC/PR).
- Runs a bootstrap on discovery to report internal AUC mean, SD, and 95% CI.

# Files in this folder

* scripts/ptdvt_covadj_bootstrap.py — end-to-end pipeline (edit config at the top).

* examples/covariates_template.csv — optional example format for covariates.

* data/ — place Excel files here locally (not committed).

* outputs/ — generated plots/CSVs (gitignored).

* Please do not commit the paper’s Excel files or generated outputs.

# Inputs you provide

- Place the following files in data/ (download them from the paper’s Supplementary Information):

- Discovery (raw intensities, 580 samples): sup data 2.xlsx, sheet “Supplementary Data 2”

- Validation (raw intensities, 100 samples): sup 11.xlsx, sheet “Supplementary Data 11”

- Optional (to mimic the paper’s normalized analyses):

- Discovery (normalized): sup data 3.xlsx, sheet “Supplementary Data 3”

- Validation (normalized, 28 features): sup 12.xlsx, sheet “Supplementary Data 12”

- Optional covariates (your own file), e.g. data/covariates.csv:
```
sample_id,age,sex,trauma_time
DISC_001,54,1,3.2
DISC_002,41,0,1.7
...
```
sex can be 0/1. The script coerces covariates to numeric and median-imputes missing values internally.


# How to run (local Python)

##  Install dependencies (add to your project env):

pandas
numpy
scipy
scikit-learn
matplotlib
openpyxl

- Example:
```
pip install pandas numpy scipy scikit-learn matplotlib openpyxl
```
## Open scripts/main.py and edit the config at the top:
```
DATA_DIR   = Path("<absolute path to this subfolder>/data")
SUP2_XLSX  = DATA_DIR / "sup data 2.xlsx"
SUP11_XLSX = DATA_DIR / "sup 11.xlsx"
DISC_SHEET = "Supplementary Data 2"
VAL_SHEET  = "Supplementary Data 11"

TOP_K        = 9          # try 9, 20, or 28
RANDOM_STATE = 42
LOG1P        = False      # True for raw non-negative intensities if you want variance stabilization
BOOTSTRAPS   = 200
TEST_SIZE    = 0.20

COV_CSV      = None       # or DATA_DIR / "covariates.csv" to enable covariate adjustment
```

## Run the script (from your Python environment):
```
python scripts/main.py
```

# Outputs
- discovery_metabolomics_all.csv, validation_metabolomics_all.csv
- (if covariates used) *_covadj.csv
- metab_all_val_roc.png (ROC curve), metab_all_val_pr.png (PR curve)
- metab_all_validation_scores.csv (per-sample probabilities)
- discovery_bootstrap_aucs.csv, discovery_bootstrap_aucs_hist.png

# Methods
- Parsing: Keep only sample columns whose Group is exactly pt-DVT or Control; use Sample name as sample_id; drop internal standards.
- Preprocessing (fit on discovery): median imputation → optional log1p → Z-scaling.
- Feature selection (discovery only): Welch’s t-test per metabolite; keep top-K.
- Model: SVM (polynomial kernel, degree 3, probability=True) trained on discovery.
- Evaluation: single pass on validation (ROC AUC, Average Precision); bootstrap on discovery with per-split fit of imputer/ranker/scaler on the training fold only.

- Example results 
* Validation: ROC AUC 0.543, PR (AP) 0.588
* Discovery bootstrap (B=200): mean AUC 0.791 ± 0.039, 95% CI [0.719, 0.864]

Interpretation: strong separation inside discovery but modest generalization on validation → likely dataset shift / preprocessing differences.


# Tips / variants
- Try the normalized sheets (sup data 3.xlsx / sup 12.xlsx) with LOG1P=False.
- Adjust TOP_K (e.g., 20 or 28).
- Swap model for robustness: SVC(kernel="linear", probability=True) or LogisticRegression.
- Provide covariates (age, sex, trauma_time, batch) via COV_CSV to residualize confounders.

# Attribution
- Data: obtain the supplementary Excel files from the published article; do not redistribute here.

- Code: licensed under the parent repository’s license.

- https://www.nature.com/articles/s41467-024-52262-0#Sec26
