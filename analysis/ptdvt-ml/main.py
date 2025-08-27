#!/usr/bin/env python3
# ptdvt_covadj_bootstrap.py
# End-to-end: parse supplements -> build CSVs -> optional covariate adjustment -> train SVM ->
# external validation -> bootstrap AUC on discovery. Saves plots and CSVs.
#
# Requires: pandas numpy scipy scikit-learn matplotlib openpyxl

from pathlib import Path
import re, math, sys, traceback
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# CONFIG (EDIT THESE) 
DATA_DIR   = Path("/Users/selina/Documents/sup data")     # folder with "sup data 2.xlsx" and "sup 11.xlsx"
SUP2_XLSX  = DATA_DIR / "sup data 2.xlsx"                 # discovery
SUP11_XLSX = DATA_DIR / "sup data 11.xlsx"                     # validation
DISC_SHEET = "Supplementary Data 2"
VAL_SHEET  = "Supplementary Data 11"

TOP_K         = 9
RANDOM_STATE  = 42
LOG1P         = False        # set True for raw, non-negative intensities
BOOTSTRAPS    = 200          # e.g., 200–1000
TEST_SIZE     = 0.20
COV_CSV       = None         # Optional: Path("/Users/you/covariates.csv") with columns: sample_id, age, sex, trauma_time, ...
                             # If provided, must include both discovery & validation sample_ids to adjust both.
# 


# ------------------------ helpers ------------------------
def _find_row(raw: pd.DataFrame, names):
    col0 = raw[0].astype(str).str.strip().str.lower()
    canon = {n.strip().lower() for n in names}
    hits = col0[col0.isin(canon)]
    return hits.index[0] if not hits.empty else None

def _norm_group_value(s):
    if s is None or (isinstance(s, float) and math.isnan(s)):
        return None
    v = str(s).strip()
    if v == "pt-DVT":
        return 1
    if v == "Control":
        return 0
    return None

def _drop_internal_standard_features(feat_cols):
    pat = re.compile(r'(?i)^(is(\.|_|$))|(?:\bis\b)|internal\s*standard|internal\b')
    return [c for c in feat_cols if not pat.search(str(c))]

def _rank_by_pvalue(X: pd.DataFrame, y: pd.Series):
    mask1 = (y == 1); mask0 = (y == 0)
    pvals = []
    for col in X.columns:
        a = X.loc[mask1, col].astype(float)
        b = X.loc[mask0, col].astype(float)
        if a.dropna().nunique() < 2 or b.dropna().nunique() < 2:
            p = 1.0
        else:
            _, p = ttest_ind(a, b, equal_var=False, nan_policy="omit")
        if math.isnan(p):
            p = 1.0
        pvals.append((col, p))
    pvals.sort(key=lambda z: z[1])
    return [c for c,_ in pvals]

def _parse_sup_sheet_with_names(xlsx_path: Path, sheet_name: str) -> pd.DataFrame:
    """Return tidy per-sample matrix with columns [sample_id, label, features...]"""
    raw = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=None)

    # locate required rows
    group_row = _find_row(raw, ["Group"])
    id_row    = _find_row(raw, ["ID"])         # header band for metabolite table
    if group_row is None or id_row is None:
        raise RuntimeError(f"Could not find Group/ID rows in {xlsx_path.name}:{sheet_name}")

    hdr = raw.loc[id_row].astype(str).str.strip()
    metabolite_col = None
    if (hdr.str.lower() == "metabolite").any():
        metabolite_col = int(hdr[hdr.str.lower() == "metabolite"].index[0])
    id_col = 0
    if (hdr.str.lower() == "id").any():
        id_col = int(hdr[hdr.str.lower() == "id"].index[0])

    # sample columns: where Group is exactly "pt-DVT" or "Control"
    group_vals = raw.loc[group_row].astype(str).str.strip()
    is_sample_col = group_vals.isin(["pt-DVT", "Control"])
    sample_cols = [c for c in raw.columns if bool(is_sample_col.get(c, False))]
    if not sample_cols:
        raise RuntimeError(f"No sample columns detected from Group row in {xlsx_path.name}:{sheet_name}")

    # sample IDs (prefer explicit sample row; fallback to the header band)
    sample_row = _find_row(raw, ["Sample name", "Sample", "Sample id", "SampleID", "samplename"])
    if sample_row is not None:
        sample_ids = raw.loc[sample_row, sample_cols].astype(str).tolist()
    else:
        sample_ids = raw.loc[id_row, sample_cols].astype(str).tolist()

    # labels
    labels_raw = raw.loc[group_row, sample_cols].astype(str).tolist()
    labels = [ _norm_group_value(v) for v in labels_raw ]
    keep_mask = np.array([v in (0,1) for v in labels], dtype=bool)
    sample_cols = [c for c,k in zip(sample_cols, keep_mask) if k]
    sample_ids  = [s for s,k in zip(sample_ids,  keep_mask) if k]
    labels      = [v for v in labels if v in (0,1)]

    # data matrix sits below the ID header row
    data_start = id_row + 1
    ids = raw.iloc[data_start:, id_col].astype(str).tolist()
    if metabolite_col is not None:
        metabolite_names = raw.iloc[data_start:, metabolite_col].astype(str).tolist()
    else:
        metabolite_names = ids

    # remove internal standards rows
    is_is = (pd.Series(ids).str.match(r"(?i)^is(\b|\.|_|-)", na=False)
             | pd.Series(metabolite_names).str.contains(r"(?i)\bis\b|internal\s*standard|internal\b", regex=True, na=False))
    keep_rows = ~is_is
    # build feature names "ID__Metabolite"
    feat_names = []
    for mid, mname, keep in zip(ids, metabolite_names, keep_rows):
        if not keep: 
            continue
        name = mname if mname and mname.lower() != "nan" else mid
        feat_names.append(f"{mid}__{name}" if mid and mid.lower() != "nan" else name)

    matrix = raw.iloc[data_start:, sample_cols]
    matrix = matrix.loc[keep_rows.values].apply(pd.to_numeric, errors="coerce")

    X = pd.DataFrame(matrix.values, index=feat_names, columns=sample_ids).T.reset_index()
    X = X.rename(columns={"index":"sample_id"})
    X.insert(1, "label", labels)
    numeric_cols = X.select_dtypes(include=["number"]).columns
    features = [c for c in numeric_cols if c not in {"label"}]
    return pd.concat([X[["sample_id","label"]], X[features]], axis=1)

def covariate_adjust(df: pd.DataFrame, cov: pd.DataFrame) -> pd.DataFrame:
    """
    Regress out covariates from each feature using LinearRegression.
    df: DataFrame with columns ['sample_id','label', features...]
    cov: DataFrame with columns ['sample_id', cov1, cov2, ...] for the SAME sample_ids.
    Returns a new df with the same columns, where feature columns are replaced by residuals.
    """
    merged = df.merge(cov, on="sample_id", how="inner", suffixes=("",""))
    if merged.shape[0] < df.shape[0]:
        missing = set(df["sample_id"]) - set(merged["sample_id"])
        print(f"[warn] covariate_adjust: dropping {len(missing)} samples with missing covariates.")
    # Identify covariate columns (non-feature, not label/sample_id)
    non_feature = {"sample_id","label"}
    feat_cols = [c for c in df.columns if c not in non_feature]
    cov_cols  = [c for c in merged.columns if c not in set(df.columns) and c not in non_feature]
    # Coerce covariates to numeric (with simple impute for any missing)
    covX = merged[cov_cols].apply(pd.to_numeric, errors="coerce")
    covX = pd.DataFrame(SimpleImputer(strategy="median").fit_transform(covX), columns=cov_cols, index=merged.index)
    out = merged[["sample_id","label"]].copy()
    # Regress out for each feature
    lr = LinearRegression()
    for col in feat_cols:
        y = pd.to_numeric(merged[col], errors="coerce").values.reshape(-1,1)
        # if all nan or constant, set residuals to nan/zero
        if np.all(np.isnan(y)) or np.nanstd(y) == 0:
            res = np.full((len(merged),), np.nan)
        else:
            y_f = np.nan_to_num(y, nan=np.nanmedian(y))
            lr.fit(covX, y_f.ravel())
            pred = lr.predict(covX)
            res = (y_f.ravel() - pred)
        out[col] = res
    return out

def bootstrap_auc_on_discovery(X: pd.DataFrame, y: pd.Series, top_k: int, B: int, test_size: float, random_state: int):
    """
    Stratified bootstrap using repeated random splits on discovery:
    - within each split, fit imputer/selector/scaler on training only
    - compute AUC on held-out test
    """
    rng = np.random.RandomState(random_state)
    sss = StratifiedShuffleSplit(n_splits=B, test_size=test_size, random_state=random_state)
    aucs = []
    for b, (tr, te) in enumerate(sss.split(X, y), 1):
        X_tr, X_te = X.iloc[tr].copy(), X.iloc[te].copy()
        y_tr, y_te = y.iloc[tr].copy(), y.iloc[te].copy()

        # impute train -> apply to test
        imp = SimpleImputer(strategy="median").fit(X_tr)
        X_tr_imp = pd.DataFrame(imp.transform(X_tr), columns=X.columns, index=X_tr.index)
        X_te_imp = pd.DataFrame(imp.transform(X_te), columns=X.columns, index=X_te.index)

        # optional log1p
        if LOG1P and (X_tr_imp.values.min() >= 0) and (X_te_imp.values.min() >= 0):
            X_tr_imp = np.log1p(X_tr_imp); X_te_imp = np.log1p(X_te_imp)

        # rank on train, select top_k
        ranked = _rank_by_pvalue(X_tr_imp, y_tr)
        feats  = ranked[:min(top_k, len(ranked))]

        # scale train -> apply to test
        sc = StandardScaler().fit(X_tr_imp[feats])
        clf = SVC(kernel="poly", degree=3, probability=True, random_state=RANDOM_STATE)
        clf.fit(sc.transform(X_tr_imp[feats]), y_tr)

        y_score = clf.predict_proba(sc.transform(X_te_imp[feats]))[:,1]
        try:
            auc = roc_auc_score(y_te, y_score)
            aucs.append(auc)
        except Exception:
            # e.g., if y_te becomes single-class (shouldn't with stratify, but just in case)
            continue
        if b % 25 == 0:
            print(f"  bootstrap {b}/{B} | AUC mean={np.mean(aucs):.3f} (n={len(aucs)})")
    return np.array(aucs)


# ------------------------ main ------------------------
def main():
    print("=== ptdvt_covadj_bootstrap ===")
    print("DATA_DIR:", DATA_DIR.resolve())
    print("Inputs:"); print(" •", SUP2_XLSX); print(" •", SUP11_XLSX)

    # 1) Parse supplements -> tidy CSVs
    disc = _parse_sup_sheet_with_names(SUP2_XLSX, DISC_SHEET)
    val  = _parse_sup_sheet_with_names(SUP11_XLSX, VAL_SHEET)

    # Save tidy CSVs
    disc_out = DATA_DIR / "discovery_metabolomics_all.csv"
    val_out  = DATA_DIR / "validation_metabolomics_all.csv"
    disc.to_csv(disc_out, index=False); print("Saved:", disc_out, "| exists:", disc_out.exists())
    val.to_csv(val_out, index=False);   print("Saved:", val_out,  "| exists:", val_out.exists())
    print("Discovery label counts:", disc["label"].value_counts().to_dict())
    print("Validation label counts:", val["label"].value_counts().to_dict())

    # 2) Optional covariate adjustment (if COV_CSV provided)
    if COV_CSV:
        cov_path = Path(COV_CSV)
        if not cov_path.exists():
            print(f"[warn] COV_CSV not found: {cov_path}. Skipping covariate adjustment.")
        else:
            cov = pd.read_csv(cov_path)
            need_cols = {"sample_id"}
            if not need_cols.issubset(set(cov.columns)):
                print(f"[warn] covariates file must include: {need_cols}. Skipping covariate adjustment.")
            else:
                # Adjust discovery and validation separately using their own covariates
                print("Applying covariate adjustment (linear residuals) ...")
                disc = covariate_adjust(disc, cov)
                val  = covariate_adjust(val, cov)
                # Overwrite CSVs after adjustment
                disc.to_csv(DATA_DIR / "discovery_metabolomics_all_covadj.csv", index=False)
                val.to_csv(DATA_DIR / "validation_metabolomics_all_covadj.csv", index=False)
                print("Saved covariate-adjusted CSVs.")

    # 3) Align features (drop IS defensively)
    non_features = {"sample_id","label"}
    feat_cols = sorted(set(disc.columns) & set(val.columns) - non_features)
    feat_cols = _drop_internal_standard_features(feat_cols)
    print("Common features (no IS):", len(feat_cols))

    # Build design matrices
    X_disc = disc[feat_cols].apply(pd.to_numeric, errors="coerce")
    X_val  = val[feat_cols].apply(pd.to_numeric, errors="coerce")
    y_disc = disc["label"].astype(int)
    y_val  = val["label"].astype(int)

    # 4) Impute on discovery -> apply to validation
    imp = SimpleImputer(strategy="median").fit(X_disc)
    X_disc_imp = pd.DataFrame(imp.transform(X_disc), columns=feat_cols, index=X_disc.index)
    X_val_imp  = pd.DataFrame(imp.transform(X_val),  columns=feat_cols, index=X_val.index)

    # optional log1p
    if LOG1P and (X_disc_imp.values.min() >= 0) and (X_val_imp.values.min() >= 0):
        X_disc_imp = np.log1p(X_disc_imp); X_val_imp = np.log1p(X_val_imp)

    # 5) Rank & select top-K on discovery
    ranked = _rank_by_pvalue(X_disc_imp, y_disc)
    top_feats = ranked[:min(TOP_K, len(ranked))]
    print("\nTop features (discovery):")
    for i, f in enumerate(top_feats, 1):
        print(f"{i:2d}. {f}")

    # 6) Train SVM on discovery & evaluate once on validation
    sc = StandardScaler().fit(X_disc_imp[top_feats])
    clf = SVC(kernel="poly", degree=3, probability=True, random_state=RANDOM_STATE)
    clf.fit(sc.transform(X_disc_imp[top_feats]), y_disc)
    y_score = clf.predict_proba(sc.transform(X_val_imp[top_feats]))[:,1]

    fpr, tpr, _ = roc_curve(y_val, y_score)
    auc = roc_auc_score(y_val, y_score)
    ap  = average_precision_score(y_val, y_score)

    # Save plots and scores
    roc_path = DATA_DIR / "metab_all_val_roc.png"
    pr_path  = DATA_DIR / "metab_all_val_pr.png"
    plt.figure(); plt.plot(fpr, tpr, lw=2, label=f"AUC={auc:.3f}")
    plt.plot([0,1],[0,1],"--", lw=1); plt.xlabel("FPR"); plt.ylabel("TPR")
    plt.title("ROC (validation)"); plt.legend(); plt.tight_layout(); plt.savefig(roc_path, dpi=200); plt.close()
    precision, recall, _ = precision_recall_curve(y_val, y_score)
    plt.figure(); plt.plot(recall, precision, lw=2, label=f"AP={ap:.3f}")
    plt.xlabel("Recall"); plt.ylabel("Precision"); plt.title("PR (validation)"); plt.legend()
    plt.tight_layout(); plt.savefig(pr_path, dpi=200); plt.close()
    scores_csv = DATA_DIR / "metab_all_validation_scores.csv"
    pd.DataFrame({"y_true": y_val.values, "y_score": y_score}).to_csv(scores_csv, index=False)

    print("\n=== Validation Results ===")
    print(f"ROC AUC: {auc:.3f}")
    print(f"Average Precision (PR AUC): {ap:.3f}")
    print("Saved:", roc_path); print("Saved:", pr_path); print("Saved:", scores_csv)

    # 7) Bootstrap AUC on discovery
    print(f"\nRunning discovery bootstrap (B={BOOTSTRAPS}, test_size={TEST_SIZE}) ...")
    aucs = bootstrap_auc_on_discovery(X_disc, y_disc, top_k=TOP_K, B=BOOTSTRAPS, test_size=TEST_SIZE, random_state=RANDOM_STATE)

    if len(aucs) == 0:
        print("[warn] No valid bootstrap AUCs computed (class issues?).")
        return

    mean_auc = float(np.mean(aucs))
    std_auc  = float(np.std(aucs))
    ci_l, ci_u = float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))

    print(f"Bootstrap AUC: mean={mean_auc:.3f} ± {std_auc:.3f} | 95% CI [{ci_l:.3f}, {ci_u:.3f}] | n={len(aucs)}")

    # Save bootstrap results
    aucs_csv = DATA_DIR / "discovery_bootstrap_aucs.csv"
    pd.DataFrame({"auc": aucs}).to_csv(aucs_csv, index=False)

    # Histogram plot
    hist_png = DATA_DIR / "discovery_bootstrap_aucs_hist.png"
    plt.figure()
    plt.hist(aucs, bins=30)
    plt.xlabel("AUC"); plt.ylabel("Count")
    plt.title(f"Discovery bootstrap AUCs (n={len(aucs)})")
    plt.tight_layout(); plt.savefig(hist_png, dpi=200); plt.close()
    print("Saved:", aucs_csv); print("Saved:", hist_png)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("\n[ERROR] Script aborted.")
        print(e)
        traceback.print_exc()
        sys.exit(1)
