#!/usr/bin/env python3
"""P2.5: conservative bulk cell-state scores and locked external model pilot.

The immune quantities are relative transcriptomic signature scores, NOT cell
fractions.  The two cohort split and predictor sets are fixed before viewing
external test performance: GSE162694 trains; GSE130970 tests.
"""
from __future__ import annotations
import os

from pathlib import Path
import gzip
import re
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss, cohen_kappa_score, mean_absolute_error
from statsmodels.miscmodels.ordinal_model import OrderedModel
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
FULL = ROOT / "data_processed" / "bulk_geo" / "full"
META = ROOT / "results" / "bulk_validation"
OUT = ROOT / "results" / "revision_p2_5" / "bulk_immune_and_models"
FIG = ROOT / "figures" / "revision_p2_5"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
SEED = 20260923

# Fixed prior marker panels. Scores are mean gene-wise z-scores within bulk dataset.
PANELS = {
    "SAMac_like": ["TREM2", "GPNMB", "CD9", "SPP1", "LGALS3", "LPL", "APOC1", "FABP5", "CTSB", "CTSD"],
    "SAMac_like_no_SPP1": ["TREM2", "GPNMB", "CD9", "LGALS3", "LPL", "APOC1", "FABP5", "CTSB", "CTSD"],
    "Macrophage_myeloid": ["LST1", "TYROBP", "AIF1", "FCER1G", "CSF1R", "C1QA", "C1QB", "C1QC"],
    "Activated_HSC_fibrogenic": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "ACTA2", "TAGLN", "PDGFRB"],
    "T_cell": ["CD3D", "CD3E", "TRBC1", "TRBC2", "LCK", "CD247"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD37", "CD74", "HLA-DRA"],
    "NK_cell": ["NKG7", "KLRD1", "TRAC", "GNLY", "PRF1", "GZMB"],
}
MODEL_SPECS = {
    "SLAMF8_only": ["SLAMF8"],
    "Locked_4_feature_transcriptomic": ["SLAMF8", "SAMac_like", "SPP1", "Activated_HSC_fibrogenic"],
}

plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":9, "axes.spines.top":False,
                     "axes.spines.right":False, "pdf.fonttype":42, "ps.fonttype":42, "savefig.dpi":600})


def gene_maps():
    import anndata as ad
    obj = ad.read_h5ad(ROOT / "data_processed" / "GSE212837_myeloid_reclustered.h5ad", backed="r")
    g = obj.var[["gene_id", "gene"]].copy()
    g["gene_id"] = g.gene_id.astype(str).str.replace(r"\..*", "", regex=True)
    g["gene"] = g.gene.astype(str).str.upper()
    gene_map = g.drop_duplicates("gene_id").set_index("gene_id")["gene"]
    info = FULL / "Homo_sapiens.gene_info.gz"
    n = pd.read_csv(info, sep="\t", compression="gzip", dtype=str, low_memory=False)
    n.columns = n.columns.str.removeprefix("#")
    entrez = n.drop_duplicates("GeneID").set_index("GeneID")["Symbol"].str.upper()
    return gene_map, entrez


def metadata(accession: str) -> pd.DataFrame:
    x = pd.read_csv(META / f"{accession}_sample_metadata.tsv", sep="\t")
    cols = [c for c in x.columns if "characteristics" in c]
    rows = []
    for _, v in x.iterrows():
        row = {"gsm": v["!Sample_geo_accession"], "title": v["!Sample_title"]}
        for c in cols:
            val = str(v[c])
            if ":" in val:
                k, value = val.split(":", 1)
                row[k.strip().lower()] = value.strip()
        rows.append(row)
    out = pd.DataFrame(rows)
    if accession == "GSE162694":
        out["sample"] = out.title.str.extract(r"(548nash\d+)", flags=re.I, expand=False)
        out["fibrosis"] = pd.to_numeric(out["fibrosis stage"], errors="coerce")
        out["normal_histology"] = out["fibrosis stage"].eq("normal liver histology")
        out.loc[out.normal_histology, "fibrosis"] = 0
    else:
        out["sample"] = out.title
        out["fibrosis"] = pd.to_numeric(out["fibrosis stage"], errors="coerce")
        out["normal_histology"] = False
    out["age_numeric"] = pd.to_numeric(out.get("age", out.get("age at biopsy")), errors="coerce")
    out["sex_female"] = out["sex"].astype(str).str.strip().str.lower().isin(["female", "f"]).astype(int)
    out["nas_numeric"] = pd.to_numeric(out.get("nas score", out.get("nafld activity score")), errors="coerce")
    return out


def expression(accession: str, gene_map: pd.Series, entrez: pd.Series) -> pd.DataFrame:
    if accession == "GSE162694":
        x = pd.read_csv(FULL / "GSE162694_raw_counts.csv.gz", index_col=0)
        x.index = gene_map.reindex(x.index.astype(str).str.split(".").str[0]).to_numpy()
        x = x.loc[x.index.notna()].groupby(level=0).sum()
        return np.log2(x.div(x.sum(axis=0), axis=1).mul(1e6) + 1)
    x = pd.read_csv(FULL / "GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv.gz", index_col=0)
    x.index = entrez.reindex(x.index.astype(str)).to_numpy()
    x = x.loc[x.index.notna()].groupby(level=0).sum()
    return np.log2(x + 1)


def score(x: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    found = [g for g in genes if g in x.index]
    if len(found) < 3:
        return pd.Series(np.nan, index=x.columns), found
    z = x.loc[found].sub(x.loc[found].mean(axis=1), axis=0).div(x.loc[found].std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), found


def build_scores(accession: str, gene_map: pd.Series, entrez: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = metadata(accession).set_index("sample")
    x = expression(accession, gene_map, entrez)
    common = meta.index.intersection(x.columns)
    meta, x = meta.loc[common].copy(), x.loc[:, common]
    out = meta[["gsm", "title", "fibrosis", "normal_histology", "age_numeric", "sex_female", "nas_numeric"]].copy()
    coverage = []
    for name, genes in PANELS.items():
        s, found = score(x, genes)
        out[name] = s.reindex(out.index)
        coverage.append({"dataset":accession, "signature":name, "genes_requested":";".join(genes),
                         "genes_found":";".join(found), "n_found":len(found)})
    for gene in ["SLAMF8", "SPP1"]:
        out[gene] = x.loc[gene] if gene in x.index else np.nan
        coverage.append({"dataset":accession, "signature":gene, "genes_requested":gene,
                         "genes_found":gene if gene in x.index else "", "n_found":int(gene in x.index)})
    return out.reset_index(names="sample"), pd.DataFrame(coverage)


def bootstrap_metrics(y: np.ndarray, p: np.ndarray, rng: np.random.Generator, nboot: int = 3000):
    values = {"auc": [], "pr_auc": [], "brier": []}
    n = len(y)
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        if np.unique(y[idx]).size < 2:
            continue
        values["auc"].append(roc_auc_score(y[idx], p[idx]))
        values["pr_auc"].append(average_precision_score(y[idx], p[idx]))
        values["brier"].append(brier_score_loss(y[idx], p[idx]))
    return {f"{k}_{q}": float(np.quantile(v, z)) for k, v in values.items() for q, z in [("ci_low", .025), ("ci_high", .975)] if v}


def calibration(y: np.ndarray, p: np.ndarray):
    eps = 1e-6
    logit = np.log(np.clip(p, eps, 1-eps) / np.clip(1-p, eps, 1-eps))
    try:
        fit = sm.GLM(y, sm.add_constant(logit), family=sm.families.Binomial()).fit()
        return float(fit.params[0]), float(fit.params[1])
    except Exception:
        return np.nan, np.nan


def model_pilot(train: pd.DataFrame, test: pd.DataFrame):
    records, predictions = [], []
    rng = np.random.default_rng(SEED)
    for name, fields in MODEL_SPECS.items():
        tr = train.dropna(subset=fields + ["fibrosis"]).copy()
        te = test.dropna(subset=fields + ["fibrosis"]).copy()
        ytr, yte = (tr.fibrosis >= 2).astype(int).to_numpy(), (te.fibrosis >= 2).astype(int).to_numpy()
        pipe = Pipeline([("scale", StandardScaler()), ("logistic", LogisticRegression(C=1.0, penalty="l2", solver="lbfgs", max_iter=500, random_state=SEED))])
        pipe.fit(tr[fields], ytr)
        for split, frame, y in [("training_apparent", tr, ytr), ("external_test", te, yte)]:
            p = pipe.predict_proba(frame[fields])[:, 1]
            rec = {"model":name, "split":split, "endpoint":"F_ge_2", "n":len(frame), "events":int(y.sum()),
                   "auc":roc_auc_score(y,p), "pr_auc":average_precision_score(y,p), "brier":brier_score_loss(y,p)}
            rec["calibration_intercept"], rec["calibration_slope"] = calibration(y,p)
            if split == "external_test": rec.update(bootstrap_metrics(y,p,rng))
            records.append(rec)
            predictions.extend(pd.DataFrame({"model":name, "split":split, "sample":frame["sample"], "fibrosis":frame.fibrosis, "F_ge_2":y, "probability":p}).to_dict("records"))

        # Proportional-odds pilot: same locked predictors, explicitly exploratory.
        # It reports discrimination agreement, not a validated clinical staging tool.
        tr_ord = tr[tr.fibrosis.between(0,4)].copy(); te_ord = te[te.fibrosis.between(0,4)].copy()
        scaler = StandardScaler().fit(tr_ord[fields])
        try:
            fit = OrderedModel(tr_ord.fibrosis.astype(int), scaler.transform(tr_ord[fields]), distr="logit").fit(method="bfgs", disp=False, maxiter=300)
            probs = np.asarray(fit.model.predict(fit.params, exog=scaler.transform(te_ord[fields])))
            classes = np.array(sorted(tr_ord.fibrosis.astype(int).unique()))
            # OrderedModel returns columns corresponding to observed endogenous categories.
            pred = classes[np.argmax(probs, axis=1)]
            expected = probs @ classes
            records.append({"model":name, "split":"external_test", "endpoint":"ordinal_F0_F4", "n":len(te_ord),
                            "events":np.nan, "auc":np.nan, "pr_auc":np.nan, "brier":np.nan,
                            "calibration_intercept":np.nan, "calibration_slope":np.nan,
                            "ordinal_mae":mean_absolute_error(te_ord.fibrosis, expected),
                            "ordinal_quadratic_weighted_kappa":cohen_kappa_score(te_ord.fibrosis.astype(int), pred, weights="quadratic")})
        except Exception as e:
            records.append({"model":name, "split":"external_test", "endpoint":"ordinal_F0_F4", "n":len(te_ord), "fit_status":f"failed: {type(e).__name__}"})
    return pd.DataFrame(records), pd.DataFrame(predictions)


def adjusted_audits(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fixed, small-number sensitivity analyses; no variable selection."""
    associations, adjusted = [], []
    for ds, raw in scores.groupby("dataset"):
        d = raw[~raw.normal_histology].copy()
        for left, right in [("SAMac_like", "SLAMF8"), ("SAMac_like", "Activated_HSC_fibrogenic"),
                            ("SLAMF8", "Activated_HSC_fibrogenic")]:
            q = d[[left, right]].dropna()
            rho, p = spearmanr(q[left], q[right])
            associations.append({"dataset":ds, "left_feature":left, "right_feature":right,
                                 "n":len(q), "rho":rho, "p_value":p})
        for feature in ["SAMac_like", "SAMac_like_no_SPP1", "Activated_HSC_fibrogenic"]:
            for covariates, label in [(["age_numeric", "sex_female"], "age_sex"),
                                      (["age_numeric", "sex_female", "nas_numeric"], "age_sex_NAS")]:
                q = d[["fibrosis", feature] + covariates].dropna()
                if len(q) < 25 or q.fibrosis.nunique() < 2 or q[feature].std() == 0:
                    continue
                design = q[covariates].copy()
                design["signature_per_sd"] = (q[feature] - q[feature].mean()) / q[feature].std()
                fit = sm.OLS(q.fibrosis, sm.add_constant(design)).fit(cov_type="HC3")
                adjusted.append({"dataset":ds, "scope":"NAFLD_only", "feature":feature,
                                 "adjustment":label, "n":len(q),
                                 "beta_fibrosis_stage_per_sd":fit.params["signature_per_sd"],
                                 "p_value":fit.pvalues["signature_per_sd"]})
    assoc = pd.DataFrame(associations)
    adj = pd.DataFrame(adjusted)
    if not assoc.empty:
        from statsmodels.stats.multitest import multipletests
        assoc["fdr_within_dataset"] = assoc.groupby("dataset").p_value.transform(lambda x:multipletests(x, method="fdr_bh")[1])
    if not adj.empty:
        from statsmodels.stats.multitest import multipletests
        adj["fdr_within_adjustment"] = adj.groupby(["dataset", "adjustment"]).p_value.transform(lambda x:multipletests(x, method="fdr_bh")[1])
    return assoc, adj


def plot(scores: pd.DataFrame, model: pd.DataFrame):
    use = scores[scores.dataset.isin(["GSE162694", "GSE130970"])].copy()
    features = list(PANELS)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), gridspec_kw={"width_ratios":[1.2,1]}, facecolor="white")
    ax = axes[0]
    d = use[use.scope.eq("NAFLD_only")].copy()
    sns.barplot(data=d, y="feature", x="rho", hue="dataset", order=features, palette={"GSE162694":"#0072B2", "GSE130970":"#D55E00"}, ax=ax, errorbar=None)
    ax.axvline(0,color="#37474F",lw=.8); ax.set_xlabel("Spearman rho with fibrosis stage"); ax.set_ylabel("")
    ax.set_title("A  Relative bulk cell-state signatures", loc="left", fontweight="bold")
    ax.legend(frameon=False, title="")
    ax.text(0,-.18,"Transcriptomic signature scores are relative measures, not estimated cell fractions.", transform=ax.transAxes, fontsize=7.4, color="#65727A")
    ax = axes[1]
    e = model[(model.split=="external_test") & (model.endpoint=="F_ge_2")].copy()
    metrics = e.melt(id_vars="model", value_vars=["auc","pr_auc","brier"], var_name="metric", value_name="value")
    sns.barplot(data=metrics, x="metric", y="value", hue="model", palette=["#7A8C99", "#B2182B"], ax=ax)
    ax.set_xlabel(""); ax.set_ylabel("External-test value")
    ax.set_title("B  Locked external transcriptomic pilot", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, title="")
    ax.text(0,-.18,"Exploratory tissue-transcriptomic model; not a clinical prediction tool.", transform=ax.transAxes, fontsize=7.4, color="#65727A")
    fig.suptitle("P2.5  Bulk cell-state signatures and external fibrosis-model pilot", x=.02, ha="left", fontweight="bold", fontsize=14)
    fig.tight_layout(rect=[0, .04, 1, .94])
    fig.savefig(FIG / "P25_bulk_immune_signature_and_locked_model_pilot.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIG / "P25_bulk_immune_signature_and_locked_model_pilot.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    gene_map, entrez = gene_maps()
    frames, cov = [], []
    for ds in ["GSE162694", "GSE130970"]:
        s,c = build_scores(ds,gene_map,entrez); s["dataset"] = ds; frames.append(s); cov.append(c)
    all_scores = pd.concat(frames, ignore_index=True)
    all_scores.to_csv(OUT / "bulk_cell_state_signature_scores.tsv", sep="\t", index=False)
    pd.concat(cov, ignore_index=True).to_csv(OUT / "bulk_cell_state_signature_coverage.tsv", sep="\t", index=False)
    rows = []
    for ds, d0 in all_scores.groupby("dataset"):
        for scope, d in [("all_biopsies", d0), ("NAFLD_only", d0[~d0.normal_histology])]:
            for feature in PANELS:
                q = d[["fibrosis", feature]].dropna(); rho,p = spearmanr(q.fibrosis,q[feature])
                rows.append({"dataset":ds,"scope":scope,"feature":feature,"n":len(q),"rho":rho,"p_value":p})
    assoc = pd.DataFrame(rows)
    from statsmodels.stats.multitest import multipletests
    assoc["fdr_within_dataset_scope"] = assoc.groupby(["dataset","scope"]).p_value.transform(lambda x:multipletests(x,method="fdr_bh")[1])
    assoc.to_csv(OUT / "bulk_cell_state_signature_fibrosis_associations.tsv", sep="\t", index=False)
    train = all_scores[(all_scores.dataset=="GSE162694") & ~all_scores.normal_histology].copy()
    test = all_scores[(all_scores.dataset=="GSE130970") & ~all_scores.normal_histology].copy()
    models,pred = model_pilot(train,test)
    models.to_csv(OUT / "locked_external_transcriptomic_model_pilot.tsv",sep="\t",index=False)
    pred.to_csv(OUT / "locked_external_transcriptomic_model_predictions.tsv",sep="\t",index=False)
    state_assoc, state_adjusted = adjusted_audits(all_scores)
    state_assoc.to_csv(OUT / "bulk_cell_state_signature_cross_feature_associations.tsv", sep="\t", index=False)
    state_adjusted.to_csv(OUT / "bulk_cell_state_signature_covariate_sensitivity.tsv", sep="\t", index=False)
    plot(assoc,models)
    summary = {"train":"GSE162694 NAFLD-only", "external_test":"GSE130970 NAFLD-only", "endpoint":"F>=2", "models":MODEL_SPECS,
               "interpretation":"relative cell-state signatures and locked exploratory tissue-transcriptomic pilot; not cell fractions or a clinical prediction tool"}
    (OUT / "README.md").write_text("# P2.5 bulk signatures and locked models\n\n"+json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(assoc[assoc.scope.eq("NAFLD_only")].to_string(index=False))
    print(models.to_string(index=False))
    print(state_assoc.to_string(index=False))
    print(state_adjusted.to_string(index=False))

if __name__ == "__main__": main()
