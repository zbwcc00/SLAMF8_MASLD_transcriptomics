#!/usr/bin/env python
import os
"""Plot a conservative, patient-aware communication audit."""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
AUDIT = ROOT / "results" / "communication" / "audit"
FIG = ROOT / "figures"
REVISION_FIG = FIG / "revision_p0"
SOURCE_DATA = ROOT / "source_data"
REVISION_FIG.mkdir(parents=True, exist_ok=True)
SOURCE_DATA.mkdir(parents=True, exist_ok=True)


def main():
    FIG.mkdir(exist_ok=True)
    sns.set_theme(style="ticks", context="paper", font_scale=1.1)
    cellchat = pd.read_csv(AUDIT / "cellchat_candidate_audit.tsv", sep="\t")
    expression = pd.read_csv(AUDIT / "patient_level_ligand_receptor_expression_tests.tsv", sep="\t")
    subtype = pd.read_csv(AUDIT / "patient_level_macrophage_subtype_expression.tsv", sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(15.8, 5.2), gridspec_kw={"width_ratios": [1.05, 1.25, 1]})
    cc = cellchat.pivot(index="ligand", columns="condition", values="n_pairs").fillna(0)
    cc = cc.reindex(["SPP1", "PDGFC", "PDGFB", "TGFB1", "CXCL16"])
    sns.heatmap(cc, annot=True, fmt=".0f", cmap="Blues", cbar=False, linewidths=1, linecolor="white", ax=axes[0])
    axes[0].set_title("CellChat candidate pairs\nGSE298719", loc="left")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")
    axes[0].tick_params(axis="x", rotation=20, pad=8)
    for tick in axes[0].get_xticklabels():
        tick.set_rotation_mode("anchor")

    expression["gene"] = expression.feature.str.extract(r"^[^_]+_(.*)$", expand=False)
    expression["compartment"] = expression.feature.str.extract(r"^(mac|hsc)_", expand=False)
    wanted = expression.gene.isin(["SLAMF8", "SPP1", "PDGFB", "PDGFC", "TGFB1", "PDGFRA", "PDGFRB", "CD44", "ITGB5"])
    plot = expression[wanted].copy()
    plot["label"] = plot.compartment.str.upper() + " " + plot.gene
    matrix = plot.set_index("label")[["masld_minus_healthy"]].sort_index()
    matrix.columns = ["MASLD − Healthy"]
    sns.heatmap(matrix, annot=True, fmt="+.2f", center=0, cmap="RdBu_r", linewidths=1, linecolor="white", ax=axes[1], cbar_kws={"label": "MASLD − Healthy median", "pad": 0.12})
    axes[1].set_title("Patient-level expression audit\n6 Healthy vs 7 MASLD", loc="left", pad=10)
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    axes[1].tick_params(axis="x", rotation=0, pad=12)
    for tick in axes[1].get_xticklabels():
        tick.set_rotation_mode("anchor")

    data = subtype[(subtype.condition == "MASLD") & (subtype.cell_group == "Mac_GPNMB_TREM2_SAMac")].dropna(subset=["SLAMF8", "SPP1"])
    rho, pvalue = spearmanr(data.SLAMF8, data.SPP1)
    axes[2].scatter(data.SLAMF8, data.SPP1, s=55, color="#b2182b", alpha=0.85)
    if len(data) >= 2:
        coeff = np.polyfit(data.SLAMF8, data.SPP1, 1)
        x = np.linspace(data.SLAMF8.min(), data.SLAMF8.max(), 50)
        axes[2].plot(x, coeff[0] * x + coeff[1], color="#333333", linewidth=1.2)
    axes[0].set_title("CellChat candidate pairs\nGSE298719", loc="left", pad=10)
    axes[2].set_title(f"GPNMB/TREM2 SAMac\nSLAMF8–SPP1: ρ={rho:.2f}, P={pvalue:.3f}, n={len(data)}", loc="left", pad=14)
    axes[2].set_xlabel("SLAMF8 mean log2(CPM+1)")
    axes[2].set_ylabel("SPP1 mean log2(CPM+1)")
    sns.despine(ax=axes[2])
    fig.suptitle("Macrophage–HSC communication audit (exploratory)", x=0.01, ha="left", y=1.04, fontsize=13)
    fig.subplots_adjust(left=0.035, right=0.985, top=0.76, bottom=0.19, wspace=0.38)
    stems = [FIG / "015_communication_audit", REVISION_FIG / "Supplementary_Figure7_GSE298719_communication_audit"]
    for stem in stems:
        fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
    source_tables = []
    for panel, table in [("cellchat_candidate_pairs", cellchat),
                         ("patient_level_expression", expression),
                         ("macrophage_subtype_expression", subtype)]:
        copy = table.copy()
        copy.insert(0, "panel", panel)
        source_tables.append(copy)
    pd.concat(source_tables, ignore_index=True, sort=False).to_csv(
        SOURCE_DATA / "Source_Data_Supplementary_Figure7_GSE298719_communication_audit.tsv.gz",
        sep="\t", index=False, compression="gzip")
    plt.close(fig)


if __name__ == "__main__":
    main()
