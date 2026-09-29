#!/usr/bin/env python3
import os
"""Create the independent GSE192741 Visium spot-level co-expression audit."""

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
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
IN_DIR = ROOT / "results" / "spatial_validation" / "GSE192741"
OUT_DIR = ROOT / "figures" / "revision_p0"
SOURCE_DIR = ROOT / "source_data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["GSM5764424", "GSM5764425", "GSM5764427", "GSM5764426", "GSM5764428"]
SAMPLE_LABELS = {
    "GSM5764424": "GSM5764424\nsteatotic",
    "GSM5764425": "GSM5764425\nsteatotic",
    "GSM5764427": "GSM5764427\nsteatotic",
    "GSM5764426": "GSM5764426\nhealthy",
    "GSM5764428": "GSM5764428\nhealthy",
}
GENES = ["SLAMF8", "SPP1", "CD68", "CD163", "ITGAV", "ITGB5"]
PANELS = ["HSC_fibrogenic", "HSC_activation"]
PANEL_LABELS = {"HSC_fibrogenic": "HSC/fibrogenic", "HSC_activation": "HSC activation"}
GENE_LABELS = {"SLAMF8": "SLAMF8", "SPP1": "SPP1", "CD68": "CD68", "CD163": "CD163", "ITGAV": "ITGAV", "ITGB5": "ITGB5"}

plt.rcParams.update({
    "font.size": 8.2,
    "axes.titlesize": 10.2,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.4,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.75,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def key_rows(corr: pd.DataFrame, gene: str) -> pd.DataFrame:
    return corr[(corr["gene"] == gene) & (corr["panel"].isin(PANELS))].copy()


def draw_dot_panel(ax: plt.Axes, corr: pd.DataFrame, panel: str) -> None:
    data = corr[corr["panel"] == panel].copy()
    data["sample"] = pd.Categorical(data["sample"], categories=SAMPLE_ORDER, ordered=True)
    data["gene"] = pd.Categorical(data["gene"], categories=GENES, ordered=True)
    data = data.sort_values(["sample", "gene"])
    x = data["gene"].cat.codes.to_numpy()
    y = data["sample"].cat.codes.to_numpy()
    pseudocount = 1e-300
    p_values = np.clip(data["p_value"].to_numpy(float), pseudocount, None)
    sizes = 28 + 28 * np.clip(-np.log10(p_values), 0, 8)
    norm = Normalize(vmin=-0.45, vmax=0.45)
    ax.scatter(x, y, c=data["rho"], s=sizes, cmap="RdBu_r", norm=norm,
               edgecolors="white", linewidths=0.45, alpha=0.92)
    ax.set_xticks(range(len(GENES)), [GENE_LABELS[g] for g in GENES], rotation=35, rotation_mode="anchor", ha="right")
    ax.set_yticks(range(len(SAMPLE_ORDER)), [SAMPLE_LABELS[s] for s in SAMPLE_ORDER])
    ax.tick_params(axis="y", pad=12)
    ax.set_xlim(-0.6, len(GENES) - 0.4)
    ax.set_ylim(len(SAMPLE_ORDER) - 0.5, -0.5)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(PANEL_LABELS[panel], loc="left", pad=7)
    ax.grid(color="#E6EAED", linewidth=0.55)
    ax.set_axisbelow(True)
    for tick in ax.get_yticklabels():
        if "steatotic" in tick.get_text():
            tick.set_color("#B33A3A")
        else:
            tick.set_color("#4B6575")


def main() -> None:
    corr = pd.read_csv(IN_DIR / "spot_level_gene_panel_correlations.tsv", sep="\t")
    corr = corr[corr["sample"].isin(SAMPLE_ORDER) & corr["gene"].isin(GENES) & corr["panel"].isin(PANELS)].copy()
    corr["sample"] = pd.Categorical(corr["sample"], categories=SAMPLE_ORDER, ordered=True)
    corr["gene"] = pd.Categorical(corr["gene"], categories=GENES, ordered=True)
    corr["panel"] = pd.Categorical(corr["panel"], categories=PANELS, ordered=True)

    summary = pd.read_csv(IN_DIR / "sample_level_spatial_summary.tsv", sep="\t")
    summary = summary[summary["sample"].isin(SAMPLE_ORDER)].copy()
    summary["sample"] = pd.Categorical(summary["sample"], categories=SAMPLE_ORDER, ordered=True)
    summary = summary.sort_values("sample")

    fig = plt.figure(figsize=(7.2, 7.5), facecolor="white")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85], hspace=0.52, wspace=0.28)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    draw_dot_panel(ax_a, corr, "HSC_fibrogenic")
    draw_dot_panel(ax_b, corr, "HSC_activation")
    ax_b.set_yticklabels([])
    ax_b.tick_params(axis="y", length=0, pad=0)
    sm = ScalarMappable(norm=Normalize(vmin=-0.45, vmax=0.45), cmap="RdBu_r")
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[ax_a, ax_b], fraction=0.025, pad=0.02, shrink=0.86)
    cbar.set_label("Spearman ρ", rotation=90, rotation_mode="anchor", labelpad=8)
    cbar.ax.tick_params(labelsize=7.5, length=2)

    ax_c = fig.add_subplot(gs[1, :])
    ax_c.axis("off")
    columns = ["Sample", "Condition", "Spots", "SPP1–fibrogenic", "SPP1–activation", "SLAMF8 vs HSC"]
    cell_text = []
    for sample in SAMPLE_ORDER:
        row = summary[summary["sample"] == sample].iloc[0]
        s_f = key_rows(corr, "SPP1")
        s_f = s_f[(s_f["sample"] == sample) & (s_f["panel"] == "HSC_fibrogenic")]["rho"].iloc[0]
        s_a = key_rows(corr, "SPP1")
        s_a = s_a[(s_a["sample"] == sample) & (s_a["panel"] == "HSC_activation")]["rho"].iloc[0]
        cell_text.append([
            sample,
            str(row["condition"]),
            f"{int(row['n_spots']):,}",
            f"{s_f:.3f}",
            f"{s_a:.3f}",
            "no stable association",
        ])
    table = ax_c.table(cellText=cell_text, colLabels=columns, cellLoc="center", colLoc="center",
                       colWidths=[0.16, 0.14, 0.10, 0.19, 0.19, 0.22], bbox=[0.01, 0.16, 0.98, 0.72],
                       edges="horizontal")
    table.auto_set_font_size(False)
    table.set_fontsize(8.0)
    for (row_idx, col_idx), cell in table.get_celld().items():
        cell.set_edgecolor("#D8DEE3")
        cell.set_linewidth(0.55)
        if row_idx == 0:
            cell.set_facecolor("#EAF0F3")
            cell.set_text_props(weight="bold", color="#263238")
            cell.set_edgecolor("white")
            cell.set_linewidth(0)
        elif row_idx % 2 == 0:
            cell.set_facecolor("#F8FAFB")
        if col_idx == 1 and row_idx > 0:
            cell.get_text().set_color("#B33A3A" if "steatotic" in cell.get_text().get_text() else "#4B6575")
    ax_c.set_title("C  Sample-level audit scope and prespecified SPP1 correlations", loc="left", pad=6)
    ax_c.text(0.01, 0.04,
              "Dots show spot-level co-expression; dot size is proportional to −log10(P). Spots are mixed measurements and are not patient-level replicates.",
              transform=ax_c.transAxes, fontsize=7.6, color="#5E6B73")

    fig.suptitle("Figure S12  Independent GSE192741 Visium spot-level co-expression audit",
                 fontsize=14.0, fontweight="bold", x=0.02, ha="left", y=0.995)
    stem = OUT_DIR / "Supplementary_Figure12_GSE192741_spatial_coexpression_audit"
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    source = corr.merge(summary[["sample", "condition", "n_spots"]], on=["sample", "condition",], how="left")
    source.to_csv(SOURCE_DIR / "Source_Data_Supplementary_Figure12_GSE192741_spatial_coexpression.tsv.gz",
                  sep="\t", index=False, compression="gzip")
    print(f"Wrote {stem}.pdf")


if __name__ == "__main__":
    main()
