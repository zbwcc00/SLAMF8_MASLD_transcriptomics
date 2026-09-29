#!/usr/bin/env python3
"""Redraw Supplementary Figure S1 as an original BioRender-inspired workflow."""
from __future__ import annotations
import os

import gzip
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, PathPatch, Polygon
from matplotlib.path import Path as MplPath
from PIL import Image


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
PACKAGE = ROOT / "submission_package_Scientific_Reports_20260929"
OUT = PACKAGE / "02_Figures" / "Supplementary"
WORK = ROOT / "figures" / "scientific_reports_revision"
SOURCE = ROOT / "source_data" / "Source_Data_Supplementary_Figure1_workflow.tsv.gz"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.8,
        "savefig.facecolor": "white",
    }
)

INK = "#263746"
MUTED = "#637381"
LINE = "#AEBBC5"
WHITE = "#FFFFFF"
BLUE = "#3C78A8"
BLUE_BG = "#EAF3F9"
TEAL = "#2D8C82"
TEAL_BG = "#E8F5F2"
PURPLE = "#8066A8"
PURPLE_BG = "#F1ECF7"
ORANGE = "#D17A3A"
ORANGE_BG = "#FBF0E6"
GREEN = "#3A8A62"
GREEN_BG = "#E9F4ED"
CORAL = "#C85F55"
CORAL_BG = "#F9ECEA"


def rounded(ax, xy, width, height, face, edge=LINE, radius=0.12, lw=0.9, z=1):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.02,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        zorder=z,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, start, end, color=MUTED, lw=1.15, style="solid", curve=0.0, z=3):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=9.5,
        linewidth=lw,
        linestyle=style,
        color=color,
        connectionstyle=f"arc3,rad={curve}",
        shrinkA=2,
        shrinkB=2,
        zorder=z,
    )
    ax.add_patch(patch)


def phase_label(ax, x, y, number, title, color):
    ax.add_patch(Circle((x, y), 0.115, facecolor=color, edgecolor="none", zorder=5))
    ax.text(x, y - 0.003, str(number), ha="center", va="center", fontsize=6.5, fontweight="bold", color=WHITE, zorder=6)
    ax.text(x + 0.17, y, title.upper(), ha="left", va="center", fontsize=6.7, fontweight="bold", color=color, zorder=6)


def icon_dataset(ax, cx, cy, color):
    for offset in (0.10, 0.05, 0.0):
        rounded(ax, (cx - 0.22 + offset, cy - 0.16 + offset), 0.33, 0.25, WHITE, color, 0.04, 0.8, 6)
        ax.plot([cx - 0.16 + offset, cx + 0.02 + offset], [cy + 0.01 + offset] * 2, color=color, lw=0.8, zorder=7)
        ax.plot([cx - 0.16 + offset, cx - 0.01 + offset], [cy - 0.05 + offset] * 2, color=color, lw=0.8, zorder=7)


def icon_filter(ax, cx, cy, color):
    verts = [(cx - 0.22, cy + 0.16), (cx + 0.22, cy + 0.16), (cx + 0.06, cy - 0.02), (cx + 0.06, cy - 0.18), (cx - 0.05, cy - 0.23), (cx - 0.05, cy - 0.02), (cx - 0.22, cy + 0.16)]
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * 5 + [MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=WHITE, edgecolor=color, lw=1.1, zorder=6))
    for dx, dy in [(-0.13, 0.08), (0.0, 0.08), (0.13, 0.08)]:
        ax.add_patch(Circle((cx + dx, cy + dy), 0.025, facecolor=color, edgecolor="none", zorder=7))


def icon_cells(ax, cx, cy, colors):
    coords = [(-0.16, 0.10), (0.0, 0.16), (0.17, 0.08), (-0.10, -0.10), (0.10, -0.10)]
    for idx, (dx, dy) in enumerate(coords):
        color = colors[idx % len(colors)]
        ax.add_patch(Circle((cx + dx, cy + dy), 0.095, facecolor=WHITE, edgecolor=color, lw=1.0, zorder=6))
        ax.add_patch(Circle((cx + dx, cy + dy), 0.030, facecolor=color, edgecolor="none", alpha=0.85, zorder=7))


def icon_trajectory(ax, cx, cy, color):
    for offset in (-0.09, 0.0, 0.09):
        ax.add_patch(Circle((cx - 0.20 + offset * 0.3, cy + offset), 0.035, facecolor=color, edgecolor="none", zorder=7))
    ax.plot([cx - 0.17, cx + 0.02, cx + 0.19], [cy, cy + 0.07, cy + 0.16], color=color, lw=1.3, zorder=6)
    ax.plot([cx + 0.02, cx + 0.20], [cy + 0.07, cy - 0.12], color=color, lw=1.3, zorder=6)
    ax.add_patch(Circle((cx + 0.20, cy + 0.16), 0.045, facecolor=WHITE, edgecolor=color, lw=1.0, zorder=7))
    ax.add_patch(Circle((cx + 0.20, cy - 0.12), 0.045, facecolor=WHITE, edgecolor=color, lw=1.0, zorder=7))


def icon_communication(ax, cx, cy, color_a, color_b):
    ax.add_patch(Circle((cx - 0.17, cy), 0.11, facecolor=WHITE, edgecolor=color_a, lw=1.2, zorder=6))
    ax.add_patch(Circle((cx + 0.17, cy), 0.11, facecolor=WHITE, edgecolor=color_b, lw=1.2, zorder=6))
    for x in (cx - 0.05, cx + 0.02, cx + 0.09):
        ax.add_patch(Circle((x, cy + 0.02), 0.025, facecolor=CORAL, edgecolor="none", zorder=7))
    arrow(ax, (cx - 0.06, cy - 0.04), (cx + 0.08, cy - 0.04), color=CORAL, lw=1.0, z=7)


def icon_network(ax, cx, cy, color):
    nodes = [(cx, cy + 0.16), (cx - 0.18, cy), (cx + 0.18, cy), (cx - 0.08, cy - 0.18), (cx + 0.12, cy - 0.17)]
    edges = [(0, 1), (0, 2), (1, 3), (2, 4), (3, 4), (1, 4)]
    for a, b in edges:
        ax.plot([nodes[a][0], nodes[b][0]], [nodes[a][1], nodes[b][1]], color=LINE, lw=0.9, zorder=5)
    for idx, (x, y) in enumerate(nodes):
        ax.add_patch(Circle((x, y), 0.052 if idx else 0.065, facecolor=color if idx == 0 else WHITE, edgecolor=color, lw=1.0, zorder=6))
    ax.plot([cx - 0.06, cx + 0.06], [cy + 0.16, cy + 0.16], color=WHITE, lw=1.2, zorder=7)


def icon_bulk(ax, cx, cy, color):
    for i, h in enumerate([0.10, 0.20, 0.29, 0.16]):
        ax.add_patch(FancyBboxPatch((cx - 0.23 + i * 0.12, cy - 0.17), 0.07, h, boxstyle="round,pad=0.01,rounding_size=0.015", facecolor=color, edgecolor="none", alpha=0.85, zorder=6))
    ax.plot([cx - 0.25, cx + 0.25], [cy - 0.18, cy - 0.18], color=MUTED, lw=0.8, zorder=6)


def icon_spatial(ax, cx, cy, color):
    for row in range(3):
        for col in range(4):
            x = cx - 0.17 + col * 0.11 + (row % 2) * 0.025
            y = cy + 0.12 - row * 0.12
            alpha = 0.35 + 0.15 * ((row + col) % 4)
            ax.add_patch(Circle((x, y), 0.035, facecolor=color, edgecolor="none", alpha=alpha, zorder=6))


def icon_liver(ax, cx, cy, color):
    verts = [
        (cx - 0.28, cy + 0.02),
        (cx - 0.18, cy + 0.20),
        (cx + 0.06, cy + 0.25),
        (cx + 0.28, cy + 0.13),
        (cx + 0.21, cy - 0.09),
        (cx - 0.03, cy - 0.20),
        (cx - 0.24, cy - 0.11),
        (cx - 0.28, cy + 0.02),
    ]
    codes = [MplPath.MOVETO] + [MplPath.CURVE3] * 6 + [MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=0.88, zorder=6))
    ax.plot([cx + 0.02, cx + 0.02], [cy + 0.18, cy - 0.12], color=WHITE, lw=1.0, alpha=0.8, zorder=7)


def method_card(ax, x, y, w, h, title, subtitle, face, color, icon):
    rounded(ax, (x, y), w, h, face, "#D5DEE5", 0.10, 0.8, 2)
    icon(ax, x + 0.38, y + h * 0.60, color)
    ax.text(x + 0.73, y + h * 0.68, title, fontsize=7.3, fontweight="bold", color=INK, ha="left", va="center", zorder=8)
    ax.text(x + 0.73, y + h * 0.38, subtitle, fontsize=5.7, color=MUTED, ha="left", va="center", linespacing=1.25, zorder=8)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.25, 5.25), dpi=600)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7.25)
    ax.axis("off")

    ax.text(0.30, 6.93, "Cross-cohort analytical workflow", fontsize=14, fontweight="bold", color=INK, ha="left", va="center")
    ax.text(0.30, 6.64, "From public human liver transcriptomes to a biopsy-resolved macrophage–HSC tissue model", fontsize=7.4, color=MUTED, ha="left", va="center")

    # Phase 1: input and QC.
    rounded(ax, (0.25, 4.70), 2.15, 1.65, BLUE_BG, "none", 0.16, 0.0, 1)
    phase_label(ax, 0.48, 6.10, 1, "Input + QC", BLUE)
    icon_liver(ax, 0.67, 5.52, BLUE)
    ax.text(1.325, 5.82, "Public liver sc/snRNA-seq", fontsize=6.4, fontweight="bold", color=INK, ha="center")
    cohort_y = [5.54, 5.30, 5.06]
    cohorts = [("GSE344087", "20 libraries"), ("GSE298719", "13 libraries"), ("GSE212837", "18 libraries")]
    for yy, (name, count) in zip(cohort_y, cohorts):
        rounded(ax, (1.10, yy - 0.09), 1.08, 0.20, WHITE, "#BFD3E1", 0.04, 0.65, 3)
        ax.text(1.18, yy, name, fontsize=5.4, fontweight="bold", color=BLUE, va="center")
        ax.text(2.08, yy, count.split()[0], fontsize=5.4, color=MUTED, ha="right", va="center", fontweight="bold")
    ax.text(0.51, 5.03, "51", fontsize=10.2, fontweight="bold", color=BLUE, ha="center")
    ax.text(0.51, 4.84, "libraries", fontsize=5.2, color=MUTED, ha="center")
    ax.text(1.13, 4.80, "QC  •  doublets  •  singlets", fontsize=5.1, color=MUTED, ha="left")

    # Phase 2: cell-state discovery.
    rounded(ax, (2.62, 4.70), 2.48, 1.65, PURPLE_BG, "none", 0.16, 0.0, 1)
    phase_label(ax, 2.85, 6.10, 2, "Cell-state discovery", PURPLE)
    icon_cells(ax, 3.08, 5.51, [BLUE, PURPLE, TEAL])
    ax.text(3.55, 5.85, "Lineage-restricted analysis", fontsize=6.1, fontweight="bold", color=INK, ha="left")
    ax.text(3.55, 5.58, "Myeloid / HSC-fibroblast", fontsize=5.7, color=MUTED, ha="left")
    ax.text(3.55, 5.31, "Per-cohort Harmony + Leiden", fontsize=5.7, color=MUTED, ha="left")
    rounded(ax, (2.87, 4.93), 0.97, 0.29, WHITE, "#CDBDDE", 0.06, 0.7, 3)
    rounded(ax, (3.95, 4.93), 0.91, 0.29, WHITE, "#CDBDDE", 0.06, 0.7, 3)
    ax.text(3.355, 5.075, "Marker-guided\nannotation", fontsize=5.5, color=PURPLE, ha="center", va="center", fontweight="bold")
    ax.text(4.405, 5.075, "Patient/sample\npseudobulk", fontsize=5.5, color=PURPLE, ha="center", va="center", fontweight="bold")

    # Phase 3: macrophage state.
    rounded(ax, (5.32, 4.70), 2.12, 1.65, TEAL_BG, "none", 0.16, 0.0, 1)
    phase_label(ax, 5.55, 6.10, 3, "Macrophage state", TEAL)
    icon_filter(ax, 5.80, 5.51, TEAL)
    ax.text(6.18, 5.84, "SLAMF8-marked state", fontsize=5.9, fontweight="bold", color=INK, ha="left")
    ax.text(6.23, 5.56, "Within-sample high / low", fontsize=5.7, color=MUTED, ha="left")
    ax.text(6.20, 5.29, "Ferro-aging + HSC scores", fontsize=5.2, color=MUTED, ha="left")
    rounded(ax, (5.46, 4.93), 1.84, 0.29, WHITE, "none", 0.06, 0.0, 3)
    ax.text(6.38, 5.075, "Cross-cohort random-effects synthesis", fontsize=5.0, color=TEAL, ha="center", va="center", fontweight="bold")

    # Phase 4: validation.
    rounded(ax, (7.66, 4.70), 2.09, 1.65, GREEN_BG, "none", 0.16, 0.0, 1)
    phase_label(ax, 7.89, 6.10, 4, "Tissue validation", GREEN)
    icon_bulk(ax, 8.05, 5.56, GREEN)
    ax.text(8.43, 5.86, "Bulk liver biopsies", fontsize=6.9, fontweight="bold", color=INK, ha="left")
    ax.text(8.43, 5.61, "GSE162694  •  n = 112", fontsize=5.7, color=MUTED, ha="left")
    ax.text(8.43, 5.39, "GSE130970  •  n = 78", fontsize=5.7, color=MUTED, ha="left")
    icon_spatial(ax, 8.08, 5.02, GREEN)
    ax.text(8.43, 5.08, "GSE192741 Visium", fontsize=5.7, fontweight="bold", color=GREEN, ha="left")
    ax.text(8.43, 4.89, "spot-level co-expression", fontsize=5.4, color=MUTED, ha="left")

    # Primary horizontal flow.
    for start, end, color in [((2.40, 5.58), (2.62, 5.58), BLUE), ((5.10, 5.58), (5.32, 5.58), PURPLE), ((7.44, 5.58), (7.66, 5.58), TEAL)]:
        arrow(ax, start, end, color=color, lw=1.25)

    # Mechanistic modules.
    ax.text(0.30, 4.48, "ORTHOGONAL ANALYSES", fontsize=6.7, fontweight="bold", color=ORANGE, ha="left")
    method_card(ax, 0.28, 3.36, 2.17, 0.92, "Trajectory", "Monocle3 + Slingshot\ntradeSeq audit", ORANGE_BG, ORANGE, icon_trajectory)
    method_card(ax, 2.63, 3.36, 2.17, 0.92, "Cell communication", "LIANA + CellChat\npatient-matched SPP1 audit", CORAL_BG, CORAL, lambda a, x, y, c: icon_communication(a, x, y, CORAL, GREEN))
    method_card(ax, 4.98, 3.36, 2.17, 0.92, "Ligand prioritization", "NicheNet\nHSC target programs", PURPLE_BG, PURPLE, icon_network)
    method_card(ax, 7.33, 3.36, 2.40, 0.92, "Counterfactual audit", "GRNBoost2 donor-aware\nSLAMF8 reduction", BLUE_BG, BLUE, icon_network)

    arrow(ax, (3.75, 4.70), (1.37, 4.28), color=PURPLE, lw=0.9, style="dashed", curve=-0.10)
    arrow(ax, (6.27, 4.70), (3.72, 4.28), color=TEAL, lw=0.9, style="dashed", curve=0.10)
    arrow(ax, (6.42, 4.70), (6.07, 4.28), color=TEAL, lw=0.9, style="dashed", curve=0.0)
    arrow(ax, (6.58, 4.70), (8.50, 4.28), color=TEAL, lw=0.9, style="dashed", curve=-0.12)

    # Convergent tissue model.
    rounded(ax, (0.28, 1.12), 9.45, 1.82, "#F8FAFB", "#D6E0E6", 0.16, 0.9, 1)
    ax.text(0.52, 2.69, "CONVERGENT TISSUE MODEL", fontsize=6.7, fontweight="bold", color=INK, ha="left")
    ax.text(0.52, 2.47, "Recurrent state identity, candidate output and fibrogenic response are tested as linked but distinct measurements", fontsize=5.8, color=MUTED, ha="left")

    nodes = [
        (0.70, 1.53, 1.42, "FCN1/VCAN-like\nmonocytes", BLUE_BG, BLUE),
        (2.55, 1.53, 1.58, "SLAMF8-marked\nSAMac-like state", TEAL_BG, TEAL),
        (4.58, 1.53, 1.18, "SPP1\noutput", CORAL_BG, CORAL),
        (6.20, 1.53, 1.46, "ITGAV/ITGB5\nHSC program", PURPLE_BG, PURPLE),
        (8.10, 1.53, 1.18, "Fibrosis-associated\ntissue state", GREEN_BG, GREEN),
    ]
    for x, y, w, label, face, color in nodes:
        rounded(ax, (x, y), w, 0.58, face, color, 0.10, 0.95, 4)
        ax.text(x + w / 2, y + 0.29, label, fontsize=6.3, fontweight="bold", color=INK, ha="center", va="center", linespacing=1.15, zorder=6)
    labels = ["state ordering", "candidate ligand", "matched association", "bulk replication"]
    for idx in range(len(nodes) - 1):
        start = (nodes[idx][0] + nodes[idx][2] + 0.03, 1.82)
        end = (nodes[idx + 1][0] - 0.03, 1.82)
        arrow(ax, start, end, color=[BLUE, TEAL, CORAL, GREEN][idx], lw=1.1)
        ax.text((start[0] + end[0]) / 2, 2.23, labels[idx], fontsize=5.0, color=MUTED, ha="center", va="center")

    # Footer and visual grammar.
    ax.plot([0.30, 9.70], [0.78, 0.78], color="#E1E7EB", lw=0.8)
    arrow(ax, (0.43, 0.48), (0.81, 0.48), color=INK, lw=1.0)
    ax.text(0.91, 0.48, "primary analytical sequence", fontsize=5.5, color=MUTED, va="center")
    arrow(ax, (2.73, 0.48), (3.11, 0.48), color=INK, lw=0.9, style="dashed")
    ax.text(3.21, 0.48, "orthogonal analysis branch", fontsize=5.5, color=MUTED, va="center")
    ax.text(9.70, 0.48, "Association and prioritization; no causal inference", fontsize=5.5, color=MUTED, ha="right", va="center", style="italic")

    fig.subplots_adjust(left=0.02, right=0.985, top=0.98, bottom=0.035)
    stem = WORK / "Figure_S1_cross_cohort_workflow"
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.03)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.03)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.03, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    tiff_path = stem.with_suffix(".tiff")
    with Image.open(tiff_path) as image:
        if image.mode != "RGB":
            image.convert("RGB").save(tiff_path, format="TIFF", compression="tiff_lzw", dpi=(600, 600))

    for ext in (".pdf", ".tiff"):
        shutil.copy2(stem.with_suffix(ext), OUT / f"Figure_S1{ext}")

    rows = [
        (1, "Input and quality control", "GSE344087; GSE298719; GSE212837", "Source-study QC, library-level scDblFinder and singlet retention"),
        (2, "Cell-state discovery", "Myeloid and HSC/fibroblast subsets", "Per-cohort normalization, Harmony, Leiden, marker-guided annotation and patient/sample pseudobulk"),
        (3, "Macrophage state definition", "SLAMF8-marked SAMac-like state", "Within-sample state contrasts, ferro-aging-related scores and cross-cohort synthesis"),
        (4, "Orthogonal analyses", "Trajectory; communication; ligand prioritization; counterfactual network", "Monocle3, Slingshot, tradeSeq, LIANA, CellChat, NicheNet and donor-aware GRNBoost2 audit"),
        (5, "Independent tissue validation", "GSE162694; GSE130970; GSE192741", "Bulk biopsy replication and Visium spot-level co-expression audit"),
        (6, "Convergent tissue model", "SLAMF8-marked SAMac-like state to SPP1-associated HSC program", "Linked state identity, candidate output and fibrosis-associated tissue response without causal inference"),
    ]
    pd.DataFrame(rows, columns=["order", "stage", "input_or_target", "analysis_or_interpretation"]).to_csv(SOURCE, sep="\t", index=False, compression="gzip")
    shutil.copy2(SOURCE, PACKAGE / "05_Source_Data" / SOURCE.name)
    print(stem)


if __name__ == "__main__":
    main()
