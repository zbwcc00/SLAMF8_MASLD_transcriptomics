#!/usr/bin/env python
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
OUT = PROJECT / "results" / "perturbation"
FIG = PROJECT / "figures"


def main() -> None:
    data = pd.read_csv(OUT / "035_SLAMF8_network_perturbation_dataset_summary.tsv", sep="\t")
    targets = ["FARG95", "Iron_homeostasis", "PUFA_ACSL4", "Lipid_peroxidation", "Senescence_SASP", "Inflammation", "ligand_TGFB1", "ligand_PDGF B", "ligand_PDGFB", "ligand_SPP1", "ligand_CXCL16", "ligand_TNF"]
    targets = [x for x in targets if x in data.target.unique()]
    plot = data[data.target.isin(targets)].pivot(index="target", columns="dataset", values="mean_delta").reindex(targets)
    FIG.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 5.8))
    sns.heatmap(plot, cmap="vlag", center=0, annot=True, fmt=".3g", linewidths=0.4, cbar_kws={"label": "Predicted Δ after SLAMF8 low-state"})
    plt.title("In silico SLAMF8 network perturbation")
    plt.xlabel("Dataset")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(FIG / "011_SLAMF8_network_perturbation_heatmap.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
