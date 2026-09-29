#!/usr/bin/env python
"""Sensitivity of counterfactual results to SLAMF8 low-state quantile."""

from __future__ import annotations
import os

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
SCRIPT = PROJECT / "scripts" / "011_virtual_SLAMF8_network_perturbation.py"
OUT = PROJECT / "results" / "perturbation"


def load_module():
    spec = importlib.util.spec_from_file_location("slamf8_perturbation", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    module = load_module()
    rows = []
    for dataset in ["GSE344087", "GSE298719", "GSE212837"]:
        print(f"sensitivity {dataset}", flush=True)
        obs, obj = module.read_sample(PROJECT, dataset, 400)
        genes = module.choose_genes(obj, 1200)
        expr = module.normalized_expression(obj, genes)
        network = pd.read_csv(OUT / f"030_{dataset}_SLAMF8_GRN_consensus.tsv", sep="\t")
        scores = module.regulon_scores(expr, network)
        features = pd.concat([expr[["SLAMF8"]], scores], axis=1)
        targets_dict = {name: obs[col].astype(float).to_numpy() for name, col in module.MODULES.items() if col in obs}
        for ligand in module.LIGANDS:
            if ligand in expr.columns:
                targets_dict[f"ligand_{ligand}"] = expr[ligand].to_numpy()
        targets = pd.DataFrame(targets_dict, index=expr.index)
        for quantile in [0.05, 0.10, 0.20]:
            slamf8_low = obs.assign(_expr=expr["SLAMF8"].to_numpy()).groupby("_subtype", observed=True)["_expr"].transform(lambda x: np.quantile(x, quantile)).to_numpy()
            pred, _ = module.ridge_counterfactual(features, targets, obs, slamf8_low)
            pred["donor"] = obs["_donor"].astype(str).to_numpy()
            for target in targets.columns:
                cell_delta = pred[f"delta_{target}"]
                donor_delta = pd.DataFrame({"donor": pred.donor, "delta": cell_delta}).groupby("donor").delta.mean()
                rows.append({"dataset": dataset, "quantile": quantile, "target": target, "donor_n": len(donor_delta), "mean_delta_patient": donor_delta.mean(), "median_delta_patient": donor_delta.median(), "fraction_down_patient": float((donor_delta < 0).mean()), "fraction_up_patient": float((donor_delta > 0).mean())})
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "037_SLAMF8_perturbation_sensitivity.tsv", sep="\t", index=False)
    focus = result[result.target.isin(["FARG95", "Iron_homeostasis", "Senescence_SASP", "Inflammation", "ligand_PDGFB", "ligand_SPP1", "ligand_CXCL16", "ligand_TNF"])]
    focus.to_csv(OUT / "038_SLAMF8_perturbation_sensitivity_focus.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
