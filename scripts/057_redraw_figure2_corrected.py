#!/usr/bin/env python
import os
"""Re-render Figure 2 using anchor-excluded FARG95 results."""
from pathlib import Path
import importlib.util
import pandas as pd
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

OUTPUT_DPI = 600

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
OUT = ROOT / "figures" / "revision_p0"
OUT.mkdir(parents=True, exist_ok=True)
source = ROOT / "scripts" / "047_redraw_figure2_top_journal.py"
spec = importlib.util.spec_from_file_location("figure2_source", source)
mod = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(mod)
mod.FIG = OUT
mod.FEATURE_LABELS["module_FARG95"] = "FARG95 (anchor-excluded)"

old_combined, old_meta = mod.read_state_tables()
old_combined = old_combined[old_combined.feature != "module_FARG95"].copy()
old_meta = old_meta[old_meta.feature != "module_FARG95"].copy()
corr_dir = ROOT / "results" / "revision_p0" / "farg95_self_containment"
corr = pd.read_csv(corr_dir / "cross_dataset_paired_tests.tsv", sep="\t")
corr = corr[corr.feature == "FARG95_anchor_excluded"].copy()
corr["feature"] = "module_FARG95"
corr["fdr"] = corr["fdr_within_dataset"]
corr["wilcoxon_statistic"] = corr["wilcoxon_statistic"]
meta = pd.read_csv(corr_dir / "cross_dataset_farg95_exclusion_meta.tsv", sep="\t")
meta = meta[meta.feature == "FARG95_anchor_excluded"].copy()
meta["feature"] = "module_FARG95"
meta = meta.rename(columns={"signed_stouffer_fdr": "signed_stouffer_fdr"})
mod.read_state_tables = lambda: (pd.concat([old_combined, corr], ignore_index=True, sort=False), pd.concat([old_meta, meta], ignore_index=True, sort=False))
mod.save_main()
mod.save_bootstrap_supplement()

source_data = ROOT / "source_data"
source_data.mkdir(parents=True, exist_ok=True)
state_combined, state_meta = mod.read_state_tables()
disease_combined, disease_meta = mod.read_disease_tables()
source_tables = [
    ("panel_A_state_paired", state_combined),
    ("panel_A_state_meta", state_meta),
    ("panel_B_disease_group", disease_combined),
    ("panel_B_disease_meta", disease_meta),
    ("panel_C_macrophage_HSC", pd.read_csv(ROOT / "results" / "033_cross_dataset_macrophage_HSC_correlations.tsv", sep="\t")),
    ("panel_D_receptor_program_meta", pd.read_csv(ROOT / "results" / "034_cross_dataset_macrophage_HSC_correlation_meta.tsv", sep="\t")),
]
source_rows = []
for panel, table in source_tables:
    copy = table.copy()
    copy.insert(0, "panel", panel)
    source_rows.append(copy)
pd.concat(source_rows, ignore_index=True, sort=False).to_csv(
    source_data / "Source_Data_Figure2_SLAMF8_state.tsv.gz", sep="\t", index=False, compression="gzip"
)
bootstrap = pd.read_csv(ROOT / "results" / "035_Figure2_state_bootstrap_ci.tsv", sep="\t")
bootstrap.to_csv(
    source_data / "Source_Data_Supplementary_Figure4_anchor_excluded_bootstrap.tsv.gz",
    sep="\t", index=False, compression="gzip"
)
for old, new in [
    (OUT / "Figure2_SLAMF8_macrophage_state_v2.pdf", OUT / "Figure2_SLAMF8_macrophage_state_v3_corrected.pdf"),
    (OUT / "Figure2_SLAMF8_macrophage_state_v2.svg", OUT / "Figure2_SLAMF8_macrophage_state_v3_corrected.svg"),
    (OUT / "Figure2_SLAMF8_macrophage_state_v2.png", OUT / "Figure2_SLAMF8_macrophage_state_v3_corrected.png"),
    (OUT / "Figure2_SLAMF8_macrophage_state_v2.tiff", OUT / "Figure2_SLAMF8_macrophage_state_v3_corrected.tiff"),
    (OUT / "Supplementary_Figure2_SLAMF8_state_bootstrap.pdf", OUT / "Supplementary_Figure4_SLAMF8_anchor_excluded_bootstrap.pdf"),
    (OUT / "Supplementary_Figure2_SLAMF8_state_bootstrap.svg", OUT / "Supplementary_Figure4_SLAMF8_anchor_excluded_bootstrap.svg"),
    (OUT / "Supplementary_Figure2_SLAMF8_state_bootstrap.png", OUT / "Supplementary_Figure4_SLAMF8_anchor_excluded_bootstrap.png"),
    (OUT / "Supplementary_Figure2_SLAMF8_state_bootstrap.tiff", OUT / "Supplementary_Figure4_SLAMF8_anchor_excluded_bootstrap.tiff"),
]:
    if old.exists():
        old.replace(new)
print("Corrected Figure 2 and state bootstrap figures written to", OUT)
