#!/usr/bin/env python3
import os
"""Create a reproducibility manifest for the P1 manuscript freeze."""
from pathlib import Path
import hashlib
import json
from datetime import datetime

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
OUT = ROOT / "results" / "revision_p0" / "p1_freeze"
OUT.mkdir(parents=True, exist_ok=True)

FILES = [
    ROOT / "018_论文结果段初稿_待导师审阅.md",
    ROOT / "030_图注与Results数字核验报告.md",
    ROOT / "031_核验后Results与Figure_legends.md",
    ROOT / "032_补充表整理与审计.md",
    ROOT / "033_Methods_manuscript_ready.md",
    ROOT / "036_P1最终冻结清单.md",
    ROOT / "figures" / "revision_p0" / "Figure2_SLAMF8_macrophage_state_v3_corrected.pdf",
    ROOT / "figures" / "revision_p0" / "Supplementary_Figure4_SLAMF8_anchor_excluded_bootstrap.pdf",
    ROOT / "figures" / "revision_p0" / "FigureS11_SLAMF8_counterfactual_audit_corrected.pdf",
    ROOT / "supplementary_tables" / "Supplementary_Tables_S1_S13.xlsx",
    ROOT / "038_P2_5结果与停止判定.md",
    ROOT / "039_Manuscript_Methods_formal_draft.md",
    ROOT / "040_Manuscript_Results_formal_draft.md",
    ROOT / "041_正式写作数字核验清单.md",
    ROOT / "042_Manuscript_Introduction_formal_draft.md",
    ROOT / "043_Introduction_citation_verification_log.md",
    ROOT / "044_审稿人导师生信综合审查_正式稿后.md",
    ROOT / "046_Manuscript_Discussion_formal_draft.md",
    ROOT / "047_Discussion_review_and_revision_log.md",
    ROOT / "048_Manuscript_Title_Abstract.md",
    ROOT / "049_Manuscript_integrated_draft.md",
    ROOT / "050_GSE192741空间审计与下一步结论.md",
    ROOT / "cache" / "discussion_citation_check_20260923" / "verified_pubmed_records.json",
    ROOT / "cache" / "discussion_citation_check_20260923" / "verified_pubmed_records.tsv",
    ROOT / "figures" / "revision_p2_5" / "P25_bulk_immune_signature_and_locked_model_pilot.pdf",
    ROOT / "figures" / "Figure1_myeloid_HSC_atlas_v2.pdf",
    ROOT / "figures" / "Figure1_myeloid_HSC_atlas_v2.svg",
    ROOT / "figures" / "Figure1_myeloid_HSC_atlas_v2.tiff",
    ROOT / "source_data" / "Source_Data_Figure1_single_cell_atlas.tsv.gz",
    ROOT / "source_data" / "README_Figure1.md",
    ROOT / "scripts" / "046_redraw_figure1_top_journal.py",
    ROOT / "results" / "revision_p2_5" / "bulk_immune_and_models" / "bulk_cell_state_signature_fibrosis_associations.tsv",
    ROOT / "results" / "revision_p2_5" / "bulk_immune_and_models" / "bulk_cell_state_signature_covariate_sensitivity.tsv",
    ROOT / "results" / "revision_p2_5" / "bulk_immune_and_models" / "locked_external_transcriptomic_model_pilot.tsv",
    ROOT / "scripts" / "061_samac_signature_specificity_audit.py",
    ROOT / "results" / "revision_p0" / "signature_specificity" / "signature_associations.tsv",
    ROOT / "results" / "revision_p0" / "signature_specificity" / "leave_one_gene_out.tsv",
    ROOT / "results" / "revision_p0" / "signature_specificity" / "matched_random_signature_null.tsv",
    ROOT / "results" / "revision_p0" / "signature_specificity" / "gene_coverage.tsv",
    ROOT / "results" / "revision_p0" / "signature_specificity" / "README.md",
    ROOT / "supplementary_tables" / "Table_S11_bulk_validation_results.tsv",
    ROOT / "scripts" / "062_samac_vs_myeloid_bootstrap_audit.py",
    ROOT / "results" / "revision_p2_5" / "samac_vs_myeloid_bootstrap_audit.tsv",
    ROOT / "scripts" / "063_build_manuscript_references.py",
    ROOT / "scripts" / "064_assemble_formal_manuscript.py",
    ROOT / "scripts" / "065_spatial_GSE192741_audit.py",
    ROOT / "results" / "spatial_validation" / "GSE192741" / "sample_level_spatial_summary.tsv",
    ROOT / "results" / "spatial_validation" / "GSE192741" / "spot_level_gene_panel_correlations.tsv",
    ROOT / "results" / "spatial_validation" / "GSE192741" / "README.md",
    ROOT / "cache" / "manuscript_references_20260923" / "references.md",
    ROOT / "cache" / "manuscript_references_20260923" / "references.ris",
]

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

rows = []
for path in FILES:
    rows.append({
        "relative_path": path.relative_to(ROOT).as_posix(),
        "exists": path.exists(),
        "bytes": path.stat().st_size if path.exists() else None,
        "sha256": sha256(path) if path.exists() else None,
    })

table_files = sorted((ROOT / "supplementary_tables").glob("Table_S[1-9]*.tsv"))
table_numbers = sorted({int(p.name.split("_", 2)[1][1:]) for p in table_files if p.name.startswith("Table_S")})
status = {
    "freeze": "P1",
    "generated_at": datetime.now().isoformat(timespec="seconds"),
    "project_root": str(ROOT),
    "key_files": rows,
    "supplementary_table_numbers": table_numbers,
    "supplementary_tables_complete": table_numbers == list(range(1, 14)),
    "revision_figures_complete": all(r["exists"] for r in rows[6:9]),
    "key_files_complete": all(r["exists"] for r in rows),
}
(OUT / "P1_freeze_manifest.json").write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8")
with (OUT / "P1_freeze_manifest.tsv").open("w", encoding="utf-8") as fh:
    fh.write("relative_path\texists\tbytes\tsha256\n")
    for row in rows:
        fh.write(f"{row['relative_path']}\t{row['exists']}\t{row['bytes']}\t{row['sha256']}\n")
print(json.dumps({k: status[k] for k in ("freeze", "supplementary_table_numbers", "supplementary_tables_complete", "revision_figures_complete", "key_files_complete")}, ensure_ascii=False))
