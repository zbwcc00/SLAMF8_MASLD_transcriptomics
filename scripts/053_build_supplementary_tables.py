#!/usr/bin/env python
import os
"""Aggregate frozen result files into manuscript supplementary tables only."""
from pathlib import Path
import pandas as pd

P = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
R, O = P / "results", P / "supplementary_tables"
O.mkdir(exist_ok=True)

def read(path):
    try:
        return pd.read_csv(path, sep="\t", compression="infer", dtype=object)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()

def combine(paths, unit, note):
    xs = []
    for path in paths:
        d = read(path)
        if not d.empty:
            # Some freeze manifests already contain a provenance field.
            for col in ("source_file", "analysis_unit", "interpretation_note"):
                if col in d.columns:
                    d = d.drop(columns=col)
            d.insert(0, "source_file", str(path.relative_to(P)).replace("\\", "/"))
            d.insert(1, "analysis_unit", unit)
            d.insert(2, "interpretation_note", note)
            xs.append(d)
    return pd.concat(xs, ignore_index=True, sort=False) if xs else pd.DataFrame()

def save(n, stem, data, description, method, sources):
    if data.empty:
        data = pd.DataFrame({"status": ["NO_ROWS_FOUND"], "source_pattern": [sources]})
    data.to_csv(O / f"Table_S{n}_{stem}.tsv", sep="\t", index=False, encoding="utf-8-sig")
    (O / f"Table_S{n}_{stem}.README.md").write_text(
        f"# Table S{n}: {stem.replace('_', ' ')}\n\n"
        f"- Description: {description}\n"
        f"- Statistical method / unit: {method}\n"
        f"- Source pattern: `{sources}`\n"
        "- Provenance: each row retains `source_file`; this is an aggregation of frozen result files.\n",
        encoding="utf-8")

# S1. Curated cohort inventory, including the role and analysis unit used here.
rows = [
 ["GSE212837","scRNA-seq","human liver","Control/NASH","patient-level single-cell","myeloid/HSC discovery","results/006_GSE212837_R_input_summary.tsv"],
 ["GSE298719","scRNA-seq","human liver biopsy","Healthy/MASLD","patient-level single-cell","communication/trajectory","results/006_GSE298719_R_input_summary.tsv"],
 ["GSE344087","scRNA-seq","human liver","no fibrosis/fibrosis","patient-level single-cell","myeloid/HSC discovery","results/006_GSE344087_R_input_summary.tsv"],
 ["GSE162694","bulk RNA-seq","human liver biopsy","NAFLD spectrum","biopsy/sample","fibrosis validation","results/bulk_validation/full_bulk_fibrosis_validation.tsv"],
 ["GSE130970","bulk RNA-seq","human liver biopsy","NAFLD spectrum","biopsy/sample","fibrosis validation","results/bulk_validation/GSE130970_full_bulk_scores.tsv"],
 ["GSE167523","bulk RNA-seq","human liver biopsy","NAFL/NASH","independent biopsy","external validation","results/bulk_validation/expanded_public/GSE167523_scores.tsv"],
 ["GSE126848","bulk RNA-seq","human liver","NAFLD","scored sample","external validation","results/bulk_validation/expanded_public/GSE126848_scores.tsv"],
 ["GSE193066","bulk RNA-seq","human liver biopsy","NAFLD longitudinal","paired visit","follow-up audit","results/bulk_validation/expanded_public/GSE193066_HUn164_paired.tsv"],
 ["GSE281797","bulk RNA-seq","human liver","MASL/MASH","sample","additional validation","results/bulk_validation/expanded_public/GSE281797_scores.tsv"],
 ["GSE83452","bulk RNA-seq","human liver","NASH intervention follow-up","paired follow-up","intervention audit","results/bulk_validation/expanded_public/GSE83452_paired_followup.tsv"],
 ["GSE192741","spatial transcriptomics (Visium)","human liver","healthy/steatotic","spot-level exploratory audit","independent spatial context","results/spatial_validation/GSE192741/sample_level_spatial_summary.tsv"]]
s1 = pd.DataFrame(rows, columns=["accession","data_type","tissue","phenotype","analysis_unit","role_in_study","manifest_or_result_source"])
s1.insert(0, "source_file", s1["manifest_or_result_source"])
s1.insert(1, "interpretation_note", "Cohort-level inventory; sample-level counts and analysis-specific denominators are reported in the linked source result files.")
save(1,"Cohort_and_sample_manifest",s1,"Cohorts used for discovery and validation.","Cohort-level manifest.","curated from results/")

ds = ["GSE212837","GSE298719","GSE344087"]
save(2,"scRNA_library_QC_and_doublet_summary",combine(sorted(R.glob("002_*_scDblFinder_summary.tsv")),"library-level QC","scDblFinder class counts; fraction is within library."),"Per-library doublet calls.","Library-level counts/fractions; cells are not patient replicates.","results/002_*")
save(3,"reclustering_cell_composition",combine(sorted(R.glob("003_*_myeloid_sample_subtype_counts.tsv"))+sorted(R.glob("004_*_hsc_sample_subtype_counts.tsv")),"sample × subtype","Descriptive cell counts; not independent patient replicates."),"Myeloid/HSC subtype counts after reclustering.","Sample × subtype cell counts.","results/003_*sample_subtype_counts; results/004_*sample_subtype_counts")
save(4,"marker_validation_summary",combine(sorted(R.glob("003_*_myeloid_cluster_signature_scores.tsv"))+sorted(R.glob("004_*_hsc_cluster_signature_scores.tsv")),"cluster-level signature score","Annotation support only; not differential-expression effect sizes."),"Cluster signature scores and curated labels.","Cluster-level scores.","results/003_*cluster_signature_scores; results/004_*cluster_signature_scores")

s5 = combine(sorted(R.glob("011_*_trajectory_summary.tsv"))+[R/"036_Figure3_trajectory_audit_summary.tsv"],"dataset-level trajectory audit","Monocle3 and Slingshot provide complementary inference.")
s5b = combine(sorted(R.glob("009_*_tradeSeq_association.tsv.gz"))+sorted(R.glob("010_*_tradeSeq_start_vs_end.tsv.gz")),"gene-level trajectory test","Endpoint significance is dataset-specific.")
save(5,"trajectory_results",pd.concat([s5,s5b],ignore_index=True,sort=False),"Trajectory audit plus tradeSeq association and endpoint tests.","Dataset audit and gene-level tradeSeq FDR.","results/009_*,010_*,011_*,036_*")
farg95_revision = R / "revision_p0" / "farg95_self_containment" / "cross_dataset_farg95_exclusion_meta.tsv"
continuous_revision = R / "revision_p0" / "continuous_meta" / "continuous_random_effects_meta.tsv"
leaveout_revision = R / "revision_p0" / "continuous_meta" / "continuous_leave_one_cohort_out.tsv"
threshold_revision = R / "revision_p0" / "continuous_meta" / "threshold_sensitivity.tsv"
s6_revision = combine([farg95_revision, continuous_revision, leaveout_revision, threshold_revision], "revision-aware cross-cohort sensitivity", "Prespecified anchor-exclusion FARG95 sensitivity, continuous SLAMF8 random-effects meta-analysis, leave-one-cohort-out and threshold sensitivity; only SLAMF8 among the five listed anchors is a member of FARG95, so the anchor-excluded and SLAMF8-excluded scores are identical.")
save(6,"SLAMF8_high_low_paired_results",pd.concat([combine(sorted(R.glob("016_*_SLAMF8_high_paired_tests.tsv"))+[R/"035_Figure2_state_bootstrap_ci.tsv",R/"022_cross_dataset_SLAMF8_high_meta.tsv"],"paired sample-level state contrast","Within-sample paired contrasts; bootstrap CIs quantify uncertainty."), s6_revision], ignore_index=True, sort=False),"SLAMF8-high/low paired tests plus revision-aware anchor-exclusion and continuous sensitivity analyses.","Paired Wilcoxon, bootstrap CI, directional Stouffer, random-effects Fisher-z meta-analysis, leave-one-cohort-out and quantile sensitivity.","results/016_*,035_*,022_*; results/revision_p0/farg95_self_containment/*; continuous_meta/*")
save(7,"macrophage_HSC_matched_associations",combine(sorted(R.glob("032_*_macrophage_HSC_correlations.tsv"))+[R/"033_cross_dataset_macrophage_HSC_correlations.tsv",R/"034_cross_dataset_macrophage_HSC_correlation_meta.tsv"],"patient-level matched association","Spearman associations; pooled rows are meta-analytic summaries."),"Matched macrophage–HSC associations.","Patient-level Spearman correlation and pooled estimates.","results/032_*,033_*,034_*")

c = R/"communication"
# Use focused macrophage→HSC, pathway and audit outputs for the submission
# table. Full LIANA consensus matrices (017) stay archived in results/.  The
# focused LIANA matrices have 85k+ rows each, so retain the top 1% by consensus
# aggregate rank (lower is stronger); the source-file field preserves the full
# matrix location and Table S8 remains spreadsheet-reviewable.
s8paths = sorted((c/"audit").glob("*.tsv"))+sorted(c.glob("019_*CellChat_macrophage_to_HSC.tsv"))+sorted(c.glob("021_*pathway_summary.tsv"))+sorted(c.glob("022_*pathway_summary.tsv"))+sorted(c.glob("023_*interactions.tsv"))
s8 = combine(s8paths,"inferred communication row or audit","Transcriptome-based inference; no spatial proximity, secretion or causality established.")
liana_focus = []
for path in sorted(c.glob("018_*LIANA_macrophage_to_HSC.tsv")):
    d = read(path)
    if not d.empty:
        d["submission_filter"] = "top 1% macrophage→HSC LIANA rows by aggregate_rank (lower is stronger)"
        d = d.loc[pd.to_numeric(d["aggregate_rank"], errors="coerce") <= pd.to_numeric(d["aggregate_rank"], errors="coerce").quantile(0.01)]
        d.insert(0,"source_file",str(path.relative_to(P)).replace("\\","/"))
        d.insert(1,"analysis_unit","inferred macrophage→HSC interaction")
        d.insert(2,"interpretation_note","Transcriptome-based inference; no spatial proximity, secretion or causality established.")
        liana_focus.append(d)
if liana_focus:
    s8 = pd.concat([s8]+liana_focus,ignore_index=True,sort=False)
save(8,"communication_LIANA_CellChat_audit",s8,"LIANA/CellChat inferred communication and audits.","Computational candidate inference only; LIANA table includes the top 1% macrophage→HSC consensus-ranked rows, while full matrices are archived.","results/communication/018_* (top 1%); 019_*,021_*,022_*,023_*; audit/*.tsv")
save(9,"cross_cohort_SPP1_patient_level_audit",combine(sorted((c/"cross_cohort_spp1").glob("*.tsv")),"patient-level ligand/receptor audit","SPP1 is context-dependent; no SLAMF8→SPP1 causality is claimed."),"SPP1 source/target expression and association audits.","Patient-level and condition-adjusted summaries.","results/communication/cross_cohort_spp1/*.tsv")
save(10,"NicheNet_ligand_activity",combine(sorted(c.glob("027_*NicheNet_ligand_activities.tsv"))+sorted(c.glob("028_*NicheNet_ligand_target_links.tsv")),"dataset-level ligand activity/link","Exploratory ligand–target prioritization, not causal validation."),"NicheNet ligand activity and target links.","Reported AUROC/AUPR and ranking.","results/communication/027_*;028_*")

b = R/"bulk_validation"
bulk = [b/"full_bulk_fibrosis_validation.tsv",b/"full_bulk_SLAMF8_module_correlations.tsv",b/"full_bulk_covariate_sensitivity.tsv"] + sorted((b/"expanded_public").glob("*.tsv"))
p25 = R / "revision_p2_5" / "bulk_immune_and_models"
bulk += [p25 / "bulk_cell_state_signature_scores.tsv", p25 / "bulk_cell_state_signature_coverage.tsv",
         p25 / "bulk_cell_state_signature_fibrosis_associations.tsv", p25 / "bulk_cell_state_signature_cross_feature_associations.tsv",
         p25 / "bulk_cell_state_signature_covariate_sensitivity.tsv", p25 / "locked_external_transcriptomic_model_pilot.tsv"]
specificity = R / "revision_p0" / "signature_specificity"
bulk += [specificity / "signature_associations.tsv", specificity / "leave_one_gene_out.tsv",
         specificity / "matched_random_signature_null.tsv", specificity / "gene_coverage.tsv",
         R / "revision_p2_5" / "samac_vs_myeloid_bootstrap_audit.tsv"]
save(11,"bulk_validation_results",combine(bulk,"bulk sample or paired biopsy","Bulk validation independent of single-cell discovery; adjusted estimates retain complete-case n."),"Independent bulk validation and sensitivity results.","Correlation, AUC, paired deltas and covariate sensitivity as reported.","results/bulk_validation/*.tsv; expanded_public/*.tsv")

p = R/"perturbation"
legacy_p12 = combine(sorted(p.glob("035_*dataset_summary.tsv"))+sorted(p.glob("036_*meta.tsv"))+sorted(p.glob("037_*sensitivity.tsv"))+sorted(p.glob("038_*focus.tsv")),"legacy archived counterfactual summary","Archived legacy output; not used for formal perturbation inference.")
revision_p12 = combine([R/"revision_p0"/"virtual_perturbation"/"donor_level_counterfactual_audit.tsv", R/"revision_p0"/"virtual_perturbation"/"cross_dataset_counterfactual_audit.tsv"], "donor-level counterfactual audit", "Model-based counterfactual prediction with donor bootstrap intervals, sign-flip audit, direction consistency and group-CV R2; not experimental knockout and no causal P value.")
save(12,"virtual_SLAMF8_perturbation_audit",pd.concat([revision_p12, legacy_p12], ignore_index=True, sort=False),"Revision-aware virtual SLAMF8 perturbation audit with archived legacy outputs.","Donor-level predicted deltas, bootstrap intervals, sign-flip permutation audit, direction consistency and group-CV R2; legacy Fisher P values are not formal inference.","results/revision_p0/virtual_perturbation/*; results/perturbation/035_*,036_*,037_*,038_*")
freeze = sorted((R/"final_freeze").glob("*.tsv"))+sorted((b/"expanded_public").glob("*manifest*.tsv"))
save(13,"data_freeze_and_unit_audit",combine(freeze,"freeze/unit audit","Reproducibility metadata rather than biological effect estimates."),"Data freeze, sample-unit checks and audit metadata.","Audit-level records; SHA256 retained where supplied.","results/final_freeze/*.tsv; expanded_public/*manifest*.tsv")
spatial = combine([
    R / "spatial_validation" / "GSE192741" / "sample_level_spatial_summary.tsv",
    R / "spatial_validation" / "GSE192741" / "spot_level_gene_panel_correlations.tsv",
], "Visium spot-level exploratory audit", "Mixed spot measurements; not patient-level independent replicates; no spatial proximity or causal signaling inferred.")
save(14, "GSE192741_spatial_coexpression_audit", spatial, "Independent human Visium spot-level co-expression audit for Figure S12.", "Within-specimen Spearman correlations between prespecified genes and HSC/fibrogenic scores.", "results/spatial_validation/GSE192741/*.tsv")

with pd.ExcelWriter(O/"Supplementary_Tables_S1_S14.xlsx",engine="openpyxl") as x:
    for n in range(1,15):
        f = next(O.glob(f"Table_S{n}_*.tsv"))
        pd.read_csv(f,sep="\t",dtype=object).to_excel(x,sheet_name=f"Table S{n}",index=False)
(O/"Supplementary_Tables_README.md").write_text("# Supplementary Tables S1–S14\n\n本目录由 `scripts/053_build_supplementary_tables.py` 从冻结的结果文件汇总而成，不新增生物学分析。\n\n- 单细胞 QC 是 library/cell 单位；样本级比较与相关分析以 patient/sample 为单位。\n- GSE192741 单独列为 Table S14，与 Figure S12 对应；spot 是混合测量单位，不作为患者级独立重复。\n- LIANA、CellChat 和 NicheNet 是转录组推断或候选排序，不能单独证明空间邻近、蛋白分泌或因果。\n- Table S12 是反事实网络预测，不是实验性 SLAMF8 敲除。\n- Table S11 的 SAMac-like signature specificity audit（leave-one-gene-out 和表达量十分位匹配随机基因集）是敏感性分析，不是独立队列。\n- 每张表均保留 `source_file`、`analysis_unit` 与 `interpretation_note` 以保证可追溯性。\n\nS1 队列；S2 QC/双细胞；S3 组成；S4 注释；S5 拟时序；S6 SLAMF8-high/low；S7 巨噬细胞–HSC 关联；S8 通讯；S9 SPP1 审计；S10 NicheNet；S11 bulk；S12 虚拟扰动；S13 冻结与单位审计；S14 GSE192741 Visium 空间共表达审计。\n",encoding="utf-8")
print(f"Wrote S1–S14 to {O}")
