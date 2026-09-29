#!/usr/bin/env python
"""Inventory PRIDE files and controlled-access NAFLD/NASH resources."""

from __future__ import annotations
import os

import json
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
OUT = ROOT / "results" / "public_cohort_search"
OUT.mkdir(parents=True, exist_ok=True)


def pride_files(accession: str) -> list[dict]:
    records = []
    for page in range(100):
        response = requests.get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files", params={"page": page}, timeout=180)
        response.raise_for_status()
        batch = response.json()
        if not batch:
            break
        records.extend(batch)
    return records


def inventory_pride() -> None:
    accessions = ["PXD052798", "PXD052784", "PXD052787", "PXD043340", "PXD073651"]
    rows = []
    for accession in accessions:
        metadata = requests.get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}", timeout=120).json()
        files = pride_files(accession)
        (OUT / f"{accession}_project.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
        for item in files:
            ftp = next((x.get("value", "") for x in item.get("publicFileLocations", []) if x.get("name") == "FTP Protocol"), "")
            rows.append({"accession": accession, "title": metadata.get("title", ""), "file_name": item.get("fileName", ""), "size_bytes": item.get("fileSizeBytes", ""), "ftp_url": ftp, "file_category": item.get("fileCategory", {}).get("name", "")})
    pd.DataFrame(rows).to_csv(OUT / "PRIDE_project_file_inventory.tsv", sep="\t", index=False)


def inventory_dbgap() -> None:
    terms = ["NAFLD", "NASH", "MASH", "MASLD", '"nonalcoholic fatty liver disease"']
    rows = []
    seen = set()
    for term in terms:
        response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params={"db": "gap", "term": term, "retmode": "json", "retmax": 200}, timeout=120)
        response.raise_for_status()
        ids = response.json()["esearchresult"]["idlist"]
        if not ids:
            continue
        summary = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", params={"db": "gap", "id": ",".join(ids), "retmode": "json"}, timeout=180).json().get("result", {})
        for uid in ids:
            item = summary.get(uid, {})
            study = item.get("d_study_results", {})
            study_id = study.get("d_study_id", "")
            if not study_id or study_id in seen:
                continue
            seen.add(study_id)
            rows.append({"query": term, "uid": uid, "study_id": study_id, "study_name": study.get("d_study_name", ""), "study_design": study.get("d_study_design", ""), "has_sra": study.get("d_study_has_sra", ""), "archive": study.get("d_study_archive", ""), "molecular_data_types": "; ".join(x.get("d_molecular_data_type_name", "") for x in study.get("d_study_molecular_data_type_list", []) or []), "participants": study.get("d_num_participants_in_subtree", ""), "url": f"https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id={study_id.split('.')[0]}"})
    pd.DataFrame(rows).to_csv(OUT / "dbGaP_NAFLD_NASH_controlled_access.tsv", sep="\t", index=False)


def write_access_links() -> None:
    links = {
        "NIDDK Central Repository": "https://repository.niddk.nih.gov/",
        "NIDDK study search": "https://repository.niddk.nih.gov/study/",
        "dbGaP": "https://www.ncbi.nlm.nih.gov/gap/",
        "EGA": "https://ega-archive.org/",
        "UK Biobank application": "https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access",
        "UK Biobank Showcase": "https://biobank.ndph.ox.ac.uk/showcase/",
        "All of Us Researcher Workbench": "https://researchallofus.org/",
        "NIDDK NAFLD Microbiome phs001837": "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001837.v1.p1",
        "NASH CRN iron deposition paper": "https://pubmed.ncbi.nlm.nih.gov/41264906/",
    }
    (OUT / "controlled_access_links.json").write_text(json.dumps(links, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    inventory_pride()
    inventory_dbgap()
    write_access_links()
    print("PRIDE and dbGaP inventories written to", OUT, flush=True)
