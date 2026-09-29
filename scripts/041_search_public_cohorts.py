#!/usr/bin/env python
import os
"""Search GEO, PRIDE, BioStudies and PubMed for public MASLD/NASH cohorts."""

from pathlib import Path
import json
import re
import requests
import pandas as pd


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
OUT = ROOT / "results" / "public_cohort_search"
OUT.mkdir(parents=True, exist_ok=True)


def ncbi_search(term: str, database: str, retmax: int = 100) -> list[str]:
    payload = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params={"db": database, "term": term, "retmode": "json", "retmax": retmax}, timeout=60).json()
    return payload["esearchresult"]["idlist"]


def ncbi_summary(ids: list[str], database: str) -> list[dict]:
    if not ids:
        return []
    result = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", params={"db": database, "id": ",".join(ids), "retmode": "json"}, timeout=60).json()["result"]
    return [result[item] for item in ids if item in result]


def search_geo() -> pd.DataFrame:
    terms = [
        "NAFLD AND liver biopsy AND fibrosis",
        "NASH AND liver transcriptome",
        "MASLD AND liver biopsy",
        "NASH Clinical Research Network",
        "NAFLD AND longitudinal",
    ]
    records = {}
    for term in terms:
        for item in ncbi_summary(ncbi_search(term, "gds", 100), "gds"):
            accession = item.get("accession")
            if accession and accession.startswith("GSE"):
                records[accession] = {"accession": accession, "title": item.get("title", ""), "summary": item.get("summary", ""), "samples": item.get("samples", ""), "platform": item.get("gpl", ""), "date": item.get("pdat", "")}
    table = pd.DataFrame(records.values())
    if len(table):
        keywords = r"NAFLD|NASH|MASLD|fatty liver|steatohepatitis|liver biopsy"
        table = table[table.title.str.contains(keywords, case=False, regex=True, na=False) | table.summary.str.contains(keywords, case=False, regex=True, na=False)].copy()
        table["relevance"] = table.title.str.contains(r"biopsy|longitudinal|follow|fibrosis|NASH|MASLD", case=False, regex=True, na=False).astype(int) + table.summary.str.contains(r"biopsy|longitudinal|follow|fibrosis|NASH|MASLD", case=False, regex=True, na=False).astype(int)
        table = table.sort_values(["relevance", "date"], ascending=[False, False])
    return table


def search_pride() -> pd.DataFrame:
    records = {}
    for keyword in ["NAFLD", "NASH", "MASH", "fatty liver", "nonalcoholic steatohepatitis"]:
        response = requests.get("https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects", params={"keyword": keyword, "page": 0, "pageSize": 100}, timeout=120)
        response.raise_for_status()
        items = response.json()
        for item in items:
            accession = item.get("accession")
            text = " ".join(str(item.get(key, "")) for key in ["title", "projectDescription", "keywords"])
            if accession and re.search(r"NAFLD|NASH|MASH|fatty liver|steatohepatitis|liver", text, flags=re.I):
                records[accession] = {"accession": accession, "title": item.get("title", ""), "description": item.get("projectDescription", ""), "submissionDate": item.get("submissionDate", ""), "numAssays": item.get("numAssays", ""), "organisms": ";".join(map(str, item.get("organisms", []) or []))}
    table = pd.DataFrame(records.values())
    if len(table):
        table["relevance"] = table.title.str.contains(r"NAFLD|NASH|MASH|fatty liver|steatohepatitis", case=False, regex=True, na=False).astype(int) * 2 + table.description.str.contains(r"NAFLD|NASH|MASH|fatty liver|steatohepatitis", case=False, regex=True, na=False).astype(int)
        table = table.sort_values(["relevance", "submissionDate"], ascending=[False, False])
    return table


def search_biostudies() -> pd.DataFrame:
    records = []
    for keyword in ["NAFLD", "NASH", "MASLD", "fatty liver"]:
        data = requests.get("https://www.ebi.ac.uk/biostudies/api/v1/search", params={"query": keyword, "pageSize": 100}, timeout=120).json()
        for item in data.get("hits", []):
            text = " ".join(str(item.get(key, "")) for key in ["accession", "title", "description"])
            if re.search(r"NAFLD|NASH|MASLD|fatty liver|steatohepatitis", text, flags=re.I):
                records.append(item)
    if not records:
        return pd.DataFrame()
    table = pd.DataFrame(records).drop_duplicates(subset=["accession"])
    return table


def search_pubmed() -> pd.DataFrame:
    terms = [
        '"NASH Clinical Research Network"',
        '"NASH CRN" AND liver biopsy',
        'NIDDK AND NAFLD AND biopsy',
        'NAFLD AND longitudinal AND biopsy',
    ]
    records = []
    for term in terms:
        records.extend(ncbi_summary(ncbi_search(term, "pubmed", 100), "pubmed"))
    if not records:
        return pd.DataFrame()
    rows = []
    for item in {str(x.get("uid")): x for x in records}.values():
        rows.append({"pmid": item.get("uid"), "title": item.get("title", ""), "pubdate": item.get("pubdate", ""), "fulljournalname": item.get("fulljournalname", ""), "elocationid": item.get("elocationid", "")})
    return pd.DataFrame(rows)


def main() -> None:
    geo = search_geo()
    pride = search_pride()
    biostudies = search_biostudies()
    pubmed = search_pubmed()
    geo.to_csv(OUT / "GEO_MASLD_NAFLD_search.tsv", sep="\t", index=False)
    pride.to_csv(OUT / "PRIDE_MASLD_NAFLD_search.tsv", sep="\t", index=False)
    biostudies.to_csv(OUT / "BioStudies_MASLD_NAFLD_search.tsv", sep="\t", index=False)
    pubmed.to_csv(OUT / "PubMed_NIDDK_NASH_CRN_search.tsv", sep="\t", index=False)
    manifest = {
        "GEO_search_url": "https://www.ncbi.nlm.nih.gov/gds/?term=NAFLD+OR+NASH+OR+MASLD",
        "PRIDE_search_url": "https://www.ebi.ac.uk/pride/archive/projects?keyword=NAFLD",
        "BioStudies_search_url": "https://www.ebi.ac.uk/biostudies/",
        "NIDDK_repository_entry": "https://repository.niddk.nih.gov/",
        "UK_Biobank_application": "https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access",
        "All_of_Us_researcher_workbench": "https://researchallofus.org/",
        "dbGaP_search": "https://www.ncbi.nlm.nih.gov/gap/",
        "EGA_search": "https://ega-archive.org/",
    }
    (OUT / "public_cohort_search_links.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("GEO", len(geo), "PRIDE", len(pride), "BioStudies", len(biostudies), "PubMed", len(pubmed), flush=True)


if __name__ == "__main__":
    main()
