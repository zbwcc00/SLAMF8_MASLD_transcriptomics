#!/usr/bin/env python
import os
from pathlib import Path
import requests

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1])) / "data_processed" / "bulk_geo" / "full"
FILES = {
    "GSE162694": ["GSE162694_raw_counts.csv.gz"],
    "GSE130970": ["GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv.gz", "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz"],
    "GSE126848": ["GSE126848_Gene_counts_raw.txt.gz"],
    "GSE167523": ["GSE167523_Raw_gene_counts_matrix.txt.gz"],
    "GSE193066": ["GSE193066_NAFLD.HUn106.gct.gz", "GSE193066_NAFLD.HUn164.gct.gz"],
}


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    for accession, files in FILES.items():
        for filename in files:
            url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/suppl/{filename}"
            response = requests.head(url, timeout=30, allow_redirects=True)
            print(accession, filename, response.status_code, response.headers.get("Content-Length"), flush=True)
            if response.status_code != 200:
                continue
            if accession not in {"GSE162694", "GSE130970", "GSE126848"}:
                continue
            path = ROOT / filename
            size = int(response.headers.get("Content-Length") or 0)
            if path.exists() and path.stat().st_size == size and size > 0:
                continue
            with requests.get(url, stream=True, timeout=180) as download:
                download.raise_for_status()
                with path.open("wb") as handle:
                    for chunk in download.iter_content(1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            print("saved", filename, path.stat().st_size, flush=True)


if __name__ == "__main__":
    main()
