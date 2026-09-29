#!/usr/bin/env python
import os
"""Download public MASLD/NAFLD expression supplements to the D-drive project."""

from pathlib import Path
import requests


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
DEST = ROOT / "data_processed" / "bulk_geo" / "expanded_public"
FILES = {
    "GSE167523_Raw_gene_counts_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167523/suppl/GSE167523_Raw_gene_counts_matrix.txt.gz",
    "GSE193066_NAFLD.HUn106.gct.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193066/suppl/GSE193066_NAFLD.HUn106.gct.gz",
    "GSE193066_NAFLD.HUn164.gct.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193066/suppl/GSE193066_NAFLD.HUn164.gct.gz",
    "GSE126848_Gene_counts_raw.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126848/suppl/GSE126848_Gene_counts_raw.txt.gz",
}


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    manifest = []
    for filename, url in FILES.items():
        path = DEST / filename
        response = requests.head(url, allow_redirects=True, timeout=60)
        response.raise_for_status()
        size = int(response.headers.get("Content-Length") or 0)
        if not path.exists() or (size and path.stat().st_size != size):
            with requests.get(url, stream=True, timeout=180) as download:
                download.raise_for_status()
                with path.open("wb") as handle:
                    for chunk in download.iter_content(1024 * 1024):
                        if chunk:
                            handle.write(chunk)
        manifest.append({"filename": filename, "url": url, "content_length": size, "local_size": path.stat().st_size})
        print(filename, path.stat().st_size, flush=True)
    import pandas as pd
    pd.DataFrame(manifest).to_csv(DEST / "download_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
