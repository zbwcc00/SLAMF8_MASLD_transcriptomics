#!/usr/bin/env python
import os
from pathlib import Path
import requests

PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
OUT = PROJECT / "data_processed" / "bulk_geo"
ACCESSIONS = ["GSE159088", "GSE163211", "GSE182060", "GSE150734", "GSE48452", "GSE63067", "GSE89632", "GSE135251", "GSE126848", "GSE167523", "GSE193066", "GSE162694", "GSE130970"]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for accession in ACCESSIONS:
        path = OUT / f"{accession}_series_matrix.txt.gz"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/matrix/{accession}_series_matrix.txt.gz"
        if path.exists() and path.stat().st_size > 1000:
            print(f"exists {accession}", flush=True)
            continue
        print(f"downloading {accession}", flush=True)
        with requests.get(url, stream=True, timeout=120) as response:
            response.raise_for_status()
            with path.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        print(f"saved {accession}: {path.stat().st_size} bytes", flush=True)


if __name__ == "__main__":
    main()
