#!/usr/bin/env python
import os
"""Download additional GEO cohorts prioritized by clinical relevance."""

from pathlib import Path
import requests
import pandas as pd


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
DEST = ROOT / "data_processed" / "bulk_geo" / "expanded_public"
FILES = {
    "GSE83452_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83452/matrix/GSE83452_series_matrix.txt.gz",
    "GSE281797_GeneCount.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281797/suppl/GSE281797_GeneCount.txt.gz",
    "GSE281797_tpm.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281797/suppl/GSE281797_tpm.txt.gz",
    "GSE319035_NAFLD_time_course.counts.20250908.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319035/suppl/GSE319035_NAFLD_time_course.counts.20250908.txt.gz",
    "GSE319035_NAFLD_time_course_sample_key.20220408.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319035/suppl/GSE319035_NAFLD_time_course_sample_key.20220408.txt.gz",
}


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    rows = []
    for filename, url in FILES.items():
        path = DEST / filename
        head = requests.head(url, allow_redirects=True, timeout=60)
        head.raise_for_status()
        expected = int(head.headers.get("Content-Length") or 0)
        if not path.exists() or (expected and path.stat().st_size != expected):
            with requests.get(url, stream=True, timeout=300) as response:
                response.raise_for_status()
                with path.open("wb") as handle:
                    for chunk in response.iter_content(1024 * 1024):
                        if chunk:
                            handle.write(chunk)
        rows.append({"filename": filename, "url": url, "content_length": expected, "local_size": path.stat().st_size})
        print(filename, path.stat().st_size, flush=True)
    pd.DataFrame(rows).to_csv(DEST / "additional_download_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
