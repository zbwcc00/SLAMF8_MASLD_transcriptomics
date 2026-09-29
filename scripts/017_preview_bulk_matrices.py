#!/usr/bin/env python
import os
from pathlib import Path
import gzip
import pandas as pd

PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
IN = PROJECT / "data_processed" / "bulk_geo"


def preview(accession: str) -> None:
    path = IN / f"{accession}_series_matrix.txt.gz"
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
        header = next(handle)
        values = [next(handle) for _ in range(3)]
    print(f"\n{accession}")
    print("header:", repr(header[:300]))
    for value in values:
        print("row:", repr(value[:300]))


def main() -> None:
    for accession in ["GSE159088", "GSE163211", "GSE182060", "GSE150734", "GSE48452", "GSE63067"]:
        preview(accession)


if __name__ == "__main__":
    main()
