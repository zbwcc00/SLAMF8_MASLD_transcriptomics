#!/usr/bin/env python
import os
from pathlib import Path
import gzip
import csv
import pandas as pd

PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
IN = PROJECT / "data_processed" / "bulk_geo"
OUT = PROJECT / "results" / "bulk_validation"
ACCESSIONS = ["GSE159088", "GSE163211", "GSE182060", "GSE150734", "GSE48452", "GSE63067", "GSE89632", "GSE135251", "GSE126848", "GSE167523", "GSE193066", "GSE162694", "GSE130970"]


def inspect(accession: str) -> pd.DataFrame:
    path = IN / f"{accession}_series_matrix.txt.gz"
    rows = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                fields = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                key, values = fields[0], fields[1:]
                rows.setdefault(key, []).append(values)
            if line.startswith("!series_matrix_table_begin"):
                break
    n = max((len(v[0]) for v in rows.values() if v), default=0)
    data = {}
    for key, value_rows in rows.items():
        for row_number, values in enumerate(value_rows):
            column = key if row_number == 0 else f"{key}_{row_number + 1}"
            data[column] = values + [None] * (n - len(values))
    table = pd.DataFrame(data)
    table.insert(0, "dataset", accession)
    return table


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_tables = []
    for accession in ACCESSIONS:
        table = inspect(accession)
        table.to_csv(OUT / f"{accession}_sample_metadata.tsv", sep="\t", index=False)
        all_tables.append(table)
        print(f"{accession}: {table.shape}")
        for col in ["!Sample_title", "!Sample_source_name_ch1", "!Sample_characteristics_ch1"]:
            if col in table:
                print(f"--- {col}")
                print(table[col].head(12).to_string(index=False))
    pd.concat(all_tables, ignore_index=True).to_csv(OUT / "bulk_geo_sample_metadata_all.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
