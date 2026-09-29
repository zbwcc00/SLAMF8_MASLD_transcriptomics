import os
from pathlib import Path
import gzip
import pandas as pd

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1])) / "data_processed" / "bulk_geo" / "full"
META = ROOT.parent.parent.parent / "results" / "bulk_validation"


def main():
    for filename in ["GSE162694_raw_counts.csv.gz", "GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv.gz", "GSE126848_Gene_counts_raw.txt.gz"]:
        accession = filename[:9]
        print("\n", filename)
        with gzip.open(ROOT / filename, "rt", errors="replace") as handle:
            for _ in range(3):
                print(repr(next(handle)[:650]))
        table = pd.read_csv(META / f"{accession}_sample_metadata.tsv", sep="\t")
        print(table[["!Sample_geo_accession", "!Sample_title"]].head(4).to_string(index=False))


if __name__ == "__main__":
    main()
