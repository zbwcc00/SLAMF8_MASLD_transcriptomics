#!/usr/bin/env python
import os
from pathlib import Path
import gzip
import pandas as pd

PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
PLAT = PROJECT / "data_processed" / "bulk_geo" / "platforms"


def main() -> None:
    for name in ["GPL11532.annot.gz", "GPL570.annot.gz"]:
        path = PLAT / name
        marker = None
        with gzip.open(path, "rt", errors="replace") as handle:
            for i, line in enumerate(handle):
                if line.startswith("!platform_table_begin"):
                    marker = i
                    break
        table = pd.read_csv(path, sep="\t", compression="gzip", skiprows=marker + 1, dtype=str)
        symbol = table["Gene symbol"]
        mapped = symbol.notna() & symbol.ne("")
        print(name, table.shape, "mapped", int(mapped.sum()))
        print(table.loc[mapped, ["ID", "Gene symbol"]].head(10).to_string(index=False))
        print("SLAMF8", table.loc[symbol.fillna("").str.contains("SLAMF8", case=False), ["ID", "Gene symbol"]].to_string(index=False))


if __name__ == "__main__":
    main()
