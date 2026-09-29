#!/usr/bin/env python
import os
from pathlib import Path
import gzip

PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
PLAT = PROJECT / "data_processed" / "bulk_geo" / "platforms"


def main() -> None:
    for path in sorted(PLAT.glob("*.gz")):
        print(f"\n{path.name}")
        with gzip.open(path, "rt", errors="replace") as handle:
            for index, line in enumerate(handle):
                if index >= 30:
                    break
                print(repr(line[:500]))


if __name__ == "__main__":
    main()
