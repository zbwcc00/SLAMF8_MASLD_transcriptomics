#!/usr/bin/env python
import os
from pathlib import Path
import requests

PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
OUT = PROJECT / "data_processed" / "bulk_geo" / "platforms"
PLATFORMS = ["GPL11532", "GPL570", "GPL21185"]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for gpl in PLATFORMS:
        candidates = [
            (f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl[:-3]}nnn/{gpl}/annot/{gpl}.annot.gz", f"{gpl}.annot.gz"),
            (f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl[:-3]}nnn/{gpl}/soft/{gpl}_family.soft.gz", f"{gpl}_family.soft.gz"),
        ]
        for url, name in candidates:
            target = OUT / name
            try:
                response = requests.get(url, stream=True, timeout=120)
                if response.status_code != 200:
                    continue
                with target.open("wb") as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
                print(f"saved {name}: {target.stat().st_size}", flush=True)
                break
            except Exception as error:
                print(f"failed {url}: {error}", flush=True)


if __name__ == "__main__":
    main()
