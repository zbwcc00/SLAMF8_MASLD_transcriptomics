#!/usr/bin/env python3
from __future__ import annotations

import ast
import hashlib
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def canonical_bytes(path: Path) -> bytes:
    data = path.read_bytes()
    if path.suffix.lower() not in {".gz", ".xlsx", ".pdf", ".tiff", ".png"}:
        data = data.replace(b"\r\n", b"\n")
    return data


def sha256(data: bytes) -> str:
    digest = hashlib.sha256()
    digest.update(data)
    return digest.hexdigest()


def main() -> int:
    failures = []
    required = [
        "README.md", "CITATION.cff", "LICENSE", "LICENSE-DATA.md", ".zenodo.json",
        "environment/requirements.txt", "environment/R_packages.tsv",
        "data/metadata/cohorts.tsv", "data/reference/FARG95_gene_set_provenance.tsv",
        "MANIFEST_SHA256.tsv",
    ]
    for relative in required:
        if not (ROOT / relative).is_file():
            failures.append(f"missing: {relative}")

    forbidden = re.compile(r"(?:(?<![A-Za-z])[A-Za-z]:[\\/]|/Users/HUAWEI|第八篇大论文)")
    for path in sorted((ROOT / "scripts").glob("*")):
        if path.suffix.lower() not in {".py", ".r"}:
            continue
        text = path.read_text(encoding="utf-8")
        if forbidden.search(text):
            failures.append(f"workstation path remains: {path.relative_to(ROOT)}")
        if path.suffix.lower() == ".py":
            try:
                ast.parse(text, filename=str(path))
            except SyntaxError as error:
                failures.append(f"python syntax: {path.name}: {error}")

    farg = ROOT / "data/reference/FARG95_gene_set_provenance.tsv"
    if farg.exists():
        genes = [line.split("\t", 1)[0] for line in farg.read_text(encoding="utf-8-sig").splitlines()[1:] if line]
        if len(genes) != 95 or len(set(genes)) != 95 or "SLAMF8" not in genes:
            failures.append("FARG95 membership check failed")

    manifest = ROOT / "MANIFEST_SHA256.tsv"
    if manifest.exists():
        for line in manifest.read_text(encoding="utf-8").splitlines()[1:]:
            relative, size, digest = line.split("\t")
            path = ROOT / relative
            data = canonical_bytes(path) if path.is_file() else b""
            if not path.is_file() or len(data) != int(size) or sha256(data) != digest:
                failures.append(f"manifest mismatch: {relative}")

    print(f"repository={ROOT}")
    print(f"failures={len(failures)}")
    for failure in failures:
        print(f"FAIL: {failure}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
