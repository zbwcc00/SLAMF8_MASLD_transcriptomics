# R environment

The frozen analysis used R 4.4.3. Exact package versions are listed in `R_packages.tsv`.

Install the matching Bioconductor release for R 4.4, then install the listed Bioconductor, CRAN and GitHub packages. Because several packages are distributed through GitHub, verify each installed version before rerunning the workflow. A version mismatch should be recorded rather than silently accepted.

Set the repository root before running R scripts:

```bash
export SLAMF8_MASLD_ROOT=/absolute/path/to/SLAMF8_MASLD_transcriptomics
```

On PowerShell:

```powershell
$env:SLAMF8_MASLD_ROOT = "D:\path\to\SLAMF8_MASLD_transcriptomics"
```
