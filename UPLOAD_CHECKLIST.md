# GitHub and Zenodo upload checklist

## GitHub

1. Create a public repository named `SLAMF8_MASLD_transcriptomics`.
2. Upload the contents of this directory, preserving hidden files such as `.github`, `.gitignore`, `.gitattributes` and `.zenodo.json`.
3. Run `python tools/validate_repository.py` locally or wait for the GitHub Actions validation workflow.
4. Create release `v1.0.0` from the exact uploaded commit.

## Zenodo

1. Sign in to Zenodo and optionally connect the GitHub repository.
2. Create a software deposit or archive the GitHub `v1.0.0` release.
3. Use the metadata in `.zenodo.json`; select open access and confirm the mixed-licence statement (MIT code; CC BY 4.0 author-generated processed data/documentation; no redistributed third-party raw data).
4. Upload `Zenodo_Deposit_SLAMF8_MASLD_v1.0.0.zip` when creating the deposit manually.
5. Reserved version DOI: `10.5281/zenodo.23042100`.
6. The DOI and GitHub release URL have been synchronized across the repository and manuscript files.
7. Replace the draft upload with the final DOI-synchronized archive, then publish the record.
8. Test the public landing page and download link outside the author account.

Do not publish until the authors approve the licensing fields and the final title/author metadata.
