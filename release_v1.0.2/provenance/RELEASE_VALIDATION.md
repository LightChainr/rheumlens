# Version 1.0.2 release validation

Validation date: 19 July 2026

## Passed checks

- Repository test suite: 14 passed under Python 3.12.
- Runtime and project version consistency: covered by `tests/test_version.py`.
- Changed Python files: Ruff checks passed.
- Source and wheel build: `uv build` produced `rheumlens-1.0.2.tar.gz` and `rheumlens-1.0.2-py3-none-any.whl`.
- Zenodo JSON: valid JSON.
- CITATION.cff: valid YAML.
- Main manuscript DOCX: ZIP/OOXML integrity passed; 27-page Letter PDF rendered successfully.
- Cover-letter DOCX: ZIP/OOXML integrity passed; one-page PDF rendered successfully.
- Supplementary DOCX: ZIP/OOXML integrity passed.
- Supplementary workbook: 25 sheets; Tables S19 and S20 present and readable.
- Main figures: eight PNG files, each at least 2,327 pixels on its shorter reported axis.
- Stable DOI substitution: release text and rendered documents reference concept DOI `10.5281/zenodo.20813922`; the previous v1.0.1 DOI is retained only as historical metadata.
- Release and repository SHA256 manifests: generated and verified.

## Non-blocking historical lint debt

A full-repository Ruff scan reports 19 pre-existing findings in historical package modules and tests. None is introduced by the v1.0.2 release files or version-consistency test; the changed-file Ruff gate passes. These style findings do not affect the release analyses or package execution.

