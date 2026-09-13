#!/usr/bin/env python3
"""Generate the cohort registry (Table S7) from the registry file and the donor tables.

The released Table S7 had been maintained by hand, had drifted out of column
alignment (donor counts sitting under `collection`, citations under `case_labels`),
covered only the five main-screen datasets while the text described seven, and left
several declared fields empty. Every field is now derived:

  * dataset identifier, citation, API donor count and download size come from
    cohorts/registry.yaml, which records the CELLxGENE Discover API response;
  * the analysis donor, case and control counts are counted from the donor label
    tables the models were actually fitted on;
  * the two datasets used only in Section 5 are included and marked as such.

The API count and the analysis count differ - donors below the minimum cell count,
or carrying a label outside the contrast, are dropped - so both are reported.
"""
from __future__ import annotations
import os, re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPO = Path(os.environ.get(
    "RHEUMLENS_REPO", ROOT.parent / "20_repo_restructure_20260907"))
REG = REPO / "cohorts" / "registry.yaml"
DONORS = REPO / "inputs" / "donor_level"
OUT = ROOT / "supplementary" / "Table_S7_cohort_registry.tsv"

# Which donor-level directories belong to which released dataset, and which of
# them the screen forms comparisons from.
DATASETS = {
    "SLE_GSE174188":    dict(role="main screen", dirs=["SLE_GSE174188_CD4",
                                                       "SLE_GSE174188_ALLCELL"]),
    "COVID_REN":        dict(role="main screen", dirs=["COVID_REN",
                                                       "COVID_REN_ASSAY_10x_5_v2"]),
    "COVID_STEPHENSON": dict(role="main screen", dirs=["COVID_STEPHENSON"]),
    "COMBAT":           dict(role="main screen", dirs=["COMBAT_COVID",
                                                       "COMBAT_INFLUENZA",
                                                       "COMBAT_CROSS"]),
    "CMV_HIHA":         dict(role="main screen", dirs=["CMV_HIHA"]),
    "SLE_GSE135779":    dict(role="Section 5 only", dirs=["SLE_GSE135779"]),
    "SLE_GSE285773":    dict(role="Section 5 only", dirs=["SLE_GSE285773_CD4"]),
}
# The two Section 5 datasets are not in the CELLxGENE registry: they were taken
# from GEO. Their accessions and citations are declared here, and their donor
# counts are still counted from the donor tables like every other row.
EXTRA = {
    "SLE_GSE135779": dict(dataset_id="GEO GSE135779",
                          citation="Nehar-Belaid et al. (2020) Nature Immunology",
                          case_labels="SLE", expected_confounding="not screened"),
    "SLE_GSE285773": dict(dataset_id="GEO GSE285773",
                          citation="NCBI GEO GSE285773 (public 2025-09-09)",
                          case_labels="SLE", expected_confounding="not screened"),
}


def parse_registry(text: str) -> dict[str, dict[str, str]]:
    """Read the flat scalar fields of each cohort entry. No YAML dependency."""
    out, cur = {}, None
    for line in text.splitlines():
        m = re.match(r"\s*-\s+name:\s*(\S+)", line)
        if m:
            cur = m.group(1)
            out[cur] = {}
            continue
        if cur is None:
            continue
        m = re.match(r"\s{4}([a-z_0-9]+):\s*(.+?)\s*$", line)
        if m and not m.group(2).startswith(">"):
            out[cur][m.group(1)] = m.group(2).strip("[]")
    return out


reg = parse_registry(REG.read_text())
verified = re.search(r"API on[\s#]+(\d{4}-\d{2}-\d{2})", REG.read_text())
verified = verified.group(1) if verified else "unrecorded"

rows = []
for name, spec in DATASETS.items():
    e = {**reg.get(name, {}), **EXTRA.get(name, {})}
    counts = []
    for d in spec["dirs"]:
        f = DONORS / d / "donor_labels.tsv"
        if not f.exists():
            raise SystemExit(f"missing donor table: {f}")
        t = pd.read_csv(f, sep="\t")
        # Two encodings are in use across the donor tables: a 0/1 `y_true` and a
        # case/control string. Normalise rather than assume either.
        col = "case_control" if "case_control" in t.columns else "y_true"
        v = t[col]
        pos = (v == 1).sum() if v.dtype != object else (v.astype(str) == "case").sum()
        neg = len(t) - pos
        counts.append((d, len(t), int(pos), int(neg)))
    from_api = name in reg
    rows.append({
        "dataset": name,
        "role": spec["role"],
        "dataset_id": e.get("dataset_id", ""),
        "citation": e.get("citation", ""),
        "disease_group": e.get("disease_group", "SLE"),
        "expected_confounding": e.get("expected_confounding", ""),
        "n_donor_released": e.get("n_donor", str(sum(n for _, n, _, _ in counts))),
        "download_gb": e.get("h5ad_gb", "n/a (GEO, not an h5ad download)"),
        "case_labels": e.get("case_labels", ""),
        "control_label": "normal" if from_api else "healthy donor",
        "analysis_units": "; ".join(
            f"{d} (n={n}, {c} case / {k} control)" for d, n, c, k in counts),
        "n_comparisons_in_screen": len(spec["dirs"]) if spec["role"] == "main screen" else 0,
        "source": ("CZ CELLxGENE Discover curation API" if from_api
                   else "NCBI GEO series metadata"),
        "verification_date": verified if from_api else "not API-verified (GEO)",
    })

df = pd.DataFrame(rows)
assert df.notna().all().all() and (df.astype(str) != "").all().all(), \
    "every field must be filled:\n" + df.astype(str).to_string()

with OUT.open("w") as fh:
    fh.write("# Cohort registry. 'n_donor_released' is the donor count reported by the\n")
    fh.write("# source for the whole dataset; 'analysis_units' gives the donor, case and\n")
    fh.write("# control counts actually modelled, after donors below the minimum cell\n")
    fh.write("# count and donors outside the contrast were dropped.\n")
    df.to_csv(fh, sep="\t", index=False)
print(f"wrote {OUT.relative_to(ROOT)}: {len(df)} datasets, "
      f"{df.n_comparisons_in_screen.sum()} screened comparisons")
