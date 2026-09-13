#!/usr/bin/env python3
"""Retrospective count of CELLxGENE Discover datasets eligible for the screen.

The five screen datasets were chosen purposively to span acquisition structures,
not sampled from every eligible public cohort. This script asks, after the fact,
how many datasets met the screen's dataset-level entry requirements on the date
the registry was verified, so that a reader can see how large the pool was and
where the five sit in it. It was not used to select anything.

Dataset-level criteria (all must hold):
  * organism Homo sapiens;
  * at least one tissue label containing "blood";
  * single-cell suspension (not nuclei-only, not spatial);
  * disease labels include "normal" and at least one other disease;
  * at least 40 distinct donor identifiers (the screen's power floor);
  * published on or before the registry verification date.

Donor counts per class are not in the dataset-level metadata, so a dataset that
passes may still fail the screen's per-class floors once donors are counted.
Datasets are collapsed to collections, because one collection often releases the
same donors as several cell-type subsets.

    python3 audit/cellxgene_eligibility_audit.py --json cxg_datasets.json \
        --registry cohorts/registry.yaml --out supplementary/Table_S9_eligibility_audit.tsv
"""
from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

import pandas as pd
import yaml

API = "https://api.cellxgene.cziscience.com/curation/v1/datasets"
CUTOFF = "2026-09-07"
MIN_DONORS = 40


def labels(field) -> list[str]:
    return [x["label"] for x in (field or [])]


def eligible(d: dict) -> bool:
    if "Homo sapiens" not in labels(d.get("organism")):
        return False
    if not any("blood" in t.lower() for t in labels(d.get("tissue"))):
        return False
    if "cell" not in (d.get("suspension_type") or []) or d.get("spatial"):
        return False
    dis = labels(d.get("disease"))
    if "normal" not in dis or len([x for x in dis if x != "normal"]) == 0:
        return False
    if len(set(d.get("donor_id") or [])) < MIN_DONORS:
        return False
    return (d.get("published_at") or "9999")[:10] <= CUTOFF


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", help="saved API response; fetched if omitted")
    ap.add_argument("--registry", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.json and Path(a.json).exists():
        data = json.loads(Path(a.json).read_text())
    else:
        with urllib.request.urlopen(API) as r:
            data = json.load(r)

    reg = yaml.safe_load(Path(a.registry).read_text())
    analysed = {c["dataset_id"]: c["name"] for c in reg["cohorts"] if c.get("dataset_id")}
    analysed_collections = {d["collection_id"] for d in data if d["dataset_id"] in analysed}

    rows = {}
    for d in data:
        if not eligible(d):
            continue
        c = rows.setdefault(d["collection_id"], {
            "collection_name": d.get("collection_name"),
            "collection_doi": d.get("collection_doi") or "",
            "n_datasets_eligible": 0, "max_donors": 0, "diseases": set(),
            "tissues": set(), "assays": set(),
            "analysed_in_screen": d["collection_id"] in analysed_collections})
        c["n_datasets_eligible"] += 1
        c["max_donors"] = max(c["max_donors"], len(set(d["donor_id"])))
        c["diseases"] |= {x for x in labels(d["disease"]) if x != "normal"}
        c["tissues"] |= set(labels(d["tissue"]))
        c["assays"] |= set(labels(d["assay"]))

    df = pd.DataFrame([{**v, "collection_id": k,
                        "diseases": "; ".join(sorted(v["diseases"])),
                        "tissues": "; ".join(sorted(v["tissues"])),
                        "assays": "; ".join(sorted(v["assays"]))}
                       for k, v in rows.items()])
    df = df.sort_values(["analysed_in_screen", "max_donors"], ascending=[False, False])
    cols = ["collection_id", "collection_name", "collection_doi", "analysed_in_screen",
            "n_datasets_eligible", "max_donors", "diseases", "tissues", "assays"]
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(f"# Table S9. Retrospective eligibility audit of CELLxGENE Discover, "
                 f"datasets published on or before {CUTOFF}.\n"
                 f"# {len(data)} datasets in the API response; "
                 f"{sum(v['n_datasets_eligible'] for v in rows.values())} datasets in "
                 f"{len(df)} collections meet the dataset-level criteria described in Section 9.1. "
                 f"Performed after the analysis; "
                 f"not used for selection.\n")
        df[cols].to_csv(fh, sep="\t", index=False)
    n_an = int(df.analysed_in_screen.sum())
    print(f"api={len(data)} eligible_datasets={sum(v['n_datasets_eligible'] for v in rows.values())} "
          f"eligible_collections={len(df)} analysed_collections_eligible={n_an} "
          f"analysed_collections_total={len(analysed_collections)}")
    print(df[cols[:6]].to_string(index=False, max_colwidth=60))
    # which analysed collections failed a criterion, and why
    for cid in sorted(analysed_collections - set(rows)):
        ds = [d for d in data if d["collection_id"] == cid and d["dataset_id"] in analysed]
        for d in ds:
            print("analysed but not eligible:", analysed[d["dataset_id"]],
                  labels(d["tissue"]), labels(d["disease"]), len(set(d["donor_id"])),
                  d.get("suspension_type"), (d.get("published_at") or "")[:10])


if __name__ == "__main__":
    main()
