"""Convert a CELLxGENE .h5ad cohort into the donor-level interface used by the audit.

Contract produced for every cohort (identical to the locked v2 inputs):

    inputs/donor_level/<COHORT>/donor_labels.tsv
    inputs/donor_level/<COHORT>/donor_log1p_cpm.parquet     index = donor_id
    inputs/design_metadata/<cohort>_donor_covariates.tsv

Why this script exists
----------------------
The v2 design block for GSE174188 was

    qc          = [log_cells_per_donor, log_mean_umi_per_cell,
                   log_mean_genes_per_cell, aggregate_pct_mito]
    demographic = [age_years, sex]

Every one of those is recomputable from the count matrix plus CELLxGENE's
guaranteed obs fields. No dataset-specific supplementary table is required.
That is what makes a multi-disease screen cheap.

Memory
------
The largest registered cohort is ~14 GB on disk. The matrix is never fully
materialised: cells are streamed in row chunks and accumulated into per-donor
gene sums (n_donor x n_gene float64, ~50 MB for 200 donors x 30k genes).

Usage
-----
    python scripts/cohorts/build_donor_level.py --cohort COVID_REN
    python scripts/cohorts/build_donor_level.py --cohort COMBAT --contrast COMBAT_CROSS
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

REGISTRY = REPO / "cohorts" / "registry.yaml"
CHUNK = 100_000

# CELLxGENE guarantees these obs columns (schema >= 3.0).
REQUIRED_OBS = ["donor_id", "disease", "sex", "development_stage", "assay", "tissue"]

# obs columns that, if present, are treated as candidate acquisition/batch variables.
BATCH_CANDIDATE_PATTERNS = [
    r"^batch", r"batch$", r"^site", r"site$", r"centre", r"center",
    r"^pool", r"pool$", r"^library", r"library_id", r"^sample_id$",
    r"^lane", r"^run", r"^chemistry", r"^processing", r"^wave", r"^timepoint",
    r"^collection", r"^institute", r"^cohort$", r"^sequencing",
]


def parse_age(stage: str) -> float:
    """Map a CELLxGENE development_stage label to an approximate numeric age.

    Handles the three shapes that actually occur in the registered cohorts:
        '25-year-old stage'          -> 25
        'fifth decade stage'         -> 45   (decade N -> 10N - 5)
        '90 year-old and over stage' -> 90
    Returns NaN when no age can be recovered (e.g. 'human adult stage').
    """
    if not isinstance(stage, str):
        return np.nan
    s = stage.lower()

    m = re.search(r"(\d+)\s*-?\s*year-?old", s)
    if m:
        return float(m.group(1))
    m = re.search(r"(\d+)\s*year-old and over", s)
    if m:
        return float(m.group(1))

    decades = {
        "first": 1, "second": 2, "third": 3, "fourth": 4, "fifth": 5,
        "sixth": 6, "seventh": 7, "eighth": 8, "ninth": 9, "tenth": 10,
    }
    m = re.search(r"(\w+) decade", s)
    if m and m.group(1) in decades:
        # "fifth decade" = ages 40-49 -> midpoint 45
        return decades[m.group(1)] * 10 - 5
    return np.nan


def pick_raw_layer(adata: ad.AnnData) -> str:
    """Decide where the raw integer counts live.

    CELLxGENE schema puts normalised values in .X and raw counts in .raw.X for
    most datasets, but some deposit raw counts directly in .X. Guess from a
    small sample: raw counts are non-negative integers.
    """
    def looks_like_counts(mat) -> bool:
        sub = mat[: min(2000, mat.shape[0])]
        sub = sub.data if sp.issparse(sub) else np.asarray(sub).ravel()
        if sub.size == 0:
            return False
        sub = sub[:100_000]
        return bool(np.all(sub >= 0) and np.allclose(sub, np.round(sub)))

    if adata.raw is not None and looks_like_counts(adata.raw.X):
        return "raw"
    if looks_like_counts(adata.X):
        return "X"
    if adata.raw is not None:
        print("  ! neither .X nor .raw.X look like integer counts; using .raw.X")
        return "raw"
    raise RuntimeError("no raw count matrix found (.X is not counts and .raw is absent)")


def detect_batch_columns(obs: pd.DataFrame, donor_col: str = "donor_id") -> list[str]:
    """Find obs columns that vary between donors and plausibly encode acquisition."""
    found = []
    for col in obs.columns:
        if col in {donor_col, "disease", "cell_type", "tissue"}:
            continue
        if not any(re.search(p, col.lower()) for p in BATCH_CANDIDATE_PATTERNS):
            continue
        nun = obs[col].nunique(dropna=True)
        if not (2 <= nun <= max(40, obs[donor_col].nunique() // 2)):
            continue
        # keep only columns that are (near-)constant within a donor
        per_donor = obs.groupby(donor_col, observed=True)[col].nunique(dropna=True)
        if (per_donor <= 1).mean() >= 0.95:
            found.append(col)
    return found


def build(cohort_cfg: dict, defaults: dict, h5ad: Path, out_name: str,
          case_labels: list[str], control_labels: list[str],
          repo_inputs: Path) -> dict:
    print(f"[{out_name}] opening {h5ad.name}")
    adata = ad.read_h5ad(h5ad, backed="r")
    obs = adata.obs.copy()
    print(f"  cells={adata.n_obs:,}  genes={adata.n_vars:,}")

    missing = [c for c in REQUIRED_OBS if c not in obs.columns]
    if missing:
        raise RuntimeError(f"obs is missing required CELLxGENE columns: {missing}")

    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["disease"] = obs["disease"].astype(str)

    # --- cell filter: tissue, then disease groups -------------------------------
    keep = pd.Series(True, index=obs.index)
    tissue_filter = defaults.get("tissue_filter")
    if tissue_filter:
        tis = obs["tissue"].astype(str).str.lower()
        if tis.str.contains(tissue_filter).any():
            keep &= tis.str.contains(tissue_filter)
            print(f"  tissue filter '{tissue_filter}': {int(keep.sum()):,} cells kept")

    wanted = set(case_labels) | set(control_labels)
    keep &= obs["disease"].isin(wanted)
    print(f"  disease filter {sorted(wanted)}: {int(keep.sum()):,} cells kept")
    if keep.sum() == 0:
        raise RuntimeError("no cells left after filtering; check case/control labels")

    obs = obs.loc[keep]
    row_idx = np.flatnonzero(keep.to_numpy())

    # --- donor label table ------------------------------------------------------
    donor_disease = obs.groupby("donor_id", observed=True)["disease"].agg(
        lambda s: s.value_counts().idxmax()
    )
    y = donor_disease.isin(case_labels).astype(int)

    # --- stream the matrix, accumulate per-donor sums ---------------------------
    layer = pick_raw_layer(adata)
    print(f"  raw counts source: .{layer}{'.X' if layer == 'raw' else ''}")
    mat = adata.raw.X if layer == "raw" else adata.X
    var_names = (adata.raw.var_names if layer == "raw" else adata.var_names).astype(str)
    n_gene = len(var_names)

    donors = donor_disease.index.tolist()
    d_pos = {d: i for i, d in enumerate(donors)}
    codes = obs["donor_id"].map(d_pos).to_numpy()

    gene_sum = np.zeros((len(donors), n_gene), dtype=np.float64)
    n_cells = np.zeros(len(donors), dtype=np.int64)
    umi_sum = np.zeros(len(donors), dtype=np.float64)
    gene_det = np.zeros(len(donors), dtype=np.float64)
    mito_sum = np.zeros(len(donors), dtype=np.float64)

    # mitochondrial genes: prefer symbols in var, fall back to feature_name
    var_df = adata.raw.var if layer == "raw" else adata.var
    sym = var_df["feature_name"].astype(str) if "feature_name" in var_df.columns \
        else pd.Series(var_names, index=var_names)
    mito_mask = sym.str.upper().str.startswith(("MT-", "MT.")).to_numpy()
    print(f"  mitochondrial genes detected: {int(mito_mask.sum())}")

    for start in range(0, len(row_idx), CHUNK):
        sel = row_idx[start:start + CHUNK]
        block = mat[sel]
        block = block if sp.issparse(block) else sp.csr_matrix(np.asarray(block))
        block = block.tocsr()
        c = codes[start:start + CHUNK]

        per_cell_umi = np.asarray(block.sum(axis=1)).ravel()
        per_cell_gene = np.asarray((block > 0).sum(axis=1)).ravel()
        per_cell_mito = np.asarray(block[:, mito_mask].sum(axis=1)).ravel() \
            if mito_mask.any() else np.zeros(block.shape[0])

        for gi in np.unique(c):
            m = c == gi
            gene_sum[gi] += np.asarray(block[m].sum(axis=0)).ravel()
            n_cells[gi] += int(m.sum())
            umi_sum[gi] += per_cell_umi[m].sum()
            gene_det[gi] += per_cell_gene[m].sum()
            mito_sum[gi] += per_cell_mito[m].sum()

        print(f"  streamed {min(start + CHUNK, len(row_idx)):,}/{len(row_idx):,} cells",
              end="\r", flush=True)
    print()

    # --- drop under-sampled donors ---------------------------------------------
    min_cells = int(defaults.get("min_cells_per_donor", 200))
    ok = n_cells >= min_cells
    if not ok.all():
        dropped = [donors[i] for i in np.flatnonzero(~ok)]
        print(f"  dropping {len(dropped)} donors with <{min_cells} cells")
    donors = [d for d, k in zip(donors, ok) if k]
    gene_sum, n_cells = gene_sum[ok], n_cells[ok]
    umi_sum, gene_det, mito_sum = umi_sum[ok], gene_det[ok], mito_sum[ok]
    y = y.loc[donors]

    if y.nunique() < 2:
        raise RuntimeError("only one class survives filtering; contrast is degenerate")

    # --- pseudobulk: log1p CPM --------------------------------------------------
    totals = gene_sum.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    pseudobulk = np.log1p(gene_sum / totals * 1e6).astype(np.float32)
    pb = pd.DataFrame(pseudobulk, index=pd.Index(donors, name="donor_id"),
                      columns=var_names)

    # --- design covariates ------------------------------------------------------
    cov = pd.DataFrame(index=pd.Index(donors, name="donor_id"))
    cov["y_true"] = y.to_numpy()
    cov["disease_label"] = donor_disease.loc[donors].to_numpy()
    cov["cells_per_donor"] = n_cells
    cov["mean_umi_per_cell"] = umi_sum / n_cells
    cov["mean_genes_per_cell"] = gene_det / n_cells
    cov["aggregate_pct_mito"] = np.where(umi_sum > 0, mito_sum / umi_sum * 100.0, np.nan)
    cov["log_cells_per_donor"] = np.log1p(cov["cells_per_donor"])
    cov["log_mean_umi_per_cell"] = np.log1p(cov["mean_umi_per_cell"])
    cov["log_mean_genes_per_cell"] = np.log1p(cov["mean_genes_per_cell"])

    first = obs.groupby("donor_id", observed=True).first()
    cov["sex"] = first.loc[donors, "sex"].astype(str).to_numpy()
    cov["assay"] = first.loc[donors, "assay"].astype(str).to_numpy()
    cov["age_years"] = [parse_age(s) for s in first.loc[donors, "development_stage"]]
    if "self_reported_ethnicity" in first.columns:
        cov["ethnicity"] = first.loc[donors, "self_reported_ethnicity"].astype(str).to_numpy()
    if "suspension_type" in first.columns:
        cov["suspension_type"] = first.loc[donors, "suspension_type"].astype(str).to_numpy()

    batch_cols = detect_batch_columns(obs)
    for col in batch_cols:
        cov[f"batch__{col}"] = first.loc[donors, col].astype(str).to_numpy()
    print(f"  batch-like obs columns registered: {batch_cols or 'NONE'}")

    n_age = int(cov["age_years"].notna().sum())
    print(f"  donors={len(donors)}  cases={int(cov.y_true.sum())} "
          f"controls={int((cov.y_true == 0).sum())}  age recovered for {n_age}")

    # --- write ------------------------------------------------------------------
    dl = repo_inputs / "donor_level" / out_name
    dm = repo_inputs / "design_metadata"
    dl.mkdir(parents=True, exist_ok=True)
    dm.mkdir(parents=True, exist_ok=True)

    pb.to_parquet(dl / "donor_log1p_cpm.parquet")
    cov[["y_true"]].reset_index().to_csv(dl / "donor_labels.tsv", sep="\t", index=False)
    cov.reset_index().to_csv(dm / f"{out_name.lower()}_donor_covariates.tsv",
                             sep="\t", index=False)

    summary = {
        "cohort": out_name,
        "dataset_id": cohort_cfg["dataset_id"],
        "citation": cohort_cfg["citation"],
        "n_donor": len(donors),
        "n_case": int(cov.y_true.sum()),
        "n_control": int((cov.y_true == 0).sum()),
        "n_cell_used": int(n_cells.sum()),
        "n_gene": int(n_gene),
        "case_labels": case_labels,
        "control_labels": control_labels,
        "raw_source": layer,
        "batch_columns": batch_cols,
        "age_recovered": n_age,
        "expected_confounding": cohort_cfg.get("expected_confounding"),
    }
    (dl / "build_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  wrote {dl}")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, help="registry cohort name")
    ap.add_argument("--contrast", default=None,
                    help="named sub-contrast, e.g. COMBAT_CROSS")
    ap.add_argument("--h5ad", default=None, help="override local h5ad path")
    ap.add_argument("--data-dir", default=None,
                    help="directory holding downloaded h5ad files")
    args = ap.parse_args()

    reg = yaml.safe_load(REGISTRY.read_text())
    defaults = reg["defaults"]
    cfg = next((c for c in reg["cohorts"] if c["name"] == args.cohort), None)
    if cfg is None:
        raise SystemExit(f"unknown cohort {args.cohort}; "
                         f"known: {[c['name'] for c in reg['cohorts']]}")

    data_dir = Path(args.data_dir) if args.data_dir else REPO / "data" / "h5ad"
    h5ad = Path(args.h5ad) if args.h5ad else data_dir / f"{cfg['name']}.h5ad"
    if not h5ad.exists():
        raise SystemExit(f"{h5ad} not found - run scripts/cohorts/fetch_cohort.py first")

    control = [defaults["control_label"]]
    case = list(cfg["case_labels"])
    out = cfg["name"]

    # COMBAT ships three contrasts out of one file.
    if args.contrast:
        if args.contrast == "COMBAT_COVID":
            case, out = ["COVID-19"], "COMBAT_COVID"
        elif args.contrast == "COMBAT_INFLUENZA":
            case, out = ["influenza"], "COMBAT_INFLUENZA"
        elif args.contrast == "COMBAT_CROSS":
            case, control, out = ["COVID-19"], ["influenza"], "COMBAT_CROSS"
        else:
            raise SystemExit(f"unknown contrast {args.contrast}")

    summary = build(cfg, defaults, h5ad, out, case, control, REPO / "inputs")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
