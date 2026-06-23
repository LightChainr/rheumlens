#!/usr/bin/env python3
"""Create a deterministic raw-count sample for a valid Geneformer rerun."""

from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd

from fold_contained_benchmark import COHORTS, gene_names


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/geneformer_rerun"
N_PER_DONOR = 500
SEED = 20260619


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    src = Path(COHORTS["GSE135779"]["h5ad"])
    marker = np.load(ROOT / "data/secondary_features/GSE135779_project_marker_lineage.npy")
    a = ad.read_h5ad(src, backed="r")
    assert len(marker) == a.n_obs
    donor = a.obs["donor_id"].astype(str).to_numpy()
    rng = np.random.default_rng(SEED)
    selected = []
    for d in sorted(pd.unique(donor)):
        idx = np.flatnonzero(donor == d)
        selected.extend(np.sort(rng.choice(idx, min(N_PER_DONOR, len(idx)), replace=False)))
    selected = np.asarray(sorted(selected), dtype=np.int64)

    # Materialize only selected raw-count rows. Geneformer must not receive the
    # normalized/log-transformed X matrix.
    counts = a.layers["counts"][selected]
    if hasattr(counts, "to_memory"):
        counts = counts.to_memory()
    obs = a.obs.iloc[selected][["donor_id", "disease"]].copy()
    names = np.array(["T", "Monocyte", "B", "NK", "Unassigned"], dtype=object)
    code = marker[selected].astype(int)
    obs["marker_lineage"] = np.where(code >= 0, names[code], "Unassigned")
    obs["source_row_index"] = selected
    var = pd.DataFrame(index=pd.Index(gene_names(a).astype(str), name="gene_symbol"))
    out = ad.AnnData(X=counts, obs=obs, var=var)
    out.obs["n_counts"] = np.asarray(out.X.sum(axis=1)).ravel().astype(np.int64)
    out.write_h5ad(OUT / "GSE135779_geneformer_500_per_donor_raw_counts.h5ad",
                   compression="gzip")
    np.savez_compressed(
        OUT / "GSE135779_matched_scgpt_embeddings.npz",
        embeddings=np.asarray(a.obsm["X_scGPT"][selected], dtype=np.float32),
        donor_id=obs["donor_id"].astype(str).to_numpy(),
        disease=obs["disease"].astype(str).to_numpy(),
        marker_lineage=obs["marker_lineage"].astype(str).to_numpy(),
        source_row_index=selected,
    )
    np.save(OUT / "GSE135779_geneformer_selected_source_rows.npy", selected)
    metadata = {
        "source": str(src), "source_n_cells": int(a.n_obs),
        "selected_n_cells": int(len(selected)), "n_per_donor": N_PER_DONOR,
        "n_donors": int(pd.unique(donor[selected]).size), "seed": SEED,
        "matrix": "raw integer counts from layers['counts']",
        "marker_labels": "project fixed-marker score; not author annotations",
    }
    (OUT / "sampling_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    a.file.close()
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
