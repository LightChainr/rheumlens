"""Recompute the archived fixed pooling operators (mean / median / 10% trimmed
mean / mean+std) from cell-level float32 embeddings, for harness validation."""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from minipq import read_column

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from path_config import CELL_EMBEDDING_ROOT, LEARNED_POOLING_INPUT_ROOT

CE = str(CELL_EMBEDDING_ROOT)
OUT = str(LEARNED_POOLING_INPUT_ROOT)
os.makedirs(OUT, exist_ok=True)

def blocks_174188():
    shards = sorted(glob.glob(CE + "/SLE_GSE174188_CD4/shards/part-*.embeddings.npy"))
    seq = []
    for sp in shards:
        seq.extend(read_column(sp.replace(".embeddings.npy", ".metadata.parquet"), "donor_id"))
    donors = sorted(set(seq))
    pos = {d: i for i, d in enumerate(donors)}
    idx = np.asarray([pos[d] for d in seq], dtype=np.int32)
    X = np.concatenate([np.load(sp) for sp in shards], axis=0)[np.argsort(idx, kind="stable")]
    counts = np.bincount(idx, minlength=len(donors))
    return donors, counts, X

def blocks_285773():
    sh = sorted(glob.glob(CE + "/SLE_GSE285773_CD4/shards/donor-*.embeddings.npy"))
    donors = [os.path.basename(s).split(".")[0].split("-", 2)[2] for s in sh]
    mats = [np.load(s) for s in sh]
    return donors, np.array([len(m) for m in mats]), np.concatenate(mats, axis=0)

for tag, fn in (("GSE174188", blocks_174188), ("GSE285773", blocks_285773)):
    donors, counts, X = fn()
    starts = np.concatenate(([0], np.cumsum(counts)[:-1]))
    mean, med, trim, std = [], [], [], []
    for s, c in zip(starts, counts):
        v = X[s:s + c].astype(np.float64)
        mean.append(v.mean(axis=0))
        med.append(np.median(v, axis=0))
        t = int(np.floor(0.10 * len(v)))
        o = np.sort(v, axis=0)
        trim.append(o[t:len(v) - t].mean(axis=0) if t else o.mean(axis=0))
        std.append(v.std(axis=0, ddof=0))
    np.save(OUT + f"/{tag}_pool_mean.npy", np.asarray(mean))
    np.save(OUT + f"/{tag}_pool_median.npy", np.asarray(med))
    np.save(OUT + f"/{tag}_pool_trimmed10.npy", np.asarray(trim))
    np.save(OUT + f"/{tag}_pool_std.npy", np.asarray(std))
    pd.DataFrame({"donor_id": donors, "n_cells": counts}).to_csv(OUT + f"/{tag}_pool_index.csv", index=False)
    print(tag, "done", len(donors), X.shape, flush=True)
print("DONE", flush=True)
