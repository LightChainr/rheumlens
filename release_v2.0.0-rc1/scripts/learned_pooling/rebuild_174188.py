import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from minipq import read_column

P = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
W = P + "/08_Applied_Sciences_special_issue_submission_20260717/server_results/extracted/patient_generalization_work"
CE = W + "/outputs/cell_embeddings_cap500/SLE_GSE174188_CD4"
OUT = P + "/14_learned_pooling_20260726/inputs"

shards = sorted(glob.glob(CE + "/shards/part-*.embeddings.npy"))
donor_seq = []
for sp in shards:
    donor_seq.extend(read_column(sp.replace(".embeddings.npy", ".metadata.parquet"), "donor_id"))
print("cells", len(donor_seq), "donors", len(set(donor_seq)), flush=True)

donors = sorted(set(donor_seq))
pos = {d: i for i, d in enumerate(donors)}
cell_idx = np.asarray([pos[d] for d in donor_seq], dtype=np.int32)
counts = np.bincount(cell_idx, minlength=len(donors))
cov = pd.read_csv(P + "/05_patient_classifier_generalization_20260716/results/01_inputs/gse174188_final_donor_covariates.tsv", sep="\t")
cov["donor_id"] = cov["donor_id"].astype(str)
exp = cov.set_index("donor_id").loc[donors, "cells_per_donor"].clip(upper=500).to_numpy()
print("counts match expected:", bool((counts == exp).all()), flush=True)

X = np.concatenate([np.load(sp) for sp in shards], axis=0)
print("loaded", X.shape, X.dtype, flush=True)
X = X[np.argsort(cell_idx, kind="stable")]
starts = np.concatenate(([0], np.cumsum(counts)[:-1]))
mean = np.add.reduceat(X.astype(np.float64), starts, axis=0) / counts[:, None]
np.save(OUT + "/GSE174188_donor_mean_f64_v2.npy", mean)
pd.DataFrame({"donor_id": donors, "n_cells": counts}).to_csv(OUT + "/GSE174188_donor_index_v2.csv", index=False)
np.save(OUT + "/GSE174188_cell_donor_idx_v2.npy", np.repeat(np.arange(len(donors)), counts).astype(np.int16))
X = X.astype(np.float16)
b = np.linspace(0, len(X), 11).astype(int)
for i in range(10):
    np.save(OUT + "/GSE174188_cells_grouped_f16_%02d.npy" % i, X[b[i]:b[i+1]])
print("chunks", b.tolist(), flush=True)
print("DONE", flush=True)
