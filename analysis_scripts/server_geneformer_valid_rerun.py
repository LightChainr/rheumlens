#!/usr/bin/env python3
"""Valid Geneformer V2 tokenization and frozen embedding extraction on GPU server."""

from pathlib import Path
import argparse
import json
import pickle

import anndata as ad
import numpy as np

BASE = Path("/autodl-fs/data/rheumlens_20260619")
REPO = BASE / "code/Geneformer"
DATA = BASE / "data/geneformer_rerun"
RESULTS = BASE / "results/geneformer_valid"
MODEL = REPO / "Geneformer-V2-104M"
SOURCE = DATA / "GSE135779_geneformer_500_per_donor_raw_counts.h5ad"
PREPARED = DATA / "official_input/GSE135779_geneformer_official_input.h5ad"
TOKENIZED = DATA / "tokenized/GSE135779_v2.dataset"


def prepare():
    RESULTS.mkdir(parents=True, exist_ok=True)
    PREPARED.parent.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(SOURCE)
    mapping = pickle.load(open(REPO / "geneformer/gene_name_id_dict_gc104M.pkl", "rb"))
    symbols = np.asarray(a.var_names.astype(str))
    clean = np.asarray([s.split("__dup", 1)[0].upper() for s in symbols])
    ensembl = np.asarray([mapping.get(s) for s in clean], dtype=object)
    mapped = np.asarray([x is not None for x in ensembl])
    # Unmapped identifiers remain unique and are filtered by the official
    # tokenizer vocabulary. They are not converted into model tokens.
    ensembl[~mapped] = [f"UNMAPPED_{i}" for i in np.flatnonzero(~mapped)]
    a.var["gene_symbol"] = symbols
    a.var["ensembl_id"] = ensembl.astype(str)
    a.obs["n_counts"] = np.asarray(a.X.sum(axis=1)).ravel().astype(np.int64)
    for c in ["donor_id", "disease", "marker_lineage"]:
        a.obs[c] = a.obs[c].astype(str)
    a.write_h5ad(PREPARED, compression="gzip")
    info = {"source": str(SOURCE), "prepared": str(PREPARED),
            "n_cells": int(a.n_obs), "n_genes_input": int(a.n_vars),
            "n_gene_symbols_mapped_to_ensembl": int(mapped.sum()),
            "mapping_fraction": float(mapped.mean()),
            "model": str(MODEL), "model_revision": "04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5",
            "tokenization": "official Geneformer V2 rank-value tokenizer; raw counts"}
    (RESULTS / "preparation_metadata.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info, indent=2), flush=True)


def tokenize():
    from geneformer import TranscriptomeTokenizer
    TOKENIZED.parent.mkdir(parents=True, exist_ok=True)
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"donor_id": "donor_id", "disease": "disease",
                               "marker_lineage": "marker_lineage", "source_row_index": "source_row_index"},
        nproc=8, chunk_size=512, model_version="V2", model_input_size=4096,
        special_token=True, collapse_gene_ids=True, use_h5ad_index=False)
    tk.tokenize_data(PREPARED.parent, TOKENIZED.parent, "GSE135779_v2",
                     file_format="h5ad", input_identifier="official_input")
    print(TOKENIZED, flush=True)


def embed(max_ncells):
    from geneformer import EmbExtractor
    RESULTS.mkdir(parents=True, exist_ok=True)
    prefix = "GSE135779_v2_104M_all" if max_ncells is None else f"GSE135779_v2_104M_pilot{max_ncells}"
    ex = EmbExtractor(model_type="Pretrained", num_classes=0, emb_mode="cls",
                      max_ncells=max_ncells, emb_layer=-1,
                      emb_label=["donor_id", "disease", "marker_lineage", "source_row_index"],
                      forward_batch_size=8, nproc=8, model_version="V2")
    embs = ex.extract_embs(MODEL, TOKENIZED, RESULTS, prefix, output_torch_embs=True)
    print(type(embs), getattr(embs, "shape", None), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=["prepare", "tokenize", "pilot", "full"])
    p.add_argument("--pilot-cells", type=int, default=512)
    args = p.parse_args()
    if args.stage == "prepare": prepare()
    elif args.stage == "tokenize": tokenize()
    elif args.stage == "pilot": embed(args.pilot_cells)
    else: embed(None)


if __name__ == "__main__":
    main()
