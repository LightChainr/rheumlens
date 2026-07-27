from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


WORKSPACE = Path(__file__).resolve().parents[1]
REPO = WORKSPACE.parent
REFRAME = REPO / "04_重构_投稿定位_20260714"
ARCHIVE = (
    REPO
    / "03_远程回传"
    / "final_gpu_figure_ready_20260709"
    / "extracted_plus"
    / "RheumLens_GPU_figure_ready_plus_20260709"
)
EXTRA = ARCHIVE / "extra"
H5AD = REPO / "05_原始数据" / "GSE174188_CELLxGENE_2025-11-08.h5ad"
GENEFORMER_METHOD = "geneformer_v2_316m_cell_sample1000_clspool_logistic_maxlen4096_seed001"
DATASETS = ("SLE_GSE135779", "SLE_GSE174188_CD4", "SLE_GSE285773_CD4")


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def label_path(dataset: str) -> Path:
    return EXTRA / "pseudobulk" / dataset / "donor_labels.tsv"


def pseudobulk_path(dataset: str) -> Path:
    return EXTRA / "pseudobulk" / dataset / "donor_log1p_cpm.parquet"


def embedding_path(dataset: str, cap: str | int = 1000) -> Path:
    if cap == "allavailable":
        method = "geneformer_v2_316m_cell_allavailable_clspool_logistic_allavailable_maxlen4096_seed001"
    else:
        method = f"geneformer_v2_316m_cell_sample{cap}_clspool_logistic_maxlen4096_seed001"
    return EXTRA / "04_models" / "Geneformer" / dataset / method / "donor_embedding.parquet"


def load_labels(dataset: str) -> pd.Series:
    labels = pd.read_csv(label_path(dataset), sep="\t", dtype={"donor_id": str})
    labels = labels.set_index("donor_id")["case_control"].astype(str).sort_index()
    if labels.index.has_duplicates or set(labels) != {"case", "control"}:
        raise RuntimeError(f"Invalid labels for {dataset}")
    return labels


def load_embedding(dataset: str, cap: str | int = 1000) -> pd.DataFrame:
    frame = pd.read_parquet(embedding_path(dataset, cap))
    frame.index = frame.index.astype(str)
    labels = load_labels(dataset)
    if frame.index.has_duplicates or set(frame.index) != set(labels.index):
        raise RuntimeError(f"Embedding donor mismatch for {dataset}, cap={cap}")
    frame = frame.loc[labels.index]
    if not np.isfinite(frame.to_numpy(dtype=np.float32)).all():
        raise RuntimeError(f"Non-finite embedding for {dataset}, cap={cap}")
    return frame


def load_pseudobulk(dataset: str) -> pd.DataFrame:
    frame = pd.read_parquet(pseudobulk_path(dataset))
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str).str.strip()
    labels = load_labels(dataset)
    if frame.index.has_duplicates or frame.columns.has_duplicates or set(frame.index) != set(labels.index):
        raise RuntimeError(f"Pseudobulk donor/feature mismatch for {dataset}")
    frame = frame.loc[labels.index]
    if not np.isfinite(frame.to_numpy(dtype=np.float32)).all():
        raise RuntimeError(f"Non-finite pseudobulk for {dataset}")
    return frame
