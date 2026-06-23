from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from rheumlens.aggregators.base import DonorAggregator
from rheumlens.types import CellDataset, DonorFeatures, FitContext
from rheumlens.utils.arrays import as_dense, grouped_indices


@dataclass(frozen=True)
class FocusQuery:
    id: str
    genes: tuple[str, ...] = ()
    compartment: str = "global"
    construction: str = "gene_set_high_low_prototype"


@dataclass
class FocusLiteAggregator(DonorAggregator):
    queries: tuple[FocusQuery, ...]
    topk_fraction: float = 0.05
    high_quantile: float = 0.9
    low_quantile: float = 0.1
    positive_quantile: float = 0.95
    random_queries: int = 0
    seed: int = 0

    def __post_init__(self) -> None:
        self.query_vectors_: dict[str, np.ndarray] = {}
        self.thresholds_: dict[str, float] = {}
        self.query_ids_: list[str] = []

    @staticmethod
    def _cosine_scores(X: np.ndarray, q: np.ndarray) -> np.ndarray:
        Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
        qn = q / (np.linalg.norm(q) + 1e-12)
        return Xn @ qn

    def _compartment_mask(self, data: CellDataset, compartment: str) -> np.ndarray:
        if compartment == "global":
            return np.ones(data.X.shape[0], dtype=bool)
        if data.cell_types is None:
            raise ValueError(f"query compartment {compartment} requires cell_types")
        return data.cell_types.astype(str) == compartment

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "FocusLiteAggregator":
        if context.expression is None:
            raise ValueError("FOCUS requires an aligned expression dataset in FitContext.expression")
        embedding = data.subset_donors(train_donors)
        expression = context.expression.subset_donors(train_donors)
        if not np.array_equal(embedding.cell_ids, expression.cell_ids):
            index = {cell: i for i, cell in enumerate(expression.cell_ids.astype(str))}
            try:
                order = np.asarray([index[x] for x in embedding.cell_ids.astype(str)], dtype=int)
            except KeyError as exc:
                raise ValueError("expression and embedding cells are not aligned") from exc
            metadata = None if expression.metadata is None else expression.metadata.iloc[order].reset_index(drop=True)
            expression = CellDataset(
                X=expression.X[order],
                cell_ids=expression.cell_ids[order],
                donor_ids=expression.donor_ids[order],
                y=expression.y[order],
                feature_names=expression.feature_names,
                cell_types=None if expression.cell_types is None else expression.cell_types[order],
                cohorts=None if expression.cohorts is None else expression.cohorts[order],
                metadata=metadata,
                name=expression.name,
            )
            if not np.array_equal(embedding.cell_ids, expression.cell_ids):
                raise ValueError("failed to reorder expression cells to embedding cells")
        emb_X = as_dense(embedding.X)
        expr_X = as_dense(expression.X)
        gene_lookup = {str(g).upper(): i for i, g in enumerate(expression.feature_names)}
        control_mask = np.asarray([context.y_by_donor[str(d)] == 0 for d in embedding.donor_ids])
        if not control_mask.any():
            raise ValueError("FOCUS requires training controls for thresholds")
        self.query_vectors_.clear()
        self.thresholds_.clear()
        self.query_ids_.clear()
        for query in self.queries:
            comp_mask = self._compartment_mask(embedding, query.compartment)
            indices = [gene_lookup[g.upper()] for g in query.genes if g.upper() in gene_lookup]
            if not indices:
                continue
            mechanism_score = expr_X[:, indices].mean(axis=1)
            eligible = np.flatnonzero(comp_mask)
            if eligible.size < 10:
                continue
            high = eligible[mechanism_score[eligible] >= np.quantile(mechanism_score[eligible], self.high_quantile)]
            low = eligible[mechanism_score[eligible] <= np.quantile(mechanism_score[eligible], self.low_quantile)]
            if len(high) == 0 or len(low) == 0:
                continue
            q = emb_X[high].mean(axis=0) - emb_X[low].mean(axis=0)
            if np.linalg.norm(q) <= 1e-12:
                continue
            key = f"{query.id}__{query.compartment}"
            self.query_vectors_[key] = q
            scores = self._cosine_scores(emb_X, q)
            threshold_pool = scores[control_mask & comp_mask]
            if threshold_pool.size == 0:
                threshold_pool = scores[control_mask]
            self.thresholds_[key] = float(np.quantile(threshold_pool, self.positive_quantile))
            self.query_ids_.append(key)
        rng = np.random.default_rng(self.seed)
        for i in range(self.random_queries):
            key = f"random_{i:03d}"
            q = rng.normal(size=emb_X.shape[1])
            self.query_vectors_[key] = q
            scores = self._cosine_scores(emb_X, q)
            self.thresholds_[key] = float(np.quantile(scores[control_mask], self.positive_quantile))
            self.query_ids_.append(key)
        if not self.query_ids_:
            raise ValueError("no FOCUS queries could be constructed")
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if not self.query_ids_:
            raise RuntimeError("aggregator is not fitted")
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        X_all = as_dense(subset.X)
        rows: list[np.ndarray] = []
        feature_names: list[str] = []
        stats = ("mean", "topk", "q95", "q99", "positive_fraction", "variance")
        for query_id in self.query_ids_:
            feature_names.extend(f"{query_id}__{stat}" for stat in stats)
        for donor in order:
            X = X_all[groups[donor]]
            values = []
            for query_id in self.query_ids_:
                scores = self._cosine_scores(X, self.query_vectors_[query_id])
                k = max(1, int(np.ceil(self.topk_fraction * len(scores))))
                top = np.partition(scores, len(scores) - k)[-k:]
                values.extend(
                    [
                        scores.mean(),
                        top.mean(),
                        np.quantile(scores, 0.95),
                        np.quantile(scores, 0.99),
                        np.mean(scores > self.thresholds_[query_id]),
                        scores.var(ddof=1) if len(scores) > 1 else 0.0,
                    ]
                )
            rows.append(np.asarray(values))
        diagnostics: dict[str, Any] = {
            "query_ids": self.query_ids_,
            "thresholds": self.thresholds_,
            "query_norms": {k: float(np.linalg.norm(v)) for k, v in self.query_vectors_.items()},
        }
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, feature_names, diagnostics=diagnostics)
