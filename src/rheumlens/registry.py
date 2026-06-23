from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from rheumlens.aggregators.classic import CompositionCLRAggregator, MeanAggregator
from rheumlens.aggregators.distance import (
    GaussianBuresReferenceAggregator,
    PrototypeAggregator,
    SlicedWassersteinReferenceAggregator,
)
from rheumlens.aggregators.focus import FocusLiteAggregator, FocusQuery
from rheumlens.aggregators.gds import GraphDonorSignatureAggregator
from rheumlens.aggregators.kme import (
    CCKMEUAggregator,
    MMDReferenceAggregator,
    MultiScaleRFFKMEAggregator,
    RobustMedianOfMeansKMEAggregator,
)
from rheumlens.aggregators.moments import (
    MeanSkewAggregator,
    MeanVarianceAggregator,
    MedianMADAggregator,
    QuantileAggregator,
    ShrinkageCovarianceAggregator,
    TailFractionAggregator,
)
from rheumlens.aggregators.red import REDAggregator
from rheumlens.aggregators.uncertainty import UDERWrapper
from rheumlens.bag_models.donorclr import DonorCLR
from rheumlens.bag_models.focus_adapter import FocusAdapterBagModel
from rheumlens.bag_models.mixmil import GaussianMixMIL
from rheumlens.bag_models.torch_models import TorchBagClassifier
from rheumlens.constants import PRIMARY_ISG_15
from rheumlens.estimators.linear import LogisticL2Estimator
from rheumlens.estimators.uncertainty import ReliabilityWeightedLogistic
from rheumlens.evaluation.engine import BagMethod, FixedMethod
from rheumlens.postprocessors.covariates import CovariateResidualizer
from rheumlens.providers.covariates import DonorCovariateProvider
from rheumlens.providers.expression import (
    CellPCAProvider,
    DonorExpressionPCAProvider,
    DonorMeanHVGProvider,
    ISGProvider,
    RawPseudobulkProvider,
)
from rheumlens.providers.frozen import FrozenCellProvider
from rheumlens.providers.multiview import CellTypePseudobulkProvider, ScFeaturesStyleProvider


@dataclass(frozen=True)
class RegisteredMethod:
    method: FixedMethod | BagMethod
    data_key: str
    expression_key: str | None = None
    learned: bool = False


def _focus_queries(config: dict[str, Any] | None) -> tuple[FocusQuery, ...]:
    if not config:
        return (FocusQuery("type_I_interferon", PRIMARY_ISG_15, "global"),)
    result: list[FocusQuery] = []
    for item in config.get("queries", []):
        compartments = item.get("compartments", ["global"])
        for compartment in compartments:
            result.append(
                FocusQuery(
                    id=str(item["id"]),
                    genes=tuple(map(str, item.get("genes", []))),
                    compartment=str(compartment),
                    construction=str(item.get("construction", "gene_set_high_low_prototype")),
                )
            )
    return tuple(result)


def _parse_method_id(method_id: str) -> tuple[str, str | None]:
    """Allow method families to be evaluated on multiple cell representations.

    Examples: `kme_multiscale@scgpt`, `kme_multiscale@geneformer`,
    `moments_mean_var@lognorm`. The suffix overrides the default data key while preserving
    the registered method name in outputs.
    """

    if "@" not in method_id:
        return method_id, None
    base, data_key = method_id.rsplit("@", 1)
    return base, data_key


def _cell_provider(defaults: dict[str, Any]) -> CellPCAProvider:
    return CellPCAProvider(
        n_components=int(defaults.get("cell_pca_dim", 50)),
        max_fit_cells=int(defaults.get("cell_pca_max_fit_cells", 200_000)),
        standardize=bool(defaults.get("cell_pca_standardize", True)),
    )


def _kme_kwargs(defaults: dict[str, Any], seed: int) -> dict[str, Any]:
    return {
        "rff_dim": int(defaults.get("kme_rff_dim", 256)),
        "scales": tuple(defaults.get("kme_scales", [0.5, 1.0, 2.0])),
        "seed": seed,
        "max_bandwidth_points": int(defaults.get("kme_max_bandwidth_points", 4000)),
        "max_diagnostic_points": int(defaults.get("kme_max_diagnostic_points", 512)),
    }


def build_method(
    method_id: str,
    defaults: dict[str, Any] | None = None,
    query_bank: dict[str, Any] | None = None,
    seed: int = 0,
) -> RegisteredMethod:
    d = defaults or {}
    base_id, data_override = _parse_method_id(method_id)
    logistic = lambda: LogisticL2Estimator(C=float(d.get("estimator_C", 1.0)), random_state=seed)
    fixed = lambda provider, aggregator, estimator=None, postprocessor=None: FixedMethod(
        id=method_id,
        provider=provider,
        aggregator=aggregator,
        estimator=estimator or logistic(),
        postprocessor=postprocessor,
        params={"base_id": base_id, "data_override": data_override},
    )

    # Family A: authoritative and classic baselines.
    if base_id == "donor_mean_hvg":
        return RegisteredMethod(fixed(DonorMeanHVGProvider(n_hvg=int(d.get("n_hvg", 2000))), None), data_override or "lognorm")
    if base_id == "donor_expression_pca":
        return RegisteredMethod(
            fixed(
                DonorExpressionPCAProvider(
                    n_hvg=int(d.get("n_hvg", 2000)),
                    n_components=int(d.get("donor_pca_dim", 25)),
                ),
                None,
            ),
            data_override or "lognorm",
        )
    if base_id in {"scgpt_mean", "geneformer_mean"}:
        key = data_override or ("scgpt" if base_id.startswith("scgpt") else "geneformer")
        return RegisteredMethod(fixed(FrozenCellProvider(), MeanAggregator()), key)
    if base_id == "raw_pseudobulk":
        return RegisteredMethod(
            fixed(RawPseudobulkProvider(n_features=int(d.get("pseudobulk_features", 2000))), None),
            data_override or "raw_counts",
        )
    if base_id == "celltype_pseudobulk":
        return RegisteredMethod(
            fixed(CellTypePseudobulkProvider(n_features_per_type=int(d.get("pseudobulk_features_per_type", 250))), None),
            data_override or "raw_counts",
        )
    if base_id == "scfeatures_multiview":
        return RegisteredMethod(
            fixed(ScFeaturesStyleProvider(n_genes_per_type=int(d.get("scfeatures_genes_per_type", 50))), None),
            data_override or "raw_counts",
        )
    if base_id == "isg_scalar":
        return RegisteredMethod(fixed(ISGProvider(PRIMARY_ISG_15, vector=False), None), data_override or "lognorm")
    if base_id == "isg_vector":
        return RegisteredMethod(fixed(ISGProvider(PRIMARY_ISG_15, vector=True), None), data_override or "lognorm")
    if base_id == "composition_clr":
        return RegisteredMethod(fixed(FrozenCellProvider(), CompositionCLRAggregator()), data_override or "scgpt")
    if base_id == "donor_table_features":
        # The NPZ contains one row/cell per donor; mean pooling is therefore identity at donor level.
        return RegisteredMethod(fixed(FrozenCellProvider(), MeanAggregator()), data_override or "donor_table")
    if base_id == "covariates_only":
        return RegisteredMethod(fixed(DonorCovariateProvider(), None), data_override or "scgpt")

    # Fold-contained residualized sensitivity methods.
    residualizer = CovariateResidualizer(alpha=float(d.get("covariate_residual_alpha", 1.0)))
    if base_id == "scgpt_mean_covres":
        return RegisteredMethod(fixed(FrozenCellProvider(), MeanAggregator(), postprocessor=residualizer), data_override or "scgpt")
    if base_id == "geneformer_mean_covres":
        return RegisteredMethod(fixed(FrozenCellProvider(), MeanAggregator(), postprocessor=residualizer), data_override or "geneformer")
    if base_id == "donor_mean_hvg_covres":
        return RegisteredMethod(fixed(DonorMeanHVGProvider(), None, postprocessor=residualizer), data_override or "lognorm")
    if base_id == "donor_expression_pca_covres":
        return RegisteredMethod(fixed(DonorExpressionPCAProvider(), None, postprocessor=residualizer), data_override or "lognorm")

    # Families B and D fixed donor aggregators are evaluated after fold-contained cell PCA by default.
    provider = _cell_provider(d)
    agg: Any | None = None
    if base_id == "moments_mean_var":
        agg = MeanVarianceAggregator()
    elif base_id == "robust_median_mad":
        agg = MedianMADAggregator()
    elif base_id == "mean_skew":
        agg = MeanSkewAggregator()
    elif base_id == "quantiles":
        agg = QuantileAggregator(quantiles=tuple(d.get("quantiles", [0.1, 0.5, 0.9, 0.95, 0.99])))
    elif base_id == "tail_fractions":
        agg = TailFractionAggregator(quantiles=tuple(d.get("tail_quantiles", [0.95, 0.99])))
    elif base_id == "shrinkage_covariance":
        agg = ShrinkageCovarianceAggregator(diagonal_only=bool(d.get("covariance_diagonal_only", False)))
    elif base_id == "kme_multiscale":
        agg = MultiScaleRFFKMEAggregator(**_kme_kwargs(d, seed), include_linear_mean=True)
    elif base_id == "robust_kme":
        agg = RobustMedianOfMeansKMEAggregator(
            **_kme_kwargs(d, seed), n_blocks=int(d.get("robust_kme_blocks", 5))
        )
    elif base_id == "mmd_control":
        agg = MMDReferenceAggregator(
            **_kme_kwargs(d, seed), cell_type_conditional=bool(d.get("mmd_cell_type_conditional", False))
        )
    elif base_id == "true_bures":
        agg = GaussianBuresReferenceAggregator(n_references=int(d.get("distance_references", 16)), seed=seed)
    elif base_id == "sliced_wasserstein":
        agg = SlicedWassersteinReferenceAggregator(
            n_directions=int(d.get("sw_directions", 64)),
            n_quantiles=int(d.get("sw_quantiles", 64)),
            n_references=int(d.get("distance_references", 16)),
            seed=seed,
        )
    elif base_id == "prototype_histogram":
        agg = PrototypeAggregator(n_prototypes=int(d.get("n_prototypes", 32)), seed=seed)
    elif base_id == "gds":
        agg = GraphDonorSignatureAggregator(seed=seed)
    elif base_id == "red":
        agg = REDAggregator(seed=seed)
    elif base_id in {"uder_meanvar", "uder_meanvar_weighted", "uder_kme", "uder_kme_weighted"}:
        if "kme" in base_id:
            base_agg = MultiScaleRFFKMEAggregator(**_kme_kwargs(d, seed), include_linear_mean=True)
        else:
            base_agg = MeanVarianceAggregator()
        agg = UDERWrapper(
            base=base_agg,
            subsample_size=int(d.get("matched_cells_per_donor", 500)),
            n_subsamples=int(d.get("uder_subsamples", 20)),
            append_log_variance=True,
            seed=seed,
        )
        estimator = (
            ReliabilityWeightedLogistic(C=float(d.get("estimator_C", 1.0)), random_state=seed)
            if base_id.endswith("weighted")
            else logistic()
        )
        return RegisteredMethod(fixed(provider, agg, estimator=estimator), data_override or "scgpt")
    elif base_id in {"cckme_u", "cckme_u_weighted"}:
        agg = CCKMEUAggregator(
            **_kme_kwargs(d, seed),
            include_linear_mean=True,
            subsample_size=int(d.get("matched_cells_per_donor", 500)),
            n_subsamples=int(d.get("cckme_subsamples", 20)),
            tau=float(d.get("cckme_tau", 1.0)),
        )
        estimator = ReliabilityWeightedLogistic(C=float(d.get("estimator_C", 1.0)), random_state=seed) if base_id.endswith("weighted") else logistic()
        return RegisteredMethod(fixed(provider, agg, estimator=estimator), data_override or "scgpt")
    elif base_id in {"focus_lite", "language_focus"}:
        agg = FocusLiteAggregator(
            queries=_focus_queries(query_bank),
            topk_fraction=float(d.get("focus_topk_fraction", 0.05)),
            high_quantile=float(d.get("focus_high_quantile", 0.90)),
            low_quantile=float(d.get("focus_low_quantile", 0.10)),
            positive_quantile=float(d.get("focus_positive_quantile", 0.95)),
            random_queries=int(query_bank.get("negative_controls", {}).get("random_queries", 0)) if query_bank else 0,
            seed=seed,
        )
        # FOCUS operates on the raw frozen embedding; queries are derived using aligned expression.
        return RegisteredMethod(
            fixed(FrozenCellProvider(), agg),
            data_override or "scgpt",
            expression_key="lognorm",
        )
    if agg is not None:
        return RegisteredMethod(fixed(provider, agg), data_override or "scgpt")

    # Family C and trainable original models.
    bag_provider = _cell_provider(d)
    if base_id == "deep_sets":
        model = TorchBagClassifier("deepsets", seed=seed, epochs=int(d.get("bag_epochs", 100)))
    elif base_id == "attention_mil":
        model = TorchBagClassifier("attention", seed=seed, epochs=int(d.get("bag_epochs", 100)))
    elif base_id == "topk_mil":
        model = TorchBagClassifier(
            "topk_attention",
            seed=seed,
            epochs=int(d.get("bag_epochs", 100)),
            topk_fraction=float(d.get("focus_topk_fraction", 0.05)),
        )
    elif base_id == "set_transformer":
        model = TorchBagClassifier("set_transformer", seed=seed, epochs=int(d.get("bag_epochs", 100)))
    elif base_id == "mixmil":
        model = GaussianMixMIL(n_components=int(d.get("mixmil_components", 16)), seed=seed)
    elif base_id == "donorclr":
        model = DonorCLR(
            hidden=int(d.get("donorclr_hidden", 64)),
            cells_per_view=int(d.get("donorclr_cells_per_view", 256)),
            epochs=int(d.get("donorclr_epochs", 100)),
            seed=seed,
        )
    elif base_id == "focus_adapter":
        path = d.get("focus_initial_queries_npy")
        if not path or not Path(path).exists():
            raise ValueError("focus_adapter requires method_defaults.focus_initial_queries_npy")
        queries = np.load(path)
        model = FocusAdapterBagModel(
            initial_queries=queries,
            rank=int(d.get("focus_adapter_rank", 4)),
            epochs=int(d.get("bag_epochs", 100)),
            seed=seed,
        )
        # Query vectors are defined in the frozen embedding space, so do not apply cell PCA.
        bag_provider = FrozenCellProvider()
    else:
        raise KeyError(f"unknown method_id: {method_id}")
    return RegisteredMethod(BagMethod(method_id, bag_provider, model), data_override or "scgpt", learned=True)
