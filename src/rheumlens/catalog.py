from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class MethodInfo:
    method_id: str
    family: str
    priority: str
    role: str
    implementation: str = "native"
    confirmatory_eligible: bool = True


METHOD_CATALOG: tuple[MethodInfo, ...] = (
    MethodInfo("donor_mean_hvg", "A", "P0", "Fold-contained donor-mean HVG expression"),
    MethodInfo("donor_expression_pca", "A", "P0", "Fold-contained donor-expression PCA"),
    MethodInfo("scgpt_mean", "A", "P0", "Frozen scGPT mean pooling"),
    MethodInfo("geneformer_mean", "A", "P0", "Frozen Geneformer mean pooling"),
    MethodInfo("raw_pseudobulk", "A", "P0", "Raw-count donor pseudobulk"),
    MethodInfo("celltype_pseudobulk", "A", "P1", "Cell-type raw-count pseudobulk"),
    MethodInfo("composition_clr", "A", "P0", "Cell composition with CLR transform"),
    MethodInfo("donor_table_features", "A", "P0", "Externally supplied donor-level feature table"),
    MethodInfo("isg_scalar", "A", "P0", "Fixed 15-gene ISG scalar"),
    MethodInfo("isg_vector", "A", "P0", "Fixed 15-gene ISG vector"),
    MethodInfo("scfeatures_multiview", "A", "P1", "scFeatures-style multiview summary"),
    MethodInfo("moments_mean_var", "B", "P0", "Mean plus diagonal variance"),
    MethodInfo("robust_median_mad", "B", "P0", "Median plus MAD"),
    MethodInfo("mean_skew", "B", "P1", "Mean plus skewness"),
    MethodInfo("quantiles", "B", "P0", "Quantile ladder"),
    MethodInfo("tail_fractions", "B", "P0", "Training-threshold tail fractions"),
    MethodInfo("shrinkage_covariance", "B", "P1", "Shrinkage covariance summary"),
    MethodInfo("kme_multiscale", "B", "P0", "Multi-scale RFF KME"),
    MethodInfo("robust_kme", "B", "P1", "Median-of-means robust KME"),
    MethodInfo("mmd_control", "B", "P1", "MMD-to-control reference"),
    MethodInfo("true_bures", "B", "P1", "Gaussian Bures distances to training references"),
    MethodInfo("sliced_wasserstein", "B", "P1", "Sliced-Wasserstein references"),
    MethodInfo("prototype_histogram", "B", "P1", "Prototype occupancy and distance"),
    MethodInfo("deep_sets", "C", "P1", "Deep Sets bag model"),
    MethodInfo("attention_mil", "C", "P0", "Gated attention MIL"),
    MethodInfo("topk_mil", "C", "P1", "Sparse top-k MIL"),
    MethodInfo("set_transformer", "C", "P2", "Set Transformer"),
    MethodInfo("mixmil", "C", "P1", "Internal mixture-MIL approximation", "approximation", False),
    MethodInfo("donorclr", "D", "P2", "Contrastive donor subsample encoder", "native", False),
    MethodInfo("cckme_u", "D", "P0", "Cell-count calibrated KME with uncertainty"),
    MethodInfo("uder_meanvar", "D", "P1", "Repeated-subsample mean/variance with representation uncertainty"),
    MethodInfo("uder_kme", "D", "P1", "Repeated-subsample KME with representation uncertainty"),
    MethodInfo("focus_lite", "D", "P0", "Fixed mechanism-query FOCUS"),
    MethodInfo("focus_adapter", "D", "P1", "Low-rank learnable FOCUS adapter", "native", False),
    MethodInfo("red", "D", "P0", "Reference-relative evidence decomposition"),
    MethodInfo("gds", "D", "P2", "Graph donor signature", "native", False),
)


def catalog_rows() -> list[dict[str, object]]:
    return [asdict(item) for item in METHOD_CATALOG]


def catalog_ids() -> tuple[str, ...]:
    return tuple(item.method_id for item in METHOD_CATALOG)
