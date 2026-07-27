"""Portable path configuration for the v2 release candidate.

All release-facing scripts default to files inside the extracted release directory.
Large optional upstream inputs can be supplied with environment variables.
"""

from __future__ import annotations

import os
from pathlib import Path


def configured_path(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    return Path(value).expanduser().resolve() if value else default.resolve()


RELEASE_ROOT = configured_path(
    "RHEUMLENS_RELEASE_ROOT",
    Path(__file__).resolve().parents[1],
)
INPUT_ROOT = configured_path("RHEUMLENS_INPUT_ROOT", RELEASE_ROOT / "inputs")
RESULTS_ROOT = configured_path("RHEUMLENS_RESULTS_ROOT", RELEASE_ROOT / "results")
FIGURES_ROOT = configured_path("RHEUMLENS_FIGURES_ROOT", RELEASE_ROOT / "figures")

DONOR_LEVEL_ROOT = INPUT_ROOT / "donor_level"
PUBLIC_METADATA_ROOT = INPUT_ROOT / "public_metadata"
DESIGN_METADATA_ROOT = INPUT_ROOT / "design_metadata"

RAW_H5AD = configured_path(
    "RHEUMLENS_GSE174188_H5AD",
    INPUT_ROOT / "optional_raw" / "GSE174188_CELLxGENE_2025-11-08.h5ad",
)
LEARNED_POOLING_INPUT_ROOT = configured_path(
    "RHEUMLENS_LEARNED_POOLING_INPUT_ROOT",
    INPUT_ROOT / "optional_learned_pooling",
)
CELL_EMBEDDING_ROOT = configured_path(
    "RHEUMLENS_CELL_EMBEDDING_ROOT",
    INPUT_ROOT / "optional_cell_embeddings_cap500",
)


def donor_file(dataset: str, filename: str) -> Path:
    return DONOR_LEVEL_ROOT / dataset / filename

