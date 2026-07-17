"""Data persistence, validation, and cross-cohort alignment utilities."""

from rheumlens.data.io import load_npz_dataset, save_npz_dataset
from rheumlens.data.operations import align_features
from rheumlens.data.validation import validate_embedding, validate_raw_counts

__all__ = [
    "align_features",
    "load_npz_dataset",
    "save_npz_dataset",
    "validate_embedding",
    "validate_raw_counts",
]
