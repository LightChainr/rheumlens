from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml


def _expand(value: Any) -> Any:
    """Recursively expand environment variables and ``~`` in YAML string values.

    This lets the server configuration use paths such as ``${RL_ROOT}/data/...`` while
    keeping numeric values, lists and mappings unchanged.
    """

    if isinstance(value, dict):
        return {key: _expand(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_expand(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_expand(item) for item in value)
    if isinstance(value, str):
        return os.path.expanduser(os.path.expandvars(value))
    return value


def load_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected mapping in {path}")
    return _expand(value)
