from __future__ import annotations

from pathlib import Path

from rheumlens.config import load_yaml


def test_yaml_expands_environment(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setenv("RL_TEST_ROOT", str(tmp_path))
    path = tmp_path / "config.yaml"
    path.write_text("paths:\n  root: ${RL_TEST_ROOT}/project\nitems:\n  - ~/example\n", encoding="utf-8")
    value = load_yaml(path)
    assert value["paths"]["root"] == str(tmp_path / "project")
    assert not value["items"][0].startswith("~")
