from __future__ import annotations

from pathlib import Path

import yaml

from rheumlens.run import main


def test_smoke_cli(tmp_path):
    config = {
        "seed": 0,
        "paths": {"root": str(tmp_path / "smoke"), "results": str(tmp_path / "smoke" / "results")},
        "evaluation": {"bootstrap_reps": 50, "estimator_C": 1.0},
        "analysis": {"reference_method": "scgpt_mean"},
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    main(["--config", str(path), "--stage", "smoke"])
    assert (tmp_path / "smoke" / "results" / "final_tables" / "method_summary.csv").exists()
