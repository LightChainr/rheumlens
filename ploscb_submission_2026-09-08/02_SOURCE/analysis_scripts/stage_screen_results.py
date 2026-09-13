#!/usr/bin/env python3
"""Turn the five per-seed screen outputs into the tables the figures and the
supplement read. Everything downstream reads these, so there is one conversion
step and no hand-copied numbers.

Writes:
  data/design_screen_all_blocks.tsv   one row per (comparison, block), seed 1
  data/FINAL_confounding_spectrum.tsv one row per comparison, the `all` block
  data/collection_block.tsv           one row per comparison, the collection block
  data/seed_stability_all.tsv         one row per (comparison, seed)
  supplementary/Table_S3_full_screen.tsv
  supplementary/Table_S11_seed_stability.tsv
"""
from __future__ import annotations
import os, re, sys
from pathlib import Path
import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from si_names import rename_columns, check_header

ROOT   = Path(__file__).resolve().parent.parent
SCREEN = Path(os.environ.get("SCREEN_DIR", ROOT / "results" / "screen"))
PRIMARY_SEED = 20260907

files = sorted(SCREEN.glob("**/design_screen.tsv"))
if not files:
    sys.exit(f"no design_screen.tsv under {SCREEN}")

def seed_of(f: Path) -> int:
    """Read the seed from the seed_NNNNNNNN directory, not from the whole path.

    Matching an 8-digit run anywhere in the path picks up the dated workspace
    directory instead, and every seed then comes back labelled the same.
    """
    for part in reversed(f.parts):
        m = re.fullmatch(r"seed_(\d{8})", part)
        if m:
            return int(m.group(1))
    raise SystemExit(f"cannot read a seed from {f}")


df = pd.concat([pd.read_csv(f, sep="\t").assign(seed=seed_of(f)) for f in files],
               ignore_index=True)

# The screen names its collection block "batch" internally; the manuscript calls
# it collection. Rename once, here, so nothing downstream has to know both.
df["block"] = df["block"].replace({"batch": "collection", "qc": "quality"})

(ROOT / "data").mkdir(exist_ok=True)
(ROOT / "supplementary").mkdir(exist_ok=True)

one = df[df.seed == PRIMARY_SEED].copy()
if one.empty:
    one = df[df.seed == df.seed.min()].copy()

one.to_csv(ROOT / "data/design_screen_all_blocks.tsv", sep="\t", index=False)
one[one.block == "all"].to_csv(ROOT / "data/FINAL_confounding_spectrum.tsv",
                               sep="\t", index=False)
one[one.block == "collection"].to_csv(ROOT / "data/collection_block.tsv",
                                      sep="\t", index=False)

seed_tab = (df[df.block == "all"]
            [["cohort", "seed", "n_design_feature", "design_auc_linear",
              "design_auc_frozen", "p_design_auc", "I_D_cv", "p_I_D",
              "disease_auc", "observed_frozen_auc", "p_standard",
              "p_collection_preserving"]]
            .sort_values(["cohort", "seed"]))
seed_tab.to_csv(ROOT / "data/seed_stability_all.tsv", sep="\t", index=False)
seed_pub = seed_tab.copy()
seed_pub.columns = rename_columns(seed_pub.columns)
bad = check_header(seed_pub.columns)
if bad:
    raise SystemExit(f"Table S11 still uses retired names: {bad}")
seed_pub.to_csv(ROOT / "supplementary/Table_S11_seed_stability.tsv",
                sep="\t", index=False)

s2 = df.rename(columns={
    "block": "variable_group", "I_D": "V_D_insample", "I_D_cv": "V_D_crossfitted",
    "I_D_null_mean": "V_D_insample_null_mean",
    "I_D_cv_null_mean": "V_D_crossfitted_null_mean",
    "I_D_cv_null_p025": "V_D_crossfitted_null_p025",
    "p_I_D": "p_V_D", "p_I_D_insample": "p_V_D_insample",
    "p_standard": "p_free", "p_collection_preserving": "p_collection_preserving",
})
cols = [c for c in [
    "cohort", "variable_group", "seed", "n_donor", "n_case", "n_design_feature",
    "V_D_insample", "p_V_D_insample", "V_D_insample_null_mean",
    "V_D_crossfitted", "p_V_D", "V_D_crossfitted_null_mean",
    "V_D_crossfitted_null_p025",
    "design_auc_linear", "design_auc_rf", "design_auc_frozen", "p_design_auc",
    "disease_auc", "observed_frozen_auc", "frozen_minus_tuned_auc",
    "p_free", "p_collection_preserving", "n_strata", "strata_definition",
    "degenerate", "degenerate_reason", "underpowered", "underpowered_reason",
] if c in s2.columns]
s2[cols].sort_values(["cohort", "variable_group", "seed"]).to_csv(
    ROOT / "supplementary/Table_S3_full_screen.tsv", sep="\t", index=False)

# ---- S10: incremental value of expression over the recorded metadata ---------
INC = Path(os.environ.get("INCREMENTAL_DIR", ROOT / "results" / "incremental"))
inc_files = sorted(INC.glob("seed_*/incremental.tsv"))
if inc_files:
    inc = pd.concat([pd.read_csv(f, sep="\t") for f in inc_files], ignore_index=True)
    inc = inc.sort_values(["cohort", "seed"])
    inc.to_csv(ROOT / "supplementary/Table_S4_incremental.tsv", sep="\t", index=False)
    inc.to_csv(ROOT / "data/incremental_all_seeds.tsv", sep="\t", index=False)
    print(f"Table S4: {len(inc)} rows over {inc.cohort.nunique()} comparisons")
else:
    print(f"WARNING: no incremental results under {INC}; Table S4 not refreshed")

print(f"seeds: {sorted(df.seed.unique())}")
print(f"comparisons: {df.cohort.nunique()}   rows: {len(df)}")
print(f"blocks present: {sorted(df.block.unique())}")
print("wrote data/*.tsv and supplementary/Table_S3_full_screen.tsv, Table_S11_seed_stability.tsv")
