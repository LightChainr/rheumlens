#!/usr/bin/env python3
"""Generate Table 1 from the screen output, across all seeds.

Table 1 was previously assembled by hand and drifted from the results it summarises:
one comparison with a minority class of 10 was labelled estimable although the
declared rule flags anything below 15 as underpowered. It is now generated, and the
verdict column is read from the screen's own gate rather than retyped.
"""
from __future__ import annotations
import os, re, sys
from pathlib import Path
import pandas as pd

SCREEN = Path(os.environ.get("SCREEN_DIR", "results/screen"))
OUT    = Path(__file__).resolve().parent.parent / "manuscript" / "table1.md"

LABEL = {
    "COMBAT_INFLUENZA": "Influenza · COMBAT",
    "COMBAT_CROSS": "Sepsis vs COVID-19 · COMBAT",
    "SLE_GSE174188_CD4": "SLE · GSE174188 (CD4)",
    "SLE_GSE174188_ALLCELL": "SLE · GSE174188 (all cells)",
    "COMBAT_COVID": "COVID-19 · COMBAT",
    "COVID_REN": "COVID-19 · Ren",
    "COVID_STEPHENSON": "COVID-19 · Stephenson",
    "COVID_REN_ASSAY_10x_5_v2": "COVID-19 · Ren (10x 5' v2)",
    "CMV_HIHA": "CMV · HIHA",
}


def seed_of(f: Path) -> int:
    """Read the seed from the seed_NNNNNNNN directory component only.

    Searching the whole path for eight digits matched the workspace directory
    name and labelled every run with the same seed.
    """
    for part in reversed(f.parts):
        m = re.fullmatch(r"seed_(\d{8})", part)
        if m:
            return int(m.group(1))
    raise SystemExit(f"cannot read a seed from {f}")


def rng(v, fmt="{:.3f}"):
    v = [x for x in v if pd.notna(x)]
    if not v:
        return "—"
    lo, hi = min(v), max(v)
    return fmt.format(lo) if lo == hi else f"{fmt.format(lo)}-{fmt.format(hi)}"


def main():
    files = sorted(SCREEN.glob("**/design_screen.tsv"))
    if not files:
        sys.exit(f"no design_screen.tsv under {SCREEN}")
    df = pd.concat([pd.read_csv(f, sep="\t").assign(seed=seed_of(f)) for f in files])

    rows = []
    for coh, g in df.groupby("cohort"):
        allb = g[g.block == "all"]
        colb = g[g.block == "batch"]          # collection-only block
        seeds = sorted(allb.seed.unique())
        minority = int(min(allb.n_case.iloc[0],
                           allb.n_donor.iloc[0] - allb.n_case.iloc[0]))
        if bool(allb.degenerate.any()):
            verdict = "**not estimable**"
        elif bool(allb.underpowered.any()):
            verdict = "**underpowered**"
        else:
            verdict = "estimable"
        rows.append(dict(
            cohort=coh,
            Comparison=LABEL.get(coh, coh),
            Donors=int(allb.n_donor.iloc[0]),
            Minority=minority,
            Seeds=len(seeds),
            all_cols=int(allb.n_design_feature.iloc[0]),
            all_auc=rng(allb.design_auc_linear),
            all_p=rng(allb.p_design_auc, "{:.4f}"),
            all_vd=rng(allb.I_D_cv),
            all_pvd=rng(allb.p_I_D, "{:.4f}"),
            all_pvd_ins=rng(allb.p_I_D_insample, "{:.4f}"),
            coll_cols=(str(int(colb.n_design_feature.iloc[0])) if len(colb) else "n/a"),
            coll_auc=(rng(colb.design_auc_linear) if len(colb) else "n/a"),
            coll_p=(rng(colb.p_design_auc, "{:.4f}") if len(colb) else "n/a"),
            frozen=rng(allb.observed_frozen_auc),
            tuned=rng(allb.disease_auc),
            p_free=rng(allb.p_standard, "{:.4f}"),
            p_coll=rng(allb.p_collection_preserving, "{:.4f}"),
            verdict=verdict,
            sort=float(allb.design_auc_linear.max()),
        ))
    t = pd.DataFrame(rows).sort_values("sort", ascending=False)

    # Thirteen columns. The tuned diagnosis AUC, the in-sample V_D pair and the
    # per-seed values are in Tables S6 and S8; putting them here as well made a
    # sixteen-column table that no reader would work through.
    head = ("| Comparison | Donors | Minority | "
            "Recorded: cols | AUC | p | V_D | p(V_D) | "
            "Collection: cols | AUC | p | "
            "Diagnosis AUC | p free | p coll.-pres. | Verdict |")
    sep = "|---|---:|---:|---:|---:|---|---:|---|---:|---:|---|---:|---|---|---|"
    lines = [head, sep]
    for _, r in t.iterrows():
        lines.append(
            f"| {r.Comparison} | {r.Donors} | {r.Minority} | "
            f"{r.all_cols} | {r.all_auc} | {r.all_p} | {r.all_vd} | {r.all_pvd} | "
            f"{r.coll_cols} | {r.coll_auc} | {r.coll_p} | "
            f"{r.frozen} | {r.p_free} | {r.p_coll} | {r.verdict} |")
    OUT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote {OUT}")

    n_all = int((t.all_p.str.split("-").str[-1].astype(float) <= 0.05).sum())
    n_all_vd = int((t.all_pvd.str.split("-").str[-1].astype(float) <= 0.05).sum())
    coll = t[t.coll_p != "n/a"]
    n_coll = int((coll.coll_p.str.split("-").str[-1].astype(float) <= 0.05).sum())
    est = t[t.verdict == "estimable"]
    est_coll = est[est.coll_p != "n/a"]
    n_coll_est = int((est_coll.coll_p.str.split("-").str[-1].astype(float) <= 0.05).sum())
    print(f"\nHEADLINE COUNTS")
    print(f"  recorded metadata significant : {n_all} of {len(t)}  (by p(V_D): {n_all_vd})")
    print(f"  collection-only significant   : {n_coll} of {len(coll)} that record any")
    print(f"  ... restricted to estimable   : {n_coll_est} of {len(est_coll)}")


if __name__ == "__main__":
    main()
