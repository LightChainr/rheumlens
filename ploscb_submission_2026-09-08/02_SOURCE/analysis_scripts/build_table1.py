#!/usr/bin/env python3
"""Generate Table 1 from the screen output, across all seeds.

Table 1 was previously assembled by hand and drifted from the results it summarises:
one comparison with a minority class of 10 was labelled estimable although the
declared rule flags anything below 15 as underpowered. It is now generated, and the
verdict column is read from the screen's own gate rather than retyped.
"""
from __future__ import annotations
import os, re, sys
import sys
from pathlib import Path
import pandas as pd

ROOT   = Path(__file__).resolve().parent.parent
SCREEN = Path(os.environ.get("SCREEN_DIR", "results/screen"))
INCREMENTAL = Path(os.environ.get("INCREMENTAL_DIR", ROOT / "results" / "incremental"))
OUT    = ROOT / "manuscript" / "table1.md"

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
    """Min-max over seeds. Signed ranges join with "to", because "-0.001-+0.001"
    reads as one malformed number."""
    v = [x for x in v if pd.notna(x)]
    if not v:
        return "\u2014"
    lo, hi = min(v), max(v)
    if lo == hi:
        return fmt.format(lo)
    joiner = " to " if "+" in fmt else "-"
    return f"{fmt.format(lo)}{joiner}{fmt.format(hi)}"


def main():
    files = sorted(SCREEN.glob("**/design_screen.tsv"))
    if not files:
        sys.exit(f"no design_screen.tsv under {SCREEN}")
    df = pd.concat([pd.read_csv(f, sep="\t").assign(seed=seed_of(f)) for f in files])

    inc_files = sorted(INCREMENTAL.glob("seed_*/incremental.tsv"))
    inc = (pd.concat([pd.read_csv(f, sep="\t") for f in inc_files], ignore_index=True)
           if inc_files else pd.DataFrame(columns=["cohort", "delta_over_metadata"]))

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
            all_auc=rng(allb.design_auc_frozen),
            all_p=rng(allb.p_design_auc, "{:.4f}"),
            all_vd=rng(allb.I_D_cv),
            all_pvd=rng(allb.p_I_D, "{:.4f}"),
            all_pvd_ins=rng(allb.p_I_D_insample, "{:.4f}"),
            coll_cols=(str(int(colb.n_design_feature.iloc[0])) if len(colb) else "n/a"),
            coll_auc=(rng(colb.design_auc_frozen) if len(colb) else "n/a"),
            coll_p=(rng(colb.p_design_auc, "{:.4f}") if len(colb) else "n/a"),
            frozen=rng(allb.observed_frozen_auc),
            tuned=rng(allb.disease_auc),
            p_free=rng(allb.p_standard, "{:.4f}"),
            p_coll=rng(allb.p_collection_preserving, "{:.4f}"),
            strata_def=str(allb.strata_definition.iloc[0]),
            delta=(rng(inc[inc.cohort == coh].delta_over_metadata, "{:+.3f}")
                   if len(inc[inc.cohort == coh]) else "n/a"),
            verdict=verdict,
            sort=float(allb.design_auc_frozen.max()),
        ))
    t = pd.DataFrame(rows).sort_values("sort", ascending=False)

    # Every AUC printed here is the frozen prespecified pipeline, which is the
    # statistic each p-value beside it was computed against. The tuned nested-CV
    # AUCs are a description, not an inferential quantity, and live in Table S4:
    # printing a tuned AUC next to a p-value computed from the frozen one invites
    # the reader to treat the p as a test of the number they can see.
    #
    # A dagger marks a comparison whose restricted permutation had no collection
    # variable to condition on and fell back to cell-count tertiles. That is a
    # sample-quality-preserving permutation, and calling it collection-preserving
    # would be wrong - those rows also show n/a in the collection columns.
    marked = []
    for _, r in t.iterrows():
        marked.append("\u2020" if not r.strata_def.startswith("batch__") else "")
    t["dag"] = marked

    # Eleven columns, down from fifteen. Pre-submission review was right that the
    # earlier version was too dense to read across a landscape page: it carried the
    # column counts of both matrices and the V_D pair as well as the effects and the
    # tests. Those are secondary and are all in Table S4, per seed. What stays is one
    # effect and one test for each of the two metadata definitions, the expression
    # AUC, both permutation tests, and the verdict - plus the increment expression
    # adds over the recorded metadata, which is the quantity a reader of a
    # patient-level classifier is actually after and which the earlier table did not
    # report anywhere.
    head = ("| Comparison | Donors | Minority | "
            "Metadata AUC | p | Collection AUC | p | "
            "Expression AUC | + over metadata | p free | p strat. | Verdict |")
    sep = "|---|---:|---:|---:|---|---:|---|---:|---:|---|---|---|"
    lines = [head, sep]
    for _, r in t.iterrows():
        lines.append(
            f"| {r.Comparison} | {r.Donors} | {r.Minority} | "
            f"{r.all_auc} | {r.all_p} | {r.coll_auc} | {r.coll_p} | "
            f"{r.frozen} | {r.delta} | {r.p_free} | {r.p_coll}{r.dag} | {r.verdict} |")
    table = "\n".join(lines)
    OUT.write_text(table + "\n")

    # Splice it into the manuscript as well. Keeping a generated table beside a
    # hand-maintained copy of the same table is how the two drift apart, and the
    # copy in the text is the one a reader sees.
    if "--no-splice" in sys.argv:
        # verify.sh regenerates the table from the released per-seed outputs and
        # then checks it against the copy in the manuscript, so it must not write
        # that copy first.
        print(table)
        return
    ms_path = OUT.parent / "manuscript_v5.md"
    ms = ms_path.read_text()
    start = ms.index("| Comparison | Donors |")
    end = ms.index("\n\n", start)
    if ms[start:end].strip() != table:
        ms_path.write_text(ms[:start] + table + ms[end:])
        print("spliced the regenerated table into manuscript_v5.md")
    print(table)
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
