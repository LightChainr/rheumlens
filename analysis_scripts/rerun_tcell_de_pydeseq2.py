#!/usr/bin/env python3
"""PyDESeq2 sensitivity analysis on project marker-defined T-cell pseudobulk.

The cell labels are project-derived fixed marker-score labels, not original
author annotations. Both unadjusted and estimable adjusted designs are saved.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "data/secondary_features/GSE135779_project_marker_T_raw_count_sums.npz"
META = ROOT / "results/evidence_package/public_metadata/structured/GSE135779_ST1b_donor_clinical.csv"
MAP = ROOT / "results/evidence_package/public_metadata/structured/GSE135779_study_name_to_donor_id.csv"
OUT = ROOT / "results/evidence_package/tcell_de_pydeseq2"


def fit(counts, metadata, design, name):
    dds = DeseqDataSet(counts=counts, metadata=metadata, design=design,
                       refit_cooks=True, n_cpus=1, quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["disease", "SLE", "HD"],
                      n_cpus=1, quiet=True)
    stat.summary()
    result = stat.results_df.sort_values(["padj", "pvalue"], na_position="last")
    result.to_csv(OUT / f"{name}_all_genes.csv", index_label="gene")
    result.loc[result["padj"].fillna(1) < .05].to_csv(
        OUT / f"{name}_fdr05.csv", index_label="gene")
    return {
        "analysis": name,
        "design": design,
        "n_donors": len(metadata),
        "n_genes_tested": int(result["pvalue"].notna().sum()),
        "n_fdr05": int((result["padj"].fillna(1) < .05).sum()),
        "min_pvalue": float(result["pvalue"].min()),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # This NPZ is generated locally by marker_lineage_secondary.py; AnnData
    # gene names are serialized by NumPy as an object string array.
    z = np.load(COUNTS, allow_pickle=True)
    donors = z["donors"].astype(str)
    genes = z["genes"].astype(str)
    count = pd.DataFrame(z["counts"], index=donors, columns=genes)
    # Duplicate symbols cannot be separate DESeq2 variables; sum them.
    count = count.T.groupby(level=0).sum().T
    count = count.loc[:, count.sum(axis=0) >= 10]

    donor_map = pd.read_csv(MAP).drop_duplicates("donor_id").set_index("donor_id")
    clinical = pd.read_csv(META).set_index("Names")
    rows = []
    for donor in donors:
        study = donor_map.loc[donor, "study_name"]
        c = clinical.loc[study]
        rows.append({"donor_id": donor,
                     "disease": "SLE" if str(c["Groups"]).endswith("SLE") else "HD",
                     "Batch": str(c["Batch"]), "Age": float(c["Age"]),
                     "Gender": str(c["Gender"])})
    metadata = pd.DataFrame(rows).set_index("donor_id")
    metadata.to_csv(OUT / "analysis_metadata.csv")
    count.to_csv(OUT / "analysis_counts.csv.gz", compression="gzip", index_label="donor_id")

    summaries = [fit(count, metadata, "~disease", "disease_only")]
    # Batch, age and gender are the pre-specified covariates available for all
    # archived pediatric donors and form a full-rank design.
    summaries.append(fit(count, metadata, "~Batch + Age + Gender + disease",
                         "batch_age_gender_adjusted"))
    pd.DataFrame(summaries).to_csv(OUT / "summary.csv", index=False)
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
