#!/usr/bin/env python3
"""One rename map from analysis column names to the names used in the manuscript.

The analysis code still calls the residual label variance `I_D` internally. Renaming
it in the code would touch every result file and every checked-in output, so the
translation happens once, here, on the way into the released supplementary tables.
Everything a reader sees uses the manuscript's vocabulary; nothing in the analysis
depends on this map.
"""
from __future__ import annotations
import re

RENAME = {
    "I_D": "V_D_insample",
    "I_D_cv": "V_D_crossfitted",
    "I_D_null_mean": "V_D_insample_null_mean",
    "I_D_cv_null_mean": "V_D_crossfitted_null_mean",
    "I_D_cv_null_p025": "V_D_crossfitted_null_p025",
    "I_D_vs_null": "V_D_insample_minus_null",
    "p_I_D": "p_V_D",
    "p_I_D_insample": "p_V_D_insample",
    "p_standard": "p_free",
    "reject_unstratified": "reject_free",
    "p_unstratified": "p_free",
    "design_c": "design_cindex",
    "disease_c": "diagnosis_cindex",
    "residualised_c": "residualised_cindex",
    "restricted_c": "restricted_cindex",
    "external_c": "external_cindex",
    "sig": "fraction_flagged",
    "variable_group": "variable_group",
}
# Retired words that must not appear in a released table header or cell label.
RETIRED = ["unstratified", "stratified", "design-preserving", "design_preserving",
           "I_D", "p_standard"]


def rename_columns(cols):
    """Apply the map to a list of column names, including `<name>_mean/_std` suffixes."""
    out = []
    for c in cols:
        if c in RENAME:
            out.append(RENAME[c]); continue
        m = re.fullmatch(r"(.+?)_(mean|std)", c)
        if m and m.group(1) in RENAME:
            out.append(f"{RENAME[m.group(1)]}_{m.group(2)}"); continue
        out.append(c)
    return out


def check_header(cols) -> list[str]:
    """Return any header entries still carrying a retired word."""
    bad = []
    for c in cols:
        for w in RETIRED:
            if w == "stratified" and "collection-preserving" in c:
                continue
            if re.search(rf"(?<![A-Za-z]){re.escape(w)}(?![A-Za-z])", c):
                bad.append(f"{c} (contains '{w}')")
    return bad
