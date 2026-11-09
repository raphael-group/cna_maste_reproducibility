#!/usr/bin/env python
"""Compare rerun CalicoST outputs with the deposit's published CalicoST outputs, per initialization.

Runpeng Luo (2026-11-02)

Dependencies: numpy, pandas, scikit-learn (calicost_env).

Usage:
    python -I calicost/compare_repro.py <deposit_calicost_dir> <repro_dir> <out_tsv> <sim> [<sim> ...]
            deposit_calicost_dir  # CalicoST_simulation_deposit/nomixing_calicost_related
            repro_dir             # calicost/repro.sh out_root
            out_tsv               # one row per (sim, rectangle)
            sim                   # simulation ids

Inputs:
    <dir>/<sim>/clone3_rectangle<r>_w1.0/{clone_labels.tsv, cnv_seglevel.tsv, cnv_genelevel.tsv,
        rdrbaf_final_nstates7_smp.npz} for both the deposit and the rerun.
Outputs:
    out_tsv: sim, rectangle, best_dep, best_rep (rectangle chosen by argmax total_llf on each side),
        total_llf_dep, total_llf_rep, n_spots, label_identical, label_ari, seg_identical,
        seg_rows_dep, seg_rows_rep, seg_cn_mismatch, gene_identical, gene_cn_mismatch.
Notes/References:
    Deposit: https://doi.org/10.5281/zenodo.14175627; best initialization as get_best_r_hmrf() in
    utils/Simulation2.ipynb.
    label_identical: same barcodes and same clone_label for every spot (label ids, not up to permutation).
    seg_identical: same (CHR, START, END) rows and same integer A/B per clone; seg_cn_mismatch counts
    differing cells when the segment rows match, else NA.
"""

import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

RECTS = range(5)


def rect_dir(root: str, sim: str, r: int) -> str:
    """Path of one initialization's output directory.

    Args:
        root: CalicoST output root holding <sim>/.
        sim: Simulation id.
        r: Rectangle initialization index.

    Returns:
        Directory path.
    """
    return f"{root}/{sim}/clone3_rectangle{r}_w1.0"


def total_llf(d: str) -> float:
    """total_llf of one initialization.

    Args:
        d: Initialization output directory.

    Returns:
        total_llf, or NaN if the npz is missing.
    """
    f = f"{d}/rdrbaf_final_nstates7_smp.npz"
    if not os.path.exists(f):
        return float("nan")
    return float(np.load(f, allow_pickle=True)["total_llf"])


def cmp_table(a: pd.DataFrame, b: pd.DataFrame) -> tuple[bool, object]:
    """Compare two integer copy-number tables.

    Args:
        a: Deposit table.
        b: Rerun table.

    Returns:
        (identical, number of differing cells or "NA" if shape/index/columns differ).
    """
    if a.shape != b.shape or not a.index.equals(b.index) or list(a) != list(b):
        return False, "NA"
    n = int((a.values != b.values).sum())
    return n == 0, n


def main(dep: str, rep: str, out_tsv: str, sims: list[str]) -> None:
    """Write one comparison row per simulation and initialization.

    Args:
        dep: Deposit nomixing_calicost_related directory.
        rep: Rerun output root.
        out_tsv: Output table.
        sims: Simulation ids.
    """
    rows = []
    for sim in sims:
        ll = {
            side: {r: total_llf(rect_dir(root, sim, r)) for r in RECTS} for side, root in (("dep", dep), ("rep", rep))
        }
        best = {s: max(v, key=lambda r: -np.inf if np.isnan(v[r]) else v[r]) for s, v in ll.items()}
        for r in RECTS:
            dd, rd = rect_dir(dep, sim, r), rect_dir(rep, sim, r)
            row = {
                "sim": sim,
                "rectangle": r,
                "best_dep": best["dep"],
                "best_rep": best["rep"],
                "total_llf_dep": ll["dep"][r],
                "total_llf_rep": ll["rep"][r],
            }
            if not os.path.exists(f"{rd}/cnv_genelevel.tsv"):
                rows.append(row)
                continue
            la = pd.read_csv(f"{dd}/clone_labels.tsv", sep="\t", index_col=0)["clone_label"]
            lb = pd.read_csv(f"{rd}/clone_labels.tsv", sep="\t", index_col=0)["clone_label"]
            common = la.index.intersection(lb.index)
            row["n_spots"] = f"{len(la)}/{len(lb)}"
            row["label_identical"] = la.index.equals(lb.index) and bool((la == lb).all())
            row["label_ari"] = round(adjusted_rand_score(la[common], lb[common]), 6)
            sa = pd.read_csv(f"{dd}/cnv_seglevel.tsv", sep="\t").set_index(["CHR", "START", "END"])
            sb = pd.read_csv(f"{rd}/cnv_seglevel.tsv", sep="\t").set_index(["CHR", "START", "END"])
            row["seg_identical"], row["seg_cn_mismatch"] = cmp_table(sa, sb)
            row["seg_rows_dep"], row["seg_rows_rep"] = len(sa), len(sb)
            ga = pd.read_csv(f"{dd}/cnv_genelevel.tsv", sep="\t", index_col=0)
            gb = pd.read_csv(f"{rd}/cnv_genelevel.tsv", sep="\t", index_col=0)
            row["gene_identical"], row["gene_cn_mismatch"] = cmp_table(ga, gb)
            rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(out_tsv, sep="\t", index=False, na_rep="NA")
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(res.to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.exit("usage: calicost/compare_repro.py <deposit_calicost_dir> <repro_dir> <out_tsv> <sim> ...")
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
