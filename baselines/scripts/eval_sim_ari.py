#!/usr/bin/env python
"""Score clone calls on the CalicoST simulations against truth_clone_labels.tsv.

Runpeng Luo (2026-10-30)

Dependencies: pandas, scikit-learn (calicost_env).

Usage:
    python -I eval_sim_ari.py <sim_dir> <bench_out> <deposit_numbat_dir> <out_tsv>
            sim_dir             # simulated_data_related/
            bench_out           # benchmark output root (numbat/, clonalscope/)
            deposit_numbat_dir  # the deposit's numbat_related/ (published numbat runs)
            out_tsv             # one row per (simulation, method)

Inputs:
    sim_dir/<sim>/truth_clone_labels.tsv: spot -> clone_<k> or normal.
    numbat: run_numbat/clone_post_<max>.tsv, clone_opt (clone 1 is numbat's normal clone).
    clonalscope_visium: <sim>/visium_tutorial/postprocess/<sim>.scRNA.cells.tsv (BARCODE, malignant, cluster);
        spots absent from it or with malignant N are scored as normal.
Outputs:
    out_tsv: sim, method, n_spots, n_pred_groups, ari, tumor_normal_acc.
Notes/References:
    ARI: Hubert & Arabie (1985), https://doi.org/10.1007/BF01908075 (sklearn adjusted_rand_score).
    Scored over every spot in truth_clone_labels.tsv, normal spots included.
"""

import glob
import os
import re
import sys

import pandas as pd
from sklearn.metrics import adjusted_rand_score


def numbat_labels(run_dir: str) -> pd.Series | None:
    """Final numbat clone per spot, or None if the run is unfinished.

    Args:
        run_dir: run_numbat() output directory.

    Returns:
        Series spot -> "normal" or "clone<k>".
    """
    posts = glob.glob(os.path.join(run_dir, "clone_post_*.tsv"))
    if not posts or not glob.glob(os.path.join(run_dir, "bulk_clones_final.tsv.gz")):
        return None
    last = max(posts, key=lambda p: int(re.search(r"clone_post_(\d+)", p).group(1)))
    df = pd.read_csv(last, sep="\t")
    return df.set_index("cell")["clone_opt"].map(lambda c: "normal" if c == 1 else f"clone{c}")


def clonalscope_labels(cells_tsv: str) -> pd.Series | None:
    """Clonalscope cluster per tumor spot, or None if missing.

    Args:
        cells_tsv: postprocessed <sim>.scRNA.cells.tsv.

    Returns:
        Series spot -> "normal" or "clone<k>" (tumor spots only).
    """
    if not os.path.exists(cells_tsv):
        return None
    df = pd.read_csv(cells_tsv, sep="\t")
    lab = df["cluster"].map(lambda c: f"clone{c}")
    lab[df["malignant"] != "T"] = "normal"
    return pd.Series(lab.values, index=df["BARCODE"])


def main(sim_dir: str, out: str, dep: str, out_tsv: str) -> None:
    """Score every simulation and method.

    Args:
        sim_dir: simulated_data_related directory.
        out: Benchmark output root.
        dep: Deposit numbat_related directory.
        out_tsv: Output table.
    """
    rows = []
    for sim in sorted(os.listdir(sim_dir)):
        truth_f = os.path.join(sim_dir, sim, "truth_clone_labels.tsv")
        if not os.path.exists(truth_f):
            continue
        truth = pd.read_csv(truth_f, sep="\t", index_col=0)["labels"]
        preds = {
            "numbat": numbat_labels(f"{out}/numbat/{sim}/max_entropy_08/run_numbat"),
            "numbat_deposit": numbat_labels(f"{dep}/{sim}/outs"),
            "clonalscope_visium": clonalscope_labels(
                f"{out}/clonalscope/{sim}/visium_tutorial/postprocess/{sim}.scRNA.cells.tsv"
            ),
        }
        for method, p in preds.items():
            if p is None:
                continue
            p = p.reindex(truth.index).fillna("normal")
            rows.append(
                {
                    "sim": sim,
                    "method": method,
                    "n_spots": len(truth),
                    "n_pred_groups": p.nunique(),
                    "ari": adjusted_rand_score(truth, p),
                    "tumor_normal_acc": ((truth == "normal") == (p == "normal")).mean(),
                }
            )
    res = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    res.to_csv(out_tsv, sep="\t", index=False)
    print(res.groupby("method")[["ari", "tumor_normal_acc", "n_pred_groups"]].agg(["count", "median", "mean"]))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit("usage: eval_sim_ari.py <sim_dir> <bench_out> <deposit_numbat_dir> <out_tsv>")
    main(*sys.argv[1:])
