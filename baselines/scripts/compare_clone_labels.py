#!/usr/bin/env python3
"""ARI between a clone labelling and a 16um reference labelling, after expanding coarse bins to 16um.

Author: Runpeng Luo (2026-10-25)

A coarse bin s_<B>um_<R>_<C>-1 covers the 16um bins (R*k + i, C*k + j), i, j in 0..k-1, k = B / 16
(the summing rule of scripts/common/aggregate_10x_matrix.R). Each covered 16um bin inherits the coarse bin's
label. ARI is computed on the 16um bins present in both labellings; bins present in only one are
counted and reported, not scored.

Dependencies: pandas, scikit-learn.

Usage:
  python compare_clone_labels.py --pred clone_labels.tsv --ref clone_labels.tsv.gz --out ari.tsv
      --pred        # prediction table; first column barcode, column `clone_label`
      --ref         # reference table; column `barcode` (or first column), column `clone_label`
      --out         # one-row summary TSV
      --contingency # optional: reference x prediction count table TSV

Outputs:
  out: pred, ref, pred_bin_um, n_ref, n_pred_16um, n_shared, n_ref_only, n_pred_only, k_pred, k_ref, ari.
  contingency: rows = reference clones, columns = predicted clones, over shared 16um bins.
"""

import argparse
import re
import sys

import pandas as pd
from sklearn.metrics import adjusted_rand_score

BC = re.compile(r"^s_(\d{3})um_(\d+)_(\d+)-1$")


def read_labels(path):
    """Read a label table and return barcode -> clone_label as str."""
    df = pd.read_csv(path, sep="\t")
    bc_col = "barcode" if "barcode" in df.columns else df.columns[0]
    if "clone_label" not in df.columns:
        sys.exit(f"{path}: no clone_label column")
    return pd.Series(df["clone_label"].astype(str).values, index=df[bc_col].astype(str).values)


def to_16um(labels):
    """Expand coarse-bin labels to 16um bins; 16um labels pass through.

    Returns:
        (labels indexed by 16um barcode, bin size in um of the input).
    """
    sizes = {int(BC.match(b).group(1)) for b in labels.index if BC.match(b)}
    if len(sizes) != 1 or len(labels) != sum(1 for b in labels.index if BC.match(b)):
        sys.exit(f"barcodes are not all s_<bin>um_<row>_<col>-1 of one bin size: {sorted(sizes)}")
    size = sizes.pop()
    if size % 16:
        sys.exit(f"bin size {size} is not a multiple of 16")
    k = size // 16
    rows = []
    for b, lab in labels.items():
        _, r, c = BC.match(b).groups()
        r, c = int(r), int(c)
        for i in range(k):
            for j in range(k):
                rows.append((f"s_016um_{r * k + i:05d}_{c * k + j:05d}-1", lab))
    out = pd.Series(dict(rows))
    return out, size


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pred", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--contingency")
    a = ap.parse_args()

    pred, size = to_16um(read_labels(a.pred))
    ref, ref_size = to_16um(read_labels(a.ref))
    if ref_size != 16:
        sys.exit(f"reference is {ref_size}um, expected 16um")
    shared = ref.index.intersection(pred.index)
    ari = adjusted_rand_score(ref[shared], pred[shared])
    row = {
        "pred": a.pred,
        "ref": a.ref,
        "pred_bin_um": size,
        "n_ref": len(ref),
        "n_pred_16um": len(pred),
        "n_shared": len(shared),
        "n_ref_only": len(ref.index.difference(pred.index)),
        "n_pred_only": len(pred.index.difference(ref.index)),
        "k_pred": pred[shared].nunique(),
        "k_ref": ref[shared].nunique(),
        "ari": round(ari, 4),
    }
    pd.DataFrame([row]).to_csv(a.out, sep="\t", index=False)
    print(pd.DataFrame([row]).drop(columns=["pred", "ref"]).to_string(index=False))
    tab = pd.crosstab(ref[shared].rename("ref"), pred[shared].rename("pred"))
    print(tab.to_string())
    if a.contingency:
        tab.to_csv(a.contingency, sep="\t")


if __name__ == "__main__":
    main()
