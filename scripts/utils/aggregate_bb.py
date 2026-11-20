#!/usr/bin/env python
# Sum one Visium HD UGP bb dir into coarser square bins, written as a bb dir
#
# Runpeng Luo (2026-11-16)
#
# Dependencies:
#   numpy, pandas, scipy, anndata (CELLTYPING_ENV, configs/envs/celltyping.yaml)
#
# Usage:
#   python aggregate_bb.py --bb-dir <bb> --assay <assay> --bin-um <X> --out-dir <out> [-v]
#           bb-dir   # UGP bb/ dir of square bins, keys s_<side>um_<row>_<col>-1_<rep>_<assay>
#           assay    # UGP assay type, e.g. VISIUM3prime
#           bin-um   # target bin side in um, a multiple of the bb dir's side above it, e.g. 32 for 16um bins
#           out-dir  # output dir; must not exist
#
# Inputs
#   <bb>/unit/<assay>/barcodes.tsv.gz: UGP keys, the allele matrices' columns
#   <bb>/multi_snp/<assay>/bb.{A,B,T}allele.npz: multi-SNP x bin allele counts (T where present)
#   <bb>/multi_snp/<assay>/bb.tsv.gz: multi-SNP blocks
#   <bb>/<assay>.h5ad: bin x gene RNA counts
# Parameters
#   bin_um: X; a bin (r, c) of side s goes to the coarse bin (r // k, c // k), k = X / s
# Outputs:
#   <out>/unit/<assay>/barcodes.tsv.gz: coarse keys s_<XXX>um_<row:05d>_<col:05d>-1_<rep>_<assay>, sorted
#   <out>/multi_snp/<assay>/bb.{A,B,T}allele.npz: multi-SNP x coarse bin, summed (csr)
#   <out>/multi_snp/<assay>/bb.tsv.gz: the blocks, copied
#   <out>/<assay>.h5ad: coarse bin x gene, summed, in the barcodes' order; obs array_row, array_col and
#       n_bins_<side>um, the bins summed into each
# Notes/References:
#   The binning rule of scripts/baselines/numbat/aggregate_bins.R and common/aggregate_10x_matrix.R, so the coarse
#   bins are those Numbat's and Clonalscope's runs read. A coarse bin holds the in-tissue bins of the bb dir only: its
#   counts are theirs, and a coarse bin is in tissue where one of them is.
#   multi_snp columns follow unit/ barcodes (UGP combine_counts_nonbulk.py); the h5ad's rows need not.
#   Only the files above are written: what plots/umi.py reads of a bb dir, and the unit barcodes naming its columns.
import argparse
import logging
import os
import re
import shutil

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger("aggregate_bb")

KEY = re.compile(r"^s_(\d{3})um_(\d+)_(\d+)-1_(.+)$")
"""A UGP key of a square bin: its side, row, column, then `<rep>_<assay>`."""


def parse_args():
    """Parse command-line arguments.

    Returns:
        argparse.Namespace with bb_dir, assay, bin_um, out_dir and verbose.
    """
    p = argparse.ArgumentParser(description="Sum one Visium HD UGP bb dir into coarser square bins")
    p.add_argument("--bb-dir", required=True)
    p.add_argument("--assay", required=True)
    p.add_argument("--bin-um", required=True, type=int)
    p.add_argument("--out-dir", required=True)
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args()


def get_coarse_bins(keys, bin_um):
    """Map each bin's UGP key to its coarse bin's.

    Args:
        keys: UGP keys s_<side>um_<row>_<col>-1_<rep>_<assay>, one side.
        bin_um: the coarse bins' side in um.

    Returns:
        (coarse, side): each key's coarse key, as a numpy string array, and the keys' side in um.

    Raises:
        ValueError: a key is not a square bin's, the keys hold several sides, or `bin_um` is not a multiple of
            their side above it.
    """
    parts = [KEY.match(k) for k in keys]
    bad = [k for k, m in zip(keys, parts) if m is None]
    if bad:
        raise ValueError(f"{len(bad)} keys not s_<side>um_<row>_<col>-1_<rep>_<assay>, e.g. {bad[0]}")
    sides = {int(m.group(1)) for m in parts}
    if len(sides) != 1:
        raise ValueError(f"keys of several bin sides: {sorted(sides)} um")
    side = sides.pop()
    if bin_um <= side or bin_um % side:
        raise ValueError(f"--bin-um {bin_um}: not a multiple of the bins' {side} um above it")
    k = bin_um // side
    coarse = [f"s_{bin_um:03d}um_{int(m.group(2)) // k:05d}_{int(m.group(3)) // k:05d}-1_{m.group(4)}" for m in parts]
    return np.array(coarse), side


def get_indicator(coarse, levels):
    """The bins x coarse bins 0/1 matrix: each bin's one coarse bin.

    Args:
        coarse: each bin's coarse key.
        levels: the coarse keys, in output order.

    Returns:
        scipy csr matrix, (bins, coarse bins), int64.
    """
    column = pd.Index(levels).get_indexer(coarse)
    return sp.csr_matrix(
        (np.ones(len(coarse), dtype=np.int64), (np.arange(len(coarse)), column)), shape=(len(coarse), len(levels))
    )


def get_coarse_obs(levels, counts, side):
    """The coarse h5ad's obs: each coarse bin's row, column and number of bins summed.

    Args:
        levels: the coarse keys.
        counts: the bins summed into each.
        side: the summed bins' side in um.

    Returns:
        pandas DataFrame indexed by the keys: array_row, array_col, n_bins_<side>um.
    """
    rows, cols = zip(*((int(m.group(2)), int(m.group(3))) for m in map(KEY.match, levels)))
    return pd.DataFrame(
        {"array_row": rows, "array_col": cols, f"n_bins_{side:03d}um": np.asarray(counts).ravel()},
        index=pd.Index(levels, dtype=object),
    )


def main():
    """Sum one UGP bb dir into coarser bins and write it as a bb dir."""
    args = parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(message)s")
    if os.path.exists(args.out_dir):
        raise FileExistsError(f"{args.out_dir} exists")
    unit = os.path.join(args.bb_dir, "unit", args.assay)
    multi = os.path.join(args.bb_dir, "multi_snp", args.assay)

    keys = pd.read_csv(os.path.join(unit, "barcodes.tsv.gz"), header=None, dtype=str)[0].to_numpy()
    coarse, side = get_coarse_bins(keys, args.bin_um)
    levels = np.sort(np.unique(coarse))
    to_coarse = get_indicator(coarse, levels)
    log.info(f"{args.bb_dir}: {len(keys):,} {side}um bins -> {len(levels):,} {args.bin_um}um bins")

    tmp = f"{args.out_dir}.tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(os.path.join(tmp, "unit", args.assay))
    os.makedirs(os.path.join(tmp, "multi_snp", args.assay))
    pd.Series(levels).to_csv(
        os.path.join(tmp, "unit", args.assay, "barcodes.tsv.gz"), header=False, index=False, compression="gzip"
    )

    for allele in ("A", "B", "T"):
        path = os.path.join(multi, f"bb.{allele}allele.npz")
        if not os.path.exists(path):
            continue
        fine = sp.load_npz(path).tocsr()
        if fine.shape[1] != len(keys):
            raise ValueError(f"{path} has {fine.shape[1]} columns for {len(keys)} barcodes")
        summed = (fine @ to_coarse.astype(fine.dtype)).tocsr()
        if summed.sum() != fine.sum():
            raise ValueError(f"{path}: the coarse counts' total differs from the bins'")
        sp.save_npz(os.path.join(tmp, "multi_snp", args.assay, f"bb.{allele}allele.npz"), summed)
        log.info(f"bb.{allele}allele.npz: {fine.shape} -> {summed.shape}")
    shutil.copy2(os.path.join(multi, "bb.tsv.gz"), os.path.join(tmp, "multi_snp", args.assay, "bb.tsv.gz"))

    adata = ad.read_h5ad(os.path.join(args.bb_dir, f"{args.assay}.h5ad"))
    if set(adata.obs_names) != set(keys):
        raise ValueError(f"{args.assay}.h5ad barcodes differ from unit/ barcodes")
    # NB: the h5ad's rows are in its own order, so they get their own indicator over the same coarse bins
    rows = get_indicator(get_coarse_bins(adata.obs_names.to_numpy(), args.bin_um)[0], levels)
    X = sp.csr_matrix(adata.X)
    summed = (rows.T.astype(X.dtype) @ X).tocsr()
    if not np.isclose(summed.sum(), X.sum()):
        raise ValueError(f"{args.assay}.h5ad: the coarse counts' total differs from the bins'")
    var = adata.var.copy()
    # NB: anndata writes pandas string arrays only on opt-in; object strings write on every version
    var.index = var.index.astype(object)
    var = var.astype({c: object for c in var.columns if pd.api.types.is_string_dtype(var[c])})
    out = ad.AnnData(X=summed, obs=get_coarse_obs(levels, rows.sum(axis=0), side), var=var)
    out.write_h5ad(os.path.join(tmp, f"{args.assay}.h5ad"))
    log.info(f"{args.assay}.h5ad: {adata.shape} -> {out.shape}")

    os.rename(tmp, args.out_dir)
    log.info(f"wrote {args.out_dir}")


if __name__ == "__main__":
    main()
