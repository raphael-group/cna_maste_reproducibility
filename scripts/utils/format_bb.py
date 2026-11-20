#!/usr/bin/env python
# Reformat one UGP bb dir into the CalicoST-style input layout
#
# Runpeng Luo (2026-10-17)
#
# Dependencies:
#   numpy, pandas, scipy, anndata, h5py, pyarrow (CELLTYPING_ENV, configs/envs/celltyping.yaml)
#
# Usage:
#   python format_bb.py --bb-dir <bb> --assay <assay> --out-dir <out> \
#                       --spatial <rep>=<dir> [...] [--barcode-mappings <parquet>] [-v]
#           bb-dir            # UGP bb/ or bb-<XXX>um/ dir
#           assay             # UGP assay type, e.g. VISIUM3prime, VISIUM
#           out-dir           # output dir; must not exist
#           spatial           # Space Ranger spatial/ dir of one slice, keyed by UGP dataset id
#           barcode-mappings  # barcode_mappings.parquet copied into spatial/ (Visium HD only)
#
# Inputs
#   <bb>/unit/<assay>/barcodes.tsv.gz: UGP keys <barcode>_<rep>_<assay>, the matrix columns
#   <bb>/multi_snp/<assay>/bb.{A,B}allele.npz: multi-SNP x barcode allele counts
#   <bb>/multi_snp/<assay>/bb.tsv.gz: multi-SNP blocks (#CHR, START, END, ...)
#   <bb>/<assay>.h5ad: barcode x gene RNA counts
#   <spatial dir>: scalefactors_json.json, images, and tissue_positions.parquet or
#       tissue_positions_list.csv
# Parameters
#   None
# Outputs:
#   <out>/barcodes.txt: one barcode per row of cell_snp_*; <barcode>, or <barcode>_<rep>
#       when the bb dir holds more than one slice
#   <out>/cell_snp_{A,B}allele.npz: barcode x multi-SNP counts (csc), transpose of bb.*
#   <out>/unique_snp_ids.npy: <chr number>_<START>_<END> per multi-SNP block
#   <out>/filtered_feature_bc_matrix.h5: 10x h5 of the h5ad counts (UGP genes)
#   <out>/spatial/: slide images, scalefactors_json.json, tissue_positions.{parquet,csv}
#       and, if given, barcode_mappings.parquet
#   Multiple slices: the h5 and spatial/ go to <out>/<rep>/ with unsuffixed barcodes.
# Notes/References:
#   CalicoST loaders (load_data, load_joint_data, combine_gene_snps in utils_IO.py) read
#   unique_snp_ids with int(id.split("_")[0]) as chromosome and int(id.split("_")[1]) as
#   position, and join cell_snp_* rows, the h5 and the positions by barcode name; joint
#   mode expects <barcode>_<sample_id> in barcodes.txt.
#   https://github.com/raphael-group/CalicoST
#   multi_snp columns follow unit/ barcodes (UGP combine_counts_nonbulk.py).
#   10x h5 format: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices
#   tissue_positions_list.csv (Space Ranger < 2.0) has no header; it is written as
#   tissue_positions.parquet with the Space Ranger 2.0+ column names and types.
#   tissue_positions.csv holds the same table with a header row, the file CalicoST reads.
import argparse
import logging
import os
import shutil

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.sparse import load_npz, save_npz

log = logging.getLogger("format_bb")

POS_COLS = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]
POS_SCHEMA = pa.schema(
    [
        ("barcode", pa.string()),
        ("in_tissue", pa.uint8()),
        ("array_row", pa.uint32()),
        ("array_col", pa.uint32()),
        ("pxl_row_in_fullres", pa.float64()),
        ("pxl_col_in_fullres", pa.float64()),
    ]
)
POS_FILES = {
    "tissue_positions.parquet",
    "tissue_positions_list.csv",
    "tissue_positions.csv",
}


def parse_args():
    """Parse command-line arguments.

    Returns:
        argparse.Namespace with bb_dir, assay, out_dir, spatial (dict rep -> dir),
        barcode_mappings and verbose.

    Raises:
        SystemExit: A --spatial entry is not <rep>=<dir>.
    """
    p = argparse.ArgumentParser(description="Reformat one UGP bb dir into the CalicoST-style layout")
    p.add_argument("--bb-dir", required=True)
    p.add_argument("--assay", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--spatial", nargs="+", required=True)
    p.add_argument("--barcode-mappings")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()
    spatial = {}
    for s in args.spatial:
        rep, sep, d = s.partition("=")
        if not sep:
            p.error(f"--spatial {s}: expected <rep>=<dir>")
        spatial[rep] = d
    args.spatial = spatial
    return args


def split_keys(keys, assay):
    """Split UGP keys <barcode>_<rep>_<assay> into barcode and rep.

    Args:
        keys: UGP barcode keys.
        assay: Expected assay suffix.

    Returns:
        (barcodes, reps) as numpy string arrays.

    Raises:
        ValueError: A key does not end with _<assay>.
    """
    parts = [k.rsplit("_", 2) for k in keys]
    bad = [k for k, x in zip(keys, parts) if len(x) != 3 or x[2] != assay]
    if bad:
        raise ValueError(f"{len(bad)} keys without suffix _<rep>_{assay}, e.g. {bad[0]}")
    return np.array([x[0] for x in parts]), np.array([x[1] for x in parts])


def snp_ids(bb_tsv):
    """Build CalicoST-parsable ids <chr number>_<START>_<END> for multi-SNP blocks.

    Args:
        bb_tsv: multi_snp bb.tsv.gz path.

    Returns:
        Object numpy array of ids, one per block in file order.

    Raises:
        ValueError: A chromosome is not numeric after removing "chr", or ids repeat.
    """
    df = pd.read_csv(bb_tsv, sep="\t", usecols=["#CHR", "START", "END"])
    chrs = df["#CHR"].astype(str).str.removeprefix("chr")
    if not chrs.str.isdigit().all():
        raise ValueError(f"{bb_tsv}: non-numeric chromosomes {sorted(set(chrs[~chrs.str.isdigit()]))}")
    ids = chrs + "_" + df["START"].astype(str) + "_" + df["END"].astype(str)
    if ids.duplicated().any():
        raise ValueError(f"{bb_tsv}: duplicated block ids")
    return ids.to_numpy(dtype=object)


def write_10x_h5(path, X, barcodes, var, attrs):
    """Write a barcode x gene count matrix as a 10x feature-barcode h5.

    Args:
        path: Output h5 path.
        X: scipy csr matrix, barcodes x genes, integer-valued.
        barcodes: Barcode strings, one per row of X.
        var: Gene table with gene_ids, gene_symbol, feature_types and genome columns.
        attrs: Root attributes (chemistry_description, software_version, ...).

    Raises:
        ValueError: X holds non-integer values.
    """
    X = X.tocsr()
    X.sort_indices()
    if not np.array_equal(X.data, np.round(X.data)):
        raise ValueError("count matrix holds non-integer values")
    with h5py.File(path, "w") as f:
        for k, v in attrs.items():
            f.attrs[k] = v
        f.attrs["filetype"] = "matrix"
        f.attrs["version"] = 2
        m = f.create_group("matrix")
        m.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        m.create_dataset("data", data=X.data.astype(np.int32), compression="gzip")
        m.create_dataset("indices", data=X.indices.astype(np.int64), compression="gzip")
        m.create_dataset("indptr", data=X.indptr.astype(np.int64), compression="gzip")
        m.create_dataset("shape", data=np.array([X.shape[1], X.shape[0]], dtype=np.int32))
        g = m.create_group("features")
        g.create_dataset("_all_tag_keys", data=np.array(["genome"], dtype="S"))
        g.create_dataset("id", data=var["gene_ids"].astype(str).to_numpy(dtype="S"))
        g.create_dataset("name", data=var["gene_symbol"].astype(str).to_numpy(dtype="S"))
        g.create_dataset("feature_type", data=var["feature_types"].astype(str).to_numpy(dtype="S"))
        g.create_dataset("genome", data=var["genome"].astype(str).to_numpy(dtype="S"))


def read_positions(sp_dir):
    """Read a Space Ranger spatial/ tissue positions table in the 2.0+ schema.

    Args:
        sp_dir: Space Ranger spatial/ dir.

    Returns:
        pyarrow.Table with POS_SCHEMA columns.

    Raises:
        FileNotFoundError: No tissue positions file in sp_dir.
    """
    pq_path = os.path.join(sp_dir, "tissue_positions.parquet")
    if os.path.exists(pq_path):
        return pq.read_table(pq_path)
    for name, header in (
        ("tissue_positions.csv", 0),
        ("tissue_positions_list.csv", None),
    ):
        path = os.path.join(sp_dir, name)
        if os.path.exists(path):
            df = pd.read_csv(path, header=header, names=POS_COLS)
            return pa.Table.from_pandas(df, schema=POS_SCHEMA, preserve_index=False)
    raise FileNotFoundError(f"{sp_dir}: no tissue positions file")


def write_positions(pos, dst):
    """Write tissue_positions.parquet and tissue_positions.csv (with header) into dst.

    Args:
        pos: pyarrow.Table with POS_SCHEMA columns.
        dst: Output spatial/ dir.
    """
    pq.write_table(pos, os.path.join(dst, "tissue_positions.parquet"))
    pos.to_pandas().to_csv(os.path.join(dst, "tissue_positions.csv"), index=False)


def write_spatial(dst, sp_dir, barcode_mappings, barcodes):
    """Write spatial/: images and scalefactors copied, positions as parquet and csv.

    Hidden files (macOS ``._*`` AppleDouble copies) are not copied.

    Args:
        dst: Output spatial/ dir.
        sp_dir: Source Space Ranger spatial/ dir.
        barcode_mappings: barcode_mappings.parquet to copy, or None.
        barcodes: Unsuffixed barcodes of this slice, checked against in-tissue positions.

    Raises:
        ValueError: A barcode is not an in-tissue position.
    """
    os.makedirs(dst)
    for name in sorted(os.listdir(sp_dir)):
        if name not in POS_FILES and not name.startswith("."):
            shutil.copyfile(os.path.join(sp_dir, name), os.path.join(dst, name))
    pos = read_positions(sp_dir)
    df = pos.select(["barcode", "in_tissue"]).to_pandas()
    missing = set(barcodes) - set(df.loc[df["in_tissue"] == 1, "barcode"])
    if missing:
        raise ValueError(f"{sp_dir}: {len(missing)} barcodes not in tissue, e.g. {next(iter(missing))}")
    write_positions(pos, dst)
    if barcode_mappings:
        shutil.copyfile(barcode_mappings, os.path.join(dst, "barcode_mappings.parquet"))
    log.info(f"{dst}: {pos.num_rows} positions, {len(os.listdir(dst))} files")


def main():
    """Reformat one bb dir; see the file header for the layout."""
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    if os.path.exists(args.out_dir):
        raise FileExistsError(f"{args.out_dir} exists")
    unit = os.path.join(args.bb_dir, "unit", args.assay)
    multi = os.path.join(args.bb_dir, "multi_snp", args.assay)
    keys = pd.read_csv(os.path.join(unit, "barcodes.tsv.gz"), header=None, dtype=str)[0].to_numpy()
    raw, reps = split_keys(keys, args.assay)
    rep_ids = [str(r) for r in dict.fromkeys(reps)]
    if set(rep_ids) != set(args.spatial):
        raise ValueError(f"slices {rep_ids} in bb vs --spatial {sorted(args.spatial)}")
    joint = len(rep_ids) > 1
    out_bc = np.char.add(np.char.add(raw, "_"), reps) if joint else raw
    log.info(f"{args.bb_dir}: {len(keys)} barcodes, slices {rep_ids}, joint={joint}")

    mats = {}
    for allele in ("A", "B"):
        M = load_npz(os.path.join(multi, f"bb.{allele}allele.npz"))
        if M.shape[1] != len(keys):
            raise ValueError(f"bb.{allele}allele.npz has {M.shape[1]} columns for {len(keys)} barcodes")
        mats[allele] = M.T.tocsc()
    ids = snp_ids(os.path.join(multi, "bb.tsv.gz"))
    if len(ids) != mats["A"].shape[1]:
        raise ValueError(f"bb.tsv.gz has {len(ids)} blocks for {mats['A'].shape[1]} matrix rows")

    adata = ad.read_h5ad(os.path.join(args.bb_dir, f"{args.assay}.h5ad"))
    if set(adata.obs_names) != set(keys):
        raise ValueError(f"{args.assay}.h5ad barcodes differ from unit/ barcodes")
    adata = adata[keys]
    lib = rep_ids[0]
    meta = adata.uns.get("spatial", {}).get(lib, {}).get("metadata", {})
    attrs = {k: meta[k] for k in ("chemistry_description", "software_version") if k in meta}

    os.makedirs(args.out_dir)
    pd.Series(out_bc).to_csv(os.path.join(args.out_dir, "barcodes.txt"), header=False, index=False)
    for allele, M in mats.items():
        save_npz(os.path.join(args.out_dir, f"cell_snp_{allele}allele.npz"), M)
    np.save(os.path.join(args.out_dir, "unique_snp_ids.npy"), ids, allow_pickle=True)
    log.info(f"cell_snp_*: {mats['A'].shape} barcodes x multi-SNP blocks")

    for rep in rep_ids:
        rdir = os.path.join(args.out_dir, rep) if joint else args.out_dir
        os.makedirs(rdir, exist_ok=True)
        sel = reps == rep
        write_10x_h5(
            os.path.join(rdir, "filtered_feature_bc_matrix.h5"),
            adata.X[sel],
            raw[sel],
            adata.var,
            attrs,
        )
        log.info(f"{rep}: filtered_feature_bc_matrix.h5 {sel.sum()} x {adata.n_vars}")
        write_spatial(
            os.path.join(rdir, "spatial"),
            args.spatial[rep],
            args.barcode_mappings,
            raw[sel],
        )
    log.info(f"wrote {args.out_dir}")


if __name__ == "__main__":
    main()
