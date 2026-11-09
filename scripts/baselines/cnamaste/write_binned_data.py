"""Write a finished cnamaste run's per-spot bins as `binned_data.npz`, read back from its inputs.

Runpeng Luo (2026-11-04)

Dependencies: a cnamaste checkout with `cnamaste.io.get_labelled_counts` (e.g. `cleanup3`), numpy, pandas, scipy.

Usage:
    uv run --project <cnamaste checkout> python write_binned_data.py <run_dir> <out_npz> [--input_dir DIR]
            run_dir    # the run's run_<id>/ directory
            out_npz    # the file written, e.g. <run_dir>/binned_data.npz

Inputs:
    run_dir: config_<id>.yaml, segmentation.tsv.gz, snp_segmentation.tsv.gz, clone_labels.tsv.gz and
        rdrbaf_final_nstates<n>_smp.npz; the config's sample sheet and references.
Parameters:
    input_dir: the directory the config's relative paths resolve against, the current directory by default.
Outputs:
    out_npz: `lengths`, bins per contig; `barcodes`, the spots in the fit's order; `single_X__rdr__*` (umis),
        `single_X__baf__*` (alt_snp_umis) and `single_total_bb_RD__*` (total_snp_umis), each the CSR parts `data`,
        `indices`, `indptr`, `shape` of its bin-major `(n_bins, n_spots)` counts, or a dense entry where the counts are
        dense; `single_base_nb_mean`, the dense `(n_bins, n_spots)` normal baseline.
Notes/References:
    cnamaste fdd0af9 (T- #107) removed `binned_data.npz`: "results.get_labelled_counts reads the run's final counts
    back from its inputs, segmentation.tsv.gz and snp_segmentation.tsv.gz, so binned_data.npz is no longer written".
    This writes the file that commit removed, in the layout of its parent's `results.write_binned_data` and
    `spatio_genomic_counts.npz_entries` (a8083b9), which the benchmark's `plots.cnamaste_page.read_clone_tracks` reads.
    `get_labelled_counts` orders the spots as `load_input_data` returns them and does not return their barcodes; the
    barcodes are taken from that call, and the run's `clone_labels.tsv.gz`, in that order, must equal the fit's
    `new_assignment`, so the spots are the fit's.
"""

import argparse
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import scipy.sparse as sp

import cnamaste.io

CHANNELS = {"single_X__rdr": "umis", "single_X__baf": "alt_snp_umis", "single_total_bb_RD": "total_snp_umis"}
"""Each `binned_data.npz` entry, by the counts field it holds."""


def get_counts_and_barcodes(run_dir: Path, input_dir: Path) -> tuple[object, np.ndarray]:
    """Get the run's labelled counts (`cnamaste.io.get_labelled_counts`) and the barcodes of their spots, in order.

    Args:
        run_dir: the run's `run_<id>` directory.
        input_dir: the directory the config's relative paths resolve against.

    Returns:
        The `SpatioGenomicLabelledMultiCounts` and its spots' barcodes.
    """
    loaded = {}
    load_input_data = cnamaste.io.load_input_data

    def record_barcodes(*args: object, **kwargs: object) -> tuple:
        result = load_input_data(*args, **kwargs)
        loaded["barcodes"] = np.asarray(result[1], dtype=str)
        return result

    with patch.object(cnamaste.io, "load_input_data", record_barcodes):
        counts = cnamaste.io.get_labelled_counts(input_dir, run_dir)
    return counts, loaded["barcodes"]


def check_spot_order(run_dir: Path, barcodes: np.ndarray) -> None:
    """Check the spots are the fit's, in its order: `clone_labels.tsv.gz`'s `clone_label` at `barcodes` equals the
    final fit's `new_assignment`.

    Args:
        run_dir: the run's `run_<id>` directory.
        barcodes: the counts' spots.

    Raises:
        ValueError: where the run holds other than one final fit, or the labels differ.
    """
    fits = sorted(run_dir.glob("rdrbaf_final_nstates*_smp.npz"))
    if len(fits) != 1:
        raise ValueError(f"{run_dir} holds {len(fits)} rdrbaf_final_nstates*_smp.npz, one wanted")
    with np.load(fits[0], allow_pickle=True) as fit:
        assignment = np.asarray(fit["new_assignment"])
    labels = pd.read_csv(run_dir / "clone_labels.tsv.gz", sep="\t", index_col="barcode")["clone_label"]
    ordered = labels.reindex(barcodes).to_numpy()
    if len(ordered) != len(assignment) or not np.array_equal(ordered, assignment):
        raise ValueError(f"{run_dir}: clone_labels.tsv.gz at the counts' spots is not {fits[0].name}'s new_assignment")


def get_entries(name: str, values: object) -> dict[str, np.ndarray]:
    """Get one bin-major `(n_bins, n_spots)` array's `binned_data.npz` entries: its CSR parts where sparse, else itself.

    Args:
        name: the entry, e.g. `single_X__rdr`.
        values: the array, CSR or dense.

    Returns:
        The entries by key.
    """
    if not sp.issparse(values):
        return {name: np.ascontiguousarray(values)}
    csr = sp.csr_matrix(values)
    csr.sum_duplicates()
    csr.sort_indices()
    return {
        f"{name}__data": csr.data,
        f"{name}__indices": csr.indices.astype(np.int64),
        f"{name}__indptr": csr.indptr.astype(np.int64),
        f"{name}__shape": np.asarray(csr.shape, dtype=np.int64),
    }


def get_binned_entries(counts: object, barcodes: np.ndarray) -> dict[str, np.ndarray]:
    """Get `binned_data.npz`'s entries of the counts, bin-major: `single_X`'s two channels, `single_total_bb_RD` and
    the dense `single_base_nb_mean`; a dense `single_X` as CalicoST's one `(n_bins, 2, n_spots)` entry.

    Args:
        counts: a `SpatioGenomicLabelledMultiCounts`.
        barcodes: its spots'.

    Returns:
        The entries by key, `lengths` and `barcodes` among them.
    """
    bin_major = {name: getattr(counts, field).get_bin_major() for name, field in CHANNELS.items()}
    entries = {}
    if sp.issparse(bin_major["single_X__rdr"]):
        for name, values in bin_major.items():
            entries |= get_entries(name, values)
    else:
        entries["single_X"] = np.stack([bin_major["single_X__rdr"], bin_major["single_X__baf"]], axis=1)
        entries |= get_entries("single_total_bb_RD", bin_major["single_total_bb_RD"])
    baseline = np.ascontiguousarray(np.asarray(counts.normal_baseline_sized.toarray(), dtype=np.float64).T)
    return {
        "lengths": np.asarray(counts.lengths, dtype=np.int64),
        "barcodes": np.asarray(barcodes, dtype=str),
        **entries,
        "single_base_nb_mean": baseline,
    }


def main() -> None:
    """Write the run's `binned_data.npz`."""
    command = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    command.add_argument("run_dir", type=Path)
    command.add_argument("out_npz", type=Path)
    command.add_argument("--input_dir", type=Path, default=Path.cwd())
    args = command.parse_args()

    counts, barcodes = get_counts_and_barcodes(args.run_dir, args.input_dir)
    if len(barcodes) != counts.n_spots or int(np.sum(counts.lengths)) != counts.n_segments:
        raise ValueError(f"{len(barcodes):_} barcodes, {counts.n_spots:_} spots, lengths {counts.lengths}")
    check_spot_order(args.run_dir, barcodes)
    entries = get_binned_entries(counts, barcodes)
    np.savez(args.out_npz, **entries)
    print(f"Wrote {args.out_npz}: {counts.n_spots:_} spots, {counts.n_segments:_} bins, {sorted(entries)}")


if __name__ == "__main__":
    main()
