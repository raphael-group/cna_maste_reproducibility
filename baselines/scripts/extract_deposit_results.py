#!/usr/bin/env python
"""Copy the clone-label and copy-number profile files of the deposit's CalicoST, Numbat and STARCH runs.

Runpeng Luo (2026-11-02)

Dependencies: numpy (calicost_env).

Usage:
    python -I extract_deposit_results.py <deposit_dir> <out_dir>
            deposit_dir  # unzipped CalicoST_simulation_deposit/
            out_dir      # e.g. <deposit_dir>/method_results

Inputs:
    nomixing_calicost_related/<sim>/clone3_rectangle<r>_w1.0/{clone_labels.tsv, cnv_genelevel.tsv,
        cnv_seglevel.tsv, rdrbaf_final_nstates*_smp.npz}
    numbat_related/<sim>/outs/{clone_post_2.tsv, bulk_clones_final.tsv.gz}
    starch_related/<sim>/{labels_STITCH_output.csv, states_STITCH_output.csv}
Outputs:
    out_dir/calicost/<sim>/{clone_labels.tsv, cnv_genelevel.tsv, cnv_seglevel.tsv}
    out_dir/calicost/best_init.tsv: sim, rectangle, total_llf (the chosen initialization).
    out_dir/numbat/<sim>/{clone_post_2.tsv, bulk_clones_final.tsv.gz} (only runs with output).
    out_dir/starch/<sim>/{labels_STITCH_output.csv, states_STITCH_output.csv}
Notes/References:
    Deposit: https://doi.org/10.5281/zenodo.14175627 (utils/Simulation2.ipynb reads these files).
    CalicoST initialization = argmax total_llf over clone3_rectangle*, as get_best_r_hmrf() in that notebook.
    Files are copied with metadata (shutil.copy2); the deposit is not modified.
"""

import glob
import os
import shutil
import sys

import numpy as np

FILES = {
    "numbat": ("outs", ["clone_post_2.tsv", "bulk_clones_final.tsv.gz"]),
    "starch": ("", ["labels_STITCH_output.csv", "states_STITCH_output.csv"]),
}
CALICOST_FILES = ["clone_labels.tsv", "cnv_genelevel.tsv", "cnv_seglevel.tsv"]


def best_init(sim_dir: str) -> tuple[int, float]:
    """Pick the CalicoST initialization with the highest total_llf.

    Args:
        sim_dir: nomixing_calicost_related/<sim> directory.

    Returns:
        (rectangle index, total_llf).

    Raises:
        FileNotFoundError: no rdrbaf_final npz is present.
    """
    ll = {}
    for f in glob.glob(f"{sim_dir}/clone3_rectangle*_w1.0/rdrbaf_final_nstates*_smp.npz"):
        r = int(f.split("rectangle")[1].split("_")[0])
        ll[r] = float(np.load(f, allow_pickle=True)["total_llf"])
    if not ll:
        raise FileNotFoundError(f"no rdrbaf_final npz under {sim_dir}")
    r = max(ll, key=ll.get)
    return r, ll[r]


def copy_files(src_dir: str, names: list[str], dst_dir: str) -> int:
    """Copy names from src_dir to dst_dir if all exist.

    Args:
        src_dir: Source directory.
        names: File names to copy.
        dst_dir: Destination directory.

    Returns:
        Number of files copied (0 if any is missing).
    """
    srcs = [os.path.join(src_dir, n) for n in names]
    if not all(os.path.exists(s) for s in srcs):
        return 0
    os.makedirs(dst_dir, exist_ok=True)
    for s in srcs:
        shutil.copy2(s, dst_dir)
    return len(srcs)


def main(dep: str, out: str) -> None:
    """Copy every simulation's result files.

    Args:
        dep: Deposit directory.
        out: Output directory.
    """
    sims = sorted(os.listdir(f"{dep}/nomixing_calicost_related"))
    rows = ["sim\trectangle\ttotal_llf"]
    n = {"calicost": 0, "numbat": 0, "starch": 0}
    for sim in sims:
        cdir = f"{dep}/nomixing_calicost_related/{sim}"
        r, llf = best_init(cdir)
        rows.append(f"{sim}\t{r}\t{llf}")
        n["calicost"] += (
            copy_files(
                f"{cdir}/clone3_rectangle{r}_w1.0",
                CALICOST_FILES,
                f"{out}/calicost/{sim}",
            )
            > 0
        )
        for m, (sub, names) in FILES.items():
            src = os.path.join(dep, f"{m}_related", sim, sub)
            n[m] += copy_files(src, names, f"{out}/{m}/{sim}") > 0
    with open(f"{out}/calicost/best_init.tsv", "w") as fh:
        fh.write("\n".join(rows) + "\n")
    print(f"{len(sims)} simulations; runs copied: {n}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: extract_deposit_results.py <deposit_dir> <out_dir>")
    main(*sys.argv[1:])
