"""Ported from cna_maste `python/cnamaste/calicost/read.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.

Readers of a CalicoST run directory tree, `<root>/nomixing_calicost_related/<sample_id>/clone*_rectangle*_w*/`.

- `results_dir` -- a sample's directory;
- `best_run` -- its rectangle initialization with the highest `total_llf`;
- `first_run` -- its first rectangle initialization, the lowest index, irrespective of `total_llf`;
- `read_run` -- one run's arrays and tables, as given to `layout.plot_combined_page`(`CalicostRun`);
- `read_estimate` -- one run's clones and copies, as scored by `metrics.metrics`;
- `seglevel` -- one run's `cnv<_ploidy>_seglevel.tsv`, at one of `CALICOST_PLOIDIES`.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scripts.cnamaste_port.qa.metrics import RunEstimate

RUN = re.compile(r"clone\d+_rectangle(\d+)_w[^/]+")
"""A run directory's name, `clone<n>_rectangle<i>_w<weight>`; group 1 is the rectangle index."""


CALICOST_PLOIDIES = ("unconstrained", "diploid", "triploid", "tetraploid")
"""CalicoST's ploidy constraints, one integer copy-number table each: `--calicost_ploidy`'s choices."""


# NB the file each of `CALICOST_PLOIDIES` names, `cnv<suffix>_seglevel.tsv`; read by `seglevel` alone. CalicoST's
#    tables only: cnamaste's own ploidy is `int_copy_num.assumed_ploidy`.
CALICOST_PLOIDY_SUFFIXES = {
    "unconstrained": "",
    "diploid": "_diploid",
    "triploid": "_triploid",
    "tetraploid": "_tetraploid",
}
"""CalicoST's integer copies, one table per ploidy constraint, as `calicost_main.py:310-370` (c1abcae) names them:
`cnv<suffix>_seglevel.tsv`, `unconstrained` its `max_medploidy` None and the rest a median ploidy at most 2, 3, 4."""


DEFAULT_CALICOST_PLOIDY = "diploid"
"""The ploidy a CalicoST run is read at by default: `run_sim_qa --calicost_ploidy`'s default."""


def seglevel(run_dir: Path, ploidy: str = DEFAULT_CALICOST_PLOIDY) -> pd.DataFrame:
    """One CalicoST run's `cnv<suffix>_seglevel.tsv` at `ploidy`.

    Args: `run_dir` a `clone*_rectangle*_w*` directory; `ploidy` one of `CALICOST_PLOIDIES`.
    Returns: the table.
    Raises: `ValueError` on another ploidy; `FileNotFoundError` naming the table where the run has none.
    """
    if ploidy not in CALICOST_PLOIDIES:
        raise ValueError(f"ploidy {ploidy!r}: one of {', '.join(CALICOST_PLOIDIES)}")
    path = run_dir / f"cnv{CALICOST_PLOIDY_SUFFIXES[ploidy]}_seglevel.tsv"
    if not path.is_file():
        raise FileNotFoundError(f"no {path.name} in {run_dir}, CalicoST's {ploidy} copies")
    return pd.read_csv(path, sep="\t")


def final_fit_path(run_dir: Path) -> Path | None:
    """The run's `rdrbaf_final_nstates<n>_smp.npz`, the first by name, or None where it has none."""
    return next(iter(sorted(run_dir.glob("rdrbaf_final_nstates*_smp.npz"))), None)


def fitted_runs(sample_dir: Path) -> list[tuple[int, Path, float]]:
    """Every `clone*_rectangle*_w*` directory holding a final fit, as `(rectangle, run directory, total_llf)`.

    Raises: `FileNotFoundError` where no run directory holds a final fit.
    """
    runs = []
    for run_dir in sorted(sample_dir.glob("clone*_rectangle*_w*")):
        found = RUN.fullmatch(run_dir.name)
        final = final_fit_path(run_dir) if found and run_dir.is_dir() else None
        if found is None or final is None:
            continue
        llf = float(np.load(final, allow_pickle=True)["total_llf"])
        runs.append((int(found.group(1)), run_dir, llf))
    if not runs:
        raise FileNotFoundError(f"no clone*_rectangle*_w*/rdrbaf_final_nstates*_smp.npz under {sample_dir}")
    return runs


def best_run(sample_dir: Path) -> tuple[int, Path, float]:
    """Pick the rectangle initialization with the highest `total_llf`.

    Reads every `clone*_rectangle*_w*` directory holding an `rdrbaf_final_nstates*_smp.npz`, of any clone count,
    rectangle index and weight. A NaN `total_llf` ranks below every number; equals go to the lowest rectangle index,
    then to the first directory by name.

    Args: `sample_dir` the sample's directory, e.g. `results_dir(root, sample_id)`.
    Returns: `(rectangle, run directory, total_llf)`.
    Raises: `FileNotFoundError` where no run directory holds a final fit.
    """
    return min(fitted_runs(sample_dir), key=lambda r: (math.inf if math.isnan(r[2]) else -r[2], r[0], r[1].name))


def read_estimate(run_dir: Path, ploidy: str = DEFAULT_CALICOST_PLOIDY) -> RunEstimate:
    """One CalicoST run's clones and copies: `clone_labels.tsv` by `BARCODES`, and its `seglevel` at `ploidy`.

    CalicoST writes no integer clones, so its clones, `clone_labels.tsv`, stand in for them: its `int_clone_ari` is
    its `clone_ari`. Needs no `binned_data.npz`.

    Args: `run_dir` a `clone*_rectangle*_w*` directory; `ploidy` one of `CALICOST_PLOIDIES`.
    Returns: the `RunEstimate`.
    """
    labels = pd.read_csv(run_dir / "clone_labels.tsv", sep="\t")
    clones = pd.Series(labels["clone_label"].astype(int).to_numpy(), index=labels["BARCODES"].astype(str).to_numpy())
    copies = seglevel(run_dir, ploidy)
    return RunEstimate(clones=clones, int_clones=clones, seglevel=copies)
