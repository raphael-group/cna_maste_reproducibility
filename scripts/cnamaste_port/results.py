"""Ported from cna_maste `python/cnamaste/results.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.

Where a run's results are written: its name (`run_id`), its files' names, and the integer-copy tables'.
"""

from __future__ import annotations

from pathlib import Path


def get_output_path(run_dir, level):
    """The run's one integer-copy table at `level` in `run_dir`, whichever ploidy names it (`get_output_path_prefix`).

    Args: `run_dir` a `run_cnamaste` output directory; `level` `seglevel`, `genelevel` or `perstate`.
    Returns: the path.
    Raises: FileNotFoundError where `run_dir` holds none, or more than one.
    """
    found = sorted(Path(run_dir).glob(f"cnv_*_{level}.tsv"))
    if len(found) != 1:
        raise FileNotFoundError(f"one cnv_<ploidy>_{level}.tsv in {run_dir}, found {[p.name for p in found]}")
    return found[0]
