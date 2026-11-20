"""The runs' joint segments, the union of their breakpoints, and each run's clone profile on them.

- `get_joint_segments` -- per chromosome, the intervals between consecutive breakpoints of any run, where a run calls;
- `get_segment_level_table` -- each run's fields per joint segment, a `<sample>.segment_level_cnp.tsv`.

A joint segment lies inside at most one segment of each run, so every run's copies are constant on it; a gene or SNP
block takes the joint segment overlapping it most (`gene_level.get_max_overlap`), the same one for every run, as
copytyping-benchmark's `joint_segments` (`analysis/scripts/cn_error.py`) takes "the union of every competing run's
breakpoints".
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scripts.cnamaste_port.config import start_time
from scripts.cnamaste_port.logger import get_logger
from scripts.metrics.gene_level import get_chromosomes, get_interval_fields
from scripts.metrics.read_method_input import MethodInput

logger = get_logger("segment_level", start_time=start_time)


def get_joint_segments(profiles: list[pd.DataFrame]) -> pd.DataFrame:
    """Get the runs' joint segments: per chromosome, each interval between consecutive breakpoints, a breakpoint any
    run's segment's `START` or `END`, kept where some run's segment holds it.

    Args:
        profiles: each run's `MethodInput.profile`, `CHR`, `START`, `END` per segment.

    Returns:
        `CHR` (no `chr`), `START`, `END` per joint segment, in genome order.
    """
    intervals = pd.concat(
        [
            pd.DataFrame(
                {
                    "CHR": get_chromosomes(profile["CHR"]),
                    "START": profile["START"].to_numpy(dtype=np.int64),
                    "END": profile["END"].to_numpy(dtype=np.int64),
                }
            )
            for profile in profiles
        ],
        ignore_index=True,
    )
    parts = []
    for chrom, group in intervals.groupby("CHR", sort=False):
        points = np.unique(group[["START", "END"]].to_numpy().ravel())
        lo, hi = points[:-1], points[1:]
        starts, ends = group["START"].to_numpy(), group["END"].to_numpy()
        held = ((starts[None, :] <= lo[:, None]) & (ends[None, :] >= hi[:, None])).any(axis=1)
        parts.append(pd.DataFrame({"CHR": chrom, "START": lo[held], "END": hi[held]}))
    segments = pd.concat(parts, ignore_index=True)
    order = pd.to_numeric(segments["CHR"], errors="coerce")
    segments = segments.assign(order=order).sort_values(["order", "CHR", "START"]).drop(columns="order")
    logger.info(f"Joint segments of {len(profiles)} runs: {len(segments):_}.")
    return segments.reset_index(drop=True)


def get_segment_level_table(segments: pd.DataFrame, runs: dict[str, MethodInput]) -> pd.DataFrame:
    """Get each run's fields per joint segment: a `<sample>.segment_level_cnp.tsv`.

    Args:
        segments: the joint segments (`get_joint_segments`).
        runs: each run by its name, e.g. `calicost`.

    Returns:
        `CHR`, `START`, `END`, then `<name> clone<k> <field>` per run, clone and field
        (`gene_level.get_interval_fields`), `uncalled` where the run has no segment there.
    """
    parts = [segments[["CHR", "START", "END"]].reset_index(drop=True)]
    for name, run in runs.items():
        parts.append(get_interval_fields(run.profile, segments).add_prefix(f"{name} ").reset_index(drop=True))
    return pd.concat(parts, axis=1)
