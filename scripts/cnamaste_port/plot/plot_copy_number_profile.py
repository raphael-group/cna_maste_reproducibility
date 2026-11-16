"""Ported from cna_maste `python/cnamaste/plot/plot_copy_number_profile.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.
"""

import numpy as np
from scripts.cnamaste_port.palette import COPY_STATES, state_colour

BAR_SHARE = 0.8
"""The share of a clone's row its bars take, centred in it; the rest is the gap to the next clone. One value for
both profiles, `plot_copy_number_profile`'s pair of allele bars and `plot_hatched_profile`'s single hatched bar, so
the two draw a clone's row alike."""


def get_intervals(pred_cnv):
    """Find contiguous intervals of equal state in `pred_cnv`, vectorized.

    Args: `pred_cnv` the state labels.
    Returns: the (start, end) index intervals and each interval's state label.
    """
    pred_cnv = np.asarray(pred_cnv)
    if len(pred_cnv) == 0:
        return [], []

    changes = np.where(pred_cnv[:-1] != pred_cnv[1:])[0] + 1
    splits = np.concatenate(([0], changes, [len(pred_cnv)]))

    intervals = [(splits[i], splits[i + 1]) for i in range(len(splits) - 1)]
    labs = pred_cnv[splits[:-1]].tolist()

    return intervals, labs


OUTLINE = 0.5
"""Points, for row outlines and an aberration's ends."""


def used_states(df_cnv, order):
    """The phased `(A, B)` the clones of `order` carry, in `COPY_STATES`' order, a pair before its flip.

    Args: `df_cnv` copy numbers; `order` clone ids.
    Returns: the `(A, B)` tuples.
    """
    seen = set()
    for cid in order:
        a = df_cnv[f"clone{cid} A"].fillna(1).to_numpy().astype(int)
        b = df_cnv[f"clone{cid} B"].fillna(1).to_numpy().astype(int)
        seen.update(zip(a.tolist(), b.tolist()))
    rank = {state: i for i, state in enumerate(COPY_STATES)}
    return sorted(seen, key=lambda s: (rank.get((max(s), min(s)), len(rank)), s[1] > s[0]))


def state_fill(a, b):
    """`(a, b)`'s profile fill, opaque over white: `palette.state_colour` composited, so a flip is a lighter fill."""
    r, g, bl, alpha = state_colour(a, b)
    return tuple(float(alpha * v + 1.0 - alpha) for v in (r, g, bl))
