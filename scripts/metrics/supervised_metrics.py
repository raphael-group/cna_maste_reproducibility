"""A run scored against a simulation's truth, its cells counted as `cnamaste.qa.metrics` counts them.

- `overlap`, `phase_free`, `decode_ab_pair`, `segment_true_ab_pairs`, `segment_called_ab_pairs`, `count_cells`,
  `get_cell_confusion`, `exact_altered`, `clone_ari` -- copied from `cnamaste.qa.metrics`, on a run's labels and
  profile in place of its `RunEstimate`; a row's true segment found by `gene_level.get_max_overlap`;
- `get_state_confusion`, `get_state_metrics` -- cells by CNA state (`cna_states`): CNA sensitivity, normal
  specificity and CNA recall;
- `compute_supervised_metrics` -- a run's metrics.

A cell is one spot the run labels on one row of the run's profile that a true segment overlaps: a segment, or a gene
where the profile is `gene_level.get_gene_profile`'s. Each row counts once, whatever its length. Its true AB pair is
its true clone's, that of the true segment overlapping the row most, its called AB pair or state its fitted clone's,
so every labelled spot is scored and no clone is matched to another. AB pairs are compared phase-free. A state is
`neu` exactly where its pair is (1, 1), so the state sensitivity and specificity count the cells
`qa.metrics.sensitivity` and `qa.metrics.specificity` count.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scripts.metrics.cna_states import get_ab_states, get_called_states
from scripts.metrics.gene_level import get_max_overlap
from scripts.metrics.read_method_input import MethodInput
from sklearn.metrics import adjusted_rand_score

from cnamaste.truth import Truth

ABPair = tuple[int, int]

SUPERVISED_METRICS = ("sensitivity", "specificity", "exact_altered", "state_altered", "clone_ari")
"""A run's metrics against the truth, in table order: `qa.metrics`' panels' names, and `state_altered`."""


def get_supervised_metrics() -> tuple[str, ...]:
    """Get the metrics against the truth: `sensitivity`, `specificity`, `exact_altered`, `state_altered`, `clone_ari`."""
    return SUPERVISED_METRICS


def overlap(true_labels: np.ndarray, fitted: np.ndarray, n_true: int, n_fitted: int) -> np.ndarray:
    """Get the contingency table of two labellings (`qa.metrics.overlap`).

    Args:
        true_labels: the true clone per kept spot.
        fitted: the fitted clone per kept spot.
        n_true: the true clones.
        n_fitted: the fitted clones.

    Returns:
        `(n_true, n_fitted)` counts, `counts[i, j]` the spots in true clone `i` and fitted clone `j`.
    """
    counts = np.zeros((n_true, n_fitted), dtype=np.int64)
    # NB `np.add.at` is unbuffered: a repeated index pair adds once per spot.
    np.add.at(counts, (true_labels, fitted), 1)
    return counts


def phase_free(codes: np.ndarray) -> np.ndarray:
    """Get `A * 1_000 + B` codes as `(minor, major)`, `(A, B)` and `(B, A)` coded alike (`qa.metrics.phase_free`).

    Args:
        codes: AB pair codes.

    Returns:
        Phase-free codes.
    """
    major, minor = codes // 1_000, codes % 1_000
    return np.asarray(np.minimum(major, minor) * 1_000 + np.maximum(major, minor))


def decode_ab_pair(code: int) -> ABPair:
    """Get `A * 1_000 + B` as `(A, B)` (`qa.metrics.decode_ab_pair`)."""
    return int(code // 1_000), int(code % 1_000)


def fitted_at_barcode(labels: pd.Series, truth: Truth) -> np.ndarray:
    """Get a run's clones at the truth's barcodes, -1 where the run labels none (`qa.metrics.fitted_at_barcode`).

    Args:
        labels: the run's clone per barcode.
        truth: the truth.

    Returns:
        One clone per truth barcode.
    """
    return labels.reindex(truth.barcodes).fillna(-1).to_numpy(dtype=np.int64)


def segment_true_ab_pairs(profile: pd.DataFrame, truth: Truth) -> np.ndarray:
    """Get each row's true AB pair per true clone: that of the true segment overlapping it most
    (`qa.metrics.segment_true_ab_pairs`).

    Args:
        profile: the run's profile, `CHR`, `START`, `END` per row.
        truth: the truth.

    Returns:
        `(rows, true clones)` phase-free codes, -1 on a row no true segment overlaps.
    """
    true = truth.profile.rename(columns={"chr": "CHR", "start": "START", "end": "END"})
    best = get_max_overlap(profile, true)
    codes = np.full((len(profile), truth.n_clones), -1, dtype=np.int64)
    found = best >= 0
    for c, clone in enumerate(truth.clones):
        a, b = truth.profile[f"{clone}_A_copy"].to_numpy(), truth.profile[f"{clone}_B_copy"].to_numpy()
        codes[found, c] = phase_free(a[best[found]] * 1_000 + b[best[found]])
    return codes


def segment_called_ab_pairs(profile: pd.DataFrame) -> np.ndarray:
    """Get each row's called AB pair per fitted clone (`qa.metrics.segment_called_ab_pairs`).

    Args:
        profile: the run's profile, `clone<k> A`, `clone<k> B` per fitted clone `k`.

    Returns:
        `(rows, fitted clones)` phase-free codes.
    """
    n = sum(1 for c in profile.columns if str(c).startswith("clone") and str(c).endswith(" A"))
    return np.stack(
        [
            phase_free(
                profile[f"clone{k} A"].to_numpy(dtype=np.int64) * 1_000
                + profile[f"clone{k} B"].to_numpy(dtype=np.int64)
            )
            for k in range(n)
        ],
        axis=1,
    )


def count_cells(true_values: np.ndarray, called_values: np.ndarray, truth: Truth, fitted: np.ndarray) -> dict:
    """Get cells by `(true, called)` value over the scored rows, each spot on its clones' rows
    (`qa.metrics.cells_by_ab_pair`, for any values).

    Args:
        true_values: `(rows, true clones)` values.
        called_values: `(rows, fitted clones)` values.
        truth: the truth.
        fitted: each truth barcode's fitted clone, -1 where the run labels none.

    Returns:
        `(true value, called value)` -> cells, sorted.
    """
    kept = fitted >= 0
    spots = overlap(truth.labels[kept], fitted[kept], truth.n_clones, called_values.shape[1])
    counts: dict[tuple[Any, Any], int] = {}
    for c, k in zip(*np.nonzero(spots), strict=True):
        pairs, n = np.unique(np.stack([true_values[:, c], called_values[:, k]], axis=1), axis=0, return_counts=True)
        for (t, f), rows in zip(pairs.tolist(), n.tolist(), strict=True):
            counts[(t, f)] = counts.get((t, f), 0) + int(spots[c, k]) * rows
    return dict(sorted(counts.items()))


def get_cell_confusion(run: MethodInput, truth: Truth) -> dict[tuple[ABPair, ABPair], int]:
    """Get cells by `(true, called)` phase-free AB pair over the rows a true segment overlaps
    (`qa.metrics.get_cell_confusion`).

    Args:
        run: the run, its profile AB pairs.
        truth: the truth.

    Returns:
        `(true AB pair, called AB pair)` -> cells, sorted.
    """
    true_ab, called_ab = segment_true_ab_pairs(run.profile, truth), segment_called_ab_pairs(run.profile)
    scored = (true_ab[:, 0] >= 0) & (called_ab[:, 0] >= 0)
    counts = count_cells(true_ab[scored], called_ab[scored], truth, fitted_at_barcode(run.labels, truth))
    return {(decode_ab_pair(t), decode_ab_pair(f)): n for (t, f), n in counts.items()}


def exact_altered(run: MethodInput, truth: Truth) -> float:
    """Get the truly altered cells called their true AB pair, phase-free: N-CNA recall (`qa.metrics.exact_altered`).

    Args:
        run: the run, its profile AB pairs.
        truth: the truth.

    Returns:
        The rate, NaN where the truth has no altered cell.
    """
    counts = get_cell_confusion(run, truth)
    altered = sum(n for (t, _), n in counts.items() if t != (1, 1))
    return sum(n for (t, f), n in counts.items() if t != (1, 1) and t == f) / altered if altered else float("nan")


def clone_ari(labels: pd.Series, truth: Truth) -> float:
    """Get the adjusted Rand index of a run's clones against the true ones, over the spots it labels
    (`qa.metrics.clone_ari`).

    Args:
        labels: the run's clone per barcode.
        truth: the truth.

    Returns:
        The ARI.
    """
    fitted = fitted_at_barcode(labels, truth)
    kept = fitted >= 0
    return float(adjusted_rand_score(truth.labels[kept], fitted[kept]))


def get_state_confusion(profile: pd.DataFrame, truth: Truth, fitted: np.ndarray) -> dict[tuple[str, str], int]:
    """Get cells by `(true state, called state)`, the diagonal included, as `get_cell_confusion` counts AB pairs.

    The scored rows are those a true segment overlaps and the run calls a state on for every clone.

    Args:
        profile: the run's profile, `CHR`, `START` and `END` per row.
        truth: the truth.
        fitted: each truth barcode's fitted clone, -1 where the run labels none.

    Returns:
        `(true state, called state)` -> cells, sorted.
    """
    called = get_called_states(profile)
    true_ab = segment_true_ab_pairs(profile, truth)
    scored = (true_ab[:, 0] >= 0) & (called != "").all(axis=1)
    return count_cells(get_ab_states(true_ab[scored]), called[scored], truth, fitted)


def get_state_metrics(confusion: dict[tuple[str, str], int]) -> dict[str, float]:
    """Get CNA sensitivity, normal specificity and CNA recall from cells by state.

    `sensitivity` is 1 - (truly altered cells called `neu`) / (truly altered cells), `specificity` 1 - (truly `neu`
    cells called altered) / (truly `neu` cells), `state_altered` (truly altered cells called their true state) / (truly
    altered cells); each NaN where its denominator is 0.

    Args:
        confusion: `get_state_confusion`'s cells.

    Returns:
        `sensitivity`, `specificity` and `state_altered`.
    """
    altered = sum(n for (t, _), n in confusion.items() if t != "neu")
    missed = sum(n for (t, f), n in confusion.items() if t != "neu" and f == "neu")
    neutral = sum(n for (t, _), n in confusion.items() if t == "neu")
    false = sum(n for (t, f), n in confusion.items() if t == "neu" and f != "neu")
    exact = sum(n for (t, f), n in confusion.items() if t != "neu" and t == f)
    return {
        "sensitivity": 1.0 - missed / altered if altered else float("nan"),
        "specificity": 1.0 - false / neutral if neutral else float("nan"),
        "state_altered": exact / altered if altered else float("nan"),
    }


def compute_supervised_metrics(run: MethodInput, truth: Truth) -> dict[str, float]:
    """Get a run's metrics against the truth (`SUPERVISED_METRICS`), on its profile's rows.

    `exact_altered` is NaN where the run calls states, not AB pairs.

    Args:
        run: the run, its profile per segment, or per gene (`gene_level.get_gene_profile`).
        truth: the truth.

    Returns:
        `SUPERVISED_METRICS` by name.

    Raises:
        ValueError: where the run labels a clone its profile has no column of.
    """
    fitted = fitted_at_barcode(run.labels, truth)
    clones = get_called_states(run.profile).shape[1]
    if fitted.max() >= clones:
        raise ValueError(f"clone {fitted.max()} labelled, but the profile has {clones} clones")
    found = get_state_metrics(get_state_confusion(run.profile, truth, fitted))
    return {
        "sensitivity": found["sensitivity"],
        "specificity": found["specificity"],
        "exact_altered": exact_altered(run, truth) if "clone0 A" in run.profile else float("nan"),
        "state_altered": found["state_altered"],
        "clone_ari": clone_ari(run.labels, truth),
    }
