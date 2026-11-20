"""A run's quality without a truth: its genome by CNA state, and its calls on the spots.

- `compute_genome_composition` -- per clone, the share of the genome, or of the genes, in each CNA state and
  `uncalled`;
- `get_coverage` -- each run's called spots, normal and tumor, and its spot coverage rate;
- `get_false_positive_tumor`, `get_false_positive_counts` -- a run's false-positive tumor spots, and each run's normal
  and tumor calls per broad cell type, its false-positive tumor spots, its called spots of a normal type, and its
  false-positive rate;
- `get_umi_comparison` -- the false-positive tumor spots' total and SNP-covering UMI against the run's other called
  spots'.

The genome is the regions' bases (`gene_level.get_autosomal_regions`), or the genes inside them
(`gene_level.get_region_genes`). A row of a run's profile weighs its bases inside the regions, or 1 per gene; a share
is the rows' weight over the regions' bases, or over the genes. Per clone, the six states and `uncalled` sum to 1.

A spot is called where the run labels it with any clone, normal (clone 0) or tumor (any other); the spot coverage rate
is #called / #in-tissue, the in-tissue spots UGP's barcodes at the spots' side (`read_method_input.read_spots`). A
labelled spot outside them counts toward none, and is logged. A run on coarser bins gives each bin's clone to every spot
inside it; its resolution is the side of those bins.

The cell types give normal truth only: a spot of a normal type, neither `EPITHELIAL` nor `LOW_CONFIDENCE`
(`read_method_input.get_normal_type`), is normal, so a run's tumor call on it, a clone other than the normal clone 0,
is a false positive, `#false-positive-tumor-spots`. Over the run's called spots of a normal type, `#normal-type-spots`,
the negatives N, that is its false-positive rate, `FALSE_POSITIVE`, FP / N (James et al., An Introduction to
Statistical Learning, 1st ed., sec. 4.4.3, Table 4.7, https://doi.org/10.1007/978-1-4614-7138-7). An epithelial,
low-confidence or untyped spot is in neither count.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scripts.cnamaste_port.config import start_time
from scripts.cnamaste_port.logger import get_logger
from scripts.metrics.cna_states import get_cna_states, get_uncalled
from scripts.metrics.read_method_input import (
    EPITHELIAL,
    LOW_CONFIDENCE,
    SpotLabels,
    get_normal_type,
    get_run_name,
)

logger = get_logger("qc_metrics", start_time=start_time)

COMPOSITION_COLUMNS = ("clone", "#spots", *get_cna_states(), get_uncalled())
"""A composition row's columns: the clone, its spots, a share per CNA state and `uncalled`."""


def compute_genome_composition(
    states: np.ndarray, weights: np.ndarray, total: float, spots: pd.Series | None = None
) -> pd.DataFrame:
    """Get, per clone, the share of the genome in each CNA state and `uncalled`.

    A row of `all` weighs each clone by its spots: the share of the cells, a labelled spot at a gene or base, in each
    state.

    Args:
        states: `(rows, clones)` states (`cna_states.get_called_states`), "" or `uncalled` where none is called.
        weights: per row, its bases inside the regions, or 1 per gene.
        total: the regions' bases, or the genes.
        spots: the run's clone per spot, or None for no `all` row.

    Returns:
        `COMPOSITION_COLUMNS`, a row per clone, `clone<k>`, then `all` where `spots` is given.

    Raises:
        ValueError: where the rows' weight exceeds `total`, segments overlapping one another.
    """
    counts = spots.value_counts() if spots is not None else pd.Series(dtype=int)
    rows = []
    for k in range(states.shape[1]):
        row: dict[str, object] = {"clone": f"clone{k}", "#spots": int(counts.get(k, 0))}
        for state in get_cna_states():
            row[state] = float(weights[states[:, k] == state].sum() / total)
        called = sum(row[state] for state in get_cna_states())
        if called > 1 + 1e-9:
            raise ValueError(f"clone{k}: {called:.6f} of the genome called, segments overlap")
        row[get_uncalled()] = 1.0 - called
        rows.append(row)
    table = pd.DataFrame(rows, columns=list(COMPOSITION_COLUMNS))
    if spots is not None and table["#spots"].sum():
        shares = [*get_cna_states(), get_uncalled()]
        weighted = table[shares].mul(table["#spots"], axis=0).sum() / table["#spots"].sum()
        table.loc[len(table)] = {"clone": "all", "#spots": int(table["#spots"].sum()), **weighted.to_dict()}
    return table


FALSE_POSITIVE = "false-positive-rate"
"""The table's column of a run's false-positive rate: its false-positive tumor spots over its called spots of a normal
type."""


UMI_MEASURES = {"total_umi": "total UMI per spot", "snp_umi": "SNP-covering UMI per spot"}
"""Each per-spot UMI measure compared (`read_method_input.get_spot_umis`), and its axis label."""


def get_false_positive_column() -> str:
    """Get the table's column of a run's false-positive rate, `FALSE_POSITIVE`."""
    return FALSE_POSITIVE


def get_umi_measures() -> dict[str, str]:
    """Get each per-spot UMI measure compared, `total_umi` and `snp_umi`, and its axis label, `UMI_MEASURES`."""
    return UMI_MEASURES


def get_coverage(labels: list[SpotLabels], spots: pd.Index) -> pd.DataFrame:
    """Get each run's coverage of the in-tissue `spots` and its normal and tumor calls on them.

    `#called` is the spots a run labels with any clone, `spot-coverage-rate` #called / #in-tissue, and of them
    `#normal-spots` those of the normal clone 0, `#tumor-spots` the rest. A labelled spot off `spots` counts toward
    none, and is logged.

    Args:
        labels: each run's `SpotLabels` (`read_method_input.get_spot_labels`).
        spots: every in-tissue spot.

    Returns:
        One row per run, in `labels`' order: `method` the run's name, `resolution` the side of the bins it ran on,
        `<side>um`, `#in-tissue`, `#called`, `spot-coverage-rate`, `#normal-spots`, `#tumor-spots`.

    Raises:
        ValueError: where `spots` is empty.
    """
    if spots.empty:
        raise ValueError("no in-tissue spots")
    rows = []
    for run, found, run_side_um in labels:
        name = get_run_name(run)
        off = int((~found.index.isin(spots)).sum())
        if off:
            logger.info(f"{name}: {off:_} labelled spots off the in-tissue spots, left out.")
        clones = found.reindex(spots).dropna()
        rows.append(
            {
                "method": name,
                "resolution": f"{run_side_um}um",
                "#in-tissue": len(spots),
                "#called": len(clones),
                "#normal-spots": int((clones == 0).sum()),
                "#tumor-spots": int((clones != 0).sum()),
            }
        )
    table = pd.DataFrame(rows)
    table.insert(4, "spot-coverage-rate", table["#called"] / table["#in-tissue"])
    return table


def get_false_positive_tumor(labels: pd.Series, normal_type: pd.Index) -> pd.Index:
    """Get a run's false-positive tumor spots: the normal-type spots it labels with a clone other than the normal 0.

    Args:
        labels: the run's clone per spot barcode.
        normal_type: the in-tissue spots of a normal type (`read_method_input.get_normal_type`).

    Returns:
        The false-positive tumor spots' barcodes, in `normal_type`'s order.
    """
    clones = labels.reindex(normal_type)
    return normal_type[(clones.notna() & (clones != 0)).to_numpy()]


def get_false_positive_counts(labels: list[SpotLabels], spots: pd.Index, broad: pd.Series) -> pd.DataFrame:
    """Get each run's calls per broad cell type and its false-positive rate.

    Of the in-tissue `spots` with a broad cell type, per type, `#normal-spots-<type>` is those the run labels with the
    normal clone 0 and `#tumor-spots-<type>` those it labels with another. Its tumor calls on a normal type are
    `#false-positive-tumor-spots` (`get_false_positive_tumor`), and its calls on a normal type, normal or tumor,
    `#normal-type-spots`. `FALSE_POSITIVE` is `#false-positive-tumor-spots` over `#normal-type-spots`, NaN where that
    is 0.

    Args:
        labels: each run's `SpotLabels` (`read_method_input.get_spot_labels`).
        spots: every in-tissue spot.
        broad: each spot's broad cell type, per barcode (`read_method_input.read_broad_cell_types`).

    Returns:
        One row per run, in `labels`' order: `method` the run's name, `resolution` the side of the bins it ran on,
        `<side>um`, the two counts per type, `EPITHELIAL` first, `LOW_CONFIDENCE` last, the rest alphabetical, a space
        in a type's name written `_`, then `#false-positive-tumor-spots`, `#normal-type-spots` and `FALSE_POSITIVE`.

    Raises:
        ValueError: where no in-tissue spot is of a normal type.
    """
    typed = broad.reindex(spots).dropna()
    normal_type = typed.index[get_normal_type(typed).to_numpy()]
    if normal_type.empty:
        raise ValueError(f"no in-tissue spot has a broad cell type other than {EPITHELIAL} or {LOW_CONFIDENCE}")
    types = sorted(set(typed), key=lambda t: (t != EPITHELIAL, t == LOW_CONFIDENCE, t))
    rows = []
    for run, found, run_side_um in labels:
        clones = found.reindex(spots)
        row: dict[str, Any] = {"method": get_run_name(run), "resolution": f"{run_side_um}um"}
        for cell_type in types:
            of = clones.reindex(typed.index[(typed == cell_type).to_numpy()])
            row[f"#normal-spots-{cell_type.replace(' ', '_')}"] = int((of == 0).sum())
            row[f"#tumor-spots-{cell_type.replace(' ', '_')}"] = int((of.notna() & (of != 0)).sum())
        false_positive = len(get_false_positive_tumor(found, normal_type))
        negatives = int(clones.reindex(normal_type).notna().sum())
        row["#false-positive-tumor-spots"] = false_positive
        row["#normal-type-spots"] = negatives
        row[FALSE_POSITIVE] = false_positive / negatives if negatives else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def get_umi_comparison(umis: pd.DataFrame, false_positive: pd.Index, called: pd.Index) -> pd.DataFrame:
    """Get each UMI measure's spots and median, of a run's false-positive tumor spots and of its other called spots.

    Args:
        umis: `total_umi` and `snp_umi` per spot (`read_method_input.get_spot_umis`).
        false_positive: the run's false-positive tumor spots (`get_false_positive_tumor`).
        called: the run's called spots; those in `umis` and not in `false_positive` are the other spots.

    Returns:
        One row per measure of `UMI_MEASURES`: `measure`, `#false-positive`, `#other`, `median-false-positive`,
        `median-other`.
    """
    in_tissue = called[called.isin(umis.index)]
    other = in_tissue[~in_tissue.isin(false_positive)]
    false_positive = false_positive[false_positive.isin(umis.index)]
    rows = []
    for measure in UMI_MEASURES:
        found, rest = umis.loc[false_positive, measure], umis.loc[other, measure]
        rows.append(
            {
                "measure": measure,
                "#false-positive": len(found),
                "#other": len(rest),
                "median-false-positive": float(np.median(found)),
                "median-other": float(np.median(rest)),
            }
        )
    return pd.DataFrame(rows)
