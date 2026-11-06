"""Per-spot UMI: each bin side's total and SNP-covering UMI, and each run's false-positive tumor spots' UMI against
its other called spots'.

- `plot_umi_page`, `plot_umi_histogram`, `plot_spot_total_umi`, `plot_spot_snp_umi` -- a histogram of per-spot UMI,
  its median marked: an AnnData's total UMI, and allele matrices' SNP-covering UMI, a row per bin side;
- `get_tick`, `plot_false_positive_umi` -- a run's name on a page, and the false-positive tumor spots' total and
  SNP-covering UMI against the run's other called spots'.

The UMI histograms' in-tissue spots are the input's spots: UGP's `bb/<assay>.h5ad` and allele matrices hold the
barcodes of Space Ranger's filtered matrix, the bins under tissue, and UGP filters genes only
(Universal-Genotyping-Pipeline `workflow/scripts/process_rna_anndata.py`). The total UMI is summed over the genes the
AnnData keeps. The SNP-covering UMI of a spot is A + B summed over its SNPs or bins: the count
`cnamaste.io.load_input_data` floors at `quality.spot_min_snp_umis`, and CalicoST's `load_data` at
`min_snpumi_perspot` (https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L72). Space
Ranger's filtered matrix:
https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from matplotlib.figure import Figure
from scripts.metrics.qc_metrics import get_umi_comparison, get_umi_measures
from scripts.metrics.read_method_input import (
    Run,
    SpotLabels,
    get_method_name,
    get_run_name,
    get_spot_sums,
)
from scripts.plots.layout import set_text_sizes
from scripts.plots.style import add_panel_letters, get_panel_rows, get_setting_label

from cnamaste.config import start_time
from cnamaste.figure_style import PAPER_WIDTH, page_style
from cnamaste.logger import get_logger
from cnamaste.plot.layout import add_axes

logger = get_logger("umi", start_time=start_time)


BAR_COLOUR = "#2a78d6"
"""The bars' fill: the dataviz palette's first categorical slot."""


INK = "#52514e"
"""The median line and its label: secondary ink."""


GRID = "#e1e0d9"
"""The gridlines: hairline grey."""


PAGE = (8.5, 11.0)
"""Inches: a US letter page, a histogram's canvas where it is drawn alone (`plot_umi_histogram` without `ax`);
`figure_style.PAPER_WIDTH` is its text column."""


PAGE_MARGIN = (PAGE[0] - PAPER_WIDTH) / 2.0
"""Inches from the page's left and top edges to its text column."""


LABEL_ROOM = (0.6, 0.4)
"""Inches inside the text column left of the axes, for the y label and ticks, and above it, for the two-line title."""


HISTOGRAM_HEIGHT = 2.0
"""Inches: a histogram's axes."""


BAR_SLOT = 0.75
"""Inches of axes per bar on a page of its own, room for its two-line label."""


BAR_WIDTH = 0.6
"""A bar's width, in bar slots."""


REFERENCE_ROOM = 0.7
"""Inches right of the axes for the reference line's label."""


TICK_ROOM = 0.6
"""Inches under the axes for a run's name, three lines: method, resolution, label (`get_tick`)."""


COMPARISON_HEADROOM = 1.6
"""A comparison histogram's y axis top, over its tallest bar: room for its key above the bars."""


COMPARISON_BINS = 60
"""Bins of a comparison histogram, from 0 to the pooled 99.5th percentile."""


def get_tick(run: Run, run_side_um: int) -> str:
    """Get a run's name under its bar: its method (`read_method_input.get_method_name`), `(<side> um)`, then its label
    (`get_setting_label`).

    Args:
        run: the run.
        run_side_um: the side in um of the bins it ran on.

    Returns:
        The name, one line each for the method, the side and the label where the run has one.
    """
    side = rf"({run_side_um} $\mu$m)"
    return "\n".join([get_method_name(run.method), side, *([get_setting_label(run.label)] if run.label else [])])


def plot_false_positive_umi(
    groups: list[tuple[SpotLabels, pd.Index]], umis: pd.DataFrame, colour: str, sample: str
) -> tuple[Figure, pd.DataFrame]:
    """Draw a page of a row per run (`get_panel_rows`), of its false-positive tumor spots' UMI against its other
    called spots'.

    Each row, lettered (a), (b), ..., its letter labelled by the run (`add_panel_letters`), its total UMI left and
    SNP-covering UMI right: the two groups' histograms as densities on shared bins (`COMPARISON_BINS`), their medians
    dashed, keyed with each group's spots and median (`metrics.qc_metrics.get_umi_comparison`), above the bars.

    Args:
        groups: each run's `SpotLabels` and its false-positive tumor spots, a row each in order.
        umis: `total_umi` and `snp_umi` per in-tissue spot (`metrics.read_method_input.get_spot_umis`).
        colour: the false-positive spots' colour, as their spatial panels draw them.
        sample: the sample's name, the table's `sample`.

    Returns:
        The page, and its numbers: `sample`, `method`, `resolution`, then `metrics.qc_metrics.get_umi_comparison`'s
        columns, a row per run and measure.
    """
    tables, titles = [], []
    with page_style() as size:
        figure, axes = get_panel_rows(len(groups))
        for row, (spot, false_positive) in zip(axes, groups, strict=True):
            table = get_umi_comparison(umis, false_positive, spot.labels.index)
            table.insert(0, "resolution", f"{spot.run_side_um}um")
            table.insert(0, "method", get_run_name(spot.run))
            tables.append(table)
            titles.append(" ".join(get_tick(spot.run, spot.run_side_um).split("\n")))
            called = spot.labels.index[spot.labels.index.isin(umis.index)]
            other = called[~called.isin(false_positive)]
            for ax, (measure, xlabel) in zip(row, get_umi_measures().items(), strict=True):
                found, rest = umis.loc[false_positive, measure], umis.loc[other, measure]
                edges = np.linspace(0.0, float(np.percentile(np.concatenate([found, rest]), 99.5)), COMPARISON_BINS + 1)
                for values, name, tint in ((rest, "other called", INK), (found, "false positive", colour)):
                    median = float(np.median(values))
                    ax.hist(
                        values,
                        bins=edges,
                        density=True,
                        histtype="step",
                        color=tint,
                        linewidth=1.0,
                        label=f"{name} ({len(values):,}), median {median:,.0f}",
                    )
                    ax.axvline(median, ymax=1.0 / COMPARISON_HEADROOM, color=tint, linestyle="--", linewidth=0.8)
                ax.set_ylim(0.0, COMPARISON_HEADROOM * ax.get_ylim()[1])
                ax.legend(
                    loc="upper right",
                    frameon=False,
                    handlelength=1.2,
                )
                ax.set_xlim(edges[0], edges[-1])
                ax.set_xlabel(xlabel)
                ax.set_ylabel("density")
                ax.grid(axis="y", color=GRID, linewidth=0.5)
                ax.set_axisbelow(True)
                ax.spines[["top", "right"]].set_visible(False)
        set_text_sizes(figure, size)
        figure.canvas.draw()
        add_panel_letters(figure, list(axes[:, 0]), size, titles)
    table = pd.concat(tables, ignore_index=True)
    table.insert(0, "sample", sample)
    return figure, table


def plot_umi_histogram(
    umi: np.ndarray, xlabel: str, sample: str, side_um: int, bins: int = 100, ax: Any = None, title: bool = True
) -> Figure:
    """Draw a histogram of per-spot UMI, its median marked, titled by sample, bin side and spots where `title`.

    Without `ax`, the figure is a whole letter `PAGE`, the axes the text column wide at its top, so a multi-page PDF
    holds one figure per page.

    Args:
        umi: one value per in-tissue spot.
        xlabel: the quantity drawn, the x axis' label.
        sample: the sample's name.
        side_um: the bins' side in um.
        bins: the histogram's bin count.
        ax: the axes to draw on, or None for a page of its own.
        title: whether the axes carry the title, or a caller titles them, e.g. a row of them.

    Returns:
        The figure `ax` is on, in the page's face (`figure_style.page_style`).
    """
    median = float(np.median(umi))
    with page_style() as size:
        if ax is None:
            figure = plt.figure(figsize=PAGE, facecolor="white")
            ax = add_axes(
                figure,
                PAGE_MARGIN + LABEL_ROOM[0],
                PAGE_MARGIN + LABEL_ROOM[1],
                PAPER_WIDTH - LABEL_ROOM[0],
                HISTOGRAM_HEIGHT,
            )
        figure, axis = ax.figure, ax
        axis.hist(umi, bins=bins, color=BAR_COLOUR, edgecolor="white", linewidth=0.3)
        axis.axvline(median, color=INK, linestyle="--", linewidth=1.0)
        axis.annotate(
            f"median {median:,.0f}",
            xy=(median, 1.0),
            xycoords=("data", "axes fraction"),
            xytext=(3, -3),
            textcoords="offset points",
            va="top",
            color=INK,
        )
        axis.set_xlabel(xlabel)
        axis.set_ylabel("spots")
        if title:
            axis.set_title(rf"{sample} ({side_um} $\mu$m)" + f"\n#in-tissue spots={umi.size:,}")
        axis.grid(axis="y", color=GRID, linewidth=0.5)
        axis.set_axisbelow(True)
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_xlim(left=0)
        set_text_sizes(figure, size)
    logger.info(
        f"Drew {xlabel}, {sample} {side_um} um: {umi.size:_} spots, {int(np.sum(umi == 0)):_} at 0, "
        f"median {median:_.1f}."
    )
    return figure


def plot_spot_total_umi(
    h5ad_file: str | Path,
    sample: str,
    side_um: int,
    bins: int = 100,
    ax: Any = None,
    title: bool = True,
) -> Figure:
    """Draw the distribution of per-spot total UMI of an AnnData.

    Args:
        h5ad_file: an AnnData of (spots, genes) UMI counts in `X`.
        sample: the sample's name.
        side_um: the bins' side in um.
        bins: the histogram's bin count.
        ax: the axes to draw on, or None for a page of its own (`plot_umi_histogram`).
        title: whether the axes carry the title (`plot_umi_histogram`).

    Returns:
        The figure.
    """
    logger.info(f"Reading {h5ad_file}.")
    adata = anndata.read_h5ad(h5ad_file)
    umi = get_spot_sums(adata.X, spot_axis=0)
    figure = plot_umi_histogram(umi, "total UMI per spot", sample, side_um, bins, ax=ax, title=title)
    return figure


def plot_spot_snp_umi(
    a_allele: Any,
    b_allele: Any,
    sample: str,
    side_um: int,
    spot_axis: int = 1,
    bins: int = 100,
    ax: Any = None,
    title: bool = True,
) -> Figure:
    """Draw the distribution of per-spot SNP-covering UMI, A + B.

    Args:
        a_allele: A-allele UMI counts, sparse or dense, one axis spots.
        b_allele: B-allele UMI counts, the shape of `a_allele`.
        sample: the sample's name.
        side_um: the bins' side in um.
        spot_axis: the axis indexing spots, 1 for UGP's (bins, spots), 0 for CalicoST's (spots, SNPs).
        bins: the histogram's bin count.
        ax: the axes to draw on, or None for a page of its own (`plot_umi_histogram`).
        title: whether the axes carry the title (`plot_umi_histogram`).

    Returns:
        The figure.

    Raises:
        ValueError: where the two matrices' shapes differ.
    """
    if a_allele.shape != b_allele.shape:
        raise ValueError(f"A and B allele shapes differ: {a_allele.shape} != {b_allele.shape}")
    umi = get_spot_sums(a_allele, spot_axis) + get_spot_sums(b_allele, spot_axis)
    figure = plot_umi_histogram(umi, "SNP-covering UMI per spot", sample, side_um, bins, ax=ax, title=title)
    return figure


def plot_umi_page(sample: str, bb_dirs: list[tuple[int, Path]], assay: str, bins: int = 100) -> Figure:
    """Draw a page of per-spot UMI histograms, a row per bin side (`style.get_panel_rows`), lettered (a),
    (b), ...

    Each row is total UMI left and SNP-covering UMI right, its letter labelled by sample, bin side and in-tissue spots
    (`style.add_panel_letters`).

    Args:
        sample: the sample's name.
        bb_dirs: each bin side in um and its UGP `bb` directory, one row each, in order.
        assay: UGP's assay type.
        bins: each histogram's bin count.

    Returns:
        The page.

    Raises:
        ValueError: where a row's AnnData and allele matrices hold different spot counts.
    """
    titles = []
    with page_style() as size:
        figure, axes = get_panel_rows(len(bb_dirs))
        for (total_ax, snp_ax), (side_um, bb_dir) in zip(axes, bb_dirs, strict=True):
            snp_dir = bb_dir / "multi_snp" / assay
            h5ad_file = bb_dir / f"{assay}.h5ad"
            plot_spot_total_umi(h5ad_file, sample, side_um, bins=bins, ax=total_ax, title=False)
            logger.info(f"Reading {snp_dir}'s A and B allele matrices.")
            a_allele, b_allele = sp.load_npz(snp_dir / "bb.Aallele.npz"), sp.load_npz(snp_dir / "bb.Ballele.npz")
            plot_spot_snp_umi(a_allele, b_allele, sample, side_um, bins=bins, ax=snp_ax, title=False)
            spots = anndata.read_h5ad(h5ad_file, backed="r")
            n_spots = spots.n_obs
            spots.file.close()
            if a_allele.shape[1] != n_spots:
                raise ValueError(f"{h5ad_file} holds {n_spots:_} spots, its allele matrices {a_allele.shape[1]:_}")
            titles.append(rf"{sample} ({side_um} $\mu$m), #in-tissue spots={n_spots:,}")
        set_text_sizes(figure, size)
        figure.canvas.draw()
        add_panel_letters(figure, list(axes[:, 0]), size, titles)
    return figure
