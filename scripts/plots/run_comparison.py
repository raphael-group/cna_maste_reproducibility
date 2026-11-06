"""A sample's runs compared: each run's proportion of cells per altered state, its state Jaccard index against the
reference run's, and the clone ARI and CNP similarity of every two runs.

- `plot_genome_composition`, `draw_genome_composition` -- each run's proportion of cells per altered state;
- `plot_cna_jaccard`, `get_confusion_rows`, `get_label_rooms` -- each run's state Jaccard index against the
  reference run's;
- `get_run_comparison_panels`, `draw_run_pair_panel`, `draw_run_pairs` -- the panels comparing the runs beside the
  methods' profiles (`layout.plot_profile_page`), the clone ARI or the CNP similarity of every two runs.
"""

from __future__ import annotations

from collections.abc import Callable
from functools import partial
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from scripts.metrics.cna_states import get_cna_states
from scripts.metrics.read_method_input import get_method_name
from scripts.metrics.unsupervised_metrics import compute_state_jaccard
from scripts.plots.cn_error import draw_skill_score_panel
from scripts.plots.layout import set_text_sizes
from scripts.plots.style import (
    add_panel_letters,
    get_group_width,
    get_legend_room,
    get_method_fill,
    get_ordered_runs,
    get_patch,
    get_run_labels,
    get_run_order,
)

from cnamaste.config import start_time
from cnamaste.figure_style import page_size, page_style
from cnamaste.logger import get_logger
from cnamaste.plot.layout import GAP, LETTER, add_axes, place

logger = get_logger("run_comparison", start_time=start_time)

PLOTTED_STATES = tuple(s for s in get_cna_states() if s not in ("neu", "bdel"))
"""The CNA states the composition draws: `del`, `amp`, `cnloh` and `bamp`."""


CONFUSION_MIN_WIDTH = 5
"""A confusion panel's least width, in states, as its row counts it (`get_confusion_rows`): room for its axis label."""


CONFUSION_ROW_STATES = 10
"""The states a row of confusion panels holds, side by side."""


COLOURBAR_WIDTH = 0.12
"""Inches across the state Jaccard page's colour bar, right of its first row (`plot_cna_jaccard`)."""


FIT_ROUNDS = 3
"""Rounds of sizing a later row of the state Jaccard page to the first row's width, its labels measured anew each
round, as a panel's labels reach past it by an amount its side sets (`plot_cna_jaccard`)."""


HEATMAP_COLOURS = ("#FFFFFF", "#2a78d6")
"""A heatmap's colours at 0 and 1, the state Jaccard indices and the metrics of every two runs."""


RUN_PAIR_METRICS = {
    "clone_ari": "Clone ARI of every two runs, on the spots both label",
    "cnp_similarity": (
        "CNP similarity of every two runs: of their most similar tumor clones, the share of the genes either alters "
        "where their CNA states match"
    ),
}
"""Each metric of every two runs `draw_run_pair_panel` draws, by its panel title."""


PAIR_CELL = 0.24
"""Inches a side of a cell of a heatmap of every two runs beside the profiles (`get_run_comparison_panels`): room for
its value, `0.00`, at the page's size."""


LEVEL_NAMES = {"ascn": "(A,B) state", "state": "CNA state", "tcn_state": "total-copy state"}
"""Each field two runs are compared at (`metrics.cna_states.get_comparison_field`), as a state Jaccard panel's title
names it."""


COMPOSITION_HEADROOM = 1.4
"""The composition panel's y axis top, over its tallest bar: room for the bars' values, upright above them."""


def draw_genome_composition(ax: Any, runs: pd.DataFrame) -> None:
    """Draw each run's share of cells in each altered state, a group per state (`PLOTTED_STATES`), a bar per run in its
    method's fill, in `runs`' order, its value to 0.001 upright above it; the axis' top `COMPOSITION_HEADROOM` times the
    tallest bar.

    Args:
        ax: the panel.
        runs: `run`, `method`, then `qc_metrics.compute_genome_composition`'s columns, each run's `all` row, its clones
            weighed by their spots.
    """
    width = get_group_width() / max(len(runs), 1)
    shares = runs[list(PLOTTED_STATES)].to_numpy(dtype=float)
    for i, method in enumerate(runs["method"]):
        offset = (i - (len(runs) - 1) / 2) * width
        values = shares[i]
        x = np.arange(len(PLOTTED_STATES)) + offset
        bars = ax.bar(x, values, width=width, color=get_method_fill(method), edgecolor="white", linewidth=0.5)
        ax.bar_label(bars, fmt="{:.3f}", padding=2, rotation=90)
    tallest = float(np.nanmax(shares, initial=0.0))
    ax.set_ylim(0.0, COMPOSITION_HEADROOM * tallest if tallest > 0 else 1.0)
    ax.set_xticks(range(len(PLOTTED_STATES)), PLOTTED_STATES)
    ax.set_xlabel("CNA state")
    ax.set_ylabel("Proportion")


def plot_genome_composition(composition: pd.DataFrame, sample: str, names: dict[str, str]) -> Figure:
    """Draw each run's proportion of cells per altered state (`draw_genome_composition`), a group per state, each run
    in its method's fill and its value above its bar, in method order (`get_ordered_runs`), keyed under it by run in that order; untitled, its
    letter labelled by the sample (`add_panel_letters`).

    Args:
        composition: `run`, `method` and `qc_metrics.compute_genome_composition`'s columns per run and clone.
        sample: the sample, in its letter's label.
        names: each run's name (`get_run_labels`).

    Returns:
        The figure, `page_size("half")`, one panel.
    """
    runs = get_ordered_runs(composition[composition["clone"] == "all"])
    with page_style() as size:
        figure, ax = plt.subplots(figsize=page_size("half"))
        draw_genome_composition(ax, runs)
        key = [get_patch(get_method_fill(m), names[r]) for r, m in runs[["run", "method"]].itertuples(index=False)]
        figure.legend(handles=key, loc="lower center", ncol=len(key), frameon=False)
        set_text_sizes(figure, size)
        figure.tight_layout(rect=(0, get_legend_room() / 2, 1, 1))
        title = f"{sample}: proportion of each run's cells per altered state, its clones weighed by their spots"
        add_panel_letters(figure, [ax], size, [title])
    return figure


def get_confusion_rows(widths: list[int]) -> list[list[int]]:
    """Get the confusion page's rows: the panels in order, a row holding up to `CONFUSION_ROW_STATES` states, a panel
    counted at least `CONFUSION_MIN_WIDTH`.

    Args:
        widths: each panel's states, in order.

    Returns:
        Each row's panels, by their index, in order.
    """
    rows: list[list[int]] = [[]]
    states = 0
    for k, width in enumerate(widths):
        width = max(width, CONFUSION_MIN_WIDTH)
        if rows[-1] and states + width > CONFUSION_ROW_STATES:
            rows.append([])
            states = 0
        rows[-1].append(k)
        states += width
    return rows


def get_label_rooms(ax: Any, renderer: Any) -> tuple[float, float, float]:
    """Get the inches an axes' labels reach past its frame, left, right and below.

    Args:
        ax: the axes, its texts at their final sizes.
        renderer: its figure's renderer.

    Returns:
        The rooms left, right and below, in inches.
    """
    frame, tight, dpi = ax.get_window_extent(renderer), ax.get_tightbbox(renderer), ax.figure.dpi
    return (frame.x0 - tight.x0) / dpi, (tight.x1 - frame.x1) / dpi, (frame.y0 - tight.y0) / dpi


def get_fill_side(rooms: list[tuple[float, float, float]], width: float) -> float:
    """Get the one side of square panels side by side, `GAP` apart past their labels, that fill `width` inches with
    their labels.

    Args:
        rooms: each panel's label rooms, left, right and below (`get_label_rooms`).
        width: the inches the panels and their labels fill.

    Returns:
        The side, in inches.
    """
    return (width - sum(left + right for left, right, _ in rooms) - GAP * (len(rooms) - 1)) / len(rooms)


def plot_cna_jaccard(
    confusions: dict[str, pd.DataFrame],
    reference: str,
    names: dict[str, str],
    levels: dict[str, str],
) -> Figure:
    """Draw each run's state Jaccard index against the reference run's (`unsupervised_metrics.compute_state_jaccard`),
    one heatmap per run, in `confusions`' order, named by the field the two are compared at (`LEVEL_NAMES`), its
    letter's label (`add_panel_letters`).

    A row is the reference's state, a column the run's, in the confusion's order, the states either run calls on a cell
    both call; a cell's value is the Jaccard index of the cells the reference calls its row's state and the run its
    column's, over the cells both call, so the diagonal is the agreement on each state; written where 0.005 or more,
    blank where neither state holds a cell. The run's states are named under it, turned a quarter. The panels fill
    rows of up to `CONFUSION_ROW_STATES` states (`get_confusion_rows`), a row's panels top-aligned `GAP` apart past
    their labels, `LETTER` over each row for its letters and `GAP` between rows; the rows' first frames share their
    left edge, and a row's first panel names the reference run on its y axis, for the row. The colour bar,
    `COLOURBAR_WIDTH` across, stands `GAP` right of the first row, as tall as its tallest panel. A first-row cell is
    `PAIR_CELL` inches a side, as the heatmaps of every two runs beside the profiles (`get_run_comparison_panels`); a
    later row's panels share one side, so the row with its labels is as wide as the first row with its colour bar
    (`get_fill_side`); a cell under `PAIR_CELL` a side is logged. The page is as wide and tall as its content; one wider
    than the text column (`page_size("full")`) is logged.

    Args:
        confusions: each compared run's `unsupervised_metrics.compute_cna_confusion`, by its key, in panel order.
        reference: the reference run's key.
        names: each run's name (`get_run_labels`).
        levels: each run's field compared, by its key, e.g. `ascn`.

    Returns:
        The figure, a panel per run.
    """
    colours = LinearSegmentedColormap.from_list("jaccard", list(HEATMAP_COLOURS))
    values, drawn = {}, {}
    for run, confusion in confusions.items():
        jaccard = compute_state_jaccard(confusion)
        called = confusion.loc[jaccard.index, jaccard.columns]
        drawn[run] = [s for s in jaccard.index if called.loc[s].sum() > 0 or called[s].sum() > 0]
        values[run] = jaccard.loc[drawn[run], drawn[run]].to_numpy(dtype=float)
    runs = list(confusions)
    sides = [PAIR_CELL * len(drawn[run]) for run in runs]
    rows = get_confusion_rows([len(drawn[run]) for run in runs])
    firsts = [row[0] for row in rows]
    with page_style() as size:
        figure = plt.figure(figsize=page_size("full"))
        axes = [add_axes(figure, 0.0, 0.0, side, side) for side in sides]
        for k, (ax, run) in enumerate(zip(axes, runs, strict=True)):
            image = ax.imshow(values[run], cmap=colours, vmin=0, vmax=1, aspect="auto")
            ax.set_xticks(range(len(drawn[run])), drawn[run])
            ax.set_yticks(range(len(drawn[run])), drawn[run])
            ax.tick_params(axis="x", labelrotation=90)
            ax.set_xlabel(names[run])
            if k in firsts:
                ax.set_ylabel(names[reference])
        bar_height = max(sides[k] for k in rows[0])
        bar = add_axes(figure, 0.0, 0.0, COLOURBAR_WIDTH, bar_height)
        figure.colorbar(image, cax=bar, label="Jaccard index of the two states' cells")
        set_text_sizes(figure, size)
        for ax, run in zip(axes, runs, strict=True):
            for (i, j), value in np.ndenumerate(values[run]):
                if value >= 0.005:
                    ink = "white" if value > 0.6 else "black"
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center", color=ink, fontsize=size)
        renderer = figure.canvas.get_renderer()
        bar_rooms = get_label_rooms(bar, renderer)
        for _ in range(FIT_ROUNDS):
            rooms = [get_label_rooms(ax, renderer) for ax in axes]
            left = max(rooms[k][0] for k in firsts)
            right = left + sum(sides[k] + rooms[k][1] + GAP for k in rows[0])
            right += sum(rooms[k][0] for k in rows[0][1:]) + bar_rooms[0] + COLOURBAR_WIDTH + bar_rooms[1]
            for row in rows[1:]:
                side = get_fill_side([rooms[k] for k in row], right - left + rooms[row[0]][0])
                for k in row:
                    sides[k] = side
                    place(axes[k], 0.0, 0.0, side, side)
        rooms = [get_label_rooms(ax, renderer) for ax in axes]
        for k in (k for row in rows[1:] for k in row):
            if sides[k] < PAIR_CELL * len(drawn[runs[k]]):
                logger.warning(f"{runs[k]}: state Jaccard cells {sides[k] / len(drawn[runs[k]]):.2f} in a side.")

        frames: list[tuple[float, float]] = [(0.0, 0.0)] * len(runs)
        left = max(rooms[k][0] for k in firsts)
        top, width = 0.0, 0.0
        for r, row in enumerate(rows):
            top += (GAP if r else 0.0) + LETTER
            x = left - rooms[row[0]][0]
            for k in row:
                frames[k] = (x + rooms[k][0], top)
                x += rooms[k][0] + sides[k] + rooms[k][1] + GAP
            if r == 0:
                bar_frame = (x + bar_rooms[0], top)
                x += bar_rooms[0] + COLOURBAR_WIDTH + bar_rooms[1] + GAP
            width = max(width, x - GAP)
            top += max(sides[k] + rooms[k][2] for k in row)
        if width > page_size("full")[0]:
            logger.warning(f"State Jaccard page {width:.2f} in wide, past {page_size('full')[0]:.2f} in.")
        figure.set_size_inches(width, top)
        for ax, (x0, y0), side in zip(axes, frames, sides, strict=True):
            place(ax, x0, y0, side, side)
        place(bar, *bar_frame, COLOURBAR_WIDTH, bar_height)
        logger.info(f"State Jaccard page {width:.2f} x {top:.2f} in, rows {rows}.")
        add_panel_letters(figure, axes, size, [LEVEL_NAMES[levels[run]] for run in runs])
    return figure


def draw_run_pairs(ax: Any, pairwise: pd.DataFrame, metric: str, order: list[str], labels: list[str]) -> Any:
    """Draw a metric of every two runs as a heatmap, a row and a column per run of `order`, each cell's value written in
    it, `-` where NaN, its colour on [0, 1] (`HEATMAP_COLOURS`).

    Args:
        ax: the panel.
        pairwise: `run_a`, `run_b` and `metric` per ordered pair of runs.
        metric: the column drawn, a key of `RUN_PAIR_METRICS`.
        order: the runs, top to bottom and left to right.
        labels: each run's tick label, in `order`.

    Returns:
        The heatmap's image, for its colour bar.
    """
    matrix = pairwise.pivot(index="run_a", columns="run_b", values=metric).reindex(index=order, columns=order)
    values = matrix.to_numpy(dtype=float)
    image = ax.imshow(values, cmap=LinearSegmentedColormap.from_list(metric, list(HEATMAP_COLOURS)), vmin=0, vmax=1)
    ax.set_xticks(range(len(order)), labels)
    ax.set_yticks(range(len(order)), labels)
    for (i, j), value in np.ndenumerate(values):
        text = "-" if np.isnan(value) else f"{value:.2f}"
        ax.text(j, i, text, ha="center", va="center", color="white" if value > 0.6 else "black")
    return image


def get_run_comparison_panels(
    runs: pd.DataFrame, pairs: dict[str, pd.DataFrame], errors: pd.DataFrame
) -> tuple[list[Callable[[Any], None]], float, Callable[[Any], list[Patch]]]:
    """Get the panels of a sample's runs compared beside the methods' profiles (`layout.plot_profile_page`): a heatmap
    per metric of every two runs (`draw_run_pairs`), a row and a column per run in method order (`get_run_order`),
    the runs named under each, upright, and left of the first alone, `PAIR_CELL` a run, each titled as
    `RUN_PAIR_METRICS` titles it; and the runs' skill scores in that order (`draw_skill_score_panel`). A run is
    named by its method where each method has one run, else as `get_run_label` names it.

    Args:
        runs: `run`, `method` and `resolution` per run, `run_benchmark.write_real_data_metrics`' wide table.
        pairs: metric -> its table of every two runs, `run_a`, `run_b` and the metric, a key of `RUN_PAIR_METRICS`.
        errors: `run`, `method` and `ss_<channel>` per run scored on the copy-number error.

    Returns:
        The heatmaps' drawers, their side in inches, and the skill scores' drawer, which returns its key's patches.
    """
    order = get_run_order(runs)
    methods = dict(zip(runs["run"], runs["method"], strict=True))
    single = not runs["method"].duplicated().any()
    names = {run: get_method_name(methods[run]) for run in order} if single else get_run_labels(runs)
    labels = [names[run] for run in order]
    squares: list[Callable[[Any], None]] = [
        partial(draw_run_pair_panel, pairwise=pairs[metric], metric=metric, order=order, labels=labels, first=k == 0)
        for k, metric in enumerate(pairs)
    ]
    return (
        squares,
        PAIR_CELL * len(order),
        partial(draw_skill_score_panel, errors=get_ordered_runs(errors), names=names),
    )


def draw_run_pair_panel(
    ax: Any, pairwise: pd.DataFrame, metric: str, order: list[str], labels: list[str], first: bool
) -> None:
    """Draw a metric of every two runs on a square panel beside the profiles (`draw_run_pairs`): the runs named under
    it, turned a quarter, and left of it where `first`; its title `RUN_PAIR_METRICS`'.

    Args:
        ax: the panel.
        pairwise: `run_a`, `run_b` and `metric` per ordered pair of runs.
        metric: the column drawn, a key of `RUN_PAIR_METRICS`.
        order: the runs, top to bottom and left to right.
        labels: each run's name, in `order`.
        first: whether the runs are named left of it.
    """
    draw_run_pairs(ax, pairwise, metric, order, labels)
    ax.tick_params(axis="x", labelrotation=90)
    ax.tick_params(axis="y", labelleft=first)
    ax.set_title(RUN_PAIR_METRICS[metric])
