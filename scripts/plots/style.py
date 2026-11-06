"""The benchmark figures' shared style: a method's colour and fill, a run's name and order, the keys, the panel
letters, and a page of rows of panels laid out tight.

- `get_method_colour`, `get_method_fill`, `get_run_fills`, `get_patch`, `get_method_order`, `get_method_key`,
  `get_setting_label`, `get_run_label`, `get_run_labels` -- a method's colour and fill, a run's fill and name, and the
  methods' key;
- `get_run_order`, `get_ordered_runs` -- a sample's runs in method order;
- `pop_titles`, `add_panel_letters` -- each panel's letter alone, its title kept for the page naming every panel;
- `get_panel_rows`, `put_key_under`, `set_axes_height`, `get_bar_positions` -- a page of rows of panels laid out
  tight (`docs/figure-rule.md` FP-5), its key moved under them, its height set by its panels', and a group's bars side
  by side.

Every box, bar and key patch of a metric figure is in its method's fill (`get_method_fill`), its colour at
`METHOD_ALPHA`, so a method is one colour on every page. A sample's run is named by its method, the side of the bins
it ran on and its setting, e.g. `CalicoST (16 um, min-snp-umi: 10)` (`get_run_label`).
"""

from __future__ import annotations

import re
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from scripts.metrics.read_method_input import get_method_name

from cnamaste.config import start_time
from cnamaste.figure_style import column_size, page_size
from cnamaste.logger import get_logger
from cnamaste.plot.layout import GAP, LETTER

logger = get_logger("style", start_time=start_time)


MIN_SNP_UMI = re.compile(r"^msu(\d+)$")
"""A CalicoST run's setting, `msu<k>`: its `min_snpumi_perspot`
(https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/arg_parse.py#L30)."""


METHOD_COLOURS = {"cnamaste": "#2a78d6", "calicost": "#eb6834", "numbat": "#1baf7a", "clonalscope": "#eda100"}
"""Each method's colour, the reference categorical palette's first four slots in `NAMES`' order, fixed per method."""


METHOD_ALPHA = 0.4
"""A method's fill opacity, its boxes', bars' and key patches' (`qa.plot.draw_groups`' box alpha)."""


GROUP_WIDTH = 0.8
"""A group's width in group positions, shared by its methods' boxes or bars."""


LEGEND_ROOM = 0.12
"""A figure's fraction under its panels for its key."""


KEY_LINE = 0.25
"""Inches of a line of a key under a page's panels (`get_panel_rows`), so the key adds to the page, not the rows."""


PANEL_ROW = 2.3
"""Inches of a row of panels (`get_panel_rows`), its letter, axes and labels: appendix page 1's per-spot UMI,
diagnosis page 2's false-positive spots' UMI and the simulation pages."""


def get_method_alpha() -> float:
    """Get a method's fill opacity, `METHOD_ALPHA`."""
    return METHOD_ALPHA


def get_group_width() -> float:
    """Get a group's width in group positions, `GROUP_WIDTH`."""
    return GROUP_WIDTH


def get_legend_room() -> float:
    """Get a figure's fraction under its panels for its key, `LEGEND_ROOM`."""
    return LEGEND_ROOM


def get_method_colour(method: str) -> str:
    """Get a method's colour: `METHOD_COLOURS`', grey for a method not in it.

    Args:
        method: e.g. `numbat`.

    Returns:
        The colour, hex.
    """
    return METHOD_COLOURS.get(method, "#808080")


def get_method_fill(method: str) -> tuple[float, float, float, float]:
    """Get a method's fill: its colour (`get_method_colour`) at `METHOD_ALPHA`, as its boxes are drawn.

    Args:
        method: e.g. `numbat`.

    Returns:
        The colour, RGBA.
    """
    return to_rgba(get_method_colour(method), METHOD_ALPHA)


def get_run_fills(runs: pd.DataFrame) -> dict[str, tuple[float, float, float, float]]:
    """Get each run's fill: its method's (`get_method_fill`) where it is its method's one run, else its method's colour
    at opacities from `METHOD_ALPHA / 2` to 1 in the runs' order, so a sweep's runs share their method's hue.

    Args:
        runs: `run` and `method` per run.

    Returns:
        The fill per run key, RGBA.
    """
    fills = {}
    for method, group in runs.groupby("method", sort=False):
        keys = list(group["run"])
        alphas = [METHOD_ALPHA] if len(keys) == 1 else np.linspace(METHOD_ALPHA / 2, 1.0, len(keys))
        fills.update({k: to_rgba(get_method_colour(str(method)), float(a)) for k, a in zip(keys, alphas, strict=True)})
    return fills


def get_patch(fill: tuple[float, float, float, float], label: str) -> Patch:
    """Get a key's patch of a fill (`get_method_fill`, `get_run_fills`).

    Args:
        fill: the patch's colour, RGBA.
        label: the patch's label.

    Returns:
        The patch.
    """
    return Patch(facecolor=fill, edgecolor=fill, label=label)


def get_method_order(table: pd.DataFrame, methods: list[str] | None) -> list[str]:
    """Get the methods to draw: `methods`, or every method of `METHOD_COLOURS` in `table`, in its order; another key,
    e.g. a cnamaste solver's run, left out.

    Args:
        table: a table with `method` per row.
        methods: the methods, in order, or None.

    Returns:
        The methods, in order.
    """
    if methods:
        return list(methods)
    found = set(table["method"])
    return [method for method in METHOD_COLOURS if method in found]


def get_method_key(methods: list[str]) -> list[Patch]:
    """Get the methods' key: a patch per method in its fill.

    Args:
        methods: the methods, in order.

    Returns:
        One patch per method.
    """
    return [get_patch(get_method_fill(m), get_method_name(m)) for m in methods]


def get_setting_label(setting: str) -> str:
    """Get a run's setting as a figure writes it: `msu<k>` (`MIN_SNP_UMI`) as `min-snp-umi: <k>`, any other as given.

    Args:
        setting: the setting, e.g. `msu10`.

    Returns:
        The label, e.g. `min-snp-umi: 10`.
    """
    found = MIN_SNP_UMI.match(setting)
    return f"min-snp-umi: {found.group(1)}" if found else setting


def get_run_label(run: str, method: str, side_um: int) -> str:
    """Get how a figure names a run: its method's name, then the side of the bins it ran on and the run key's
    `msu<k>` setting (`MIN_SNP_UMI`), as `get_setting_label` writes it; every other part of the key, e.g. a run's id,
    is left out, so a method is named alike whatever its key.

    Args:
        run: the run key, e.g. `calicost_016um_msu10`.
        method: its method, e.g. `calicost`.
        side_um: the side in um of the bins it ran on.

    Returns:
        The name, e.g. `CalicoST (16 um, min-snp-umi: 10)`, `cn-maste (16 um)` of `cnamaste-ac88605f`.
    """
    side = rf"{side_um} $\mu$m"
    parts = re.split(r"[_-]", run.removeprefix(method))
    settings = [get_setting_label(part) for part in parts if MIN_SNP_UMI.match(part)]
    return f"{get_method_name(method)} ({', '.join([side, *settings])})"


def get_run_labels(table: pd.DataFrame) -> dict[str, str]:
    """Get each run's name (`get_run_label`).

    Args:
        table: `run`, `method` and `resolution`, `<side>um`, per run, e.g. `compare_real_data.py`'s false discovery.

    Returns:
        The name per run key.
    """
    return {
        run: get_run_label(run, method, int(str(side).removesuffix("um")))
        for run, method, side in table[["run", "method", "resolution"]].itertuples(index=False)
    }


def get_run_order(runs: pd.DataFrame) -> list[str]:
    """Get the runs in method order (`get_method_order`), a method's runs together, in the table's order.

    Args:
        runs: `run` and `method` per run.

    Returns:
        The run keys, in order.
    """
    return [run for method in get_method_order(runs, None) for run in runs.loc[runs["method"] == method, "run"]]


def get_ordered_runs(runs: pd.DataFrame) -> pd.DataFrame:
    """Get `runs`' rows in method order (`get_run_order`), so bars and keys follow the methods' order.

    Args:
        runs: `run` and `method` per run, a row each.

    Returns:
        The rows, in order, indexed from 0.
    """
    rank = {run: k for k, run in enumerate(get_run_order(runs))}
    return runs.iloc[np.argsort(runs["run"].map(rank).to_numpy(), kind="stable")].reset_index(drop=True)


def pop_titles(axes: list[Any]) -> list[str]:
    """Clear each panel's title, so its page lays out without it, and get them, for `add_panel_letters`.

    Args:
        axes: the panels, in letter order.

    Returns:
        Each panel's title, empty for none.
    """
    titles = [ax.get_title() for ax in axes]
    for ax in axes:
        ax.set_title("")
    return titles


def add_panel_letters(figure: Figure, axes: list[Any], size: float, titles: list[str] | None = None) -> None:
    """Letter each panel (a), (b), ... at its top left, left of its tick labels, as `qa.plot.plot_validation` does;
    a panel's title is its letter's label, named on the page of panel titles (`layout.plot_panel_titles`).

    Args:
        figure: the laid-out figure.
        axes: its panels, in letter order.
        size: the letters' font size.
        titles: each panel's title, in letter order, or None for none.
    """
    renderer = figure.canvas.get_renderer()
    width = figure.get_size_inches()[0] * figure.dpi
    for ax, name, title in zip(axes, "abcdefghijkl", titles or [""] * len(axes)):
        x0 = ax.get_tightbbox(renderer).x0 / width
        figure.text(x0, ax.get_position().y1, f"({name})", label=title, va="bottom", ha="left", fontsize=size)


def get_panel_rows(rows: int, columns: int = 2, sharex: bool = False, key_lines: int = 0) -> tuple[Figure, np.ndarray]:
    """Get a page of `rows` rows of `columns` panels, the text column wide, `PANEL_ROW` a row and `KEY_LINE` a line of
    the key under them, laid out by
    matplotlib's constrained layout, each panel padded `GAP` left and right and `LETTER` above and below, past its
    labels, so a panel's letter (`add_panel_letters`) fits over it (`docs/figure-rule.md` FP-5); a page taller than the
    text block less its caption (`page_size("full")`) is logged.

    Args:
        rows: the rows.
        columns: the panels a row.
        sharex: whether the panels share their x axis, its tick labels under the last row alone.
        key_lines: the lines of the key under the panels, `figure.legend(loc="outside lower center")`.

    Returns:
        The page and its axes, `(rows, columns)`.
    """
    height = PANEL_ROW * rows + KEY_LINE * key_lines
    if height > page_size("full")[1]:
        logger.warning(f"{rows} rows of panels take {height:.2f} in, past {page_size('full')[1]:.2f} in.")
    figure, axes = plt.subplots(
        rows, columns, figsize=column_size(height), layout="constrained", squeeze=False, sharex=sharex
    )
    figure.get_layout_engine().set(w_pad=GAP, h_pad=LETTER)
    return figure, axes


def put_key_under(figure: Figure, legend: Any) -> None:
    """Fix a laid-out page's layout and move its key `GAP` under its panels' lowest labels, so the page cropped to its
    content (`pages.SAVE_OPTIONS`) keeps none of the room the layout engine left round the key (FP-5).

    Args:
        figure: the page, drawn once, its key `legend` placed under the panels.
        legend: the key.
    """
    renderer = figure.canvas.get_renderer()
    height = figure.get_figheight()
    bottom = min(ax.get_tightbbox(renderer).y0 for ax in figure.axes) / (figure.dpi * height)
    figure.set_layout_engine("none")
    legend.set_loc("upper center")
    legend.set_bbox_to_anchor((0.5, bottom - GAP / height), transform=figure.transFigure)


def set_axes_height(figure: Figure, ax: Any, height: float) -> None:
    """Set the height of a page laid out by constrained layout, `hspace` 0, so `ax` is `height` inches tall: its
    labels and pads are fixed in inches, so the page changes by `ax`'s change times its grid's rows; drawn twice.

    Args:
        figure: the page, its texts at their final sizes.
        ax: a panel of its grid, its rows as tall as each other.
        height: the panel's height, inches.
    """
    figure.canvas.draw()
    drawn = ax.get_position().height * figure.get_figheight()
    figure.set_figheight(figure.get_figheight() + ax.get_gridspec().nrows * (height - drawn))
    figure.canvas.draw()


def get_bar_positions(drawn: np.ndarray, centre: float, width: float) -> np.ndarray:
    """Get the positions of a group's drawn bars, side by side in the runs' order, centred on the group's tick.

    Args:
        drawn: True per run with a value in the group.
        centre: the group's tick.
        width: a bar's width.

    Returns:
        One position per drawn run.
    """
    count = int(drawn.sum())
    return centre + (np.arange(count) - (count - 1) / 2) * width
