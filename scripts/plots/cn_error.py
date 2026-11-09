"""A sample's runs' copy-number error (`metrics.cn_error`): their RDR and BAF mean squared error, the diploid
prediction's, and their skill scores.

- `get_copy_number_rows` -- the runs of the real-data wide table scored on it;
- `plot_cn_error`, `draw_mse`, `draw_skill_scores` -- its page;
- `draw_skill_score_panel` -- the skill scores alone, and their key.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, LogLocator, MultipleLocator, NullFormatter
from scripts.plots.layout import set_text_sizes
from scripts.plots.style import (
    add_panel_letters,
    get_bar_positions,
    get_group_width,
    get_ordered_runs,
    get_panel_rows,
    get_patch,
    get_run_fills,
    pop_titles,
    put_key_under,
    set_axes_height,
)

from cnamaste.figure_style import page_style

CN_ERROR_CHANNELS = (("rdr", "RDR"), ("baf", "BAF"))
"""The copy-number error's channels, as `(suffix, label)`: `mse_<suffix>` and `ss_<suffix>` of `cn_error`'s table."""


SKILL_SCORE_LABEL = r"$1 - \mathrm{MSE} / \mathrm{MSE}_{(1,1)}$"
"""The skill score's axis label, `cn_error.get_skill_score`'s `1 - mse / null`, the null the diploid, (1,1),
prediction's."""


SKILL_SCORE_TOP = 0.75
"""The skill-score axis' top, unless a score passes it."""


SKILL_SCORE_STEP = 0.25
"""The skill-score axis' tick step, so its top, `SKILL_SCORE_TOP`, is a tick."""


CN_ERROR_AXES_HEIGHT = 1.8
"""Inches of a copy-number error panel's axes, taller than `SKILL_SCORE_LABEL` at the axis-label size
(`layout.get_axis_label_size`), the page's height set to it (`style.set_axes_height`)."""


def get_copy_number_rows(wide: pd.DataFrame) -> pd.DataFrame:
    """Get the runs of `analysis.real_data`'s wide table scored on the copy-number error: those whose profile
    gives copies, a count of RDR bins scored.

    Args:
        wide: the wide table, a row per run.

    Returns:
        Its rows of those runs.
    """
    return wide[wide["#rdr-bins"].notna()].reset_index(drop=True)


def draw_mse(ax: Any, errors: pd.DataFrame, fills: dict[str, tuple[float, float, float, float]]) -> None:
    """Draw each run's mean squared error on a log axis ticked per decade, as the channels' errors differ by a decade: a
    group per channel (`CN_ERROR_CHANNELS`) named by the bins scored, a bar per run with an error in its fill from the
    axis' foot, centred on the tick (`get_bar_positions`), its value upright inside a bar of a decade or more, else
    above it, the diploid prediction's error on the run's pseudobulk a black tick across it.

    Args:
        ax: the panel.
        errors: `run`, `method`, then `#<channel>-bins`, `mse_<channel>` and `null_<channel>` per run.
        fills: each run's fill (`get_run_fills`).
    """
    width = get_group_width() / max(len(errors), 1)
    values = errors[[f"{kind}_{c}" for c, _ in CN_ERROR_CHANNELS for kind in ("mse", "null")]].to_numpy(dtype=float)
    drawn = values[np.isfinite(values) & (values > 0)]
    low = float(10.0 ** np.floor(np.log10(drawn.min()))) if drawn.size else 1e-3
    runs = errors["run"].to_numpy()
    for j, (channel, _) in enumerate(CN_ERROR_CHANNELS):
        mse, null = errors[f"mse_{channel}"].to_numpy(dtype=float), errors[f"null_{channel}"].to_numpy(dtype=float)
        found = np.isfinite(mse)
        for i, x in zip(np.flatnonzero(found), get_bar_positions(found, j, width), strict=True):
            ax.bar(x, mse[i] - low, bottom=low, width=width, color=fills[runs[i]])
            ax.hlines(null[i], x - width / 2, x + width / 2, color="black", linewidth=1.5)
            if mse[i] >= 10.0 * low:
                ax.text(x, math.sqrt(low * mse[i]), f"{mse[i]:.3f}", ha="center", va="center", rotation=90)
            else:
                ax.text(x, max(mse[i], null[i]) * 1.1, f"{mse[i]:.3f}", ha="center", va="bottom", rotation=90)
    ax.set_yscale("log")
    ax.set_ylim(low, 2.0 * float(drawn.max()) if drawn.size else 1.0)
    ax.yaxis.set_major_locator(LogLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    bins = [f"{label}\n{int(errors[f'#{c}-bins'].max()):,} bins" for c, label in CN_ERROR_CHANNELS]
    ax.set_xticks(range(len(CN_ERROR_CHANNELS)), bins)
    ax.set_xlim(-0.5, len(CN_ERROR_CHANNELS) - 0.5)
    ax.set_title("Mean squared error")
    ax.set_ylabel("MSE")


def draw_skill_scores(ax: Any, errors: pd.DataFrame, fills: dict[str, tuple[float, float, float, float]]) -> None:
    """Draw each run's copy-number skill score against the diploid prediction, a group per channel
    (`CN_ERROR_CHANNELS`), a bar per run with a score in its fill, centred on the tick (`get_bar_positions`), its value
    upright above it, or above 0 for a bar below 0; the axis' top `SKILL_SCORE_TOP`, or the highest score past it, a
    tick every `SKILL_SCORE_STEP`.

    Args:
        ax: the panel.
        errors: `run`, `method`, then `ss_<channel>` per run (`cn_error.get_cn_error_table`).
        fills: each run's fill (`get_run_fills`).
    """
    width = get_group_width() / max(len(errors), 1)
    scores = errors[[f"ss_{c}" for c, _ in CN_ERROR_CHANNELS]].to_numpy(dtype=float)
    lowest = float(np.nanmin(scores)) if np.isfinite(scores).any() else 0.0
    highest = float(np.nanmax(scores)) if np.isfinite(scores).any() else 0.0
    runs = errors["run"].to_numpy()
    for j in range(len(CN_ERROR_CHANNELS)):
        found = np.isfinite(scores[:, j])
        for i, x in zip(np.flatnonzero(found), get_bar_positions(found, j, width), strict=True):
            ax.bar(x, scores[i, j], width=width, color=fills[runs[i]])
            ax.text(x, max(scores[i, j], 0.0) + 0.02, f"{scores[i, j]:.2f}", ha="center", va="bottom", rotation=90)
    ax.axhline(0.0, color="black", linewidth=0.5)
    ax.set_ylim(min(lowest - 0.05, 0.0), max(SKILL_SCORE_TOP, highest))
    ax.yaxis.set_major_locator(MultipleLocator(SKILL_SCORE_STEP))
    ax.set_xticks(range(len(CN_ERROR_CHANNELS)), [label for _, label in CN_ERROR_CHANNELS])
    ax.set_xlim(-0.5, len(CN_ERROR_CHANNELS) - 0.5)
    ax.set_title("Skill score")
    ax.set_ylabel(SKILL_SCORE_LABEL)


def draw_skill_score_panel(ax: Any, errors: pd.DataFrame, names: dict[str, str]) -> list[Patch]:
    """Draw each run's skill scores (`draw_skill_scores`) in its fill (`get_run_fills`), and get the key naming them.

    Args:
        ax: the panel.
        errors: `run`, `method` and `ss_<channel>` per run.
        names: each run's name (`get_run_labels`).

    Returns:
        A patch per run, in `errors`' order.
    """
    fills = get_run_fills(errors)
    draw_skill_scores(ax, errors, fills)
    return [get_patch(fills[run], names[run]) for run in errors["run"]]


def plot_cn_error(errors: pd.DataFrame, title: str, names: dict[str, str]) -> Figure:
    """Draw each run's copy-number error side by side: its mean squared error (`draw_mse`) and skill scores
    (`draw_skill_scores`), a group per channel (`CN_ERROR_CHANNELS`), a bar per run in its fill (`get_run_fills`), in
    method order (`get_ordered_runs`), the key's order; every panel's runs named by the key alone, three to a row,
    under the panels (`style.get_panel_rows`, `style.put_key_under`); the page as tall as `CN_ERROR_AXES_HEIGHT` a
    panel makes it.

    Args:
        errors: `run`, `method` and `cn_error.get_cn_error_columns()` per run.
        title: the figure's title.
        names: each run's name (`get_run_labels`).

    Returns:
        The figure, the text column wide, one row of panels and the key.
    """
    errors = get_ordered_runs(errors)
    fills = get_run_fills(errors)
    key = [
        *(get_patch(fills[r], names[r]) for r in errors["run"]),
        Line2D([], [], color="black", linewidth=1.5, label="Diploid prediction"),
    ]
    with page_style() as size:
        figure, panels = get_panel_rows(1, key_lines=math.ceil(len(key) / 3))
        axes = list(panels.flat)
        draw_mse(axes[0], errors, fills)
        draw_skill_scores(axes[1], errors, fills)
        legend = figure.legend(handles=key, loc="outside lower center", ncol=3, frameon=False, borderaxespad=0)
        figure.suptitle(title)
        titles = pop_titles(axes)
        set_text_sizes(figure, size)
        set_axes_height(figure, axes[0], CN_ERROR_AXES_HEIGHT)
        put_key_under(figure, legend)
        add_panel_letters(figure, axes, size, titles)
    return figure
