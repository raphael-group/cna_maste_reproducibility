"""The simulation pages: the simulations' metrics per simulated class, their runtime and memory, and cnamaste's
solvers.

- `sim_groups`, `group_order`, `draw_groups`, `not_computed`, `compact_size` -- copied from `cnamaste.qa.plot`;
- `drop_simulated_classes` -- a simulation table less some simulated classes;
- `get_class_medians` -- each method's median per simulated class;
- `plot_simulation_metrics`, `plot_simulation_resources` -- per simulated class, a box per method;
- `plot_solver_ablation` -- per simulated class, a box per cnamaste solver.

The pages draw `metrics.simulation_metrics.get_simulation_metrics`' table as `qa.plot.plot_validation` draws
`run_sim_qa`'s, for any number of methods: within each simulated class (`sim_groups`, ordered by `group_order`) a box
per method side by side, its samples overlaid (`draw_groups`), each method in its colour
(`style.get_method_colour`). A panel no method has a value on reads "not computed". Rates are drawn on [-0.05, 1.05],
runtime and memory on log axes ticked at 1, 2 and 5 per decade; the solver page's normal specificity as
1 - specificity on a log axis, a run at specificity 1 drawn at a floor ticked 0.00 (`get_complement`).
The methods are `style.METHOD_COLOURS`' alone; a cnamaste solver's run beside the main one, its rows' `solver` set, is
drawn on its own page with the main solver (`plot_solver_ablation`), as cnamaste `cleanup3`'s `run_sim_qa` draws its
solvers (`qa.plot.plot_validation`, `ddc060b`), one page for all of them.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter, NullLocator
from matplotlib.transforms import offset_copy
from scripts.cnamaste_port.config import start_time
from scripts.cnamaste_port.figure_style import page_style
from scripts.cnamaste_port.logger import get_logger
from scripts.cnamaste_port.plot.layout import GAP, LETTER
from scripts.metrics.read_method_input import get_method_name
from scripts.metrics.simulation_metrics import (
    get_level_column,
    get_level_metrics,
    get_levels,
    get_runtime_columns,
)
from scripts.plots.layout import set_text_sizes
from scripts.plots.style import (
    add_panel_letters,
    get_group_width,
    get_method_alpha,
    get_method_colour,
    get_method_key,
    get_method_order,
    get_panel_rows,
    get_patch,
    put_key_under,
    set_axes_height,
)

logger = get_logger("simulation", start_time=start_time)


SAMPLE_ID = re.compile(
    r"^numcnas(?P<numcnas>\d+(?:\.\d+)*)_cnasize(?P<cnasize>[\d.]+(?:[eE][+-]?\d+)?)"
    r"_ploidy(?P<ploidy>\d+(?:\.\d+)?)_random(?P<random>\d+)$"
)
"""A simulated sample's id, `numcnas<a.b>_cnasize<s>_ploidy<p>_random<r>` (`qa.plot.SAMPLE_ID`)."""


JITTER = 0.1
"""Half-width of the points' horizontal jitter, in group positions (`qa.plot.JITTER`)."""


METRIC_PANELS = (
    ("specificity", "Normal specificity"),
    ("spot-coverage-rate", "Spot coverage rate"),
    ("sensitivity", "CNA sensitivity"),
    ("clone_ari", "Clone ARI"),
    ("state_altered", "CNA recall"),
    ("exact_altered", r"$\mathbb{N}$-CNA recall"),
)
"""The metrics' panels, row by row, as `(column, label)`: rates, drawn on [-0.05, 1.05]; `spot-coverage-rate` is
one value at every level (`simulation_metrics.get_level_metrics` leaves it out)."""


METRIC_AXES_HEIGHT = 1.43
"""Inches of a metric panel's axes (`plot_simulation_metrics`), the page's height set to it (`set_axes_height`)."""


COMPLEMENT = "1 - specificity"
"""The solver page's column of 1 - normal specificity, drawn on a log axis (`get_complement`)."""


SOLVER_PANELS = (
    ("sensitivity", "CNA sensitivity"),
    (COMPLEMENT, r"1 $-$ normal specificity"),
    ("exact_altered", r"$\mathbb{N}$-CNA recall"),
    ("clone_ari", "Clone ARI"),
)
"""The solver page's panels, row by row, as `(column, label)`: cnamaste `cleanup3`'s `qa.plot.PANELS` (`ddc060b`),
its normal specificity as 1 - specificity on a log axis, its N-clone ARI the clone ARI of the integer clones,
`clone_labels_int.tsv`, every cnamaste run's."""


FLOOR_TICKS = {1e-3: r"$10^{-3}$", 1e-2: r"$10^{-2}$", 1e-1: r"$10^{-1}$", 1.0: "1.00"}
"""A log rate panel's ticks above its floor, by value (`set_floor_axis`)."""


RESOURCE_PANELS = (("Runtime [mins]", "Runtime (min)"), ("Memory [GB]", "Memory (GB)"))
"""The resources' panels, as `(column, label)`: drawn on log axes."""


LOG_TICKS = 12
"""The most ticks a log axis's locator places, so a short axis over decades keeps its 1, 2 and 5 ticks: matplotlib's
`LogLocator` thins its ticks to the axis length where `numticks` is `"auto"`, to none at 0.7 in over 2.5 decades."""


SOLVER_COLOURS = {
    "expansion-lattice-lattice": get_method_colour("cnamaste"),
    "expansion-calicost-calicost": "#8c564b",
    "opengm_icm-lattice-lattice": "#332288",
    "opengm_icm-calicost-calicost": "#88ccee",
}
"""Each cnamaste solver's colour, in drawing order: the main solver, `batch_zenodo_run.sh`'s `QA_SOLVER`, cnamaste's
own; the others outside `METHOD_COLOURS`, every two at least 20 apart in OKLab x 100 under normal, protan, deutan and
tritan vision (Machado et al. 2009, doi:10.1109/TVCG.2009.113)."""


def compact_size(size: float) -> str:
    """Get a CNA size in bases as `10M`, `50M`, `1.5M`, `500K` (`qa.plot.compact_size`).

    Args:
        size: the size in bases.

    Returns:
        The size with a K, M or G suffix where at least a thousand.
    """
    for scale, unit in ((1e9, "G"), (1e6, "M"), (1e3, "K")):
        if size >= scale:
            return f"{size / scale:g}{unit}"
    return f"{size:g}"


def sim_groups(sample_ids: Iterable[str]) -> pd.DataFrame:
    """Get simulated sample ids' parameters and group; an id that does not parse is dropped, logged
    (`qa.plot.sim_groups`).

    Args:
        sample_ids: ids such as `numcnas1.2_cnasize5e7_ploidy2_random0`, repeats allowed.

    Returns:
        One row per parsed id, in input order: `sample_id`; `numcnas` the string (`"1.2"`); `cnasize`, `ploidy`
        floats; `random` int; `group` the id less `_random<r>`; `label` `(1,2, 50M)`.
    """
    rows = []
    for sample_id in sample_ids:
        found = SAMPLE_ID.match(str(sample_id))
        if found is None:
            logger.warning(f"{sample_id}: not numcnas<a.b>_cnasize<s>_ploidy<p>_random<r>, left out of the figure")
            continue
        numcnas, cnasize = found["numcnas"], float(found["cnasize"])
        rows.append(
            {
                "sample_id": str(sample_id),
                "numcnas": numcnas,
                "cnasize": cnasize,
                "ploidy": float(found["ploidy"]),
                "random": int(found["random"]),
                "group": str(sample_id)[: found.start("random") - len("_random")],
                "label": f"({numcnas.replace('.', ',')}, {compact_size(cnasize)})",
            }
        )
    return pd.DataFrame(rows, columns=["sample_id", "numcnas", "cnasize", "ploidy", "random", "group", "label"])


def drop_simulated_classes(table: pd.DataFrame, classes: tuple[tuple[str, float], ...]) -> pd.DataFrame:
    """Get a simulation table's rows outside `classes`.

    Args:
        table: the simulation table, `sample_id` per row.
        classes: the simulated classes left out, `(numcnas, cnasize)` as `sim_groups` reads them, e.g. `("6.3", 5e7)`.

    Returns:
        The other rows, in order.
    """
    groups = sim_groups(table["sample_id"].unique())
    dropped = [(numcnas, cnasize) in classes for numcnas, cnasize in zip(groups["numcnas"], groups["cnasize"])]
    return table[~table["sample_id"].isin(groups.loc[dropped, "sample_id"])]


def group_order(groups: pd.DataFrame) -> list[str]:
    """Get the groups in drawing order: `cnasize` largest first, then the sum of `numcnas`' components largest first,
    then by name (`qa.plot.group_order`).

    Args:
        groups: `sim_groups`' result.

    Returns:
        Each group once.
    """
    unique = groups.drop_duplicates("group").assign(
        total=lambda d: d["numcnas"].map(lambda n: sum(int(p) for p in n.split(".")))
    )
    return list(unique.sort_values(["cnasize", "total", "group"], ascending=[False, False, True])["group"])


def get_class_medians(table: pd.DataFrame) -> pd.DataFrame:
    """Get each method's median per simulated class over its replicates, NaN skipped: the line each box of the
    simulation pages is drawn at.

    Args:
        table: the simulations' wide table (`simulation_metrics.merge_levels`).

    Returns:
        One row per simulated class and method, classes in `group_order`, methods in `get_method_order`'s order then
        the table's: `sim_id`, the class (`sim_groups`' `group`); `method`; `#replicates`, its runs in the class; then
        the medians of `spot-coverage-rate`, each level's metrics (`simulation_metrics.get_level_column`) and the
        runtime and memory.
    """
    groups = sim_groups(table["sample_id"].unique())
    columns = [
        "spot-coverage-rate",
        *(get_level_column(metric, level) for level in get_levels() for metric in get_level_metrics()),
        *get_runtime_columns(),
    ]
    rows = table.merge(groups[["sample_id", "group"]], on="sample_id").rename(columns={"group": "sim_id"})
    by = rows.groupby(["sim_id", "method"])
    medians = by[columns].median().join(by.size().rename("#replicates")).reset_index()
    methods = get_method_order(table, None)
    methods += [method for method in dict.fromkeys(table["method"]) if method not in methods]
    ranks = {"sim_id": group_order(groups), "method": methods}
    medians = medians.sort_values(["sim_id", "method"], key=lambda s: s.map(ranks[str(s.name)].index))
    logger.info(f"Medians of {len(methods)} methods over {len(ranks['sim_id'])} simulated classes.")
    return medians[["sim_id", "method", "#replicates", *columns]].reset_index(drop=True)


def get_complement(values: pd.Series) -> tuple[pd.Series, float]:
    """Get 1 - value of rates, and the floor a rate of 1 is drawn at on a log axis, where 0 has no place: the power of
    10 a decade under the lowest positive complement's decade, 0.1 where none is positive.

    Args:
        values: the rates, at most 1, NaN where not computed.

    Returns:
        The complements, 0 at the floor, NaN kept; and the floor.
    """
    complement = 1 - pd.to_numeric(values, errors="coerce")
    positive = complement[complement > 0]
    floor = float(10.0 ** (np.floor(np.log10(positive.min())) - 1)) if len(positive) else 0.1
    return complement.mask(complement == 0, floor), floor


def set_floor_axis(ax: Any, floor: float) -> None:
    """Set a rate panel's y axis to a log axis from `floor` to 1, as rates are drawn on [0, 1]: from half the floor to
    2, the floor ticked `0.00` and `FLOOR_TICKS` above it, no minor ticks.

    Args:
        ax: the panel, drawn.
        floor: its floor (`get_complement`).
    """
    ax.set_yscale("log")
    ax.set_ylim(floor / 2, 2)
    ticks = {floor: "0.00", **{value: label for value, label in FLOOR_TICKS.items() if value > floor}}
    ax.set_yticks(list(ticks), list(ticks.values()))
    ax.yaxis.set_minor_locator(NullLocator())


def not_computed(ax: Any) -> None:
    """Write "not computed" across the middle of `ax` (`qa.plot.not_computed`)."""
    ax.text(0.5, 0.5, "not computed", transform=ax.transAxes, ha="center", va="center", color="0.4")
    ax.set_yticks([])


def draw_groups(
    ax: Any, values: list[np.ndarray], seed: int, offset: float, width: float, box: str, whisker: str
) -> None:
    """Draw a box per non-empty group at its position `0, 1, ...` plus `offset`, its samples as jittered points
    (`qa.plot.draw_groups`).

    Args:
        ax: the panel.
        values: each group's samples, NaN dropped.
        seed: the jitter's generator seed.
        offset: the boxes' offset in group positions.
        width: the boxes' width in group positions.
        box: the box colour.
        whisker: the whisker colour.
    """
    rng = np.random.default_rng(seed)
    drawn = [k for k, v in enumerate(values) if v.size]
    if not drawn:
        return
    ax.boxplot(
        [values[k] for k in drawn],
        positions=[k + offset for k in drawn],
        widths=width,
        showfliers=False,
        patch_artist=True,
        boxprops={"facecolor": box, "alpha": get_method_alpha(), "edgecolor": whisker, "linewidth": 0.5},
        whiskerprops={"color": whisker, "linewidth": 0.5},
        capprops={"color": whisker, "linewidth": 0.5},
        medianprops={"color": "black", "linewidth": 0.8},
        manage_ticks=False,
    )
    for k in drawn:
        x = k + offset + rng.uniform(-JITTER * width / 0.6, JITTER * width / 0.6, size=values[k].size)
        ax.scatter(x, values[k], s=3, color="black", alpha=0.7, edgecolors="black", linewidths=0.2, zorder=3)


def draw_panels(
    axes: Any,
    table: pd.DataFrame,
    methods: list[str],
    panels: tuple[tuple[str, str], ...],
    series: str = "method",
    colours: dict[str, str] | None = None,
    floors: dict[str, float] | None = None,
) -> None:
    """Draw one column of `table` per panel, a box per simulated class and series, the classes along x.

    Args:
        axes: the panels, as many as `panels`.
        table: the simulation table, `sample_id`, `series` and `panels`' columns.
        methods: the series, in order, side by side in each class.
        panels: `(column, label)` per panel.
        series: the column naming each row's series, `method` by default.
        colours: each series' colour, or None for its method's (`get_method_colour`).
        floors: a log column's floor (`get_complement`), its panel's axis set by `set_floor_axis`, or None.
    """
    groups = sim_groups(table["sample_id"].unique())
    order, labels = group_order(groups), dict(zip(groups["group"], groups["label"], strict=True))
    table = table.merge(groups[["sample_id", "group"]], on="sample_id")
    width = get_group_width() / len(methods)
    for k, (ax, (column, label)) in enumerate(zip(axes, panels, strict=True)):
        ax.set_ylabel(label)
        ax.set_xlim(-0.5, max(len(order), 1) - 0.5)
        floor = (floors or {}).get(column)
        drawn = False
        for i, method in enumerate(methods):
            rows = table[table[series] == method]
            found = pd.to_numeric(rows[column], errors="coerce") if column in rows else pd.Series(dtype=float)
            values = [found[rows["group"] == g].dropna().to_numpy(float) for g in order]
            if any(v.size for v in values):
                colour = get_method_colour(method) if colours is None else colours[method]
                offset = (i - (len(methods) - 1) / 2) * width
                draw_groups(ax, values, seed=k, offset=offset, width=0.9 * width, box=colour, whisker=colour)
                drawn = True
        if not drawn:
            not_computed(ax)
        elif floor is not None:
            set_floor_axis(ax, floor)
        elif column in dict(RESOURCE_PANELS):
            ax.set_yscale("log")
            ax.yaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0), numticks=LOG_TICKS))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
            ax.yaxis.set_minor_formatter(NullFormatter())
        else:
            ax.set_ylim(-0.05, 1.05)
        ax.set_xticks(range(len(order)), [labels[g] for g in order], rotation=45, ha="right")
        ax.set_xlabel("Simulated class")


def plot_simulation_metrics(
    table: pd.DataFrame, methods: list[str] | None = None, columns: tuple[str, ...] | None = None
) -> Figure:
    """Draw each method's metrics per simulated class (`METRIC_PANELS`), those of `columns` alone where given, two a
    row (`get_panel_rows`), the classes named under the last row alone, no x-axis label; keyed `GAP` right of the
    last panel, one method a line, centred on its axes; untitled, each panel's letter labelled by its metric and, for a
    metric scored per level (`simulation_metrics.get_level_metrics`), the level the runs are scored at
    (`add_panel_letters`). Each panel is padded `LETTER / 2` above and below, so two rows are `LETTER` apart past
    their labels, room for a letter; the page is as tall as `METRIC_AXES_HEIGHT` a panel makes it
    (`set_axes_height`).

    Args:
        table: `simulation_metrics.get_simulation_metrics`' rows.
        methods: the methods to draw, in order, or None (`get_method_order`).
        columns: the metrics drawn, an even number of `METRIC_PANELS`' columns, or None for all of them.

    Returns:
        The figure, the text column wide, panels lettered (a), (b), ... in `METRIC_PANELS`' order.
    """
    methods = get_method_order(table, methods)
    levels = ", ".join(dict.fromkeys(table["level"].astype(str)))
    panels = tuple(panel for panel in METRIC_PANELS if columns is None or panel[0] in columns)
    with page_style() as size:
        figure, axes = get_panel_rows(len(panels) // 2, sharex=True)
        figure.get_layout_engine().set(h_pad=LETTER / 2, hspace=0.0)
        draw_panels(list(axes.flat), table, methods, panels)
        for ax in axes.flat:
            ax.set_xlabel("")
        last = axes[-1, -1]
        last.legend(
            handles=get_method_key(methods),
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
            bbox_transform=offset_copy(last.transAxes, figure, x=GAP, units="inches"),
            frameon=False,
            borderaxespad=0,
        )
        set_text_sizes(figure, size)
        set_axes_height(figure, axes[0, 0], METRIC_AXES_HEIGHT)
        titles = [
            f"{label}, simulations scored per {levels}" if column in get_level_metrics() else f"{label}, simulations"
            for column, label in panels
        ]
        add_panel_letters(figure, list(axes.flat), size, titles)
    return figure


def plot_simulation_resources(table: pd.DataFrame, methods: list[str] | None = None) -> Figure:
    """Draw each method's runtime and memory per simulated class (`RESOURCE_PANELS`), side by side (`get_panel_rows`),
    keyed under them; a method with neither left out; each panel's letter labelled by its label (`add_panel_letters`).

    Args:
        table: `simulation_metrics.get_simulation_metrics`' rows.
        methods: the methods to draw, in order, or None (`get_method_order`).

    Returns:
        The figure, the text column wide, one `PANEL_ROW` row, panels lettered (a)-(b).
    """
    columns = [column for column, _ in RESOURCE_PANELS]
    measured = set(table.loc[table[columns].notna().any(axis=1), "method"])
    methods = [m for m in get_method_order(table, methods) if m in measured]
    with page_style() as size:
        figure, axes = get_panel_rows(1, columns=len(RESOURCE_PANELS), key_lines=1)
        draw_panels(list(axes.flat), table, methods, RESOURCE_PANELS)
        key = figure.legend(
            handles=get_method_key(methods),
            loc="outside lower center",
            ncol=len(methods),
            frameon=False,
            borderaxespad=0,
        )
        set_text_sizes(figure, size)
        figure.canvas.draw()
        put_key_under(figure, key)
        add_panel_letters(figure, list(axes.flat), size, [f"{label}, simulations" for _, label in RESOURCE_PANELS])
    return figure


def plot_solver_ablation(table: pd.DataFrame) -> Figure:
    """Draw cnamaste's solvers' metrics per simulated class (`SOLVER_PANELS`), 2 x 2 (`get_panel_rows`), a box per
    solver side by side in `SOLVER_COLOURS`' order and colour, keyed under them, untitled, each panel's letter labelled
    by its metric and the level the runs are scored at (`add_panel_letters`). Normal specificity is drawn as
    1 - specificity on a log axis, a run at specificity 1 at the panel's floor, ticked 0.00 (`get_complement`).

    Args:
        table: `simulation_metrics.get_simulation_metrics`' rows, a cnamaste run's `solver` per row, NaN for another
            method's.

    Returns:
        The figure, the text column wide, `PANEL_ROW` a row, panels lettered (a)-(d).

    Raises:
        ValueError: where a solver has no colour in `SOLVER_COLOURS`.
    """
    rows = table[table["solver"].notna()]
    unknown = sorted(set(rows["solver"]) - set(SOLVER_COLOURS))
    if unknown:
        raise ValueError(f"solvers {unknown} have no colour in SOLVER_COLOURS")
    solvers = [solver for solver in SOLVER_COLOURS if solver in set(rows["solver"])]
    levels = ", ".join(dict.fromkeys(rows["level"].astype(str)))
    name = get_method_name("cnamaste")
    complement, floor = get_complement(rows["specificity"])
    rows = rows.assign(**{COMPLEMENT: complement})
    with page_style() as size:
        figure, axes = get_panel_rows(2, sharex=True, key_lines=2)
        draw_panels(
            list(axes.flat),
            rows,
            solvers,
            SOLVER_PANELS,
            series="solver",
            colours=SOLVER_COLOURS,
            floors={COMPLEMENT: floor},
        )
        for ax in axes[:-1].flat:
            ax.set_xlabel("")
        key = [get_patch(to_rgba(SOLVER_COLOURS[s], get_method_alpha()), f"{name} {s}") for s in solvers]
        legend = figure.legend(handles=key, loc="outside lower center", ncol=2, frameon=False, borderaxespad=0)
        set_text_sizes(figure, size)
        figure.canvas.draw()
        put_key_under(figure, legend)
        titles = [f"{label}, cnamaste's solvers, simulations scored per {levels}" for _, label in SOLVER_PANELS]
        add_panel_letters(figure, list(axes.flat), size, titles)
    return figure
