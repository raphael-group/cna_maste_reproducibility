"""`plot_annotation_page`: one Visium HD sample's H&E slide, cell types and each method's clones, and each method's clone
profile on cnplot's `GenomeAxis`, together or a page each; a method not run is a placeholder. Integer copies are
drawn in cnamaste's profile colours (`cnamaste.palette`, `plot_copy_number_profile.state_fill`), (1,1) white
(`NEUTRAL_COPY_COLOUR`). `plot_profile_page`: each method's clone profile on one set of bins as cnamaste's own page
draws its profile (`plot_copy_number_profile.plot_state_profile`; `draw_state_profile`, `draw_copy_profile`), beside
panels it is given.
"""

from __future__ import annotations

import re
import string
from collections.abc import Callable
from fractions import Fraction
from typing import Any, NamedTuple

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.transforms import Affine2D, offset_copy
from scipy.spatial import cKDTree
from scripts.cnamaste_port.config import start_time
from scripts.cnamaste_port.figure_style import CLONE_PALETTE, DPI, page_size, page_style, uniform
from scripts.cnamaste_port.logger import get_logger
from scripts.cnamaste_port.plot.layout import (
    GAP,
    MARGIN,
    TICK_LINE,
    GenomicAxis,
    add_axes,
    draw_slide,
    frame_extent,
    framed,
    place,
    square_handle,
)
from scripts.cnamaste_port.plot.plot_copy_number_profile import (
    BAR_SHARE,
    OUTLINE,
    get_intervals,
    state_fill,
    used_states,
)
from scripts.cnamaste_port.utils import cast_clone_label

# NB cnplot writes its own rcParams on import (`cnplot_style.apply_style`); the context restores cnamaste's.
with mpl.rc_context():
    from cnplot import GenomeAxis, decorate_genome_axis, plot_heatmap

logger = get_logger(__name__, start_time=start_time)


def get_lattice(xy: np.ndarray) -> tuple[float, float]:
    """The spots' square lattice: its pitch, the median distance to a spot's nearest neighbour (the spacing `slide_grid`
    reads a tissue pixel by), and its turn, the mean direction of those neighbours, modulo 90 degrees.

    Args: `xy` (n_spots, 2) the spots' places.
    Returns: the pitch, in `xy`'s units, and the turn in degrees.
    """
    distance, nearest = cKDTree(xy).query(xy, k=2)
    pitch = float(np.median(distance[:, 1]))
    step = xy[nearest[:, 1]] - xy
    # NB the mean direction modulo 90 degrees, on the circle: 4 theta, so 0 and 90 degrees are one direction.
    # NB rounded before the modulo, so a direction a rounding error below 0 reads 0 rather than 90.
    turn = float(
        np.round(np.degrees(np.angle(np.mean(np.exp(4j * np.arctan2(step[:, 1], step[:, 0]))))) / 4.0, 6) % 90.0
    )
    return pitch, turn


TITLE_LINE = 0.16
"""Inches above a panel of `plot_annotation_page` for its title."""

AXIS_LABEL_SCALE = 1.25
"""An axis label's font size over the page's other text: axis labels the largest text, keys and tick labels the page's
size."""

AXIS_LABEL_LINE = 0.16
"""Inches outside a spatial grid's tick labels for its axis label."""

CLONE_ROW = 0.2
"""Inches of profile per tumour clone on `plot_annotation_page`."""

PROFILE_GAP = 0.05
"""Inches between two methods' profiles on `plot_annotation_page`."""

PROFILE_LETTER_GAP = 0.06
"""Inches between a profile's letter and the column of clone names right of it."""

CHR_OVERHANG = 0.2
"""Inches right of the profiles, where the last chromosome's name, at its chromosome's middle, may run past the genome's
end."""

CHR_NAME_GAP = 0.04
"""Inches between two chromosome names on one row under the profiles (`get_chromosome_name_places`)."""

CHR_NAME_ROWS = 2
"""The rows of chromosome names under the profiles (`get_chromosome_name_places`)."""

CHR_NAME_PAD = 0.04
"""Inches between the last profile and the first row of chromosome names under it."""

KEY_SWATCH_HEIGHT = 0.08
"""Inches: a key swatch's height."""

KEY_ROW = 0.24
"""Inches: a row of a key (`draw_key`), its swatches and their names under them."""

KEY_GROUP_GAP = 0.25
"""The least gap between a key's groups, in swatch widths (`get_key_rows`)."""

KEY_LABEL_ROOM = 1.2
"""A key swatch's width over its widest name (`get_key_swatch`)."""

KEY_NAME = "CNA"
"""The name beside the profiles' first key, the state tables' states, else the integer profiles' `(A,B)` states; the
`(A,B)` key under the states' is unnamed, the one name covering both."""

TITLE_LIST_LINE = 0.18
"""Inches per line of `plot_panel_titles`."""

TITLE_WORD = re.compile(r"(?:\$[^$]*\$|[^\s$])+")
"""A word of a panel title: a run of characters without a space, a mathtext span `$...$` whole, spaces and all."""

PANEL_LETTER = re.compile(r"\([a-z]\)")
"""A panel's letter as drawn, `(a)`, `(b)`, ...; its label, where set, the panel's title (`get_panel_titles`)."""

METHOD_NAME_LINE = 0.14
"""Inches over a profile of `plot_profile_page` for its method's name."""

GENOME_NAME_ROOM = 0.3
"""Inches under `plot_profile_page`'s last profile for the chromosome names, two lines, as
`cnamaste.plot.layout.GenomicAxis.frame` staggers them."""

ROW_NAME_GAP = 0.06
"""Inches between `plot_profile_page`'s column of clone names and the profiles."""

FOOT_PANEL_HEIGHT = 1.2
"""Inches: the least height of `plot_profile_page`'s panel under its square panels; it grows to the left column's
height."""

NEUTRAL_COPY_STATE = "(1,1)"
"""The integer profiles' neutral `(A,B)` state: no CNA."""

NEUTRAL_COPY_COLOUR = "#ffffff"
"""`NEUTRAL_COPY_STATE`'s fill, white, as Numbat's clone profile fills `neu` (`cnv_heatmap`, R/vis.R), where cnamaste's
profile fills (1,1) light grey (`palette.COPY_STATES`)."""

KEY_COLUMN_SPACING = 0.6
"""Font sizes between the entries of a row of `plot_profile_page`'s foot panel's key."""

KEY_PANEL_GAP = 0.04
"""Inches between `plot_profile_page`'s foot panel's tick labels and its key under them."""

CLONE_TICK_LENGTH = 4.0
"""Points: a profile row's tick, left of the axes."""

CLONE_TICK_PAD = 4.0
"""Points between a profile row's tick and its clone name."""

CLONE_COLUMN = re.compile(r"^clone(\d+)$")
"""A state table's column of one clone's state, `clone<k>`."""

FALSE_POSITIVE = "False-positive tumor"
"""A false-positive tumor spot's name in a key: a tumor call on a spot of a normal cell type."""

SPOT_IMAGE_REPEAT = 5
"""Pixels a side a spot's cell is stored at in a spatial panel's image (`draw_spot_image`), so the image's pixels an inch
are this many times the panel's spots an inch, toward the 600 the lab asks of a bitmap of many points
(lab figure rule FX-3)."""

FALSE_POSITIVE_COLOUR = "#2a78d6"
"""A false-positive tumor spot's colour: OKLab dE >= 20 from clone I and II's and light grey under protan, deutan and
tritan simulation (Machado et al. 2009, doi:10.1109/TVCG.2009.113)."""


def get_genome_axis(region_bed: str, chrom_sizes: str) -> GenomeAxis:
    """cnplot's genome axis on `region_bed`'s regions of `chrom_sizes`' chromosomes, in its order, the stretches between
    regions collapsed, X, Y and M left out (cnplot's `DEFAULT_EXCLUDED_CHROMS`).

    Args: `region_bed` the regions, e.g. `data/hg38.regions.bed`, a chromosome's
        arms; `chrom_sizes` a headerless `<chromosome> <length>` file, e.g.
        `data/hg38.chrom.sizes`.
    Returns: the axis.
    """
    return GenomeAxis(region_bed, chrom_sizes)


def get_cnplot_bins(table: pd.DataFrame) -> pd.DataFrame:
    """`table` with its chromosomes named as cnplot names them, `#CHR`, `chr`-prefixed, in place of `CHR`.

    Args: `table` `CHR` and any columns.
    Returns: `table`'s other columns, then `#CHR`.
    """
    return table.drop(columns="CHR").assign(**{"#CHR": "chr" + table["CHR"].astype(str).str.removeprefix("chr")})


def get_cnplot_segments(table: pd.DataFrame, genome_axis: GenomeAxis) -> pd.DataFrame:
    """`table`'s segments on `genome_axis`: `#CHR` named as its chromosomes are, `chr`-prefixed, and each segment cut at
    its regions, one row per piece, the rest of the segment's columns kept.

    A segment through a centromere is so drawn on both arms; cnplot would drop it (`contain=True`) or keep one arm.

    Args: `table` `CHR`, `START`, `END` per segment, and any columns; `genome_axis` cnplot's axis (`get_genome_axis`).
    Returns: `#CHR`, `START`, `END` per piece, then `table`'s other columns.
    """
    regions = genome_axis.regions[["#CHR", "START", "END"]].rename(
        columns={"START": "region_start", "END": "region_end"}
    )
    pieces = get_cnplot_bins(table).merge(regions, on="#CHR")
    pieces = pieces[(pieces["START"] < pieces["region_end"]) & (pieces["END"] > pieces["region_start"])]
    pieces = pieces.assign(
        START=np.maximum(pieces["START"], pieces["region_start"]), END=np.minimum(pieces["END"], pieces["region_end"])
    )
    columns = ["#CHR", "START", "END"]
    return pieces[columns + [c for c in table.columns if c not in ("CHR", "START", "END")]].reset_index(drop=True)


def get_copy_states(segments: pd.DataFrame) -> pd.DataFrame:
    """`segments`' integer copies as one state per clone, `clone<k>` `"(A,B)"`, as cnamaste's profile key writes them
    (`plot_copy_number_profile.plot_state_key`).

    Args: `segments` `get_cnplot_segments` of a `cnv_<ploidy>_seglevel.tsv`, `clone<k> A`, `clone<k> B` per clone.
    Returns: `#CHR`, `START`, `END`, then `clone<k>` per clone, in clone order.
    Raises: `ValueError` where a clone's copy number is missing.
    """
    states = segments[["#CHR", "START", "END"]].astype({"START": int, "END": int})
    for k in sorted(int(c.split()[0].removeprefix("clone")) for c in segments.columns if c.endswith(" A")):
        a, b = segments[f"clone{k} A"], segments[f"clone{k} B"]
        if a.isna().any() or b.isna().any():
            raise ValueError(f"clone {k} has segments with no copy number")
        states[f"clone{k}"] = "(" + a.astype(int).astype(str) + "," + b.astype(int).astype(str) + ")"
    return states


def get_profile_clones(states: pd.DataFrame, integer: bool) -> list[int]:
    """The clones a profile draws: every `clone<k>` but clone 0, the normal clone, which an integer profile draws where
    it holds a state other than `(1,1)`.

    Args: `states` a state table, `clone<k>` per clone, e.g. `get_copy_states`'; `integer` whether its states are
        integer copies.
    Returns: the clones, in clone order.
    """
    clones = sorted(int(m.group(1)) for m in map(CLONE_COLUMN.match, states.columns) if m)
    return [k for k in clones if k != 0 or (integer and not (states["clone0"] == NEUTRAL_COPY_STATE).all())]


def get_copy_state_colours(segments: list[pd.DataFrame], clones: list[list[int]]) -> dict[str, str]:
    """Each `"(A,B)"` the drawn clones hold, in cnamaste's key order (`plot_copy_number_profile.used_states`), in its
    cnamaste profile colour (`plot_copy_number_profile.state_fill`): `palette.COPY_STATES`, a `B > A` state lighter;
    `NEUTRAL_COPY_STATE` in `NEUTRAL_COPY_COLOUR`.

    Args: `segments` the integer profiles, `clone<k> A`, `clone<k> B` per clone; `clones` each one's drawn clones.
    Returns: state -> hex colour.
    """
    pairs = [
        table[[f"clone{k} A", f"clone{k} B"]].set_axis(["clone0 A", "clone0 B"], axis=1)
        for table, drawn in zip(segments, clones, strict=True)
        for k in drawn
    ]
    if not pairs:
        return {}
    states = used_states(pd.concat(pairs, ignore_index=True), [0])
    colours = {f"({a},{b})": mcolors.to_hex(state_fill(a, b)) for a, b in states}
    return {state: NEUTRAL_COPY_COLOUR if state == NEUTRAL_COPY_STATE else c for state, c in colours.items()}


def get_copy_state_groups(colours: dict[str, str]) -> dict[Fraction | None, dict[str, str]]:
    """Group `(A,B)` states by BAF, min(A, B) / (A + B), ascending, so a state and its phase flip share a group, each
    group in `colours`' order; `(0,0)`, with no allele, a group of its own, first.

    Args: `colours` `"(A,B)"` -> colour, e.g. `get_copy_state_colours`'.
    Returns: BAF, None for `(0,0)`, -> its states' colours, by BAF.
    """
    groups: dict[Fraction | None, dict[str, str]] = {}
    for state, colour in colours.items():
        a, b = (int(n) for n in state.strip("()").split(","))
        groups.setdefault(Fraction(min(a, b), a + b) if a + b else None, {})[state] = colour
    return dict(sorted(groups.items(), key=lambda group: (group[0] is not None, group[0] or 0)))


def draw_spot_image(ax: Any, xy: np.ndarray, rgba: np.ndarray) -> None:
    """Draw each spot as one cell of an image on the spots' square lattice (`get_lattice`), in its colour, the cells
    with no spot clear.

    A cell is one pitch a side and cells tile without overlap, so no spot covers a neighbour, and the order the spots
    are given in draws nothing; squares drawn one by one at a panel's scale, a bin less than a device pixel across,
    cover their neighbours and leave the last drawn on top. A cell is stored as `SPOT_IMAGE_REPEAT` pixels a side, as
    `interpolation="none"` embeds the image unsampled in a PDF at its own pixels.

    Args: `ax` the axes; `xy` (n_spots, 2) the spots' places; `rgba` (n_spots, 4) their colours.
    Raises: `ValueError` where two spots fall in one cell, the spots on no one square lattice.
    """
    pitch, turn = get_lattice(xy)
    to_lattice = Affine2D().rotate_deg(-turn).scale(1.0 / pitch)
    uv = to_lattice.transform(xy)
    start = uv.min(axis=0)
    cells = np.round(uv - start).astype(int)
    if len(np.unique(cells, axis=0)) != len(cells):
        raise ValueError("spots share a lattice cell: they are on no one square lattice")
    image = np.zeros((cells[:, 1].max() + 1, cells[:, 0].max() + 1, 4))
    image[cells[:, 1], cells[:, 0]] = rgba
    to_axes = Affine2D().translate(*start).scale(pitch).rotate_deg(turn)
    height, width = image.shape[:2]
    stored = np.repeat(np.repeat(image, SPOT_IMAGE_REPEAT, axis=0), SPOT_IMAGE_REPEAT, axis=1)
    ax.imshow(
        stored, extent=(-0.5, width - 0.5, height - 0.5, -0.5), interpolation="none", transform=to_axes + ax.transData
    )


def draw_labels(ax: Any, physical: pd.DataFrame, labels: pd.Series, colours: dict[Any, str]) -> None:
    """Draw each spot of `labels` at its `physical` place in its label's colour (`draw_spot_image`).

    Args: `ax` the axes; `physical` `x, y` um per barcode; `labels` a label per barcode; `colours` label -> colour.
    """
    rgba = {label: mcolors.to_rgba(colour) for label, colour in colours.items()}
    xy = physical.loc[labels.index, ["x", "y"]].to_numpy(float)
    draw_spot_image(ax, xy, np.array([rgba[label] for label in labels.to_numpy()]))


def dot_centromeres(ax: Any, genome_axis: GenomeAxis) -> None:
    """Dot the lines cnplot dashes where a chromosome's regions part (`draw_segment_boundaries`): on
    `data/hg38.regions.bed`'s arms, its centromere.

    Args: `ax` a profile cnplot drew on `genome_axis`; `genome_axis` its axis.
    """
    at = {gap.start for gap in genome_axis.interior_gaps}
    for line in ax.lines:
        xs = np.unique(line.get_xdata())
        if len(xs) == 1 and float(xs[0]) in at:
            line.set(linestyle=":", linewidth=0.8, alpha=1.0)


def draw_clone_states(
    ax: Any, segments: pd.DataFrame, genome_axis: GenomeAxis, colours: dict[str, str], names: dict[int, str]
) -> None:
    """Draw the categorical states of `names`' clones on `genome_axis`, one row each, the first on top, each state in
    its colour, by cnplot's `plot_heatmap`, each row named by `names`.

    Args: `ax` the axes; `segments` `get_cnplot_segments` of a state table, `clone<k>` the state per clone;
        `genome_axis` cnplot's axis; `colours` state -> colour, every state drawn; `names` each clone drawn -> its
        row's name, top to bottom.
    """
    order = list(colours)
    clones = list(names)
    codes = np.array([[order.index(s) for s in segments[f"clone{k}"]] for k in reversed(clones)])
    rows = [names[k] for k in reversed(clones)]
    plot_heatmap(
        ax,
        codes,
        segments[["#CHR", "START", "END"]],
        genome_axis,
        row_labels=np.array(rows),
        plot_chrname=False,
        cmap=mcolors.ListedColormap([colours[s] for s in order]),
        norm=mcolors.BoundaryNorm(np.arange(len(order) + 1) - 0.5, len(order)),
        show_block_labels=False,
        rasterized=False,
    )
    ax.set_yticks((np.arange(len(clones)) + 0.5) / len(clones), rows)
    ax.tick_params(axis="y", left=True, length=CLONE_TICK_LENGTH, pad=CLONE_TICK_PAD)
    for spine in ax.spines.values():
        spine.set_visible(False)


def get_chromosome_names(genome_axis: GenomeAxis) -> list[str]:
    """Get each chromosome's name under the profiles, its `chr` dropped, e.g. `1`, aligned to `xtick_chrs`.

    Args: `genome_axis` cnplot's axis.
    Returns: the names, in genome order.
    """
    return [str(name).removeprefix("chr") for name in genome_axis.xlab_chrs]


def get_chromosome_name_places(centres: list[float], widths: list[float]) -> list[tuple[int, float]]:
    """Get each chromosome name's place under the profiles, in genome order: the first of `CHR_NAME_ROWS` rows where it
    starts `CHR_NAME_GAP` past the last name there, at its chromosome's middle; where no row has room, the row whose
    last name ends first, the name moved right to clear it.

    Args: `centres` each name's middle, inches from the axes' left edge, in genome order; `widths` each name's width,
        inches.
    Returns: each name's row, 0 the first, and its middle, inches from the axes' left edge.
    """
    ends = [-np.inf] * CHR_NAME_ROWS
    places = []
    for centre, width in zip(centres, widths, strict=True):
        row = next((r for r, end in enumerate(ends) if centre - width / 2.0 >= end + CHR_NAME_GAP), None)
        if row is None:
            row = int(np.argmin(ends))
            centre = ends[row] + CHR_NAME_GAP + width / 2.0
        ends[row] = centre + width / 2.0
        places.append((row, centre))
    return places


class ChromosomeNames(NamedTuple):
    """The chromosome names under a profile (`get_chromosome_name_layout`): each name (`get_chromosome_names`), its
    chromosome's middle and its row and middle (`get_chromosome_name_places`), inches from the axes' left edge, a
    row's height, and the inches the rows take under the axes, `CHR_NAME_PAD` with them."""

    names: list[str]
    centres: list[float]
    places: list[tuple[int, float]]
    line: float
    height: float


def get_chromosome_name_layout(
    figure: Any, renderer: Any, genome_axis: GenomeAxis, width: float, size: float
) -> ChromosomeNames:
    """Get the chromosome names under a profile `width` inches wide, upright at their chromosomes' middles, measured at
    `size`.

    Args: `figure`, `renderer` the page and its renderer; `genome_axis` cnplot's axis; `width` the profile's width,
        inches; `size` the names' font size.
    Returns: the names' layout.
    """
    names = get_chromosome_names(genome_axis)
    extents = [get_text_extent(figure, renderer, name, size) for name in names]
    start, end = genome_axis.ch_coords[0], genome_axis.chr_end
    centres = [(x - start) / (end - start) * width for x in genome_axis.xtick_chrs]
    places = get_chromosome_name_places(centres, [w for w, _ in extents])
    line = max((h for _, h in extents), default=0.0)
    height = CHR_NAME_PAD + line * (max((row for row, _ in places), default=-1) + 1)
    return ChromosomeNames(names, centres, places, line, height)


class ProfileStates(NamedTuple):
    """Each method's profile as a page draws it (`get_profile_states`), per title: its states on the axis' segments,
    their colours and its tumour clones drawn; the keys, name -> its groups of state -> colour; and per title the
    labelled tumour clones its table leaves out, named as its rows name them."""

    states: dict[str, pd.DataFrame]
    colours: dict[str, dict[str, str]]
    tumour: dict[str, list[int]]
    keys: dict[str, list[dict[str, str]]]
    missing: dict[str, list[str]]


def get_profile_states(
    profiles: dict[str, pd.DataFrame | None],
    genome_axis: GenomeAxis,
    state_colours: dict[str, str] | None,
    clones: dict[str, pd.Series | None],
) -> ProfileStates:
    """Get each profile's states on `genome_axis`, segments cut at its regions (`get_cnplot_segments`): a table of
    states in `state_colours`, every such table keyed together; integer copies as `(A,B)` states
    (`get_copy_states`) in cnamaste's profile colours (`get_copy_state_colours`), keyed together, grouped by BAF
    (`get_copy_state_groups`); the first key named `KEY_NAME`, the other unnamed. Each key names the states its
    profiles draw.

    Args: `profiles` title -> `CHR`, `START`, `END` and either `clone<k> A`, `clone<k> B` or `clone<k>` state, or None
        for a placeholder; `genome_axis` cnplot's axis; `state_colours` state -> colour of every state table, in key
        order, or None; `clones` title -> clone per barcode, or None.
    Returns: the profiles' states, colours, tumour clones, keys and labelled clones left out; a placeholder in none.
    """
    segments = {
        title: None if profile is None else get_cnplot_segments(profile, genome_axis)
        for title, profile in profiles.items()
    }
    return get_state_tables(segments, state_colours, clones)


def get_state_tables(
    tables: dict[str, pd.DataFrame | None], state_colours: dict[str, str] | None, clones: dict[str, pd.Series | None]
) -> ProfileStates:
    """Get each profile's states per row of its table, as `get_profile_states` gets them on cnplot's segments: integer
    copies as `(A,B)` states (`get_copy_states`); a table of states; the first key named `KEY_NAME`, the other
    unnamed.

    Args: `tables` title -> `#CHR`, `START`, `END` and either `clone<k> A`, `clone<k> B` or `clone<k>` state per row, or
        None for a placeholder; `state_colours` state -> colour of every state table, in key order, or None; `clones`
        title -> clone per barcode, or None.
    Returns: the profiles' states, colours, tumour clones, keys and labelled clones left out; a placeholder in none.
    """
    segments = {title: table for title, table in tables.items() if table is not None}
    integer = [title for title, s in segments.items() if any(c.endswith(" A") for c in s.columns)]
    states = {title: get_copy_states(s) if title in integer else s for title, s in segments.items()}
    tumour = {title: get_profile_clones(s, title in integer) for title, s in states.items()}
    copy_colours = get_copy_state_colours([segments[t] for t in integer], [tumour[t] for t in integer])
    colours = {title: copy_colours if title in integer else (state_colours or {}) for title in states}
    # NB one key for the state tables, one group; one for the integer profiles, a group per BAF.
    drawn = {
        state for title in states if title not in integer for k in tumour[title] for state in states[title][f"clone{k}"]
    }
    keys: dict[str, list[dict[str, str]]] = {}
    if drawn:
        keys[KEY_NAME] = [{state: colour for state, colour in (state_colours or {}).items() if state in drawn}]
    if copy_colours:
        keys["" if keys else KEY_NAME] = list(get_copy_state_groups(copy_colours).values())
    # NB a labelled tumour clone the profile's table leaves out, named as its rows name it.
    labelled = {title: labels for title, labels in clones.items() if labels is not None}
    missing = {
        title: [
            cast_clone_label(f"clone {k}")
            for k in sorted(set(labelled[title].astype(int)) - {0})
            if k not in tumour[title]
        ]
        for title in tumour
        if title in labelled
    }
    return ProfileStates(states, colours, tumour, keys, missing)


def draw_profile(
    ax: Any,
    states: pd.DataFrame | None,
    genome_axis: GenomeAxis,
    colours: dict[str, str],
    names: dict[int, str],
    chromosomes: ChromosomeNames | None,
) -> None:
    """Draw one method's profile on `ax`: each named clone's states a row (`draw_clone_states`), centromeres dotted
    (`dot_centromeres`); with `chromosomes`, the genome axis' frame and the chromosome names under it, each a row lower
    where crowded. A profile with no clone named is an empty frame, a placeholder.

    Args: `ax` the profile's axes; `states` its states on the axis' segments (`ProfileStates.states`), or None;
        `genome_axis` cnplot's axis; `colours` state -> colour; `names` each clone drawn -> its row's name, top to
        bottom, empty for none; `chromosomes` the names' layout (`get_chromosome_name_layout`), or None for none.
    """
    if not names or states is None:
        ax.set_xticks([])
        ax.set_yticks([])
        return
    draw_clone_states(ax, states, genome_axis, colours, names)
    if chromosomes is not None:
        decorate_genome_axis(ax, genome_axis, plot_chrname=False, hide_spines=True)
        for x, name, centre, (row, placed) in zip(
            genome_axis.xtick_chrs, chromosomes.names, chromosomes.centres, chromosomes.places, strict=True
        ):
            ax.annotate(
                name,
                xy=(x, 0.0),
                xycoords=ax.get_xaxis_transform(),
                xytext=(72.0 * (placed - centre), -72.0 * (CHR_NAME_PAD + row * chromosomes.line)),
                textcoords="offset points",
                ha="center",
                va="top",
                annotation_clip=False,
            )
    dot_centromeres(ax, genome_axis)


def get_key_layout(
    keys: dict[str, list[dict[str, str]]], figure: Any, renderer: Any, size: float, width: float
) -> tuple[float, dict[str, list[list[int]]]]:
    """Get the profiles' keys' one swatch width (`get_key_swatch`), so their swatches line up, and each key's groups
    packed into rows `width` inches wide (`get_key_rows`).

    Args: `keys` name -> its groups, state -> colour each; `figure`, `renderer` the page and its renderer; `size` the
        names' font size; `width` inches of a row.
    Returns: the swatch's width, inches, and each key's rows.
    """
    swatch = max((get_key_swatch(groups, figure, renderer, size) for groups in keys.values()), default=0.0)
    return swatch, {name: get_key_rows(groups, swatch, width) for name, groups in keys.items()}


def draw_keys(
    figure: Any,
    keys: dict[str, list[dict[str, str]]],
    rows: dict[str, list[list[int]]],
    swatch: float,
    left: float,
    top: float,
    width: float,
) -> None:
    """Draw the profiles' keys one under another (`draw_key`), the first `top` inches from the page's top, their
    swatches from `left`, `KEY_ROW` a row.

    Args: `figure` the page; `keys` name -> its groups; `rows` each key's rows, `swatch` the swatch's width
        (`get_key_layout`); `left`, `top`, `width` inches.
    """
    for name, groups in keys.items():
        height = KEY_ROW * len(rows[name])
        draw_key(add_axes(figure, left, top, width, height), groups, rows[name], name, swatch)
        top += height


def get_axis_label_size(size: float) -> float:
    """Get an axis label's font size on a page whose other text is `size`: `AXIS_LABEL_SCALE` times it.

    Args: `size` the page's font size, `figure_style.page_style`'s.
    Returns: the size in points.
    """
    return AXIS_LABEL_SCALE * size


def set_text_sizes(figure: Any, size: float) -> None:
    """Set every text on `figure` to `size` (`figure_style.uniform`), then each axis label, a colorbar's among them, to
    `get_axis_label_size`, so no key, tick label or panel letter is larger than an axis label.

    Args: `figure` the page; `size` the page's font size.
    """
    uniform(figure, size)
    for ax in figure.axes:
        for axis in (ax.xaxis, ax.yaxis):
            axis.label.set_fontsize(get_axis_label_size(size))


def get_text_extent(figure: Any, renderer: Any, text: str, size: float) -> tuple[float, float]:
    """Inches of `text` on `figure` at `size`, across and up.

    Args: `figure`, `renderer` the page and its renderer; `text` the text; `size` its font size.
    Returns: the width and the height.
    """
    drawn = figure.text(0.0, 0.0, text, fontsize=size)
    extent = drawn.get_window_extent(renderer)
    drawn.remove()
    return float(extent.width / figure.dpi), float(extent.height / figure.dpi)


def get_key_swatch(groups: list[dict[str, str]], figure: Any, renderer: Any, size: float) -> float:
    """A key swatch's width in inches: `KEY_LABEL_ROOM` times its widest state name at `size`.

    Args: `groups` the key's groups, state -> colour each; `figure`, `renderer` the page and its renderer; `size` the
        names' font size.
    Returns: the width.
    """
    widths = [get_text_extent(figure, renderer, state, size)[0] for group in groups for state in group]
    return KEY_LABEL_ROOM * max(widths, default=0.0)


def get_key_rows(groups: list[dict[str, str]], swatch: float, room: float) -> list[list[int]]:
    """Pack a key's groups into rows of `room`, in order, each row taking the next group while it fits `KEY_GROUP_GAP`
    swatches after the last; a group wider than `room` is a row of its own.

    Args: `groups` the key's groups; `swatch` a swatch's width in inches; `room` inches of a row.
    Returns: the groups' indices per row.
    """
    rows: list[list[int]] = []
    filled = 0.0
    for i, group in enumerate(groups):
        width = swatch * len(group)
        if rows and filled + KEY_GROUP_GAP * swatch + width <= room:
            rows[-1].append(i)
            filled += KEY_GROUP_GAP * swatch + width
        else:
            rows.append([i])
            filled = width
    return rows


def draw_key(ax: Any, groups: list[dict[str, str]], rows: list[list[int]], name: str, swatch: float) -> None:
    """Draw a copy-number key, laid out as cnplot's (`_plot_cnv_legend`): `name` in the left margin, right-aligned
    with the profiles' clone names (`draw_clone_states`), then each row's groups, a swatch per state, its name under
    it, from the axes' left edge; a row of two or more groups justified to the key's widest row, its groups at least
    `KEY_GROUP_GAP` swatches apart, the gaps between them equal.

    Args: `ax` the key's axes, `KEY_ROW` a row, its left edge the profiles'; `groups` state -> colour per group, in
        order; `rows` the groups' indices per row (`get_key_rows`); `name` the key's name; `swatch` a swatch's width in
        inches (`get_key_swatch`).
    """
    width, height = (
        ax.get_position().width * ax.figure.get_figwidth(),
        ax.get_position().height * ax.figure.get_figheight(),
    )
    ax.set_xlim(0.0, width)
    ax.set_ylim(0.0, height)
    ax.axis("off")
    fills = [swatch * sum(len(groups[i]) for i in row) for row in rows]
    span = max(fill + KEY_GROUP_GAP * swatch * (len(row) - 1) for fill, row in zip(fills, rows, strict=True))
    for r, row in enumerate(rows):
        y = height - r * KEY_ROW - KEY_SWATCH_HEIGHT
        gap = (span - fills[r]) / (len(row) - 1) if len(row) > 1 else 0.0
        x = 0.0
        for i in row:
            for state, colour in groups[i].items():
                ax.add_patch(Rectangle((x, y), swatch, KEY_SWATCH_HEIGHT, facecolor=colour, edgecolor="black"))
                ax.text(x + swatch / 2.0, y - 0.03, state, ha="center", va="top")
                x += swatch
            x += gap
    name_x = -(CLONE_TICK_LENGTH + CLONE_TICK_PAD) / 72.0
    ax.text(name_x, height - KEY_SWATCH_HEIGHT / 2.0, name, ha="right", va="center", clip_on=False)


def get_alternate_tick_labels(values: list[float], labels: list[str]) -> list[str]:
    """Every other tick's label, counting from the tick at 0 where the axis holds one, else from the outermost, the
    unit `coordinates` writes once, `K`, on the outermost tick labelled; `frame_extent`'s round steps run their labels
    together on a panel of `plot_annotation_page`.

    Args: `values` the ticks; `labels` their labels, `coordinates`', or empty for an axis left unlabelled.
    Returns: the labels, every other one empty.
    """
    if not values:
        return labels
    unit = "K" if any(label.endswith("K") for label in labels) else ""
    zeros = np.flatnonzero(np.isclose(values, 0.0))
    start = int(zeros[0]) if len(zeros) else int(np.argmax(np.abs(values)))
    shown = [i for i in range(len(values)) if (i - start) % 2 == 0]
    outer = max(shown, key=lambda i: abs(values[i]))
    return [
        label.removesuffix(unit) + (unit if i == outer else "") if i in shown else "" for i, label in enumerate(labels)
    ]


def get_key_count(labels: pd.Series, label: Any, in_tissue: pd.Index | None) -> str:
    """A spatial key's entry: `label`, then its spots in brackets, or, given `in_tissue`, its in-tissue spots' share of
    them, e.g. `Clone I (12%)`.

    Args: `labels` a label per barcode; `label` one of them; `in_tissue` the in-tissue spots, or None.
    Returns: the entry.
    """
    spots = labels.index[(labels == label).to_numpy()]
    if in_tissue is None:
        return f"{label} ({len(spots):,})"
    return f"{label} ({spots.isin(in_tissue).sum() / len(in_tissue):.0%})"


def get_clone_colours(clones: list[pd.Series], clone_ids: list[int] | None = None) -> dict[int, str]:
    """Get each clone's colour on a page: clone 0 light grey, as in `plot_clones_spatial`, the others `CLONE_PALETTE`
    spread over clone I and every clone of `clones` and `clone_ids`, so pages given the same ids colour a clone alike.

    Args: `clones` a clone per barcode, per method; `clone_ids` clones the palette spreads over beside theirs, or None.
    Returns: clone -> hex colour, in clone order.
    """
    ids = sorted(set().union({0, 1}, *(set(c.astype(int)) for c in clones), clone_ids or ()))
    tumour = [k for k in ids if k != 0]
    palette = sns.color_palette(CLONE_PALETTE, len(tumour)).as_hex()
    return {k: "lightgrey" if k == 0 else palette[tumour.index(k)] for k in ids}


def get_clone_names(labels: pd.Series, clone_colours: dict[int, str]) -> tuple[pd.Series, dict[str, str]]:
    """Each spot's clone by its name (`cast_clone_label`), and the colour of each name `labels` holds, in clone order.

    Args: `labels` a clone per barcode; `clone_colours` clone -> colour.
    Returns: the name per barcode, and name -> colour.
    """
    found = sorted(set(labels.astype(int)))
    names = {k: cast_clone_label(f"clone {k}") for k in found}
    return labels.astype(int).map(names), {names[k]: clone_colours[k] for k in found}


def get_profile_row_names(names: pd.Series | None, clones: list[int], in_tissue: pd.Index | None) -> dict[int, str]:
    """Name each of a profile's `clones` as the page's keys name it (`get_key_count`), e.g. `Clone I (12%)`, or by its
    name alone (`cast_clone_label`) where the method's spots are not labelled.

    Args: `names` a clone name per barcode (`get_clone_names`), or None; `clones` the clones drawn; `in_tissue` the
        in-tissue spots, or None.
    Returns: clone -> row name, in `clones`' order.
    """
    rows = {k: cast_clone_label(f"clone {k}") for k in clones}
    if names is None:
        return rows
    return {k: get_key_count(names, row, in_tissue) for k, row in rows.items()}


def get_label_left(ax: Any, renderer: Any) -> float:
    """Get the inches from the page's left edge to the left of a panel's axis labels, tick labels and ticks, or of its
    frame where it has none: where its letter aligns.

    Args: `ax` the panel, placed; `renderer` the page's renderer.
    Returns: the inches.
    """
    boxes = [ax.get_window_extent(renderer), *(axis.get_tightbbox(renderer) for axis in (ax.xaxis, ax.yaxis))]
    return float(min(box.x0 for box in boxes if box is not None) / ax.figure.dpi)


def draw_inside_key(
    ax: Any, labels: pd.Series, colours: dict[Any, str], in_tissue: pd.Index | None, size: float, ncol: int = 1
) -> None:
    """Key a spatial panel inside its top left, untitled, `ncol` labels a row, read row by row, each with its spots or
    its share of the in-tissue spots (`get_key_count`), at `size`.

    Args: `ax` the panel; `labels` a label per barcode; `colours` label -> colour, in key order; `in_tissue` the
        in-tissue spots, or None; `size` the key's font size; `ncol` the labels a row.
    """
    entries = [(square_handle(c), get_key_count(labels, label, in_tissue)) for label, c in colours.items()]
    # NB matplotlib fills a key's columns first; this order reads it row by row.
    handles, texts = zip(*[entries[i] for c in range(ncol) for i in range(c, len(entries), ncol)], strict=True)
    ax.legend(
        list(handles),
        list(texts),
        ncol=ncol,
        loc="upper left",
        frameon=False,
        handlelength=0.8,
        handletextpad=0.4,
        labelspacing=0.3,
        borderpad=0.0,
        prop={"size": size},
    )


def count_covered_spots(ax: Any, xy: np.ndarray) -> int:
    """Count the spots a spatial panel's key covers: those whose place falls inside the key's extent as drawn.

    Args: `ax` the panel, placed, its key drawn; `xy` (n_spots, 2) the spots' places, the panel's data coordinates.
    Returns: the count.
    """
    extent = ax.get_legend().get_window_extent(ax.figure.canvas.get_renderer())
    (x0, y0), (x1, y1) = ax.transData.inverted().transform([[extent.x0, extent.y0], [extent.x1, extent.y1]])
    return int(np.sum((xy[:, 0] >= x0) & (xy[:, 0] <= x1) & (xy[:, 1] >= y0) & (xy[:, 1] <= y1)))


def plot_annotation_page(
    physical: pd.DataFrame,
    epithelial: pd.Series | None,
    clones: dict[str, pd.Series | None],
    profiles: dict[str, pd.DataFrame | None],
    genome_axis: GenomeAxis,
    state_colours: dict[str, str] | None = None,
    he_frame: pd.DataFrame | None = None,
    clone_ids: list[int] | None = None,
    false_positive: dict[str, pd.Index] | None = None,
    panels_per_row: int = 2,
    in_tissue: pd.Index | None = None,
    keys_inside: bool = False,
    fit_page: bool = False,
) -> Any:
    """Draw a page the text column wide: (a) a grid of spatial panels `panels_per_row` wide, the H&E slide where
    `he_frame` is given, each spot's cell type where `epithelial` is, then each method's clones, each keyed right, or,
    with `keys_inside`, untitled inside its top left, a label a row (`draw_inside_key`), the spots a key covers logged
    (`count_covered_spots`); then each method's clone profile, keyed below. A key counts each label's spots, or, with
    `in_tissue`, their share of the in-tissue spots, the spots of `in_tissue` alone counted. Without `profiles` the
    page is the spatial panels alone. Every panel, spatial or profile, is lettered (a), (b), ... in page order, its
    letter alone, its title the letter's label, named on a page of its own (`plot_panel_titles`); a spatial panel's
    letter is over it, from the left of its axis label, or of its ticks where it has none (`get_label_left`).

    Spots are cells of an image on their lattice, one pitch a side, at their `physical` place (`draw_spot_image`),
    every panel framed on the extent of them all and the slide's origin, so each axis reads from 0, every other tick
    labelled (`get_alternate_tick_labels`). A panel is as wide as the column leaves beside its column of keys,
    at the spots' framed aspect, and the page as tall as its panels need, so may run past `figure_style.TEXT_HEIGHT`;
    with `fit_page`, the panels are narrowed, the grid centred, so the page is at most `page_size("full")` tall, the
    text block less its caption, and prints unscaled.
    A clone is one colour in every panel, clone 0 light grey as in `plot_clones_spatial`; each panel's key names its
    own clones, each with its spots in the panel. A cell type is not epithelial, in clone 0's colour, or epithelial, in
    clone I's, keyed in that order.
    A title of `false_positive` draws its false-positive tumor spots in `FALSE_POSITIVE_COLOUR` over its clones, keyed
    `FALSE_POSITIVE`. A method whose clones are None is a placeholder: its title over an empty grey frame.

    A profile stacks the drawn clones (`get_profile_clones`), `CLONE_ROW` each, each row named as the page's keys name
    its clone (`get_profile_row_names`), with `in_tissue` its share of the in-tissue spots, the profiles' left edge
    right of the margin where their longest name needs it, on `genome_axis` (`draw_profile`), its states and keys as
    `get_profile_states` gives them. The profiles stack `PROFILE_GAP` apart, each lettered level with its first row,
    the letters in one column `PROFILE_LETTER_GAP` left of the clone names, so a letter takes no line of its own; the
    chromosome names are under the last profile drawn (`get_chromosome_name_layout`). Each key is laid out as cnplot's
    copy-number key (`draw_keys`), its swatches from the profiles' left edge, every key's swatch one width, its groups
    packed into rows (`get_key_layout`), each row justified, the first key `GAP` under the chromosome names. Every
    colour is opaque. A tumour clone of `clones` that its profile's table leaves out is named beside the profile's
    title, its letter's label. A method whose profile is None is a placeholder, and a profile with no tumour clone is
    drawn as one: an empty row of `CLONE_ROW`.

    Args: `physical` `x, y` um, y up, per barcode, every spot drawn; `epithelial` whether a spot's cell type is
        epithelial, per barcode, the spots drawn, or None for no cell-type panel; `clones` title -> clone per barcode,
        or None for a placeholder, in panel order; `profiles` title -> `CHR`, `START`, `END` and either `clone<k> A`,
        `clone<k> B` or `clone<k>` state, or None for a placeholder, in row order; `genome_axis` cnplot's axis
        (`get_genome_axis`); `state_colours` state -> colour of every state table, in key order, its key writing each
        state as named; `he_frame` the slide, `x, y` um (`draw_slide`), or None for no H&E panel; `clone_ids` clones
        the palette spreads over beside `clones`' own, so pages given the same ids colour a clone alike, or None;
        `false_positive` a title of `clones` -> its false-positive tumor spots, drawn in place of their clone, or None;
        `panels_per_row` the spatial panels per row, each followed by its key; `in_tissue` the in-tissue spots, or
        None; `keys_inside` whether each spatial panel's key is inside its top left; `fit_page` whether the page is at
        most `page_size("full")` tall.
    Returns: the figure.
    Raises: `ValueError` where `fit_page` leaves the spatial panels no height.
    """
    drawn_profiles = get_profile_states(profiles, genome_axis, state_colours, clones)
    tumour, keys, missing = drawn_profiles.tumour, drawn_profiles.keys, drawn_profiles.missing
    labelled = {title: labels for title, labels in clones.items() if labels is not None}

    clone_colours = get_clone_colours(list(labelled.values()), clone_ids)
    # NB each labelled method's spots by clone name, its false-positive tumor spots as `FALSE_POSITIVE`.
    named: dict[str, tuple[pd.Series, dict[str, str]]] = {}
    for title, labels in labelled.items():
        names, colours = get_clone_names(labels, clone_colours)
        if false_positive is not None and title in false_positive:
            names = names.where(~names.index.isin(false_positive[title]), FALSE_POSITIVE)
            colours = {**colours, FALSE_POSITIVE: FALSE_POSITIVE_COLOUR}
        named[title] = (names, colours)
    row_names = {
        title: get_profile_row_names(named[title][0] if title in named else None, drawn_clones, in_tissue)
        for title, drawn_clones in tumour.items()
    }
    grid: list[tuple[str, Any, dict[Any, str], str | None]] = []
    if he_frame is not None:
        grid.append(("H&E", None, {}, None))
    if epithelial is not None:
        kinds = {True: "Epithelial", False: "Non-epithelial"}
        grid.append(
            (
                "Cell type",
                epithelial.astype(bool).map(kinds),
                {kinds[False]: clone_colours[0], kinds[True]: clone_colours[1]},
                "Cell type",
            )
        )
    for title, labels in clones.items():
        if labels is None:
            grid.append((title, None, {}, None))
            continue
        names, colours = named[title]
        grid.append((title, names, colours, "Clone"))
    x, y = physical["x"], physical["y"]
    box = (min(x.min(), 0.0), max(x.max(), 0.0), min(y.min(), 0.0), max(y.max(), 0.0))
    start_x, end_x, start_y, end_y, _ = framed(box)
    ratio = (end_y - start_y) / (end_x - start_x)
    letters = iter(string.ascii_lowercase)

    with page_style() as size:
        figure = plt.figure(figsize=page_size("full"), dpi=DPI, facecolor="white")
        page_width = figure.get_figwidth()
        width = page_width - MARGIN["left"] - MARGIN["right"]

        panels = [add_axes(figure, MARGIN["left"], MARGIN["top"], 1.0, 1.0) for _ in grid]
        if he_frame is not None:
            draw_slide(panels[0], he_frame, physical[["x", "y"]].to_numpy())
        options = {
            "loc": "upper left",
            "frameon": False,
            "handlelength": 0.8,
            "handletextpad": 0.4,
            "labelspacing": 0.3,
            "borderaxespad": 0.0,
            "alignment": "left",
        }
        renderer = figure.canvas.get_renderer()
        for ax, (title, labels, colours, name) in zip(panels, grid, strict=True):
            if labels is None:
                continue
            draw_labels(ax, physical, labels, colours)
            if keys_inside:
                draw_inside_key(ax, labels, colours, in_tissue, size)
                continue
            entries = [(square_handle(c), get_key_count(labels, label, in_tissue)) for label, c in colours.items()]
            handles, texts = zip(*entries, strict=True)
            ax.legend(list(handles), list(texts), title=name, **options)
        # NB each column of keys right of its panels as wide as its widest; none with the keys inside.
        columns = [
            0.0
            if keys_inside
            else GAP
            + max(
                (
                    ax.get_legend().get_window_extent(renderer).width / figure.dpi
                    for ax in panels[col::panels_per_row]
                    if ax.get_legend() is not None
                ),
                default=0.0,
            )
            for col in range(panels_per_row)
        ]

        # NB under a panel its x ticks, labelled on each column's last panel; below the grid the profiles and their
        #    keys. None depends on the panels' size.
        rows = -(-len(grid) // panels_per_row)
        labelled_x = [k + panels_per_row >= len(grid) for k in range(len(grid))]
        tick_room = [TICK_LINE + AXIS_LABEL_LINE if x else 0.0 for x in labelled_x]
        row_panels = [range(r * panels_per_row, min((r + 1) * panels_per_row, len(grid))) for r in range(rows)]
        under = [max(tick_room[k] for k in ks) for ks in row_panels]
        # NB a placeholder or a profile with no tumour clone is one empty row, so the rows keep their method order.
        heights = {title: CLONE_ROW * max(len(tumour.get(title, [])), 1) for title in profiles}
        last = next((title for title in reversed(profiles) if tumour.get(title)), None)
        # NB the profiles start where the page's margin, or their letters and longest row name with its tick, leave room.
        name_width = max(
            (get_text_extent(figure, renderer, row, size)[0] for rows in row_names.values() for row in rows.values()),
            default=0.0,
        )
        letter_width = max(
            (
                get_text_extent(figure, renderer, f"({string.ascii_lowercase[len(grid) + k]})", size)[0]
                for k in range(len(profiles))
            ),
            default=0.0,
        )
        name_room = name_width + (CLONE_TICK_LENGTH + CLONE_TICK_PAD) / 72.0
        profile_left = max(MARGIN["left"], GAP + letter_width + PROFILE_LETTER_GAP + name_room)
        profile_width = width - (profile_left - MARGIN["left"]) - CHR_OVERHANG
        chromosomes = get_chromosome_name_layout(figure, renderer, genome_axis, profile_width, size)
        chr_line = chromosomes.height
        rows_profiles = (
            TITLE_LINE
            + sum(heights.values())
            + PROFILE_GAP * (len(heights) - 1)
            + (chr_line if last is not None else 0.0)
            + GAP
        )
        swatch, key_rows = get_key_layout(keys, figure, renderer, size, profile_width)
        key_band = KEY_ROW * sum(len(r) for r in key_rows.values())
        height_profiles = rows_profiles + key_band if profiles else 0.0

        # NB a panel as wide as the column leaves beside its column of keys, at the spots' framed aspect; with
        #    `fit_page`, no wider than the page's height leaves it, the grid centred.
        panel_w = (width - sum(columns) - (panels_per_row - 1) * GAP) / panels_per_row
        shift = 0.0
        if fit_page and grid:
            fixed = 2.0 * MARGIN["top"] + rows * (TITLE_LINE + GAP) + sum(under) + height_profiles
            fitted = (page_size("full")[1] - fixed) / (rows * ratio)
            if fitted <= 0.0:
                raise ValueError(
                    f"fit_page: {fixed:.2f} in of the page is not panels, past {page_size('full')[1]:.2f} in"
                )
            shift = max(panel_w - fitted, 0.0) * panels_per_row / 2.0
            panel_w = min(panel_w, fitted)
        panel_h = panel_w * ratio
        row_top = [sum(TITLE_LINE + panel_h + under[i] + GAP for i in range(r)) for r in range(rows)]
        height_a = row_top[-1] + TITLE_LINE + panel_h + under[-1] + GAP if grid else 0.0
        page_height = MARGIN["top"] + height_a + height_profiles + MARGIN["top"]
        figure.set_size_inches(page_width, page_height)
        logger.info(
            f"Annotation page {page_width:.2f} x {page_height:.2f} in, spatial panels {panel_w:.2f} x {panel_h:.2f} in."
        )

        top = MARGIN["top"]
        spatial_letters = []
        for k, (ax, (title, _, _, _)) in enumerate(zip(panels, grid, strict=True)):
            row, col = divmod(k, panels_per_row)
            frame_extent(ax, box, labels=col == 0)
            left = MARGIN["left"] + shift + sum(panel_w + columns[c] + GAP for c in range(col))
            place(ax, left, top + row_top[row] + TITLE_LINE, panel_w, panel_h)
            if ax.get_legend() is not None and not keys_inside:
                ax.get_legend().set_bbox_to_anchor((1.0 + GAP / panel_w, 1.0), transform=ax.transAxes)
            for axis in (ax.xaxis, ax.yaxis):
                found = [t.get_text() for t in axis.get_ticklabels()]
                axis.set_ticklabels(get_alternate_tick_labels(list(axis.get_ticklocs()), found))
            if title in clones and clones[title] is None:
                for spine in ax.spines.values():
                    spine.set_visible(True)
                    spine.set_color("0.6")
            # NB x labelled on each column's last panel, which an odd count leaves a row short on the right.
            if not labelled_x[k]:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel(r"$x$ ($\mu$m)", fontsize=get_axis_label_size(size))
            if col == 0:
                ax.set_ylabel(r"$y$ ($\mu$m)", fontsize=get_axis_label_size(size))
            spatial_letters.append((top + row_top[row], next(letters), title, ax))
        top += height_a

        if profiles:
            top += TITLE_LINE
            # NB each profile's letter in one column left of the clone names, level with its first row.
            letter_x = (profile_left - name_room - PROFILE_LETTER_GAP) / page_width
            for k, title in enumerate(profiles):
                top += PROFILE_GAP if k else 0.0
                note = f": {', '.join(missing[title])} not in its copy-number table" if missing.get(title) else ""
                y = 1.0 - (top + CLONE_ROW / 2.0) / page_height
                figure.text(letter_x, y, f"({next(letters)})", label=f"{title}{note}", va="center", ha="right")
                draw_profile(
                    add_axes(figure, profile_left, top, profile_width, heights[title]),
                    drawn_profiles.states.get(title),
                    genome_axis,
                    drawn_profiles.colours.get(title, {}),
                    row_names.get(title, {}) if tumour.get(title) else {},
                    chromosomes if title == last else None,
                )
                top += heights[title] + (chr_line if title == last else 0.0)
            top += GAP
            # NB the keys under the chromosome names, `GAP` between.
            draw_keys(figure, keys, key_rows, swatch, profile_left, top, profile_width)

        set_text_sizes(figure, size)
        # NB a spatial panel's letter from the left of its labels, measured at the page's sizes.
        for letter_top, name, title, ax in spatial_letters:
            draw_panel_letter(figure, letter_top, name, title, get_label_left(ax, renderer))
        set_text_sizes(figure, size)
        for ax, (title, labels, colours, _) in zip(panels, grid, strict=True):
            if not keys_inside or ax.get_legend() is None:
                continue
            # NB a key over spots is laid a column wider at a time while it fits its panel; the fewest covered kept.
            xy = physical.loc[labels.index, ["x", "y"]].to_numpy(float)
            covered = {1: count_covered_spots(ax, xy)}
            for ncol in range(2, len(colours) + 1):
                if not min(covered.values()):
                    break
                draw_inside_key(ax, labels, colours, in_tissue, size, ncol)
                if ax.get_legend().get_window_extent(renderer).width > ax.get_window_extent(renderer).width:
                    break
                covered[ncol] = count_covered_spots(ax, xy)
            best = min(covered, key=lambda ncol: covered[ncol])
            draw_inside_key(ax, labels, colours, in_tissue, size, best)
            if covered[best]:
                logger.warning(f"{title}: its key covers {covered[best]:_} spots, {best} labels a row.")
            else:
                logger.info(f"{title}: its key covers no spot, {best} labels a row.")
    return figure


def draw_state_profile(
    ax: Any, states: pd.DataFrame, genome: GenomicAxis, colours: dict[str, str], neutral: str, order: list[int]
) -> None:
    """Draw a table of categorical states as cnamaste's `plot_state_profile` draws integer copies: a row per clone of
    `order`, top to bottom, on `genome`'s bins, each run of one state within a chromosome one rectangle in its colour,
    a state other than `neutral` outlined at its ends, an uncalled bin left blank, each row framed.

    Args: `ax` the axes; `states` `clone<k>` the state per bin, NaN where uncalled; `genome` the bins' axis; `colours`
        state -> colour, every state drawn; `neutral` the state left unoutlined; `order` the clones, top to bottom.
    """
    contig = np.repeat(np.arange(len(genome.lengths)), genome.lengths)
    found = list(colours)
    row = 1.0 / len(order)
    for k, clone in enumerate(order):
        height = BAR_SHARE * row
        y0 = 1.0 - (k + 1) * row + (row - height) / 2.0
        values = states[f"clone{clone}"].to_numpy(dtype=object)
        codes = np.array([found.index(v) if v in colours else -1 for v in values])
        for start, end in get_intervals(contig * 1_000 + codes + 1)[0]:
            if codes[start] < 0:
                continue
            fill = colours[values[start]]
            ax.add_patch(Rectangle((start, y0), end - start, height, facecolor=fill, edgecolor="none", linewidth=0))
            if values[start] != neutral:
                ax.vlines([start, end], ymin=y0, ymax=y0 + height, linewidth=OUTLINE, colors="black", zorder=3)
        ax.add_patch(
            Rectangle((0, y0), genome.n_bins, height, facecolor="none", edgecolor="black", linewidth=OUTLINE, zorder=3)
        )
    genome.frame(ax, labels=False)
    ax.set_ylim(0.0, 1.0)
    ax.set_yticks([1.0 - (k + 0.5) * row for k in range(len(order))])
    ax.set_yticklabels([cast_clone_label(f"clone {clone}") for clone in order])
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_copy_profile(ax: Any, table: pd.DataFrame, genome: GenomicAxis, order: list[int]) -> None:
    """Draw integer copies as cnamaste's `plot_state_profile` draws them, a missing copy number 1, each `(A,B)` in its
    `get_copy_state_colours` colour, (1,1) `NEUTRAL_COPY_COLOUR` (`draw_state_profile`).

    Args: `ax` the axes; `table` `CHR` or `#CHR`, `START`, `END` and `clone<k> A`, `clone<k> B` per bin of `genome`;
        `genome` the bins' axis; `order` the clones, top to bottom.
    """
    copies = table.fillna({c: 1 for c in table.columns if c.endswith((" A", " B"))}).rename(columns={"CHR": "#CHR"})
    colours = get_copy_state_colours([copies], [order])
    draw_state_profile(ax, get_copy_states(copies), genome, colours, NEUTRAL_COPY_STATE, order)


def plot_profile_page(
    clones: dict[str, pd.Series | None],
    profiles: dict[str, pd.DataFrame | None],
    methods: dict[str, str],
    neutral: str,
    state_colours: dict[str, str] | None = None,
    clone_ids: list[int] | None = None,
    in_tissue: pd.Index | None = None,
    square_panels: list[Callable[[Any], None]] | None = None,
    square_side: float = 0.0,
    foot_panel: Callable[[Any], list[Any]] | None = None,
) -> Any:
    """Draw a page of the methods' clone profiles on one set of bins, the text column wide, in two columns. Left: (a)
    each method's profile, top to bottom in `profiles`' order, its method's name over it, integer copies and a table of
    states as cnamaste's own `plot_state_profile` draws copies (`draw_copy_profile`, `draw_state_profile`), on the
    bins' axis (`cnamaste.plot.layout.GenomicAxis`), the chromosomes named under the last drawn as it names them; the
    keys under them (`draw_keys`), the states as `get_state_tables` keys them. Right, with `square_panels`: a row of
    square panels `square_side` inches a side, `GAP` apart past each one's labels; under them, with `foot_panel`, a
    panel as wide as the row past its labels, its key's handles keyed `KEY_PANEL_GAP` under its tick labels, in one row
    where they fit its width, else one a row.

    The bins drawn are those of the chromosomes where a drawn clone (`get_profile_clones`) of any profile holds a CNA,
    a called state other than `neutral` or `NEUTRAL_COPY_STATE`. An uncalled bin of integer copies is drawn (1,1), as
    `plot_state_profile` draws it. Each clone's row is named in two lines, its name over its share of the in-tissue
    spots (`get_profile_row_names`), the names left-aligned in one column for every profile. The profiles' rows fill
    the left column to the right column's height, each at least `CLONE_ROW` and a name's height; where the left column
    is taller, the foot panel, `FOOT_PANEL_HEIGHT` tall at least, grows to it. Each panel is lettered over it, from
    the left of its labels (`get_label_left`), its letter alone; a right panel's axes title, cleared, is its letter's
    label, and (a)'s names each profile's title, top to bottom, a tumour clone its table leaves out noted, named on the
    page of panel titles (`plot_panel_titles`).

    Args: `clones` title -> clone per barcode, or None, the spots the profiles' rows are named by; `profiles` title ->
        `CHR`, `START`, `END` per bin, the same bins in each, in genome order, and either `clone<k> A`, `clone<k> B` or
        `clone<k>` state, NaN where uncalled, or None for a placeholder, top to bottom; `methods` title -> the name over
        its profile; `neutral` a state table's neutral state, e.g. `NEU`; `state_colours` state -> colour of every
        state table, in key order; `clone_ids` clones the palette spreads over beside `clones`', or None; `in_tissue`
        the in-tissue spots, or None; `square_panels` each draws a square panel on the axes it is given, in letter
        order, or None for the profiles alone; `square_side` their side, inches; `foot_panel` draws the panel under them
        on the axes it is given and returns its key's handles, or None.
    Returns: the figure.
    """
    found = {title: table for title, table in profiles.items() if table is not None}
    bins = next(iter(found.values()))[["CHR", "START", "END"]].reset_index(drop=True)
    integer = {title for title, table in found.items() if any(c.endswith(" A") for c in table.columns)}
    tables = {
        title: table.fillna({c: 1 for c in table.columns if c.endswith((" A", " B"))}).rename(columns={"CHR": "#CHR"})
        for title, table in found.items()
    }
    drawn = get_state_tables({title: tables.get(title) for title in profiles}, state_colours, clones)
    altered = np.zeros(len(bins), dtype=bool)
    for title, states in drawn.states.items():
        values = states[[f"clone{k}" for k in drawn.tumour[title]]]
        normal = NEUTRAL_COPY_STATE if title in integer else neutral
        altered |= (values.notna() & (values != normal)).any(axis=1).to_numpy()
    chromosomes = bins["CHR"].to_numpy()
    kept = np.isin(chromosomes, chromosomes[altered]) if altered.any() else np.ones(len(bins), dtype=bool)
    genome = GenomicAxis.from_bins(chromosomes[kept])
    logger.info(f"Profiles drawn on {kept.sum():_} of {len(bins):_} bins, chromosomes {', '.join(genome.names)}.")

    labelled = {title: labels for title, labels in clones.items() if labels is not None}
    clone_colours = get_clone_colours(list(labelled.values()), clone_ids)
    row_names = {
        title: {
            k: name.replace(" (", "\n(")
            for k, name in get_profile_row_names(
                get_clone_names(labelled[title], clone_colours)[0] if title in labelled else None, tumour, in_tissue
            ).items()
        }
        for title, tumour in drawn.tumour.items()
    }
    rows = {title: max(len(drawn.tumour.get(title, [])), 1) for title in profiles}
    last = next((title for title in reversed(profiles) if drawn.tumour.get(title)), None)
    notes = [
        f"{title} ({', '.join(drawn.missing[title])} not in its copy-number table)"
        if drawn.missing.get(title)
        else title
        for title in profiles
    ]
    letters = iter(string.ascii_lowercase)

    with page_style() as size:
        figure = plt.figure(figsize=page_size("full"), dpi=DPI, facecolor="white")
        page_width = figure.get_figwidth()
        renderer = figure.canvas.get_renderer()
        # NB the right column drawn first, its text at the page's sizes, so its labels' room is measured.
        squares = [add_axes(figure, 0.0, 0.0, square_side, square_side) for _ in square_panels or []]
        for ax, draw in zip(squares, square_panels or [], strict=True):
            draw(ax)
        foot_ax = add_axes(figure, 0.0, 0.0, 1.0, FOOT_PANEL_HEIGHT) if foot_panel is not None and squares else None
        handles = foot_panel(foot_ax) if foot_panel is not None and foot_ax is not None else []
        right_axes = [*squares, *([] if foot_ax is None else [foot_ax])]
        titles = [ax.get_title() for ax in right_axes]
        for ax in right_axes:
            ax.set_title("")
        set_text_sizes(figure, size)
        frames = [(ax.get_window_extent(renderer), ax.get_tightbbox(renderer)) for ax in right_axes]
        lefts = [(frame.x0 - tight.x0) / figure.dpi for frame, tight in frames]
        belows = [(frame.y0 - tight.y0) / figure.dpi for frame, tight in frames]
        square_below = max(belows[: len(squares)], default=0.0)
        right_width = sum(lefts[: len(squares)]) + len(squares) * square_side + GAP * max(len(squares) - 1, 0)
        foot_width = right_width - (lefts[-1] if foot_ax is not None else 0.0)
        key_ax, key_height = None, 0.0
        if handles:
            key_ax = add_axes(figure, 0.0, 0.0, foot_width, 1.0)
            key_ax.axis("off")
            for ncol in (len(handles), 1):
                key = key_ax.legend(
                    handles=handles,
                    ncol=ncol,
                    loc="upper left",
                    frameon=False,
                    handlelength=0.8,
                    handletextpad=0.4,
                    labelspacing=0.3,
                    columnspacing=KEY_COLUMN_SPACING,
                    borderaxespad=0.0,
                    borderpad=0.0,
                    prop={"size": size},
                )
                if key.get_window_extent(renderer).width / figure.dpi <= foot_width:
                    break
            key_height = key.get_window_extent(renderer).height / figure.dpi
        right_height = 0.0
        if squares:
            right_height = TITLE_LINE + square_side + square_below
        if foot_ax is not None:
            right_height += GAP + TITLE_LINE + FOOT_PANEL_HEIGHT + belows[-1] + KEY_PANEL_GAP + key_height

        # NB the left column: (a) and the clone names from `GAP`, the profiles past the widest name.
        extents = [
            get_text_extent(figure, renderer, row, size) for names in row_names.values() for row in names.values()
        ]
        name_width = max((w for w, _ in extents), default=0.0)
        name_height = max((h for _, h in extents), default=0.0)
        profile_left = GAP + name_width + ROW_NAME_GAP
        right_left = page_width - MARGIN["right"] - right_width
        profile_right = right_left - GAP if squares else page_width - MARGIN["right"]
        profile_width = profile_right - profile_left
        swatch, key_rows = get_key_layout(drawn.keys, figure, renderer, size, profile_width)
        left_fixed = (
            TITLE_LINE
            + METHOD_NAME_LINE * len(profiles)
            + PROFILE_GAP * (len(profiles) - 1)
            + (GENOME_NAME_ROOM if last is not None else 0.0)
            + GAP
            + KEY_ROW * sum(len(r) for r in key_rows.values())
        )
        count = max(sum(rows.values()), 1)
        row_height = max(CLONE_ROW, name_height, (right_height - left_fixed) / count)
        left_height = left_fixed + row_height * count
        foot_height = FOOT_PANEL_HEIGHT + max(left_height - right_height, 0.0)
        page_height = 2.0 * MARGIN["top"] + max(left_height, right_height)
        figure.set_size_inches(page_width, page_height)
        logger.info(
            f"Profile page {page_width:.2f} x {page_height:.2f} in, profiles {profile_width:.2f} in wide, "
            f"{row_height:.2f} in a clone."
        )

        top = MARGIN["top"]
        draw_panel_letter(figure, top, next(letters), f"Clone profiles, top to bottom: {'; '.join(notes)}", GAP)
        top += TITLE_LINE
        for k, title in enumerate(profiles):
            top += (PROFILE_GAP if k else 0.0) + METHOD_NAME_LINE
            height = row_height * rows[title]
            ax = add_axes(figure, profile_left, top, profile_width, height)
            tumour = drawn.tumour.get(title, [])
            if not tumour:
                ax.set_xticks([])
                ax.set_yticks([])
            elif title in integer:
                draw_copy_profile(ax, tables[title].iloc[kept].reset_index(drop=True), genome, tumour)
            else:
                states = drawn.states[title].iloc[kept].reset_index(drop=True)
                draw_state_profile(ax, states, genome, drawn.colours[title], neutral, tumour)
            if tumour:
                ax.set_yticklabels([row_names[title][c] for c in tumour], ha="left")
                ax.tick_params(axis="y", pad=72.0 * (name_width + ROW_NAME_GAP))
            at = offset_copy(ax.transAxes, figure, y=2.0, units="points")
            ax.text(0.0, 1.0, methods[title], transform=at, ha="left", va="bottom")
            if title == last:
                genome.frame(ax, labels=True)
            top += height + (GENOME_NAME_ROOM if title == last else 0.0)
        draw_keys(figure, drawn.keys, key_rows, swatch, profile_left, top + GAP, profile_width)

        top, left = MARGIN["top"], right_left
        for ax, room, title in zip(squares, lefts, titles, strict=False):
            place(ax, left + room, top + TITLE_LINE, square_side, square_side)
            draw_panel_letter(figure, top, next(letters), title, get_label_left(ax, renderer))
            left += room + square_side + GAP
        if foot_ax is not None:
            top += TITLE_LINE + square_side + square_below + GAP
            place(foot_ax, right_left + lefts[-1], top + TITLE_LINE, foot_width, foot_height)
            draw_panel_letter(figure, top, next(letters), titles[-1], right_left)
            if key_ax is not None:
                key_top = top + TITLE_LINE + foot_height + belows[-1] + KEY_PANEL_GAP
                place(key_ax, right_left + lefts[-1], key_top, foot_width, key_height)
        set_text_sizes(figure, size)
    return figure


def draw_panel_letter(figure: Any, top: float, name: str, title: str, x0: float) -> None:
    """Place a panel's letter `(name)` `x0` inches from the page's left edge and `top` inches from its top, its title
    the letter's label, named on the page of panel titles (`plot_panel_titles`).

    Args: `figure` the page; `top`, `x0` in inches; `name` the letter; `title` the panel's title.
    """
    width, height = figure.get_size_inches()
    figure.text(x0 / width, 1.0 - top / height, f"({name})", label=title, va="top", ha="left")


def get_panel_titles(figure: Any) -> list[str]:
    """Get a page's panels as named on the page of panel titles: each letter `(x)` drawn with its panel's title as its
    label, `(x) <title>`, a title's lines joined by commas, in letter order.

    Args: `figure` the page.
    Returns: a line per titled panel.
    """
    letters = [t for t in figure.findobj(Text) if PANEL_LETTER.fullmatch(t.get_text()) and t.get_label()]
    return [f"{t.get_text()} {', '.join(t.get_label().splitlines())}" for t in sorted(letters, key=Text.get_text)]


def get_wrapped_lines(figure: Any, renderer: Any, text: str, size: float, width: float) -> list[str]:
    """Split `text` into lines no wider than `width` inches at `size`, at its spaces between words (`TITLE_WORD`); a
    word wider than `width` is a line of its own.

    Args: `figure`, `renderer` the page and its renderer; `text` the text; `size` its font size; `width` inches.
    Returns: the lines, in order.
    """
    words = TITLE_WORD.findall(text)
    lines = words[:1] or [""]
    for word in words[1:]:
        joined = f"{lines[-1]} {word}"
        if get_text_extent(figure, renderer, joined, size)[0] <= width:
            lines[-1] = joined
        else:
            lines.append(word)
    return lines


def plot_panel_titles(figures: list[Any]) -> Any:
    """Draw the page naming the panels of `figures`, each panel drawn with its letter alone: per page with a titled
    panel a heading, `Page <k>`, then the page's own title where it has one, then per panel its letter and title
    (`get_panel_titles`), wrapped to the text column (`get_wrapped_lines`), a line after the first from past its letter.

    Args: `figures` the pages before it, in order.
    Returns: the figure.
    """
    pages = {}
    for k, page in enumerate(figures, 1):
        titles = get_panel_titles(page)
        if titles:
            pages[f"Page {k}: {page.get_suptitle()}" if page.get_suptitle() else f"Page {k}"] = titles
    lines = [(line, line in pages) for heading, titles in pages.items() for line in ("", heading, *titles)][1:]
    with page_style() as size:
        figure = plt.figure(figsize=page_size("full"), dpi=DPI, facecolor="white")
        page_width = figure.get_figwidth()
        renderer = figure.canvas.get_renderer()
        width = page_width - MARGIN["left"] - MARGIN["right"]
        rows = []
        for line, heading in lines:
            letter = None if heading else PANEL_LETTER.match(line)
            indent = get_text_extent(figure, renderer, f"{letter.group(0)} ", size)[0] if letter else 0.0
            parts = get_wrapped_lines(figure, renderer, line, size, width - indent)
            rows += [(part, heading, indent if j else 0.0) for j, part in enumerate(parts)]
        height = 2.0 * MARGIN["top"] + TITLE_LIST_LINE * len(rows)
        figure.set_size_inches(page_width, height)
        for i, (row, heading, indent) in enumerate(rows):
            y = 1.0 - (MARGIN["top"] + TITLE_LIST_LINE * (i + 0.5)) / height
            figure.text((MARGIN["left"] + indent) / page_width, y, row, va="center", weight="bold" if heading else None)
        uniform(figure, size)
    return figure
