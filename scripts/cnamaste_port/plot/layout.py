"""Ported from cna_maste `python/cnamaste/plot/layout.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.

Pages for run_cnamaste's figures: a fixed page, panels placed in inches, one genomic axis.

Pages a third, half or full text block tall.
`plot_combined_page`: (a) H&E slide and clones, (b) copy-number profile, (c) each clone's RDR and BAF tracks.
`plot_he_page`: (a) H&E slide and (b) each spot's H&E class.
`plot_truth_page`: a simulated sample's (a) planted clones and (b) planted copy-number profile on a run's bins;
`plot_truth_spatial_page` the clones alone. Genomic panels share one `GenomicAxis`, spatial panels
one extent; every text is `[tool.cnamaste.figures]`'s `font_size`, clones named as the run names them.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D, Transform, offset_copy
from scipy.spatial import cKDTree
from scripts.cnamaste_port.figure_style import stagger, style

MARGIN = {"left": 0.62, "right": 0.05, "top": 0.05, "bottom": 0.35}
"""Inches around the panels: tick labels left, staggered chromosome names below."""


LETTER = 0.16
"""Inches of a panel letter's line."""


TICK_LINE = 0.25
"""Inches under the spatial panels for their tick labels."""


GAP = 0.1
"""Inches between panels."""


PROFILE_ROW = 0.165
"""Inches of profile per clone. The page's height is fixed, so what the profile does not take goes to the clone
tracks below it."""


NAME_LINE = 0.14
"""Inches above a clone's tracks for its name and state key."""


TRACK_GAP = 0.05
"""Inches between a clone's RDR and BAF tracks, so their end ticks do not meet."""


@dataclass(frozen=True)
class GenomicAxis:
    """Bins in genome order, one run per chromosome: the x axis every genomic panel draws on.

    Args: `lengths` bins per chromosome; `names` the chromosomes' names.
    """

    lengths: np.ndarray
    names: tuple[str, ...]

    @classmethod
    def from_bins(cls, chromosomes: Any) -> GenomicAxis:
        """Build the axis from each bin's chromosome, in bin order (`cnv_<ploidy>_seglevel.tsv`'s `CHR`).

        Args: `chromosomes` per-bin chromosome names.
        Returns: a `GenomicAxis`.
        """
        chromosomes = np.asarray(chromosomes)
        starts = np.flatnonzero(np.r_[True, chromosomes[1:] != chromosomes[:-1]])
        lengths = np.diff(np.r_[starts, chromosomes.size])
        names = tuple(str(c) for c in chromosomes[starts])
        assert len(set(names)) == len(names), "a chromosome's bins are not contiguous"
        return cls(lengths, names)

    @property
    def n_bins(self) -> int:
        return int(np.sum(self.lengths))

    @property
    def starts(self) -> np.ndarray:
        return np.r_[0, np.cumsum(self.lengths)[:-1]]

    def frame(self, ax: Any, labels: bool = False) -> None:
        """Set `ax` on bins `[0, n_bins)` with a line at each chromosome's start.

        With `labels`, names each chromosome under its middle, a clashing name on a second line, and `chr` at the left.

        Args: `ax` the axes; `labels` whether to name chromosomes.
        """
        ax.set_xlim(0, self.n_bins)
        ax.set_xticks([])
        for x in self.starts[1:]:
            ax.axvline(x, color="black", linewidth=0.5, zorder=3)
        if not labels:
            return
        size = style()["font_size"]
        points = ax.get_position().width * ax.figure.get_size_inches()[0] * 72.0 / self.n_bins
        middles = [(start + length / 2.0) * points for start, length in zip(self.starts, self.lengths)]
        for start, length, name, line in zip(self.starts, self.lengths, self.names, stagger(middles, self.names, size)):
            at = offset_copy(ax.get_xaxis_transform(), ax.figure, y=-2.0 - line * 1.1 * size, units="points")
            ax.text(start + length / 2.0, 0.0, name, transform=at, ha="center", va="top")
        at = offset_copy(ax.transAxes, ax.figure, x=-4.0, y=-2.0, units="points")
        ax.text(0.0, 0.0, "chr", transform=at, ha="right", va="top")


def add_axes(figure: Any, x0: float, top: float, width: float, height: float) -> Any:
    """Add an axes `width` by `height` inches, its left edge `x0` and its top `top` inches from the page's left and top.

    Args: `figure` the page; `x0`, `top`, `width`, `height` in inches.
    Returns: the axes.
    """
    w, h = figure.get_size_inches()
    return figure.add_axes((x0 / w, 1.0 - (top + height) / h, width / w, height / h))


def place(ax: Any, x0: float, top: float, width: float, height: float) -> None:
    """Move `ax` to `width` by `height` inches at `x0`, `top` inches down.

    Args: `ax` the axes; `x0`, `top`, `width`, `height` in inches.
    """
    w, h = ax.figure.get_size_inches()
    ax.set_position((x0 / w, 1.0 - (top + height) / h, width / w, height / h))


def slide_grid(he_frame: pd.DataFrame, physical_coords: np.ndarray) -> tuple[Any, np.ndarray, np.ndarray, np.ndarray]:
    """The slide's pixels against the spots: where they are drawn, the spot nearest each, and which are on the tissue.

    A pixel is on the tissue where a spot lies within one spot spacing (the spots' median nearest-neighbour distance)
    of it; the slide's own frame places each pixel at its `x, y`, in the spots' physical coordinates
    (`he.slide_image_and_meta`, taken to um by `spatial.physical_affine`).

    Args: `he_frame` the per-pixel frame, `x, y` um; `physical_coords` (n_spots, 2) um, y up.
    Returns: the transform taking pixel `(column, row)` to the spots' `(x, y)` as drawn; each pixel's `(row, column)`;
        its nearest spot; whether it is on the tissue.
    """
    rows, cols = he_frame["array_row"].to_numpy(), he_frame["array_col"].to_numpy()
    xy = he_frame[["x", "y"]].to_numpy(np.float64)
    # NB exact on any three pixels not in a line; a sample of them keeps the fit cheap on a large image.
    every = max(1, len(rows) // 10_000)
    design = np.column_stack([cols[::every], rows[::every], np.ones(len(rows[::every]))]).astype(np.float64)
    fit = np.linalg.lstsq(design, xy[::every], rcond=None)[0].T
    to_axes = Affine2D(np.vstack([fit[0], fit[1], [0.0, 0.0, 1.0]]))
    tree = cKDTree(physical_coords)
    spacing = float(np.median(tree.query(physical_coords, k=2)[0][:, 1]))
    distance, nearest = tree.query(xy)
    return to_axes, np.column_stack([rows, cols]), nearest, distance <= spacing


def draw_image(ax: Any, to_axes: Any, at: np.ndarray, colours: np.ndarray, shape: tuple[int, int]) -> None:
    """Draw `colours` at pixels `at` on a white image of `shape`, through `to_axes` onto the spots' coordinates."""
    image = np.ones((*shape, 3))
    image[at[:, 0], at[:, 1]] = colours
    height, width = shape
    ax.imshow(
        image, extent=(-0.5, width - 0.5, height - 0.5, -0.5), interpolation="none", transform=to_axes + ax.transData
    )


def draw_slide(ax: Any, he_frame: pd.DataFrame, physical_coords: np.ndarray) -> None:
    """Draw the H&E slide `he.slide_image_and_meta` returns under the spots, drawn at their physical `(x, y)`, white
    off the tissue.

    Args: `ax` the axes; `he_frame` the per-pixel frame, `x, y` um; `physical_coords` (n_spots, 2) um, y up.
    """
    to_axes, at, _, on = slide_grid(he_frame, physical_coords)
    rgb = he_frame[["red", "green", "blue"]].to_numpy(np.float64)
    rgb = np.clip(rgb / (255.0 if rgb.max() > 1.0 else 1.0), 0.0, 1.0)
    shape = (int(at[:, 0].max()) + 1, int(at[:, 1].max()) + 1)
    draw_image(ax, to_axes, at[on], rgb[on], shape)


def square_handle(colour: Any) -> Line2D:
    return Line2D([0], [0], marker="s", color="w", markerfacecolor=colour, markeredgewidth=0, markersize=6)


def coordinates(values: Any, thousands: bool | None = None) -> list[str]:
    """Format spatial ticks: in thousands where `thousands` or any reaches 1,000, else one decimal.

    In thousands the unit is written once, on the outermost tick (`0, 10, 20, 30K`); the inner ticks are bare numbers
    on the same scale. Adds a decimal to every tick while two would read alike.

    Args: `values` the tick values; `thousands` to force thousands.
    Returns: the tick labels.
    """
    values = [float(v) for v in values]
    if not (max(abs(v) for v in values) >= 1_000 if thousands is None else thousands):
        # NB whole numbers plain: `frame_extent` puts every tick on a round step, so a decimal says nothing.
        return [f"{v:.0f}" if all(float(x).is_integer() for x in values) else f"{v:.1f}" for v in values]
    outer = int(np.argmax([abs(v) for v in values]))
    for digits in range(4):
        out = ["0" if v == 0 else f"{v / 1_000:.{digits if abs(v) >= 1_000 else max(digits, 1)}f}" for v in values]
        if len(set(out)) == len(set(values)):
            break
    out[outer] += "K" if values[outer] != 0 else ""
    return out


def rounded_end(reach: float) -> tuple[float, float]:
    """The tick step and the end of a spatial axis reaching `reach`: a round number at or above it, from zero.

    1K up to 10K, 5K up to 50K, 10K beyond, so a panel carries three or four ticks whatever the slide's units:
    1,000 ends at 1K, 26K at 30K, 48K at 50K.

    Below 1,000 the coordinates are a lattice's own rather than a slide's pixels -- a simulation's grid runs to tens
    or hundreds -- and the step is the smallest of 1, 2 and 5 times a power of ten that ends the axis within four
    steps: 59 ends at 60 by 20, 140 at 150 by 50, 500 at 600 by 200. A 1K step there draws the spots into a corner.

    Args: `reach` the larger of the panel's two extents, in the coordinates' own units.
    Returns: the step and the end.
    """
    if reach < 1_000.0:
        reach = max(reach, 1e-9)
        power = 10.0 ** np.floor(np.log10(reach / 4.0))
        step = next(float(m * power) for m in (1.0, 2.0, 5.0, 10.0) if np.ceil(reach / (m * power)) <= 4)
    elif reach <= 10_000.0:
        step = 1_000.0
    else:
        step = 5_000.0 if reach <= 50_000.0 else 10_000.0
    return step, float(np.ceil(reach / step) * step)


def framed(box: tuple[float, float, float, float]) -> tuple[float, float, float, float, float]:
    """The frame `frame_extent` draws on `box`: each axis from the round step at or below its least coordinate to the
    one at or above its greatest, one step for both, from the smaller span.

    Args: `box` the spots' physical `(x0, x1, y0, y1)`.
    Returns: `(start_x, end_x, start_y, end_y, step)`.
    """
    x0, x1, y0, y1 = box
    # NB one step for both axes, from the smaller span: the panel is drawn at equal aspect, so two steps would put two
    #    different grids on one square, and a step from the larger span ends the shorter axis far past its spots.
    step, _ = rounded_end(min(x1 - x0, y1 - y0))
    ends = (np.floor(x0 / step), np.ceil(x1 / step), np.floor(y0 / step), np.ceil(y1 / step))
    return (*(float(e * step) for e in ends), step)


def frame_extent(ax: Any, box: tuple[float, float, float, float], labels: bool, pad: float = 0.02) -> None:
    """Frame `ax` on `box` padded by `pad` per side, aspect equal, ticks on round steps from each axis's start to end.

    `box` is the spots' physical `(x0, x1, y0, y1)` as drawn, y up (`platforms.visium.physical_coords`); each axis runs from the
    round step at or below its least coordinate to the one at or above its greatest, and its ticks read the
    coordinates themselves.

    Args: `ax` the axes; `box` `(x0, x1, y0, y1)`; `labels` whether to label the y ticks, which panels side by side
        share (every panel labels its x ticks); `pad` the fraction of each side.
    """
    start_x, end_x, start_y, end_y, step = framed(box)
    dx, dy = pad * (end_x - start_x), pad * (end_y - start_y)
    ax.set_xlim(start_x - dx, end_x + dx)
    ax.set_ylim(start_y - dy, end_y + dy)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("on")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_bounds(start_x, end_x)
    ax.spines["left"].set_bounds(start_y, end_y)
    xs = np.arange(start_x, end_x + step / 2.0, step)
    ys = np.arange(start_y, end_y + step / 2.0, step)
    ax.set_xticks(xs)
    ax.set_yticks(ys)
    thousands = max(abs(start_x), abs(end_x), abs(start_y), abs(end_y)) >= 1_000
    ax.set_xticklabels(coordinates(xs, thousands))
    ax.set_yticklabels(coordinates(ys, thousands) if labels else [])


def track_furniture(ax: Any, name: str) -> None:
    """Keep only the end y ticks, and name the track `name` beside them.

    The bottom label's digits rest on its tick (baseline-aligned: digits have no descender) and the top label's
    hang from its tick (baseline lowered by the digits' height, read from the glyphs), so both stay within the
    track's height and one track's 0.0 does not meet the top label of the track below.

    Args: `ax` the track; `name` its label.
    """
    low, high = ax.get_yticks()[[0, -1]]
    ax.set_yticks([low, high])
    bottom, top = ax.set_yticklabels([f"{low:.1f}", f"{high:.1f}"])
    for label in (bottom, top):
        label.set_verticalalignment("baseline")
    top.set_transform(top.get_transform() + _Hang(top))
    ax.set_ylabel(name)


class _Hang(Transform):
    """Lower `label` by its digits' height, so they hang from its anchor; read at every draw, as the page's fonts are
    set after layout (`uniform`) and again on the way out (`restyle`).

    Args: `label` the text.
    """

    input_dims = output_dims = 2
    is_affine = False

    def __init__(self, label: Any) -> None:
        super().__init__()
        self.label = label

    def transform_non_affine(self, values: Any) -> Any:
        height = TextPath((0, 0), self.label.get_text(), prop=self.label.get_fontproperties()).get_extents().y1
        return np.asarray(values, dtype=float) - [0.0, height * self.label.figure.dpi / 72.0]
