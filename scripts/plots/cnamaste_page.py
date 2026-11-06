"""cnamaste's own page, `cnamaste.plot.layout.plot_combined_page` broken into panels: the H&E slide and its clones
left, its copy-number profile and each clone's pseudobulk RDR and BAF per bin right, as `plot_clones_genomic` draws
them, on the run's integer clones.

- `get_cnamaste_input_mismatch`, `get_integer_fit`, `read_clone_tracks`, `read_normal_fractions` -- a run's inputs:
  its fit on its integer clones, its per-spot bins, and its decoder's normal fraction per clone;
- `CloneTracks`, `plot_cnamaste_page`, `draw_spatial_frame`, `draw_clone_tracks`, `draw_name_key`,
  `recolour_neutral_bins`, `get_track_label_room` -- the page.
"""

from __future__ import annotations

import string
from pathlib import Path
from typing import Any, NamedTuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PathCollection
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure
from scripts.plots.layout import (
    AXIS_LABEL_LINE,
    count_covered_spots,
    draw_copy_profile,
    draw_inside_key,
    draw_labels,
    draw_panel_letter,
    get_alternate_tick_labels,
    get_axis_label_size,
    get_clone_colours,
    get_clone_names,
    get_label_left,
    get_text_extent,
    set_text_sizes,
)

from cnamaste.config import start_time
from cnamaste.figure_style import DPI, page_size, page_style
from cnamaste.logger import get_logger
from cnamaste.palette import TRACK_NORMAL
from cnamaste.plot.layout import (
    GAP,
    LETTER,
    MARGIN,
    NAME_LINE,
    PROFILE_ROW,
    TICK_LINE,
    TRACK_GAP,
    GenomicAxis,
    add_axes,
    coordinates,
    draw_slide,
    frame_extent,
    framed,
    place,
    track_furniture,
)
from cnamaste.plot.plot_genomic import plot_clones_genomic
from cnamaste.spatio_genomic_counts import SpatioGenomicCounts, from_npz
from cnamaste.utils import cast_clone_label

logger = get_logger("cnamaste_page", start_time=start_time)


NORMAL_CLONE_FRACTION = 1.0
"""The normal clone's normal fraction: the decoder holds clone 0 at (1, 1) and f = 1
(`hmm_ab_decoder._hmm_ab_ploidy_decoder`'s `normal_clone`)."""


FIT_KEYS = ("new_assignment", "pred_cnv", "new_log_mu", "new_p_binom", "new_log_mu_shift")
"""The final fit's entries `plot_clones_genomic` reads, `new_log_mu_shift` where the fit has one."""


NORMAL_FRACTION_LOG = "normal fractions per clone: "
"""run_cnamaste's log of its decoder's normal fraction per final clone, a list in clone order, its final clones 0 to
N-1 (`scripts/run_cnamaste.py`: `np.unique(res_combine["new_assignment"])`, then
`logger.info(f"{decoder} normal fractions per clone: {hmm_ab_ploidy_decoded.normal_fractions.tolist()}.")`)."""


SPATIAL_PANEL_HEIGHT = 2.0
"""Inches: a spatial panel's height on `plot_cnamaste_page`, its left column; its width follows the frame's aspect."""


TRACK_HEIGHT = 0.4
"""Inches: a track's least height on `plot_cnamaste_page`; past it the page grows rather than the tracks shrink."""


CHR_NAME_ROOM = 0.3
"""Inches under `plot_cnamaste_page`'s profile for its chromosome names, two lines
(`cnamaste.plot.layout.GenomicAxis.frame`)."""


TRACK_KEY_MARKER = 3.0
"""Points: a dot's size in a clone's key over its tracks, as `cnamaste.plot.layout.name_line` sizes them."""


NEUTRAL_DOT = "grey"
"""A (1,1) bin's dot on `plot_cnamaste_page`'s tracks: neutral grey, the profile's (1,1) white, dark enough to show on
white, where cnamaste's tracks draw steel blue (`palette.TRACK_NORMAL`) beside the LOH blues."""


TRACK_LABEL_PAD = 1.0
"""Points between a track's channel name, `RDR` or `BAF`, and its tick labels (`draw_clone_tracks`)."""


TRACK_RDR_YLIM = 4.0
"""An RDR track's top, ticked 0 to it by 1 (`plot_genomic.plot_clones_genomic`'s `rdr_ylim`); a bin above it is
not drawn."""


class CloneTracks(NamedTuple):
    """A cnamaste run's per-spot bins and its final fit, as `cnamaste.plot.plot_genomic.plot_clones_genomic` takes them.

    Args:
        lengths: bins per chromosome, `binned_data.npz`'s.
        single_X: (bins, 2, spots) read depth and B count, `binned_data.npz`'s.
        single_base_nb_mean: (bins, spots) normal baseline, `binned_data.npz`'s.
        single_total_bb_RD: (bins, spots) SNP-covering depth, `binned_data.npz`'s.
        fit: `rdrbaf_final_nstates<n>_smp.npz`'s `pred_cnv`, `new_log_mu`, `new_p_binom`, and `new_assignment` and
            `new_log_mu_shift` on the integer clones (`get_integer_fit`).
    """

    lengths: np.ndarray
    single_X: Any
    single_base_nb_mean: Any
    single_total_bb_RD: Any
    fit: dict[str, np.ndarray]


def get_cnamaste_input_mismatch(labels: Path, seglevel: Path) -> str:
    """Get why a cnamaste run's inputs disagree, or "": the clones are its integer clones, `clone_labels_int.tsv[.gz]`,
    and its copies, the fit (`rdrbaf_final_nstates<n>_smp.npz`) and its bins (`binned_data.npz`) are read from the
    same run directory as they are.

    Args:
        labels: `--cnamaste_labels`.
        seglevel: `--cnamaste_seglevel`.

    Returns:
        The mismatch, or "" where there is none.
    """
    if not labels.name.startswith("clone_labels_int.tsv"):
        return f"--cnamaste_labels {labels}: a run's integer clones, clone_labels_int.tsv[.gz], wanted"
    if seglevel.resolve().parent != labels.resolve().parent:
        return f"--cnamaste_seglevel {seglevel} is not in the run directory of --cnamaste_labels, {labels.parent}"
    return ""


def get_integer_fit(fit: dict[str, np.ndarray], barcodes: np.ndarray, labels: pd.Series) -> dict[str, np.ndarray]:
    """Get a cnamaste fit on its integer clones, as `cnamaste.plot.layout.plot_combined_page` draws it given
    `int_groups`: each spot's `new_assignment` its integer clone, `clone_labels_int.tsv`'s, and `new_log_mu_shift`
    each integer clone's first clone's, in the integer clones' order. run_cnamaste merges clones whose `(A, B)` agree
    into integer clones, each named by a member clone, "whose copies it calls" (`qa.metrics.merged_at_barcode`).

    Args:
        fit: the final fit, `new_assignment` the fitted clone per spot in `barcodes`' order.
        barcodes: the fit's spots.
        labels: the integer clone per barcode.

    Returns:
        The fit, its `new_assignment` and `new_log_mu_shift` the integer clones'.

    Raises:
        ValueError: where a spot has no integer clone, or the integer clones are not a merge of the fitted clones, each
            named by one of its own.
    """
    integer = labels.reindex(barcodes)
    if integer.isna().any():
        raise ValueError(f"{int(integer.isna().sum()):_} of the fit's spots have no integer clone")
    integer = integer.to_numpy().astype(np.int64)
    pairs = pd.DataFrame({"fitted": np.asarray(fit["new_assignment"]), "integer": integer}).drop_duplicates()
    if pairs["fitted"].duplicated().any() or not all(
        k in set(pairs.loc[pairs["integer"] == k, "fitted"]) for k in pairs["integer"]
    ):
        raise ValueError(f"the integer clones are not a merge of the fitted clones: {pairs.to_dict('records')}")
    ids = np.unique(integer)
    logger.info(f"Fitted clones by integer clone: {pairs.groupby('integer')['fitted'].apply(sorted).to_dict()}.")
    integer_fit = fit | {"new_assignment": integer}
    if "new_log_mu_shift" in fit:
        integer_fit["new_log_mu_shift"] = np.asarray(fit["new_log_mu_shift"]).reshape(-1)[ids]
    return integer_fit


def get_normal_fractions(run_dir: Path, clones: list[int], given: dict[int, float]) -> dict[int, float]:
    """Get the decoder's normal fraction of each of a cnamaste run's integer clones, clone 0 held at
    `NORMAL_CLONE_FRACTION`, as the decoder holds it: `given`, else its log's, `cnamaste_<id>.log`
    (`NORMAL_FRACTION_LOG`), an integer clone's that of the fitted clone naming it (`get_integer_fit`).

    Args:
        run_dir: the run's directory.
        clones: its integer clones, each the id of the fitted clone naming it.
        given: clone -> its normal fraction, `--cnamaste_normal_fraction`'s, or empty to read the log.

    Returns:
        Each clone's normal fraction, clone 0 first.

    Raises:
        ValueError: where `given` is empty and the directory holds other than one `cnamaste_*.log`, the log other than
            one line of the fractions, or the line no fraction of a clone.
    """
    if given:
        return {0: NORMAL_CLONE_FRACTION, **given}
    logs = sorted(run_dir.glob("cnamaste_*.log"))
    if len(logs) != 1:
        raise ValueError(f"{run_dir} holds {len(logs)} cnamaste_*.log, one wanted")
    lines = [line for line in logs[0].read_text().splitlines() if NORMAL_FRACTION_LOG in line]
    if len(lines) != 1:
        raise ValueError(f"{logs[0]} logs {len(lines)} lines of {NORMAL_FRACTION_LOG!r}, one wanted")
    listed = lines[0].split(NORMAL_FRACTION_LOG, 1)[1].split("]", 1)[0].lstrip("[")
    logged = [float(value) for value in listed.split(",")]
    missing = [clone for clone in clones if clone >= len(logged)]
    if missing:
        raise ValueError(f"{logs[0]} logs {len(logged)} clones' fractions, none of clone(s) {missing}")
    fractions = {0: NORMAL_CLONE_FRACTION, **{clone: logged[clone] for clone in sorted(clones) if clone != 0}}
    logger.info(f"{logs[0]}: normal fraction per integer clone {fractions}.")
    return fractions


def read_clone_tracks(run_dir: Path, labels: pd.Series) -> CloneTracks | None:
    """Read a cnamaste run's per-spot bins and final fit, as `plot_cnamaste_page` draws them: its
    `binned_data.npz`, written where the run's config sets `run.write_binned_data` (`utils.write_binned_data`), and
    its `rdrbaf_final_nstates<n>_smp.npz`'s `FIT_KEYS`, on its integer clones (`get_integer_fit`).

    Args:
        run_dir: the run's directory.
        labels: the run's integer clone per barcode, its `clone_labels_int.tsv`.

    Returns:
        The bins and the fit, or None where the run wrote no `binned_data.npz`.

    Raises:
        ValueError: where the directory holds other than one final fit, or the fit labels other than the bins' spots.
    """
    binned_file = run_dir / "binned_data.npz"
    if not binned_file.exists():
        return None
    fits = sorted(run_dir.glob("rdrbaf_final_nstates*_smp.npz"))
    if len(fits) != 1:
        raise ValueError(f"{run_dir} holds {len(fits)} rdrbaf_final_nstates*_smp.npz, one wanted")
    with np.load(fits[0]) as loaded:
        fit = {key: np.asarray(loaded[key]) for key in FIT_KEYS if key in loaded.files}
    logger.info(f"Reading {binned_file}'s per-spot bins.")
    with np.load(binned_file) as binned:
        barcodes = np.asarray(binned["barcodes"]).astype(str)
        if len(fit["new_assignment"]) != len(barcodes):
            raise ValueError(f"{fits[0]} labels {len(fit['new_assignment']):_} spots, {binned_file} {len(barcodes):_}")
        tracks = CloneTracks(
            np.asarray(binned["lengths"]),
            from_npz(binned, "single_X"),
            from_npz(binned, "single_base_nb_mean"),
            from_npz(binned, "single_total_bb_RD"),
            get_integer_fit(fit, barcodes, labels),
        )
    return tracks


def draw_spatial_frame(ax: Any, box: tuple[float, float, float, float], xlabel: bool) -> float:
    """Frame a spatial panel as `layout.plot_annotation_page` frames its own (`cnamaste.plot.layout.frame_extent`), its
    x axis ending at the spots' right edge, `box`'s x1, its last tick the round step at or below it: every other tick
    labelled (`layout.get_alternate_tick_labels`), y labelled, and x labelled where `xlabel`, else its tick labels left
    off.

    Args:
        ax: the panel.
        box: the spots' and the slide origin's extent, `(x0, x1, y0, y1)` um.
        xlabel: whether the panel labels its x axis.

    Returns:
        The panel's width over its height, its axes at equal aspect.
    """
    frame_extent(ax, box, labels=True)
    start_x, end_x, start_y, end_y, _ = framed(box)
    low, high = ax.get_xlim()
    xs = [float(x) for x in ax.get_xticks() if x <= box[1]]
    ax.set_xlim(low, box[1] + high - end_x)
    ax.set_xticks(xs)
    ax.set_xticklabels(coordinates(xs, max(abs(start_x), abs(end_x), abs(start_y), abs(end_y)) >= 1_000))
    ax.spines["bottom"].set_bounds(xs[0], xs[-1])
    for axis in (ax.xaxis, ax.yaxis):
        found = [t.get_text() for t in axis.get_ticklabels()]
        axis.set_ticklabels(get_alternate_tick_labels(list(axis.get_ticklocs()), found))
    ax.set_ylabel(r"$y$ ($\mu$m)")
    if xlabel:
        ax.set_xlabel(r"$x$ ($\mu$m)")
    else:
        ax.tick_params(labelbottom=False)
    return float(np.diff(ax.get_xlim())[0] / np.diff(ax.get_ylim())[0])


def recolour_neutral_bins(ax: Any) -> None:
    """Recolour the dots and error bars `plot_clones_genomic` draws in `palette.TRACK_NORMAL`, the (1,1) bins', to
    `NEUTRAL_DOT`, each alpha kept.

    Args:
        ax: a track.
    """
    normal, grey = to_rgb(TRACK_NORMAL), to_rgb(NEUTRAL_DOT)
    for collection in ax.collections:
        dots = isinstance(collection, PathCollection)
        colours = np.array(collection.get_facecolors() if dots else collection.get_edgecolors())
        hit = np.isclose(colours[:, :3], normal).all(axis=1) if colours.size else np.zeros(0, dtype=bool)
        if hit.any():
            colours[hit, :3] = grey
            (collection.set_facecolors if dots else collection.set_edgecolors)(colours)


def draw_name_key(ax: Any, name: str, entries: list[tuple[Any, str]]) -> None:
    """Write a clone's name over its first track, left, and its key on the same line, right, as
    `cnamaste.plot.layout.name_line` does: one row, its columns closer where it would meet the name, then as few rows
    as clear the name, read row by row.

    Args:
        ax: the clone's first track.
        name: the clone's name.
        entries: its key's dots and labels, in reading order; empty for no key.
    """
    end = ax.text(0.0, 1.02, name, transform=ax.transAxes, ha="left", va="bottom")
    if not entries:
        return
    renderer = ax.figure.canvas.get_renderer()
    count = len(entries)
    for rows, spacing in ((1, 1.0), (1, 0.6), (1, 0.3), *((r, 1.0) for r in range(2, count + 1))):
        ncol = -(-count // rows)
        # NB matplotlib fills a key's columns first; this order reads it row by row.
        handles, labels = zip(*[entries[i] for c in range(ncol) for i in range(c, count, ncol)], strict=True)
        key = ax.legend(
            list(handles),
            list(labels),
            loc="lower right",
            bbox_to_anchor=(1.0, 1.0),
            ncol=ncol,
            frameon=False,
            handlelength=0.6,
            handletextpad=0.3,
            columnspacing=spacing,
            borderaxespad=0.0,
            borderpad=0.0,
        )
        if key.get_window_extent(renderer).x0 > end.get_window_extent(renderer).x1:
            break


def draw_clone_tracks(
    axes: dict[tuple[str, int], Any],
    tracks: CloneTracks,
    seglevel: pd.DataFrame,
    genome: GenomicAxis,
    normal_fractions: dict[int, float],
) -> dict[tuple[str, int], Any]:
    """Draw each clone's pseudobulk RDR and BAF per bin on its tracks (`plot_clones_genomic`), then furnish them as
    `cnamaste.plot.layout.plot_combined_page` does: each track's end ticks and channel (`track_furniture`), each
    clone's name once, over its first track, a tumour clone's key of its states' shares of its bins right of its name
    (`draw_name_key`), `(A,B)` written as the profiles page writes it. The normal clone's RDR track is removed, its RDR
    1 by construction, and its key with it; the (1,1) bins' dots are grey (`recolour_neutral_bins`); the chromosome
    names `plot_clones_genomic` writes under the last track are removed, the profile under the tracks naming them.

    The black lines are `plot_clones_genomic`'s own expected lines, as run_cnamaste's page draws them: each clone's
    decoded `(A, B)` per bin mixed with normal (1, 1) cells at the clone's normal fraction `f_c`
    (`hmm_ab_decoder.mixture`), RDR `f_c + (1 - f_c)(A + B) / 2` over the clone's `Z_c = exp(new_log_mu_shift_c)`, on
    the dots' scale, their count over the normal baseline; BAF `(f_c + (1 - f_c) B) / (2 f_c + (1 - f_c)(A + B))` and
    its mirror. The normal clone, (1, 1) throughout, reads BAF 1/2 at any `f_c`.

    Args:
        axes: (channel, clone) -> its track, `RDR` then `BAF` per clone in clone order, `BAF` alone where the bins
            carry no read depth; clone 0 the normal clone.
        tracks: the run's bins and fit.
        seglevel: the run's `cnv_<ploidy>_seglevel.tsv`, a row per bin.
        genome: the bins' axis.
        normal_fractions: clone -> the decoder's normal fraction `f_c`, every clone of `axes`.

    Returns:
        The tracks drawn, in `axes`' order, the normal clone's RDR left out.
    """
    plot_clones_genomic(
        tracks.lengths,
        tracks.single_X,
        tracks.single_base_nb_mean,
        tracks.single_total_bb_RD,
        df_cnv=seglevel,
        res_combine=tracks.fit,
        axes=list(axes.values()),
        genome=genome,
        normal_fractions={f"clone {clone}": normal_fractions[clone] for _, clone in axes},
        rdr_ylim=TRACK_RDR_YLIM,
    )
    keys: dict[int, list[tuple[Any, str]]] = {}
    for (channel, clone), ax in axes.items():
        recolour_neutral_bins(ax)
        track_furniture(ax, channel)
        ax.yaxis.labelpad = TRACK_LABEL_PAD
        legend = ax.get_legend()
        if legend is not None:
            for handle in legend.legend_handles:
                handle.set_markersize(TRACK_KEY_MARKER)
            texts = [t.get_text().replace(", ", ",") for t in legend.get_texts()]
            keys[clone] = list(zip(legend.legend_handles, texts, strict=True))
            legend.remove()
        # NB `plot_clones_genomic`'s clone name, upright left of the y axis, and its counts over the track.
        for text in [t for t in ax.texts if t.get_rotation() == 90.0 or t.get_position() == (0.0, 1.02)]:
            text.remove()
    last = list(axes.values())[-1]
    for text in [t for t in last.texts if t.get_text() in {*genome.names, "chr"}]:
        text.remove()
    drawn = {key: ax for key, ax in axes.items() if key != ("RDR", 0)}
    if ("RDR", 0) in axes:
        axes[("RDR", 0)].remove()
    for clone in dict.fromkeys(clone for _, clone in drawn):
        first = next(ax for (_, k), ax in drawn.items() if k == clone)
        draw_name_key(first, cast_clone_label(f"clone {clone}"), [] if clone == 0 else keys.get(clone, []))
    return drawn


def get_track_label_room(figure: Figure, names: list[str], size: float) -> float:
    """Get the inches left of `plot_cnamaste_page`'s right column its axes' labels take, at matplotlib's ticks and
    pads: a track's end tick labels and its channel's name, rotated, at the axis-label size (`track_furniture`), the
    profile's clone names (`layout.draw_copy_profile`), and `chr`, 4 points left of the profile
    (`cnamaste.plot.layout.GenomicAxis.frame`).

    Args:
        figure: the page.
        names: the profile's clone names.
        size: the page's font size.

    Returns:
        The widest, in inches.
    """
    renderer = figure.canvas.get_renderer()
    pad = float(plt.rcParams["ytick.major.pad"]) / 72.0
    ticks = max(get_text_extent(figure, renderer, f"{v:.1f}", size)[0] for v in (0.0, 1.0, TRACK_RDR_YLIM))
    channel = max(get_text_extent(figure, renderer, c, get_axis_label_size(size))[1] for c in ("RDR", "BAF"))
    tick = float(plt.rcParams["ytick.major.size"]) / 72.0
    track = ticks + tick + pad + TRACK_LABEL_PAD / 72.0 + channel
    profile = max(get_text_extent(figure, renderer, name, size)[0] for name in names) + pad
    return max(track, profile, get_text_extent(figure, renderer, "chr", size)[0] + 4.0 / 72.0)


def plot_cnamaste_page(
    physical: pd.DataFrame,
    he_frame: pd.DataFrame | None,
    labels: pd.Series,
    seglevel: pd.DataFrame,
    tracks: CloneTracks,
    normal_fractions: dict[int, float],
    title: str,
    clone_ids: list[int] | None = None,
    in_tissue: pd.Index | None = None,
) -> Figure:
    """Draw cnamaste's own page, `cnamaste.plot.layout.plot_combined_page`, in two columns. Left: (a) the H&E slide,
    (b) the run's clones, keyed inside its top left, one clone a row, each with its share of the in-tissue spots
    (`layout.draw_inside_key`), the spots the key covers logged (`layout.count_covered_spots`). Right: from (c), a
    letter per clone in clone order, its pseudobulk RDR then BAF per bin, the normal clone's BAF alone, each clone
    named once over its first track, its letter on that line, a tumour clone's states keyed right of its name
    (`draw_clone_tracks`); under them, the next letter, its copy-number profile, a row per integer clone, on its bins,
    (1,1) white (`layout.draw_copy_profile`), the chromosomes named under it. (b), the tracks and the profile draw the
    integer clones, as run_cnamaste's page does given `int_groups`: each named by and
    drawn at its first clone, its tracks pooling its spots.

    The spatial panels are `SPATIAL_PANEL_HEIGHT` tall, framed on the run's spots and the slide's origin, the x axis
    ending at the spots' right edge (`draw_spatial_frame`), the clones coloured as `layout.plot_annotation_page`
    colours them (`layout.get_clone_colours`), so a clone is one colour on every page given the same `clone_ids`. The
    right column starts `GAP` right of them past its axes' labels as measured (`get_track_label_room`), a track's
    channel `TRACK_LABEL_PAD` from its tick labels. The tracks fill the left column's height, each at least
    `TRACK_HEIGHT`; the page is as tall as its taller column. Each panel shows its letter alone, its title the letter's
    label (`layout.draw_panel_letter`), each column's letters from the left of its axis labels (`layout.get_label_left`).

    The black lines are each clone's decoded `(A, B)` mixed at its normal fraction, as run_cnamaste's page draws them
    (`draw_clone_tracks`).

    Args:
        physical: `x, y` um, y up, per barcode, the run's spots among them.
        he_frame: the slide, `x, y` um (`cnamaste.plot.layout.draw_slide`), or None for an empty frame.
        labels: the run's clone per barcode, its `clone_labels_int.tsv` at the spots drawn.
        seglevel: the run's `cnv_<ploidy>_seglevel.tsv`, a row per bin of `tracks`, `clone<k> A` and `clone<k> B` per
            clone of its fit, the integer clones' among them.
        tracks: the run's bins and fit.
        normal_fractions: clone -> the decoder's normal fraction, every integer clone's, as run_cnamaste's log
            writes them (`normal fractions per clone`).
        title: the run's name, in its panels' titles.
        clone_ids: clones the palette spreads over beside `labels`', or None.
        in_tissue: the in-tissue spots, or None to count the spots.

    Returns:
        The page, the text column wide.

    Raises:
        ValueError: where `seglevel`'s chromosomes do not run as the bins' do.
    """
    order = [int(k) for k in np.unique(tracks.fit["new_assignment"])]
    genome = GenomicAxis.from_bins(seglevel["CHR"])
    if not np.array_equal(genome.lengths, tracks.lengths):
        raise ValueError(f"seglevel's bins per chromosome {genome.lengths} are not the bins' {tracks.lengths}")
    rdr = SpatioGenomicCounts.of(tracks.single_base_nb_mean).count_positive() > 0
    channels = ("RDR", "BAF") if rdr else ("BAF",)
    xy = physical.loc[labels.index, ["x", "y"]]
    x, y = xy["x"], xy["y"]
    box = (min(x.min(), 0.0), max(x.max(), 0.0), min(y.min(), 0.0), max(y.max(), 0.0))
    names, colours = get_clone_names(labels, get_clone_colours([labels], clone_ids))

    with page_style() as size:
        figure = plt.figure(figsize=page_size("full"), dpi=DPI, facecolor="white")
        page_width = figure.get_figwidth()
        left, panel_h = MARGIN["left"], SPATIAL_PANEL_HEIGHT
        # NB every panel drawn first; its size and place are set once both columns' heights are known.
        he_ax, clone_ax = (add_axes(figure, left, 0.0, 1.0, 1.0) for _ in range(2))
        if he_frame is not None:
            draw_slide(he_ax, he_frame, xy.to_numpy())
        draw_labels(clone_ax, physical, names, colours)
        panel_w = panel_h * [draw_spatial_frame(ax, box, xlabel=ax is clone_ax) for ax in (he_ax, clone_ax)][-1]
        draw_inside_key(clone_ax, names, colours, in_tissue, size)

        column_x = left + panel_w + GAP
        right = column_x + get_track_label_room(figure, [cast_clone_label(f"clone {k}") for k in order], size)
        width = page_width - MARGIN["right"] - right
        profile_ax = add_axes(figure, right, 0.0, width, PROFILE_ROW * len(order))
        draw_copy_profile(profile_ax, seglevel, genome, order)
        track_axes = {(channel, k): add_axes(figure, right, 0.0, width, 1.0) for k in order for channel in channels}
        logger.info(f"Drawing {len(order)} clones' pseudobulk {', '.join(channels)} on {genome.n_bins:_} bins.")
        # NB the keys `name_line` fits beside each clone's name measured at the page's size.
        with plt.rc_context({"font.size": size}):
            track_axes = draw_clone_tracks(track_axes, tracks, seglevel, genome, normal_fractions)
        per_clone = {k: [key for key in track_axes if key[1] == k] for k in order}

        # NB the columns' heights in inches: the left column's sets the tracks', each at least TRACK_HEIGHT.
        renderer = figure.canvas.get_renderer()
        name_lines = {}
        for k in order:
            key = track_axes[per_clone[k][0]].get_legend()
            height = 0.0 if key is None else key.get_window_extent(renderer).height / figure.dpi
            name_lines[k] = max(NAME_LINE, height)
        b_top = MARGIN["top"] + LETTER + panel_h + GAP
        left_height = b_top + LETTER + panel_h + TICK_LINE + AXIS_LABEL_LINE + MARGIN["top"]
        fixed = (
            2.0 * MARGIN["top"]
            + sum(name_lines.values())
            + GAP / 2.0 * (len(order) - 1)
            + TRACK_GAP * sum(len(keys) - 1 for keys in per_clone.values())
            + GAP
            + LETTER
            + PROFILE_ROW * len(order)
            + CHR_NAME_ROOM
        )
        track = max(TRACK_HEIGHT, (left_height - fixed) / len(track_axes))
        page_height = max(left_height, fixed + track * len(track_axes))
        figure.set_size_inches(page_width, page_height)
        logger.info(f"cnamaste page {page_width:.2f} x {page_height:.2f} in, tracks {track:.2f} in.")
        if page_height > page_size("full")[1]:
            logger.warning(f"cnamaste page {page_height:.2f} in tall, past {page_size('full')[1]:.2f} in.")

        place(he_ax, left, MARGIN["top"] + LETTER, panel_w, panel_h)
        place(clone_ax, left, b_top + LETTER, panel_w, panel_h)
        top, clone_tops = MARGIN["top"], {}
        for j, k in enumerate(order):
            top += GAP / 2.0 if j else 0.0
            clone_tops[k] = top
            top += name_lines[k]
            for i, key in enumerate(per_clone[k]):
                top += TRACK_GAP if i else 0.0
                place(track_axes[key], right, top, width, track)
                top += track
        d_top = top + GAP
        place(profile_ax, right, d_top + LETTER, width, PROFILE_ROW * len(order))
        genome.frame(profile_ax, labels=True)
        set_text_sizes(figure, size)
        # NB each column's letters from the left of its axis labels, measured at the page's sizes.
        left_x = min(get_label_left(ax, renderer) for ax in (he_ax, clone_ax))
        right_x = min(get_label_left(ax, renderer) for ax in (profile_ax, *track_axes.values()))
        draw_panel_letter(figure, MARGIN["top"], "a", "H&E", left_x)
        draw_panel_letter(figure, b_top, "b", f"{title}, clones", left_x)
        letters = iter(string.ascii_lowercase[2:])
        for k in order:
            drawn = " and ".join(channel for channel, _ in per_clone[k])
            normal = ", the normal clone" if k == 0 else ""
            clone_title = f"{title}, {cast_clone_label(f'clone {k}')}'s pseudobulk {drawn} per bin{normal}"
            draw_panel_letter(figure, clone_tops[k], next(letters), clone_title, right_x)
        draw_panel_letter(figure, d_top, next(letters), f"{title}, copy-number profile", right_x)
        set_text_sizes(figure, size)
        covered = count_covered_spots(clone_ax, xy.to_numpy())
        if covered:
            logger.warning(f"(b)'s key covers {covered:_} spots.")
        else:
            logger.info("(b)'s key covers no spot.")
    return figure
